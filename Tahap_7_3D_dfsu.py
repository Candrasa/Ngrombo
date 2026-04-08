# =============================================================================
# DFSU Processing Pipeline v8.4 - FIXED
# Perbaikan:
#   [FIX-1] curdir: masking NODATA sebelum arctan2 agar nilai arah valid
#   [FIX-2] dropna: filter berbasis NODATAVALUE bukan NaN untuk NODATA=-9999
#   [FIX-3] Rename field: curspdms->cspdms, curspdkt->cspdkt (DBF-safe 6 char)
#   [FIX-4] Rename field: curU->cur_U, curV->cur_V konsisten (DBF-safe 5 char)
#           Termasuk: waveU->wave_U, waveV->wave_V, windU->wind_U, windV->wind_V
#   [FIX-5] Unicode logging: safe_log() + UTF-8 stream handler
#   [FIX-6] ROLEALIAS, ITEMROLES, ROLES3D, TIFFONLY3D diselaraskan nama baru
# =============================================================================

import mikeio
from mikeio import Grid2D
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.transform import from_bounds
import os, gc, logging, time, threading, queue, io, sys, codecs
from shapely.geometry import Point, box
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# LOGGING SETUP - UTF-8 safe (FIX-5)
# =============================================================================

#logging.basicConfig(level=logging.INFO,
#                    format="%(asctime)s %(levelname)-8s %(message)s")
#log = logging.getLogger(__name__)

_UNICODE_MAP = {
    "\u2713": "[OK]", "\u2717": "[FAIL]", "\u2718": "[FAIL]",
    "\u2192": "->",   "\u00b0": "deg",    "\u2714": "[OK]",
    "\u26a0": "[!]",  "\u274c": "[X]",
}

def safe_log(msg: str) -> str:
    for ch, rep in _UNICODE_MAP.items():
        msg = msg.replace(ch, rep)
    return msg.encode("ascii", errors="replace").decode("ascii")

def setup_logging_safe(log_file="dfsu_pipeline_v84.log"):
    global log
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    try:
        stream = io.TextIOWrapper(sys.stdout.buffer,
                                  encoding="utf-8", errors="replace",
                                  line_buffering=True) \
                 if hasattr(sys.stdout, "buffer") else sys.stdout
    except Exception:
        stream = sys.stdout
    logger.addHandler(logging.StreamHandler(stream)._class_(stream)
                      if False else logging.StreamHandler(stream))
    logger.handlers[-1].setFormatter(fmt)
    if log_file:
        try:
            fh = logging.FileHandler(log_file, encoding="utf-8")
            fh.setFormatter(fmt)
            logger.addHandler(fh)
        except Exception:
            pass
    log = logger   # update global
    return logger

# ============================================================
# STEP 4: Panggil setup — log sudah pasti terdefinisi
# ============================================================
log = setup_logging_safe()

# ============================================================
# STEP 5: Semua kode selanjutnya aman memanggil log.*
# ============================================================
log.info("Pipeline DFSU v8.4 dimulai")  # ← tidak akan NameError

# =============================================================================
# KONFIGURASI
# =============================================================================
DFSUFILES = [
    r"E:\Dukungan Data Osemet\2026\Pulau Gundul 2026 dfsu\3D_Pulau_Gundul.dfsu",
    r"E:\Dukungan Data Osemet\2026\Pulau Gundul 2026 dfsu\SW_Pulau_Gundul.dfsu",
]
BASEOUTPUT    = r"E:\Dukungan Data Osemet\2026\Pulau Gundul 2026 dfsu\IHDC"
OUTPUTGDB     = os.path.join(BASEOUTPUT, "Output.gdb")
OUTPUTTIF     = os.path.join(BASEOUTPUT, "GeoTiff")
BBOX          = (109.4, -6.5, 111.2, -5.5)
DX, DY        = 0.083, 0.083              # ~9 km
DAILYDATES    = ["2026-04-20","2026-04-21","2026-04-22","2026-04-24","2026-04-25","2026-04-26"]
HOURLYDATE    = "2026-04-23"
PROCESS3D     = True
TARGETDEPTHS  = [11.4, 20.4, 29.5]
# Peta label kedalaman: actual -> label singkat (maks 3 char agar nama file aman)
DEPTHLABELMAP = {11.4: "10", 20.4: "20", 29.5: "30"}
DEPTHTOLERANCE= 2.0
MAXWORKERS    = 16
RAMSEMAPHORE  = threading.Semaphore(8)
MS_TO_KNOT    = 1.0 / 0.514444
NODATAVALUE   = -9999.0

# Output flags
OUTPUTMODEGDB     = False
OUTPUTGEOTIFF     = True
OUTPUTPOINTGDB    = False
OUTPUTPOLYGONGDB  = False
OUTPUTRASTERGDB   = False
OUTPUTPOINTSHP    = True
OUTPUTPOLYGONSHP  = False

for d in [OUTPUTTIF]:
    os.makedirs(d, exist_ok=True)

GDBWRITEQUEUE = queue.Queue()
GDBWRITESTOP  = threading.Event()

# =============================================================================
# ROLE MAPPING
# [FIX-4] curU -> cur_U, curV -> cur_V, waveU -> wave_U, dsb.
# =============================================================================
ITEMROLES = {
    "surfelev":  ["Surface elevation", "Water level", "Zeta"],
    "curspeed":  ["Current speed"],
    "curdir":    ["Current direction"],
    "cur_U":     ["U velocity", "u-velocity"],          # FIX-4
    "cur_V":     ["V velocity", "v-velocity"],          # FIX-4
    "wavehs":    ["Sign. Wave Height"],
    "wavedir":   ["Peak Wave Direction", "Mean Wave Direction"],
    "waveper":   ["Peak Wave Period"],
    "windspeed": ["Wind speed"],
    "winddir":   ["Wind direction"],
}

POINTROLES       = ["curspeed","curdir","wavedir","winddir","waveper"]
POLYGONROLES     = ["surfelev","curspeed","wavehs","windspeed"]
TIFFONLY         = ["curspeed","wavehs"]
TIFFONLY3D       = ["cur_U","cur_V"]                    # FIX-4
ROLES3D          = ["curspeed","curdir","cur_U","cur_V"]# FIX-4
ROLESSURFACEONLY = ["surfelev","wavehs","wavedir","waveper","windspeed","winddir"]

# [FIX-3] [FIX-4] ROLEALIAS: nama field DBF maks 10 char, preferred <= 8
ROLEALIAS = {
    "surfelev":  "surfelv",
    "curspeed":  "curspd",
    "curdir":    "curdir",
    "cur_U":     "cur_U",    # 5 char, aman
    "cur_V":     "cur_V",    # 5 char, aman
    "wavehs":    "wavehs",
    "wavedir":   "wavedir",
    "waveper":   "waveper",
    "windspeed": "wndspd",
    "winddir":   "wnddir",
    "wave_U":    "wave_U",   # FIX-4 konsisten
    "wave_V":    "wave_V",
    "wind_U":    "wind_U",
    "wind_V":    "wind_V",
    # [FIX-3] nama kecepatan arus: cspdkt/cspdms (6 char, jauh dari batas 10)
    "cspdms":    "cspdms",   # kecepatan arus m/s
    "cspdkt":    "cspdkt",   # kecepatan arus knot
    "wndspdms":  "wndspdms", # kecepatan angin m/s
    "wndspdkt":  "wndspdkt", # kecepatan angin knot
}

FCPREFIX = {"point":"pt", "polygon":"py", "raster":"rs"}
DIRECTIONCONVENTION = {"curdir":"goingto", "wavedir":"comingfrom", "winddir":"comingfrom"}
CONVENTIONLABEL     = {"goingto":"Menuju ke", "comingfrom":"Dari"}

# =============================================================================
# KELAS KLASIFIKASI
# =============================================================================
CURSPEEDCLASSES = [
    (0.000, 0.051,  "Sgt Lemah"),
    (0.051, 0.257,  "Lemah"),
    (0.257, 0.514,  "Sedang"),
    (0.514, 0.772,  "Kuat"),
    (0.772, 9999.,  "Sgt Kuat"),
]
WAVEHSCLASSES = [
    (0.00, 0.10,   "Tenang"),
    (0.10, 0.50,   "Halus"),
    (0.50, 1.25,   "Sedang"),
    (1.25, 2.50,   "Agak Kasar"),
    (2.50, 4.00,   "Kasar"),
    (4.00, 6.00,   "Sgt Kasar"),
    (6.00, 9999.,  "Tinggi"),
]
SURFELEVCLASSES = [
    (-9999., 0.00, "Bwh MSL"),
    (0.00,  0.10,  "SS1-Glassy"),
    (0.10,  0.50,  "SS2-Rippled"),
    (0.50,  1.25,  "SS3-Wavelets"),
    (1.25,  2.50,  "SS4-Slight"),
    (2.50,  9999., "SS5-Moderate"),
]
WINDSPEEDCLASSES = [
    (0.0,  0.3,   "Calm"),
    (0.3,  1.5,   "Light Air"),
    (1.5,  3.3,   "Lt Breeze"),
    (3.3,  5.5,   "Gntl Breeze"),
    (5.5,  8.0,   "Mod Breeze"),
    (8.0,  10.8,  "Fresh Breeze"),
    (10.8, 13.9,  "Str Breeze"),
    (13.9, 17.2,  "Near Gale"),
    (17.2, 20.7,  "Gale"),
    (20.7, 24.5,  "Str Gale"),
    (24.5, 28.4,  "Storm"),
    (28.4, 32.7,  "Vlt Storm"),
    (32.7, 9999., "Hurricane"),
]
ROLECLASSES = {
    "curspeed":  (CURSPEEDCLASSES,  "klsbmkg", "valuems"),
    "wavehs":    (WAVEHSCLASSES,    "klsbmkg", "valuem"),
    "surfelev":  (SURFELEVCLASSES,  "klswmo",  "valuem"),
    "windspeed": (WINDSPEEDCLASSES, "klsbft",  "valuems"),
}
COMPASS16 = [
    (0.00,   11.25,  "N"),  (11.25,  33.75,  "NNE"),
    (33.75,  56.25,  "NE"), (56.25,  78.75,  "ENE"),
    (78.75,  101.25, "E"),  (101.25, 123.75, "ESE"),
    (123.75, 146.25, "SE"), (146.25, 168.75, "SSE"),
    (168.75, 191.25, "S"),  (191.25, 213.75, "SSW"),
    (213.75, 236.25, "SW"), (236.25, 258.75, "WSW"),
    (258.75, 281.25, "W"),  (281.25, 303.75, "WNW"),
    (303.75, 326.25, "NW"), (326.25, 348.75, "NNW"),
    (348.75, 360.00, "N"),
]

# =============================================================================
# UTILITAS
# =============================================================================
def classify_array(arr, classtable):
    result = np.full(arr.shape, "Unknown", dtype=object)
    for lo, hi, label in classtable:
        result[(arr >= lo) & (arr < hi)] = label
    return result

def direction_to_compass(dir_arr):
    """[FIX-1] Default "?" bukan "N"; NaN dan NODATA dilewati."""
    result  = np.full(dir_arr.shape, "?", dtype=object)
    # Mask: abaikan NODATA (-9999) dan NaN/inf
    valid   = np.isfinite(dir_arr) & (dir_arr != NODATAVALUE)
    dir_mod = np.mod(dir_arr, 360.0)
    for lo, hi, label in COMPASS16:
        mask = valid & (dir_mod >= lo) & (dir_mod < hi)
        result[mask] = label
    return result

def clean_array(arr):
    return np.where(np.isfinite(arr), arr, NODATAVALUE).astype(np.float32)

def safe_label(label):
    return label.replace("-","").replace(",","").replace(".","").replace(" ","")

def build_fc_name(fctype, role, label, depthm=None):
    prefix    = FCPREFIX.get(fctype, "fc")
    roles     = ROLEALIAS.get(role, role)[:8]
    if depthm is not None:
        ldepth  = DEPTHLABELMAP.get(float(depthm), str(round(depthm)))
        depthsfx = f"d{int(ldepth)}m"
    else:
        depthsfx = ""
    lsafe = safe_label(label)
    name  = f"{prefix}{roles}{lsafe}{depthsfx}"
    if len(name) > 64:
        name = name[:64]
        log.warning(safe_log(f"[!] Nama FC dipotong: {name}"))
    if name and name[0].isdigit():
        name = f"fc{name}"[:64]
    return name

def build_grid_manual(bbox=BBOX, dx=DX, dy=DY):
    grid = Grid2D(bbox=bbox, dx=dx, dy=dy, projection="LONG/LAT")
    log.info(safe_log(f"[OK] Grid2D nx={grid.nx}, ny={grid.ny} "
                      f"dx={dx*111000.0:.0f}m dy={dy*110500.0:.0f}m"))
    return grid

def build_polygons(grid):
    xx, yy  = np.meshgrid(grid.x, grid.y)
    dx2, dy2 = grid.dx/2, grid.dy/2
    polys   = [box(xi-dx2, yi-dy2, xi+dx2, yi+dy2)
               for xi, yi in zip(xx.flatten(), yy.flatten())]
    log.info(safe_log(f"[OK] Total poligon: {len(polys)}"))
    return polys, xx, yy

# =============================================================================
# DETEKSI DFSU
# =============================================================================
def detect_dfsu_type(dsopen):
    try:
        nlayers = getattr(dsopen.geometry, "n_layers", None)
        if nlayers is None or nlayers <= 1:
            log.info("Tipe DFSU: 2D")
            return "2D", 1, None, None
        zvals     = dsopen.geometry.element_coordinates[:, 2]
        layerids  = dsopen.geometry.layer_ids
        unique_lids = np.unique(layerids)
        zunique   = np.array([zvals[layerids == lid].mean() for lid in unique_lids])
        log.info(safe_log(f"[OK] Tipe DFSU 3D {nlayers} layer"))
        log.info(f"Layer IDs unik: {unique_lids.tolist()}")
        log.info(f"Kedalaman: {np.abs(zunique).round(1).tolist()} m")
        log.info(f"Total elemen 3D: {len(layerids)}, "
                 f"Elemen 2D/layer: {(layerids == unique_lids[0]).sum()}")
        return "3D", nlayers, zunique, unique_lids
    except AttributeError as e:
        log.info(f"Tipe DFSU 2D (fallback): {e}")
        return "2D", 1, None, None

def find_nearest_layer(zunique, unique_lids, target_depthm):
    targetz = -abs(target_depthm)
    diffs   = np.abs(zunique - targetz)
    idx     = int(np.argmin(diffs))
    actualz = zunique[idx]
    layerid = unique_lids[idx]
    diffm   = abs(actualz - targetz)
    if diffm > DEPTHTOLERANCE:
        log.warning(safe_log(
            f"[!] Layer -{target_depthm}m -> z={actualz:.1f}m "
            f"selisih={diffm:.1f}m > toleransi={DEPTHTOLERANCE}m"))
    log.debug(f"Layer target -{target_depthm}m -> layerid={layerid}, z={actualz:.2f}m")
    return int(layerid), actualz

def detect_roles(dfsufile):
    dsopen    = mikeio.open(dfsufile)
    available = {item.name for item in dsopen.items}
    detected  = {}
    for role, candidates in ITEMROLES.items():
        for candidate in candidates:
            if candidate in available:
                detected[role] = candidate
                break
    log.info(f"FILE: {os.path.basename(dfsufile)} | Terdeteksi: {len(detected)} role")
    for role, name in detected.items():
        log.info(f"  {role:12s} -> {name}")
    return detected, available

# =============================================================================
# EKSTRAKSI LAYER 3D (perbaikan IndexError v8.3 dipertahankan)
# =============================================================================
def extract_layer_from_3d(ds3dbbox, dsopengeom, bboxmask, layerid, target_depthm):
    try:
        all_layerids   = dsopengeom.layer_ids
        unique_lids_all = np.unique(all_layerids)
        if layerid not in unique_lids_all:
            log.error(f"layerid={layerid} tidak ada. Available: {unique_lids_all.tolist()}")
            return None, None, None
        layermask_full = (all_layerids == layerid)
        combined_mask  = bboxmask & layermask_full
        n_layer_in_bbox = combined_mask.sum()
        if n_layer_in_bbox == 0:
            log.warning(f"Tidak ada elemen layer {layerid} dalam bbox. Skip.")
            return None, None, None
        all_ex = dsopengeom.element_coordinates[:, 0]
        all_ey = dsopengeom.element_coordinates[:, 1]
        elemxl = all_ex[combined_mask]
        elemyl = all_ey[combined_mask]
        bbox_indices      = np.where(bboxmask)[0]
        layerinbboxpos    = np.where(combined_mask[bboxmask])[0]
        nbboxelem         = len(bbox_indices)
        if len(layerinbboxpos) == 0:
            log.warning(f"layerinbboxpos kosong untuk layer {layerid}.")
            return None, None, None
        maxpos = layerinbboxpos.max()
        if maxpos >= nbboxelem:
            log.error(safe_log(
                f"[X] IndexError sebelum isel: maxpos={maxpos} >= nbboxelem={nbboxelem}"))
            return None, None, None
        dslayer2d = ds3dbbox.isel(element=layerinbboxpos.tolist())
        log.debug(f"dslayer2d shape: element={len(layerinbboxpos)}")
        return dslayer2d, elemxl, elemyl
    except IndexError as e:
        log.error(f"IndexError extract layer {layerid}: {e}")
        return None, None, None
    except Exception as e:
        log.error(f"Error extract layer {layerid}: {e}", exc_info=True)
        return None, None, None

# =============================================================================
# INTERPOLASI & VEKTOR
# =============================================================================
def fast_interp_to_grid(elemx, elemy, elemvalues, gridxx, gridyy):
    if len(elemx) == 0 or len(elemvalues) == 0:
        log.warning("fast_interp_to_grid: input kosong -> NODATA")
        return np.full(gridxx.shape, NODATAVALUE).flatten()
    points  = np.column_stack((elemx, elemy))
    gridpts = np.column_stack((gridxx.flatten(), gridyy.flatten()))
    result  = LinearNDInterpolator(points, elemvalues)(gridpts)
    nanmask = np.isnan(result)
    if nanmask.any():
        result[nanmask] = NearestNDInterpolator(points, elemvalues)(gridpts[nanmask])
    return result

def compute_uv_components(direction_arr, magnitude_arr, convention):
    rad = np.deg2rad(direction_arr)
    if convention == "goingto":
        return magnitude_arr * np.sin(rad), magnitude_arr * np.cos(rad)
    else:
        return -magnitude_arr * np.sin(rad), -magnitude_arr * np.cos(rad)

def ds_to_role_arrays(dsslice, elemx, elemy, gridxx, gridyy, detected_roles):
    """
    Interpolasi semua role ke grid reguler.
    [FIX-1] Masking NODATA sebelum arctan2 agar curdir valid.
    [FIX-3] Tambah cspdms, cspdkt (menggantikan curspdms/curspdkt)
    [FIX-4] Key rolearrays menggunakan cur_U, cur_V (konsisten dengan ROLEALIAS)
    """
    rolearrays = {}
    nelem = len(elemx)
    log.debug(f"ds_to_role_arrays: {nelem} elemen input")

    for role, itemname in detected_roles.items():
        try:
            raw  = dsslice[itemname].to_numpy()
            vals = raw.flatten()
            if len(vals) != nelem:
                log.warning(f"Dimensi mismatch: role={role} vals={len(vals)}, elemx={nelem}. Skip.")
                continue
            if not np.any(np.isfinite(vals)):
                log.warning(f"{role}: semua nilai NaN/inf. Skip.")
                continue
            rolearrays[role] = fast_interp_to_grid(elemx, elemy, vals, gridxx, gridyy)
        except KeyError:
            log.warning(f"Item {itemname!r} tidak ada di dsslice. Skip role {role}.")
        except Exception as e:
            log.warning(f"Skip {role}: {e}")

    # --- Hitung U/V dari curdir jika tersedia (goingto) ---
    if "curdir" in rolearrays:
        mag = rolearrays.get("curspeed", np.ones_like(rolearrays["curdir"]))
        rolearrays["cur_U"], rolearrays["cur_V"] =             compute_uv_components(rolearrays["curdir"], mag, "goingto")

    # --- [FIX-1] Hitung curdir dari cur_U/cur_V jika curdir tidak ada ---
    if "curdir" not in rolearrays and "cur_U" in rolearrays and "cur_V" in rolearrays:
        u = rolearrays["cur_U"].copy()
        v = rolearrays["cur_V"].copy()
        # Masking titik NODATA: set NaN sebelum arctan2
        nodata_mask = (np.abs(u - NODATAVALUE) < 1.0) | (np.abs(v - NODATAVALUE) < 1.0)
        u[nodata_mask] = np.nan
        v[nodata_mask] = np.nan
        dir_rad = np.arctan2(u, v)           # atan2(east, north) -> goingto
        dir_deg = np.degrees(dir_rad) % 360  # normalisasi 0-360
        dir_deg[nodata_mask] = np.nan        # tetap NaN di NODATA
        rolearrays["curdir"] = dir_deg
        log.info("curdir dihitung dari arctan2(cur_U, cur_V) [going-to]")
        # Diagnostik
        valid_count = np.sum(np.isfinite(dir_deg))
        log.info(safe_log(f"[OK] curdir valid: {valid_count}/{len(dir_deg)} titik"))

    # --- U/V gelombang ---
    if "wavedir" in rolearrays:
        mag = rolearrays.get("wavehs", np.ones_like(rolearrays["wavedir"]))
        rolearrays["wave_U"], rolearrays["wave_V"] =             compute_uv_components(rolearrays["wavedir"], mag, "comingfrom")

    # --- U/V angin ---
    if "winddir" in rolearrays:
        mag = rolearrays.get("windspeed", np.ones_like(rolearrays["winddir"]))
        rolearrays["wind_U"], rolearrays["wind_V"] =             compute_uv_components(rolearrays["winddir"], mag, "comingfrom")

    # --- [FIX-3] Kecepatan dalam satuan baru ---
    if "curspeed" in rolearrays:
        rolearrays["cspdkt"] = rolearrays["curspeed"] * MS_TO_KNOT  # knot
        rolearrays["cspdms"] = rolearrays["curspeed"]               # m/s alias

    if "windspeed" in rolearrays:
        rolearrays["wndspdkt"] = rolearrays["windspeed"] * MS_TO_KNOT

    log.debug(f"rolearrays keys: {list(rolearrays.keys())}")
    return rolearrays

# =============================================================================
# SIMPAN SHAPEFILE ARUS (COMBINED) - v8.4
# Kolom: curdir, dirlabel, konvensi, cspdms, cspdkt, cur_U, cur_V, klsbmkg,
#        depthm, lon, lat, datetime
# [FIX-1] curdir valid (NaN di titik NODATA)
# [FIX-2] dropna berbasis NODATAVALUE threshold, bukan hanya NaN
# [FIX-3] cspdms, cspdkt (6 char, aman DBF)
# [FIX-4] cur_U, cur_V (5 char, aman DBF)
# =============================================================================
def save_current_shp_combined(rolearrays, gridxx, gridyy, label, depthm=None):
    """
    Schema atribut SHP arus:
    Field      | Tipe   | Deskripsi
    -----------|--------|-----------------------------
    curdir     | DOUBLE | Arah arus derajat (1 desimal) goingto
    dirlabel   | TEXT   | Label kompas 16 arah (N, NNE, ... )
    konvensi   | TEXT   | "Menuju ke"
    cspdms     | DOUBLE | Kecepatan arus m/s (4 desimal)   [FIX-3]
    cspdkt     | DOUBLE | Kecepatan arus knot (1 desimal)  [FIX-3]
    cur_U      | DOUBLE | Komponen U vektor arus           [FIX-4]
    cur_V      | DOUBLE | Komponen V vektor arus           [FIX-4]
    klsbmkg    | TEXT   | Klasifikasi BMKG
    depthm     | DOUBLE | Kedalaman layer (0=surface)
    lon        | DOUBLE | Bujur
    lat        | DOUBLE | Lintang
    datetime   | TEXT   | Label waktu
    """
    curkeys = ["curdir","curspeed","cur_U","cur_V"]
    if not any(k in rolearrays for k in curkeys):
        return
    try:
        depthval = float(depthm) if depthm is not None else 0.0
        pts = [Point(float(x), float(y))
               for x, y in zip(gridxx.flatten(), gridyy.flatten())]
        n = len(pts)

        # --- Kolom dasar ---
        data = {
            "lon":      gridxx.flatten(),
            "lat":      gridyy.flatten(),
            "datetime": np.full(n, label, dtype=object),
            "depthm":   np.full(n, depthval),
        }

        # --- [FIX-1] Arah arus (sudah NaN di NODATA) ---
        if "curdir" in rolearrays:
            dirarr = rolearrays["curdir"].copy()
            data["curdir"]   = np.round(dirarr, 1)    # 1 desimal derajat
            data["dirlabel"] = direction_to_compass(dirarr)
            data["konvensi"] = np.full(n, "Menuju ke", dtype=object)
            # Diagnostik: cek apakah dirlabel konstan (indikasi bug)
            unique_lbl = set(data["dirlabel"][data["dirlabel"] != "?"])
            if len(unique_lbl) == 1:
                log.warning(safe_log(
                    f"[!] dirlabel mungkin konstan: {unique_lbl} | "
                    f"label={label} depth={depthm}"))

        # --- [FIX-3] Kecepatan: cspdms & cspdkt ---
        if "curspeed" in rolearrays:
            spd = rolearrays["curspeed"]
            spd_masked = np.where(np.abs(spd - NODATAVALUE) < 1.0, np.nan, spd)
            data["cspdms"]  = np.round(spd_masked, 4)
            data["cspdkt"]  = np.round(spd_masked * MS_TO_KNOT, 1)
            data["klsbmkg"] = classify_array(spd, CURSPEEDCLASSES)

        # --- [FIX-4] Komponen vektor: cur_U, cur_V ---
        if "cur_U" in rolearrays:
            u = rolearrays["cur_U"].copy()
            u[np.abs(u - NODATAVALUE) < 1.0] = np.nan
            data["cur_U"] = np.round(u, 4)
        if "cur_V" in rolearrays:
            v = rolearrays["cur_V"].copy()
            v[np.abs(v - NODATAVALUE) < 1.0] = np.nan
            data["cur_V"] = np.round(v, 4)

        # --- Nama file ---
        lsafe   = safe_label(label)
        ldepth  = DEPTHLABELMAP.get(float(depthm), str(round(depthm))) if depthm else None
        depthsfx = f"d{int(ldepth)}m" if ldepth else ""
        fcname  = f"ptcur{lsafe}{depthsfx}"

        gdf = gpd.GeoDataFrame(data, geometry=pts, crs="EPSG:4326")

        # [FIX-2] Filter: buang baris dimana cspdms NaN DAN curdir NaN
        # (cara bersih: nilai -9999 sudah dikonversi ke NaN sebelumnya)
        filter_cols = [c for c in ["curdir","cspdms"] if c in data]
        if filter_cols:
            gdf.dropna(subset=filter_cols, how="all", inplace=True)

        if not gdf.empty:
            shppath = os.path.join(FALLBACKPNT, f"{fcname}.shp")
            gdf.to_file(shppath)
            log.info(safe_log(
                f"[OK] Current SHP {fcname}.shp {len(gdf)} records | "
                f"curdir_valid={gdf['curdir'].notna().sum() if 'curdir' in gdf else 0}"))
        del gdf; gc.collect()
    except Exception as e:
        log.error(f"ERROR save_current_shp_combined {label}: {e}", exc_info=True)

# =============================================================================
# SIMPAN SHAPEFILE GELOMBANG (COMBINED)
# [FIX-4] wave_U, wave_V (menggantikan waveU, waveV)
# =============================================================================
def save_wave_shp_combined(rolearrays, gridxx, gridyy, label, depthm=None):
    wavekeys = ["wavedir", "wavehs", "waveper"]
    if not any(k in rolearrays for k in wavekeys):
        return
    try:
        depthval = float(depthm) if depthm is not None else 0.0
        pts = [Point(float(x), float(y))
               for x, y in zip(gridxx.flatten(), gridyy.flatten())]

        # -------------------------------------------------------
        # PERBAIKAN: cek key eksplisit, bukan pakai 'or' pada array
        # -------------------------------------------------------
        if "wavedir" in rolearrays:
            anchor_arr = rolearrays["wavedir"]
        elif "wavehs" in rolearrays:
            anchor_arr = rolearrays["wavehs"]
        elif "waveper" in rolearrays:
            anchor_arr = rolearrays["waveper"]
        else:
            log.warning(f"[WARN] save_wave_shp: tidak ada anchor array. Skip {label}")
            return

        n = len(anchor_arr)
        data = {
            "lon":      gridxx.flatten(),
            "lat":      gridyy.flatten(),
            "datetime": np.full(n, label, dtype=object),
            "depthm":   np.full(n, depthval),
        }

        if "wavedir" in rolearrays:
            dirarr = rolearrays["wavedir"]
            data["wavedir"]  = np.round(dirarr, 1)
            data["dirlabel"] = direction_to_compass(dirarr)
            data["konvensi"] = np.full(n, "Dari", dtype=object)

        if "wavehs" in rolearrays:
            hs = rolearrays["wavehs"]
            hs_masked = np.where(np.abs(hs - NODATAVALUE) < 1.0, np.nan, hs)
            data["wavehs"]  = np.round(hs_masked, 2)
            data["klsbmkg"] = classify_array(hs, WAVEHSCLASSES)

        if "waveper" in rolearrays:
            data["waveper"] = np.round(rolearrays["waveper"], 1)

        if "wave_U" in rolearrays:
            data["wave_U"]  = np.round(rolearrays["wave_U"], 4)

        if "wave_V" in rolearrays:
            data["wave_V"]  = np.round(rolearrays["wave_V"], 4)

        lsafe    = safe_label(label)
        ldepth   = DEPTHLABELMAP.get(float(depthm), str(round(depthm))) \
                   if depthm is not None else None
        depthsfx = f"d{int(ldepth)}m" if ldepth else ""
        fcname   = f"ptwave{lsafe}{depthsfx}"

        gdf = gpd.GeoDataFrame(data, geometry=pts, crs="EPSG:4326")

        wavecols = [c for c in ["wavedir", "wavehs", "waveper"] if c in data]
        if wavecols:
            gdf.dropna(subset=wavecols, how="all", inplace=True)

        if not gdf.empty:
            shppath = os.path.join(FALLBACKPNT, f"{fcname}.shp")
            gdf.to_file(shppath)
            log.info(safe_log(f"[OK] Wave SHP {fcname}.shp {len(gdf)} records"))
        else:
            log.warning(f"[WARN] Wave SHP kosong setelah dropna: {label}")

        del gdf; gc.collect()

    except Exception as e:
        log.error(f"ERROR save_wave_shp_combined {label}: {e}", exc_info=True)

# =============================================================================
# RASTER (GeoTiff / GDB)
# =============================================================================
def write_raster(rolearrays, grid, role, label, gdbpath, depthm):
    if role not in rolearrays:
        return
    arr2d   = rolearrays[role].reshape(grid.ny, grid.nx).astype(np.float32)
    arrclean = clean_array(np.flipud(arr2d))
    transform = from_bounds(BBOX[0], BBOX[1], BBOX[2], BBOX[3], grid.nx, grid.ny)
    fcname  = build_fc_name("raster", role, label, depthm)
    if OUTPUTGEOTIFF:
        tifpath = os.path.join(OUTPUTTIF, f"{fcname}.tif")
        try:
            with rasterio.open(
                tifpath, "w", driver="GTiff",
                height=grid.ny, width=grid.nx, count=1, dtype="float32",
                crs="EPSG:4326", transform=transform,
                compress="lzw", nodata=NODATAVALUE
            ) as dst:
                dst.write(arrclean, 1)
                dst.update_tags(role=role, label=label,
                               depthm=str(depthm) if depthm else "surface",
                               nodata=str(NODATAVALUE), source="MIKE 21 FM")
            log.info(safe_log(f"[OK] GeoTiff {fcname}.tif"))
        except Exception as e:
            log.error(f"ERROR GeoTiff {fcname}: {e}")

# =============================================================================
# OUTPUT PATH
# =============================================================================
FALLBACKPNT  = os.path.join(BASEOUTPUT, "FallbackSHP", "point")
FALLBACKPOLY = os.path.join(BASEOUTPUT, "FallbackSHP", "poly")
for d in [FALLBACKPNT, FALLBACKPOLY]:
    os.makedirs(d, exist_ok=True)

# =============================================================================
# ENSURE GDB
# =============================================================================
def ensure_gdb(gdbpath):
    if not OUTPUTMODEGDB:
        log.info("OUTPUTMODEGDB=False, GDB dilewati.")
        return False
    try:
        import arcpy
        arcpy.env.overwriteOutput = True
        gdbfolder = os.path.dirname(gdbpath)
        gdbname   = os.path.basename(gdbpath).replace(".gdb","")
        os.makedirs(gdbfolder, exist_ok=True)
        if not arcpy.Exists(gdbpath):
            arcpy.management.CreateFileGDB(gdbfolder, gdbname)
            log.info(safe_log(f"[OK] GDB dibuat: {gdbpath}"))
        sr = arcpy.SpatialReference(4326)
        for fdname in ["HD","Raster"]:
            fdpath = os.path.join(gdbpath, fdname)
            if not arcpy.Exists(fdpath):
                arcpy.management.CreateFeatureDataset(gdbpath, fdname, sr)
        log.info("GDB siap.")
        return True
    except ImportError:
        log.warning("ArcPy tidak tersedia.")
        return False
    except Exception as e:
        log.error(f"ERROR ensure_gdb: {e}")
        return False

# =============================================================================
# WRITE POINT FC (ArcPy GDB)
# [FIX-3] cspdms, cspdkt  [FIX-4] cur_U, cur_V
# =============================================================================
def write_point_fc(rolearrays, gridxx, gridyy, pointrole, label, gdbpath, depthm):
    if not OUTPUTPOINTGDB or pointrole not in rolearrays:
        return
    try:
        import arcpy
        arcpy.env.overwriteOutput = True
        dirarr   = rolearrays[pointrole]
        convention = DIRECTIONCONVENTION.get(pointrole, "goingto")
        convlabel  = CONVENTIONLABEL[convention]
        compass    = direction_to_compass(dirarr)
        depthval   = float(depthm) if depthm is not None else 0.0
        fcname     = build_fc_name("point", pointrole, label, depthm)
        fcpath     = os.path.join(gdbpath, "HD", fcname)
        sr         = arcpy.SpatialReference(4326)
        if arcpy.Exists(fcpath):
            arcpy.management.Delete(fcpath)
        arcpy.management.CreateFeatureclass(
            os.path.join(gdbpath,"HD"), fcname, "POINT", spatial_reference=sr)

        if pointrole == "curdir":
            # [FIX-3] cspdms/cspdkt  [FIX-4] cur_U/cur_V
            fieldsdef = {
                "curdir":   ("DOUBLE", rolearrays.get("curdir")),
                "dirlabel": ("TEXT",   compass),
                "konvensi": ("TEXT",   np.full(dirarr.shape, convlabel, dtype=object)),
                "cspdms":   ("DOUBLE", rolearrays.get("cspdms")),   # FIX-3
                "cspdkt":   ("DOUBLE", rolearrays.get("cspdkt")),   # FIX-3
                "cur_U":    ("DOUBLE", rolearrays.get("cur_U")),    # FIX-4
                "cur_V":    ("DOUBLE", rolearrays.get("cur_V")),    # FIX-4
                "klsbmkg":  ("TEXT",   classify_array(
                    rolearrays.get("curspeed", np.zeros_like(dirarr)), CURSPEEDCLASSES)),
                "depthm":   ("DOUBLE", np.full(dirarr.shape, depthval)),
            }
        elif pointrole == "wavedir":
            hs = rolearrays.get("wavehs", np.full(dirarr.shape, np.nan))
            fieldsdef = {
                "wavedir":  ("DOUBLE", rolearrays.get("wavedir")),
                "dirlabel": ("TEXT",   compass),
                "konvensi": ("TEXT",   np.full(dirarr.shape, convlabel, dtype=object)),
                "wave_U":   ("DOUBLE", rolearrays.get("wave_U")),   # FIX-4
                "wave_V":   ("DOUBLE", rolearrays.get("wave_V")),   # FIX-4
                "wavehs":   ("DOUBLE", hs),
                "waveper":  ("DOUBLE", rolearrays.get("waveper")),
                "klsbmkg":  ("TEXT",   classify_array(hs, WAVEHSCLASSES)),
                "depthm":   ("DOUBLE", np.full(dirarr.shape, depthval)),
            }
        elif pointrole == "winddir":
            ws = rolearrays.get("windspeed", np.full(dirarr.shape, np.nan))
            fieldsdef = {
                "winddir":  ("DOUBLE", rolearrays.get("winddir")),
                "dirlabel": ("TEXT",   compass),
                "konvensi": ("TEXT",   np.full(dirarr.shape, convlabel, dtype=object)),
                "wind_U":   ("DOUBLE", rolearrays.get("wind_U")),   # FIX-4
                "wind_V":   ("DOUBLE", rolearrays.get("wind_V")),   # FIX-4
                "wndspdms": ("DOUBLE", ws),
                "wndspdkt": ("DOUBLE", rolearrays.get("wndspdkt")),
                "klsbft":   ("TEXT",   classify_array(ws, WINDSPEEDCLASSES)),
                "depthm":   ("DOUBLE", np.full(dirarr.shape, depthval)),
            }
        else:
            return

        for fname, (ftype, _) in fieldsdef.items():
            fn10 = fname[:10]
            if ftype == "TEXT":
                arcpy.management.AddField(fcpath, fn10, "TEXT", field_length=20)
            else:
                arcpy.management.AddField(fcpath, fn10, "DOUBLE")
        arcpy.management.AddField(fcpath, "lon",      "DOUBLE")
        arcpy.management.AddField(fcpath, "lat",      "DOUBLE")
        arcpy.management.AddField(fcpath, "datetime", "TEXT", field_length=30)

        fieldnames = ["SHAPE@XY","lon","lat","datetime"] +                      [f[:10] for f in fieldsdef]
        arrays_list = list(fieldsdef.values())
        ptsx = gridxx.flatten(); ptsy = gridyy.flatten()

        with arcpy.da.InsertCursor(fcpath, fieldnames) as cursor:
            for i in range(len(dirarr)):
                if not np.isfinite(dirarr[i]):
                    continue
                lon = float(ptsx[i]); lat = float(ptsy[i])
                row = [(lon, lat), lon, lat, label]
                for _, arr in arrays_list:
                    if arr is None:
                        row.append(None)
                    elif isinstance(arr[i], (str, np.str_)):
                        row.append(str(arr[i]))
                    else:
                        v = float(arr[i])
                        row.append(None if not np.isfinite(v) else v)
                cursor.insertRow(row)
        log.info(safe_log(f"[OK] Point GDB {convlabel} {fcname}"))
    except ImportError:
        if OUTPUTPOINTSHP:
            save_current_shp_combined(rolearrays, gridxx, gridyy, label, depthm)
    except Exception as e:
        log.error(f"ERROR write_point_fc {pointrole}/{label}: {e}", exc_info=True)

# =============================================================================
# GDB WRITER THREAD
# =============================================================================
def gdb_writer_thread(gdbpath):
    log.info("GDB Writer Thread dimulai.")
    while not GDBWRITESTOP.is_set() or not GDBWRITEQUEUE.empty():
        try:
            task = GDBWRITEQUEUE.get(timeout=2.0)
        except queue.Empty:
            continue
        try:
            t = task["type"]
            if t == "point":
                write_point_fc(task["rolearrays"], task["gridxx"], task["gridyy"],
                               task["pointrole"], task["label"], gdbpath, task["depthm"])
            elif t == "raster":
                write_raster(task["rolearrays"], task["grid"],
                             task["role"], task["label"], gdbpath, task["depthm"])
        except Exception as e:
            log.error(f"GDB Writer ERROR: {e}", exc_info=True)
        finally:
            GDBWRITEQUEUE.task_done()
    log.info("GDB Writer Thread selesai.")

# =============================================================================
# PROCESS ONE SLICE
# =============================================================================
def process_one_slice(dsslice, elemx, elemy, gridxx, gridyy, grid, polys,
                      detected_roles, label, depthm=None, activeroles=None):
    try:
        if activeroles is None:
            activeroles = detected_roles
        if len(elemx) == 0:
            log.warning(f"process_one_slice {label}/{depthm}: elemx kosong. Skip.")
            return
        rolearrays = ds_to_role_arrays(dsslice, elemx, elemy, gridxx, gridyy, activeroles)
        if not rolearrays:
            log.warning(f"process_one_slice {label}/{depthm}: rolearrays kosong. Skip.")
            return
        racopy = {k: v.copy() for k, v in rolearrays.items()}

        # --- Point SHP (output utama) ---
        if OUTPUTMODEGDB and OUTPUTPOINTGDB:
            for pr in POINTROLES:
                if pr in activeroles:
                    GDBWRITEQUEUE.put({"type":"point","rolearrays":racopy,
                                       "gridxx":gridxx,"gridyy":gridyy,
                                       "pointrole":pr,"label":label,"depthm":depthm})
        if OUTPUTPOINTSHP:
            save_current_shp_combined(racopy, gridxx, gridyy, label, depthm)
            save_wave_shp_combined(racopy, gridxx, gridyy, label, depthm)

        # --- Raster/GeoTiff ---
        raster_roles = (
            [r for r in POLYGONROLES  if r in activeroles] +
            [r for r in TIFFONLY      if r in activeroles] +
            [r for r in TIFFONLY3D    if r in activeroles]  # cur_U, cur_V
        )
        if OUTPUTGEOTIFF or (OUTPUTMODEGDB and OUTPUTRASTERGDB):
            for role in raster_roles:
                GDBWRITEQUEUE.put({"type":"raster","rolearrays":racopy,
                                   "grid":grid,"role":role,"label":label,"depthm":depthm})
        del rolearrays, racopy; gc.collect()
    except Exception as e:
        log.error(f"ERROR process_one_slice {label}/{depthm}: {e}", exc_info=True)

# =============================================================================
# PROCESS FULL SLICE (Surface + 3D layers)
# =============================================================================
def process_full_slice(dssliceall, dsopengeom, elemxsurf, elemysurf,
                       gridxx, gridyy, grid, polys, detected_roles, label,
                       dfsutype, zunique, unique_lids, bboxmask):
    log.info(f"Surface {label}")
    process_one_slice(dssliceall, elemxsurf, elemysurf, gridxx, gridyy,
                      grid, polys, detected_roles, label,
                      depthm=None, activeroles=detected_roles)

    if dfsutype == "3D" and PROCESS3D and zunique is not None:
        roles3d = {role: item for role, item in detected_roles.items()
                   if role in ROLES3D}
        if not roles3d:
            return
        for target_depth in TARGETDEPTHS:
            log.info(f"Layer -{target_depth}m {label}")
            try:
                layerid, actualz = find_nearest_layer(zunique, unique_lids, target_depth)
                dslayer, elemxl, elemyl = extract_layer_from_3d(
                    dssliceall, dsopengeom, bboxmask, layerid, target_depth)
                if dslayer is None:
                    log.error(safe_log(
                        f"[X] Gagal extract layer -{target_depth}m. Skip."))
                    continue
                n_elem_layer = len(elemxl)
                log.debug(f"Layer -{target_depth}m: {n_elem_layer} elemen, "
                          f"z={actualz:.2f}m, layerid={layerid}")
                process_one_slice(dslayer, elemxl, elemyl, gridxx, gridyy,
                                  grid, polys, roles3d, label,
                                  depthm=float(target_depth), activeroles=roles3d)
                del dslayer; gc.collect()
            except Exception as e:
                log.error(safe_log(
                    f"[X] Gagal layer -{target_depth}m {label}: {e}"), exc_info=True)
                continue

# =============================================================================
# WORKER HARIAN / JAM
# =============================================================================
def process_daily_worker(datestr, dsall, dsopengeom, elemxsurf, elemysurf,
                         gridxx, gridyy, grid, polys, detected_roles,
                         dfsutype, zunique, unique_lids, bboxmask):
    with RAMSEMAPHORE:
        try:
            log.info(f"Daily {datestr}")
            dsday = dsall.sel(time=datestr)
            if len(dsday.time) == 0:
                log.warning(f"Skip {datestr}: tidak ada data.")
                return
            dsavg = dsday.mean(axis=0)
            label = datestr.replace("-","")
            process_full_slice(dsavg, dsopengeom, elemxsurf, elemysurf,
                               gridxx, gridyy, grid, polys, detected_roles, label,
                               dfsutype, zunique, unique_lids, bboxmask)
            del dsday, dsavg; gc.collect()
        except Exception as e:
            log.error(f"ERROR Daily {datestr}: {e}", exc_info=True)

def process_hourly_worker(timeidx, dsday, dsopengeom, elemxsurf, elemysurf,
                          gridxx, gridyy, grid, polys, detected_roles,
                          dfsutype, zunique, unique_lids, bboxmask):
    try:
        dshour = dsday.isel(time=timeidx)
        ts     = pd.Timestamp(dsday.time[timeidx])
        label  = ts.strftime("%Y-%m-%d %H%M%S")
        process_full_slice(dshour, dsopengeom, elemxsurf, elemysurf,
                           gridxx, gridyy, grid, polys, detected_roles, label,
                           dfsutype, zunique, unique_lids, bboxmask)
        del dshour; gc.collect()
    except Exception as e:
        log.error(f"ERROR Hourly idx={timeidx}: {e}", exc_info=True)

def run_daily_parallel(dsall, dsopengeom, elemxsurf, elemysurf,
                       gridxx, gridyy, grid, polys, detected_roles,
                       dfsutype, zunique, unique_lids, bboxmask):
    log.info("=" * 60)
    log.info(f"PROSES HARIAN Workers={MAXWORKERS} Tipe={dfsutype}")
    log.info("=" * 60)
    wfn = partial(process_daily_worker,
                  dsall=dsall, dsopengeom=dsopengeom,
                  elemxsurf=elemxsurf, elemysurf=elemysurf,
                  gridxx=gridxx, gridyy=gridyy, grid=grid, polys=polys,
                  detected_roles=detected_roles, dfsutype=dfsutype,
                  zunique=zunique, unique_lids=unique_lids, bboxmask=bboxmask)
    with ThreadPoolExecutor(max_workers=MAXWORKERS) as ex:
        futures = {ex.submit(wfn, d): d for d in DAILYDATES}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as e:
                log.error(f"GAGAL Daily {futures[fut]}: {e}")
    log.info("Proses harian selesai.")

def run_hourly_parallel(dshourly, dsopengeom, elemxsurf, elemysurf,
                        gridxx, gridyy, grid, polys, detected_roles,
                        dfsutype, zunique, unique_lids, bboxmask):
    log.info("=" * 60)
    log.info(f"PROSES PER JAM {HOURLYDATE} Workers={MAXWORKERS}")
    log.info("=" * 60)
    nsteps = len(dshourly.time)
    log.info(f"Total timestep: {nsteps}")
    wfn = partial(process_hourly_worker,
                  dsday=dshourly, dsopengeom=dsopengeom,
                  elemxsurf=elemxsurf, elemysurf=elemysurf,
                  gridxx=gridxx, gridyy=gridyy, grid=grid, polys=polys,
                  detected_roles=detected_roles, dfsutype=dfsutype,
                  zunique=zunique, unique_lids=unique_lids, bboxmask=bboxmask)
    with ThreadPoolExecutor(max_workers=MAXWORKERS) as ex:
        list(ex.map(wfn, range(nsteps)))
    log.info("Proses per jam selesai.")

# =============================================================================
# PIPELINE UTAMA
# =============================================================================
def run_pipeline_for_file(dfsufile, grid, polys, gridxx, gridyy):
    if not os.path.exists(dfsufile):
        log.warning(f"Skip: tidak ditemukan {dfsufile}")
        return
    log.info("=" * 60)
    log.info(f"FILE: {os.path.basename(dfsufile)}")
    log.info("=" * 60)
    detected_roles, _ = detect_roles(dfsufile)
    if not detected_roles:
        log.warning("Tidak ada role. Skip.")
        return
    dsopen = mikeio.open(dfsufile)
    dfsutype, nlayers, zunique, unique_lids = detect_dfsu_type(dsopen)
    items_to_read = list(detected_roles.values())
    log.info(f"Membaca {len(items_to_read)} item ke RAM...")
    t_read = time.time()
    all_dates   = sorted(set(DAILYDATES + [HOURLYDATE]))
    time_slice  = slice(min(all_dates), max(all_dates))
    dsall       = dsopen.read(items=items_to_read, time=time_slice)
    log.info(f"Loaded {len(dsall.time)} timesteps {time.time()-t_read:.1f}s")
    ex  = dsopen.geometry.element_coordinates[:, 0]
    ey  = dsopen.geometry.element_coordinates[:, 1]
    bboxmask   = (ex >= BBOX[0]) & (ex <= BBOX[2]) & (ey >= BBOX[1]) & (ey <= BBOX[3])
    nelembbox  = bboxmask.sum()
    nelemtotal = len(bboxmask)
    log.info(f"Elemen BBox: {nelembbox}/{nelemtotal}")
    dsbbox     = dsall.isel(element=np.where(bboxmask)[0].tolist())
    log.info(f"dsbbox dims: element={nelembbox}")
    # Surface coords (layer teratas = semua elemen 2D dalam bbox untuk 2D file,
    # atau elemen surface untuk 3D)
    elemxsurf  = ex[bboxmask]
    elemysurf  = ey[bboxmask]
    dshourly   = dsbbox.sel(time=HOURLYDATE)
    run_daily_parallel(dsbbox, dsopen.geometry, elemxsurf, elemysurf,
                       gridxx, gridyy, grid, polys, detected_roles,
                       dfsutype, zunique, unique_lids, bboxmask)
    run_hourly_parallel(dshourly, dsopen.geometry, elemxsurf, elemysurf,
                        gridxx, gridyy, grid, polys, detected_roles,
                        dfsutype, zunique, unique_lids, bboxmask)
    del dsall, dsbbox, dshourly; gc.collect()

# =============================================================================
# ENTRY POINT
# =============================================================================
if __name__ == "__main__":
    t0 = time.time()
    log.info("Pipeline DFSU v8.4 dimulai")
    log.info(f"BBox: {BBOX}")
    log.info(f"Resolusi: DX={DX} DY={DY} ~{DX*111000:.0f}m")
    log.info(f"Workers: {MAXWORKERS}")
    log.info(f"Mode 3D: {PROCESS3D} | Depths: {TARGETDEPTHS}m")
    log.info(f"NoData: {NODATAVALUE}")
    log.info("-" * 50)
    log.info("[FIX-3] Nama field kecepatan: cspdms, cspdkt (6 char)")
    log.info("[FIX-4] Nama field vektor  : cur_U, cur_V, wave_U, wave_V (5 char)")
    log.info("[FIX-1] curdir: arctan2 dengan masking NODATA")
    log.info("[FIX-5] Logging: UTF-8 safe, karakter Unicode dikonversi")
    log.info("-" * 50)

    ensure_gdb(OUTPUTGDB) if OUTPUTMODEGDB else None
    writer = threading.Thread(target=gdb_writer_thread, args=(OUTPUTGDB,),
                               name="GDB-Writer", daemon=True)
    writer.start()
    log.info("GDB Writer thread aktif.")

    grid        = build_grid_manual()
    polys, xx, yy = build_polygons(grid)

    for dfsufile in DFSUFILES:
        run_pipeline_for_file(dfsufile, grid, polys, xx, yy)

    log.info("Menunggu GDB Writer menyelesaikan antrian...")
    GDBWRITEQUEUE.join()
    GDBWRITESTOP.set()
    writer.join(timeout=30)
    log.info("-" * 50)
    log.info(safe_log(f"[OK] SELESAI TOTAL {(time.time()-t0)/60:.1f} menit"))
    log.info(f"GDB    : {OUTPUTGDB}")
    log.info(f"GeoTiff: {OUTPUTTIF}")
    if OUTPUTPOINTSHP or OUTPUTPOLYGONSHP:
        log.info(f"SHP    : {os.path.join(BASEOUTPUT, 'FallbackSHP')}")
