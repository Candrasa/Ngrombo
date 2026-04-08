# DFSU Pipeline MIKE DHI → Shapefile & GeoTiff
## Pipeline Konversi Data Oseanografi MIKE 21/3 FM ke Format Spasial

![Python](https://img.shields.io/badge/Python-3.9%2B-blue)
![MikeIO](https://img.shields.io/badge/MikeIO-2.0%2B-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)
![License](https://img.shields.io/badge/License-MIT-yellow)

---

## Deskripsi

Script Python ini digunakan untuk mengonversi file **DFSU** dari model **MIKE 21/3 Flexible Mesh (DHI)** menjadi data spasial siap pakai dalam bentuk:

- **Shapefile Point** untuk arus laut dan gelombang.
- **GeoTiff Raster** untuk visualisasi spasial.
- **Output per tanggal, per jam, dan per kedalaman**.
- **Output grid reguler 9 km** dengan resolusi `DX = DY = 0.083°`.

Script ini disusun untuk kebutuhan analisis oseanografi operasional, termasuk wilayah pesisir dan laut lepas, dengan dukungan komputasi paralel dan validasi dimensi agar lebih aman saat memproses file DFSU berukuran besar.

---

## Tujuan Script

Tujuan utama script ini adalah:

1. Membaca data DFSU 2D dan 3D dari model MIKE.
2. Mendeteksi item bawaan seperti arus, gelombang, dan angin.
3. Melakukan interpolasi data unstructured mesh ke grid reguler.
4. Menghasilkan shapefile point dan raster GeoTiff.
5. Membuat output harian, hourly, dan layer kedalaman tertentu.
6. Menyediakan atribut spasial yang siap dipakai di ArcGIS, QGIS, dan analisis lanjutan.

---

## Alur Proses

Script ini bekerja melalui beberapa tahap utama:

### 1. Inisialisasi library
Script memanggil library utama seperti:

- `mikeio`
- `numpy`
- `pandas`
- `geopandas`
- `rasterio`
- `scipy.interpolate`
- `shapely`
- `threading`
- `concurrent.futures`

Library tersebut dipakai untuk membaca DFSU, mengolah array, membangun shapefile, dan membuat raster.

### 2. Setup logging
Script menyiapkan sistem logging agar semua proses tercatat ke console dan file log.  
Logging dibuat aman untuk environment seperti Jupyter, Spyder, atau ArcGIS Pro dengan fallback jika `sys.stdout.buffer` tidak tersedia.

### 3. Konfigurasi input dan output
Bagian konfigurasi berisi:

- daftar file DFSU,
- folder output,
- bounding box,
- resolusi grid,
- tanggal harian,
- tanggal hourly,
- kedalaman target,
- jumlah worker paralel,
- flag output.

### 4. Deteksi tipe DFSU
Script memeriksa apakah file DFSU termasuk:

- **2D**, atau
- **3D**

Untuk 3D, script juga membaca layer depth dan layer ID.

### 5. Deteksi role/variabel
Script mencari item bawaan seperti:

- `curspeed`
- `curdir`
- `cur_U`
- `cur_V`
- `wavehs`
- `wavedir`
- `waveper`
- `windspeed`
- `winddir`

Jika ada nama item yang cocok, variabel itu dipakai dalam proses berikutnya.

### 6. Pembuatan grid reguler
Script membangun grid reguler berdasarkan bounding box dan resolusi `DX`, `DY`.  
Grid ini dipakai sebagai target interpolasi dari mesh fleksibel DFSU.

### 7. Ekstraksi layer kedalaman 3D
Jika file adalah DFSU 3D, script memilih layer terdekat untuk kedalaman:

- 11.4 m
- 20.4 m
- 29.5 m

Script juga memberi validasi toleransi kedalaman agar layer yang dipilih benar-benar mendekati target.

### 8. Interpolasi data ke grid
Data dari mesh unstructured diinterpolasi ke grid reguler dengan:

- `LinearNDInterpolator`
- `NearestNDInterpolator` sebagai fallback

Ini memastikan nilai tetap tersedia walaupun ada titik yang tidak terjangkau interpolasi linear.

### 9. Konversi arah dan kecepatan
Untuk arus:

- `curdir` dihitung dari `cur_U` dan `cur_V` bila diperlukan.
- `cspdms` disimpan dalam meter per detik.
- `cspdkt` disimpan dalam knot.

Untuk gelombang:

- `wavedir` disimpan sebagai arah datang gelombang.
- `wavehs` disimpan sebagai tinggi gelombang signifikan.
- `waveper` disimpan sebagai periode gelombang.

### 10. Penulisan output shapefile
Script menyimpan output point shapefile:

- `ptcur...shp` untuk arus
- `ptwave...shp` untuk gelombang

Setiap file memiliki atribut:

- arah,
- label kompas,
- konvensi arah,
- kecepatan/tinggi,
- komponen vektor,
- kelas BMKG,
- koordinat,
- waktu,
- kedalaman.

### 11. Penulisan output GeoTiff
Selain shapefile, script juga dapat menulis raster GeoTiff untuk:

- kecepatan arus,
- tinggi gelombang,
- variabel lain yang dipilih.

### 12. Proses harian dan hourly
Script memproses data dalam dua mode:

- **Harian**: rata-rata harian.
- **Per jam**: setiap timestep hourly.

### 13. Paralelisasi
Script memakai `ThreadPoolExecutor` agar proses lebih cepat pada komputer dengan banyak core/thread.  
Jumlah worker dapat diatur sesuai kapasitas CPU dan RAM.

---

## Struktur Data Input

### File arus 3D
File DFSU arus dipakai untuk menghasilkan:

- `curdir`
- `cspdms`
- `cspdkt`
- `cur_U`
- `cur_V`
- layer kedalaman

### File gelombang
File DFSU gelombang dipakai untuk menghasilkan:

- `wavedir`
- `wavehs`
- `waveper`
- `wave_U`
- `wave_V`

---

## Struktur Output

### Shapefile Point Arus
Contoh nama file:

- `ptcur20260420.shp`
- `ptcur20260420d10m.shp`
- `ptcur20260423000000.shp`

Kolom utama:

- `curdir`
- `dirlabel`
- `konvensi`
- `cspdms`
- `cspdkt`
- `cur_U`
- `cur_V`
- `klsbmkg`
- `depthm`
- `lon`
- `lat`
- `datetime`

### Shapefile Point Gelombang
Contoh nama file:

- `ptwave20260420.shp`
- `ptwave20260423000000.shp`

Kolom utama:

- `wavedir`
- `dirlabel`
- `konvensi`
- `wavehs`
- `waveper`
- `wave_U`
- `wave_V`
- `klsbmkg`
- `depthm`
- `lon`
- `lat`
- `datetime`

### GeoTiff Raster
Contoh nama file:

- `rscurspd20260420.tif`
- `rswavehs20260420.tif`

---

## Konfigurasi Utama

Bagian ini biasanya diubah sesuai kebutuhan:

```python
DFSUFILES = [
    r"path/to/3D_Current.dfsu",
    r"path/to/Wave.dfsu",
]

BASEOUTPUT = r"path/to/output"
BBOX = (109.4, -6.5, 111.2, -5.5)
DX, DY = 0.083, 0.083
DAILYDATES = ["2026-04-20", "2026-04-21"]
HOURLYDATE = "2026-04-23"
TARGETDEPTHS = [11.4, 20.4, 29.5]
MAXWORKERS = 16
OUTPUTGEOTIFF = True
OUTPUTPOINTSHP = True
```

---

## Cara Menjalankan

```bash
python dfsu_pipeline_v84.py
```

Jika memakai virtual environment:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python dfsu_pipeline_v84.py
```

---

## Validasi dan Keamanan Data

Script ini dilengkapi validasi untuk:

- dimensi data input,
- keberadaan item DFSU,
- kesesuaian layer kedalaman,
- potensi `IndexError`,
- potensi `ValueError` akibat `or` pada numpy array,
- pengisian `NODATA`,
- konsistensi penamaan field.

---

## Catatan Penting

Beberapa penyesuaian yang diterapkan:

- `curU` dan `curV` dibuat konsisten menjadi `cur_U` dan `cur_V`.
- Nama field kecepatan arus dipendekkan menjadi `cspdms` dan `cspdkt`.
- Arah arus dihitung dengan masking `NODATA` agar tidak menghasilkan nilai palsu.
- Logging dibuat aman untuk environment notebook dan IDE.

---

## Riwayat Versi

| Versi | Tanggal | Perubahan |
|---|---|---|
| v8.1 | 2026-04-06 | Pipeline dasar DFSU 2D/3D |
| v8.2 | 2026-04-06 | Perbaikan IndexError layer 3D |
| v8.3 | 2026-04-07 | Perbaikan ekstraksi layer dan bboxmask |
| v8.4 | 2026-04-08 | Perbaikan arah arus, rename field, logging, konsistensi role |

---

## Lisensi

Proyek ini menggunakan lisensi MIT.
