# Urban RF-EMF Voronoi Mapping in Tirana

This repository contains the Python workflow used to generate Voronoi-based maps of outdoor RF-EMF exposure in two urban areas of Tirana, Albania.

The script creates publication-ready maps for total RF-EMF exposure and selected mobile communication frequency bands.

## Study areas


The figures are generated for:

- Sheshi Karl Topia
- Sheshi Mustafa Qemal Ataturk

## Main output

The script generates maps for:

- `E_total`
- `LTE800`
- `GSM900`
- `LTE1800`
- `UMTS2100`
- `LTE2600`
- `NR3500`

Each figure is exported as:

- PNG
- TIFF
- PDF

## Repository structure

```text
urban-rf-emf-voronoi-tirana/
├── scripts/
│   └── generate_voronoi_maps.py
├── data/
│   └── urban_rf_emf_dataset.csv
├── results/
│   └── figures/
│       └──Sheshi_Karl_Topia_E_total_Voronoi_Final.png
│       └──Sheshi_Karl_Topia_LTE800_Voronoi_Final.png    
│       └──Sheshi_Mustafa_Qemal_Ataturk_NR3500_Voronoi_Final.png   
│         
├── requirements.txt
├── LICENSE
└── README.md
```

## Input data

Place the anonymized dataset here:

```text
data/urban_rf_emf_dataset.csv
```

The required column names are described in `data_dictionary.md`.

## Installation

Create and activate a virtual environment, then install dependencies:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## Run

From the repository root:

```bash
python scripts/generate_voronoi_maps.py
```

Generated figures will be saved in:

```text
result/figures/paper_voronoi_final_renamed/
```

## Notes

- The script use Voronoi polygons to preserve the spatial discreteness of measured RF-EMF values.
- The method avoids continuous interpolation and represents each polygon according to the nearest measurement point.
- Basemap titles are loaded through `contextily`, so an internet connection is required when running the script.

## Citation

If this code is used in a publication, citate dhe associated article or repository DOI when available. 




























