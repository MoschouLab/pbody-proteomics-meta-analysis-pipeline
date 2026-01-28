# P-body Proteomics Meta-Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](https://doi.org/)

**A reproducible bioinformatics pipeline for comprehensive proteomics meta-analysis of P-body dynamics under heat stress conditions.**

---

## 📋 Overview

This pipeline integrates proteomics data from different experimental approaches (Affinity Purification and Proximity-Dependent Labeling) to characterize P-body protein composition across non stress and heat stress conditions in *Arabidopsis thaliana*.

### Key Features

- ✅ **Automated data integration** from previously published proteomics datasets
- ✅ **Comprehensive GO enrichment analysis** with statistical rigor (FDR correction)
- ✅ **Network-aware analysis** using STRING PPI database
- ✅ **Statistical overlap testing** with bootstrap validation
- ✅ **Professional visualizations** (Venn diagrams, UpSet plots, composite figures)
- ✅ **CLI parametrizable** for flexibility and reproducibility
- ✅ **FAIR principles** compliant (Findable, Accessible, Interoperable, Reusable)

---

## 🎯 Scope and Limitations

### What This Pipeline Does ✅

This pipeline provides a **fully reproducible analysis** of the specific proteomics dataset from:

> Liu C, et al. (2023). An actin remodeling role for Arabidopsis processing bodies revealed by their proximity interactome. *EMBO J*. 42(9):e111885.

It implements:
- Data loading and quality control
- UniProt annotation enrichment
- GO enrichment analysis with statistical rigor
- Network-aware overlap analysis
- Publication-ready visualizations

### What This Pipeline Does NOT Do ❌

This is **not a general-purpose proteomics analysis tool**. It is:
- Designed for the specific Excel format from Liu et al. (2023)
- Hardcoded for column names: `protein`, `AP_NS`, `AP_HS`, `PDL_NS`, `PDL_HS`
- Optimized for P-body proteomics analysis specifically

### Can You Use It for Your Data? 🤔

**If your data has the same structure** (same column names, Excel format):
- ✅ Yes, place your file in `data/raw/` and run the pipeline

**If your data has different structure**:
- ⚠️ You will need to modify `scripts/01_load_pbody_proteomics.py` to match your column names
- ⚠️ Requires basic Python knowledge

**For general proteomics analysis**:
- 🔧 Consider using established tools like MaxQuant, Perseus, or Proteome Discoverer
- 🔧 This pipeline is best used as a **reference implementation** to adapt for similar analyses

---

## 🧬 Scientific Background

**P-bodies** (Processing bodies) are biomolecular condensates involved in mRNA regulation and stress response. This pipeline addresses:

**Research Question**: How does P-body protein composition change under non stress vs heat stress, and what functional pathways are involved?

**Approach**: Meta-analysis of two complementary proteomics methods:
- **AP (Affinity Purification)**: Captures core P-body proteins
- **PDL (Proximity-Dependent Labeling)**: Captures proximal proteome

**Key Finding**: Only 1 gene (**AT1G08370**, mRNA-decapping enzyme) is consistently enriched across all 4 experimental conditions, establishing it as the universal P-body marker.

---

## 📊 Pipeline Architecture
```
Input Data (Excel files)
    ↓
[01] Load Proteomics Data
    ↓ (4774 proteins)
[02] Fetch UniProt Annotations
    ↓ (99.79% annotated, +21.7% improvement)
[03] Cell Wall Enrichment Analysis
    ↓ (Fisher's exact test, 225 CW genes)
[04] Comprehensive GO Enrichment
    ↓ (1277 GO terms tested, FDR corrected)
[05] Venn/UpSet Overlap Analysis
    ↓ (Bootstrap validation, PPI network analysis)
Output: Tables, Figures, Statistical Reports
```

---

## 🚀 Quick Start

### Prerequisites

- Python 3.10
- Conda/Miniconda
- 8 GB RAM minimum
- ~2 GB disk space

### Installation
```bash
# Clone repository
git clone https://github.com/MoschouLab/pbody-proteomics-pipeline.git
cd pbody-proteomics-pipeline

# Create conda environment
conda env create -f environment.yml
conda activate pbody_proteomics

# Verify installation
python scripts/01_load_pbody_proteomics.py --help
```

### Basic Usage
```bash
# Run complete pipeline (default parameters)
./run_pipeline.sh

# Or run scripts individually
python scripts/01_load_pbody_proteomics.py
python scripts/02_fetch_protein_annotations.py
python scripts/03_cell_wall_enrichment_analysis.py
python scripts/04_go_enrichment_comprehensive.py --threshold 0.8
python scripts/05_venn_upset_analysis.py --threshold 0.8 --bootstraps 1000
```

### Expected Runtime

| Script | Runtime | Output |
|--------|---------|--------|
| Script 01 | ~2s | 4774 proteins loaded |
| Script 02 | 4-5 min (first run), <1s (cached) | 99.79% annotated |
| Script 03 | ~1s (cached) | 225 cell wall genes |
| Script 04 | ~15s | 104 significant GO terms |
| Script 05 | ~9s | 12 intersections, 8 Venn diagrams |
| **Total** | **~6 minutes** | **50+ output files** |

---

## 📂 Repository Structure
```
pbody-proteomics-pipeline/
├── data/
│   ├── raw/                          # Input datasets (2 Excel files)
│   └── processed/                    # Intermediate and final outputs
│       ├── pbody_proteomics_annotated.csv
│       ├── go_enrichment_*.csv
│       ├── ppi_edges.csv (155K interactions)
│       └── intersections/ (12 gene lists)
├── results/
│   └── figures/
│       ├── venn_diagrams/ (9 SVG files)
│       └── composite_overlap_analysis.png
├── scripts/
│   ├── 01_load_pbody_proteomics.py
│   ├── 02_fetch_protein_annotations.py
│   ├── 03_cell_wall_enrichment_analysis.py
│   ├── 04_go_enrichment_comprehensive.py
│   ├── 05_venn_upset_analysis.py
│   └── process_string_ppi.py
├── logs/                             # Execution logs
├── environment.yml                   # Conda environment
├── requirements.txt                  # Python dependencies
├── run_pipeline.sh                   # Master execution script
├── LICENSE                           # MIT License
├── MIGRATION_LOG.md                  # Development history
└── README.md                         # This file
```

---

## 🔬 Input Data

### Required Files

Place the following files in `data/raw/`:

1. **`EMBOJ-2022-111885_SourceFile1.xlsx`** (871 KB)
   - **Source**: [DOI: 10.15252/embj.2022111885]
   - **Citation**: Liu C, Mentzelopoulou A, Muhammad A, Volkov A, Weijers D, Gutierrez-Beltran E, Moschou PN. An actin remodeling role for Arabidopsis processing bodies revealed by their proximity interactome. *EMBO J*. 2023 May 2;42(9):e111885. PMID: 36741000; PMCID: PMC10152145.
   - **Description**: AP and PDL proteomics datasets (No Stress -NS- and Heat Stress -HS- conditions). 

### Optional: PPI Network Data

For network-aware analysis (Script 05):
- Automatically downloads from STRING database (v12.0)
- Or provide custom `ppi_edges.csv` in `data/processed/`

---

## 📈 Output Files

### Core Outputs

| File | Description | Used For |
|------|-------------|----------|
| `pbody_proteomics_annotated.csv` | Annotated proteomics table (4774 proteins) | All downstream analyses |
| `cell_wall_enrichment_summary.csv` | Fisher's exact test results | Cell wall analysis |
| `go_enrichment_*_*.csv` | GO enrichment results (12 files) | Functional analysis |
| `overlap_statistics_FDR.csv` | Statistical overlap testing | Overlap validation |
| `venn_4way_all_conditions.svg` | 4-way Venn diagram | Visualization |
| `upset_all_conditions.svg` | UpSet plot | Intersection analysis |

### Supplementary Outputs

- 12 intersection gene lists (`intersection_*.csv`)
- 8 additional Venn diagrams (2-way, 3-way)
- Bootstrap confidence intervals
- Directionality analysis (stress response)
- Core P-body gene GO analysis

---

## ⚙️ Configuration

### CLI Parameters

**Script 03 & 04: GO Enrichment**
```bash
python scripts/04_go_enrichment_comprehensive.py \
    --threshold 0.8 \          # log2FC cutoff (default: 0.8)
    --min-genes 5 \            # Minimum genes per GO term
    --top-n 15 \               # Top N terms to plot
    --no-gprofiler             # Skip g:Profiler validation
```

**Script 05: Overlap Analysis**
```bash
python scripts/05_venn_upset_analysis.py \
    --threshold 0.8 \          # log2FC cutoff
    --bootstraps 1000 \        # Bootstrap iterations
    --no-network \             # Skip PPI analysis
    --no-bootstrap             # Skip bootstrap (faster)
```

### Key Thresholds

- **Enrichment threshold**: log₂FC > 0.8 (default)
- **Statistical significance**: FDR < 0.05 (Benjamini-Hochberg)
- **PPI confidence**: ≥400 (medium confidence, STRING scale 0-1000)
- **Bootstrap iterations**: 1000 (95% CI)

---

## 🧪 Reproducibility

### Computational Environment

- **OS**: Linux (Ubuntu 24.04, tested on ThinkStation P5)
- **Python**: 3.10.x
- **Key dependencies**:
  - pandas 2.0.3
  - numpy 1.24.3
  - scipy 1.11.1
  - matplotlib 3.10.8
  - statsmodels 0.14.0
  - networkx 3.4.2

### Random Seed

All stochastic analyses use fixed seed (`RANDOM_SEED = 42`) for reproducibility.

### Caching

- UniProt annotations cached (2.7 MB) - speeds up re-runs
- GO annotations cached (2.8 MB) - instant re-execution
- PPI network cached (155K edges) - one-time download

---

## 📊 Key Results

### Proteome-wide Enrichment

| Condition | Enriched Genes | Cell Wall Genes | Top GO Term | FDR |
|-----------|----------------|-----------------|-------------|-----|
| AP_NS | 299 | 35 (2.44×) | Protein phosphorylation | 4.6×10⁻⁷ |
| AP_HS | 157 | 12 (1.59×) | Ribosome biogenesis | 0.073 |
| PDL_NS | 582 | 53 (1.88×) | Apoplast (3.80×) | 2.1×10⁻⁶ |
| PDL_HS | 196 | 4 (0.44×) | mRNA binding (2.29×) | <0.05 |

### Overlap Analysis

- **Total unique genes**: 1,095
- **Core gene (all 4 conditions)**: 1 (AT1G08370: mRNA-decapping enzyme)
- **Significant overlaps**: 2 (AP_NS∩AP_HS: 3.76×, PDL_NS∩PDL_HS: 1.93×)
- **Stress response**: Net depletion (799 lost vs 270 gained)

### Network Analysis

- **PPI network**: 4,063 proteins, 155,125 interactions
- **Coherence**: 2.4-3.0 (highly interconnected)
- **Coverage**: 85.6% of proteomics dataset

---

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Reporting Issues

Please use [GitHub Issues](https://github.com/MoschouLab/pbody-proteomics-pipeline/issues) to report bugs or request features.

---

## 📄 Citation

If you use this pipeline in your research, please cite:
```bibtex
@software{pbody_proteomics_pipeline,
  author = {Moya-Cuevas, J and Moschou, PN.},
  title = {P-body Proteomics Meta-Analysis Pipeline},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/MoschouLab/pbody-proteomics-pipeline},
  note = {Version 1.0.0}
}
```

**Related publication** (under review):
> Moschou, P.N., et al. (2026). "Cell Wall Dynamics and Biomolecular Condensates in Plants." *Plant Communications*.

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Third-Party Data

- **UniProt**: [CC BY 4.0](https://www.uniprot.org/help/license)
- **STRING**: [CC BY 4.0](https://string-db.org/cgi/access?footer_active_subpage=licensing)
- **GO Consortium**: [CC BY 4.0](http://geneontology.org/docs/go-citation-policy/)

---

## 👥 Authors

**Lead Bioinformatics Scientist**:
- **Jose Moya-Cuevas** - Postdoctoral Researcher
  - IMBB-FORTH, University of Crete
  - ORCID: [0000-0001-9537-8556]
  - Email: jose_moya@imbb.forth.gr

**Principal Investigator**:
- **Panagiotis N. Moschou** - Associate Professor
  - IMBB-FORTH, University of Crete
  - ERC Consolidator Grant PLANTEX
  - ORCID: [0000-0001-7212-0595]
  - Email: panagiotis.moschou@imbb.forth.gr
  - Lab: [Moschou Lab](https://moschoulab.com)

---

## 🙏 Acknowledgments

- **ERC Consolidator Grant PLANTEX** for funding
- **UniProt**, **STRING**, and **GO Consortium** for databases
- Original proteomics datasets from [DOI: 10.15252/embj.2022111885]

---

## 📚 Related Resources

- [Moschou Lab Website](https://pmoschoulab.org/)

---

**Version**: 1.0.0  
**Last Updated**: 2026-01-23  
**Status**: Production-ready
