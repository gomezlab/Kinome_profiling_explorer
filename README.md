# Explorer for Kinome Profiling Data

https://gomezlab.github.io/Kinome_profiling_explorer/

An interactive dashboard for exploring large-scale kinome profiling datasets. This tool provides a unified view of compound-kinase interaction potencies and inhibition Kds across ~1750 compounds and ~500 kinases, totalling ~900,000 combinations.

## Data Overview

This explorer integrates high-quality kinome-wide profiling data, specifically focusing on:
- **Percent Inhibition**: Measured at a 1µM drug concentration.
- **Apparent Kd (nM)**: Quantitative binding affinity derived from dose-response profiling.
- **Extended Metadata**: Includes compound clinical trial phases, PubChem identifiers, and kinase-specific links (Entrez/UniProt).

## How to Use the Viewer

1.  **Toggle Mode**: Choose between **Compound Explorer** (to see all targets of a specific drug) or **Kinase Explorer** (to see all drugs targeting a specific kinase).
2.  **Selection**: Search for a compound or kinase using the searchable dropdown menus.
3.  **View Modes**: 
    - **Inhibition View**: Displays sensitivity as % inhibition at 1µM. Higher bars indicate stronger inhibition.
    - **Kd View**: Displays potency as apparent dissociation constants (nM). Lower values (and longer bars) indicate stronger binding affinity. Only potent binders (< 500 nM) are shown in the summary plot.
4.  **Info Panels**: View metadata for your selection, including clinical status, selectivity scores, and external database links.
5.  **Data Table**: Filter the raw data using the table search box and export your current filtered view as a CSV for further analysis.

## Data Mining and Processing (Methods)

The data provided in this viewer was synthesized from three major kinome profiling initiatives. All chemical entities were standardized using PubChem CIDs, and all protein targets were mapped to Entrez Gene IDs to ensure cross-dataset compatibility.

### Data Sources
- **Kuster Dataset**: Integrated from mass-spectrometry based kinobeads profiling. Global proteome-wide data was retrieved from the MASSIVE repository ([MSV000092248](https://massive.ucsd.edu/v05/MSV000092248/)) and processed as described in [Reinecke et al. (2024)](https://pubmed.ncbi.nlm.nih.gov/37904048/).
- **LINCS Dataset**: Derived from the Library of Integrated Network-based Cellular Signatures, including imputed kinome-wide profiles as described in [Joisa et al (2023)](https://pubmed.ncbi.nlm.nih.gov/38025707/).
- **Klaeger Dataset**: Mined from the 2107 foundational target landscape of clinical kinase inhibitors as detailed in [Berginski and Joisa et al. (2023)](https://pubmed.ncbi.nlm.nih.gov/36809237/).

### Processing Pipeline
1.  **Kuster Processing**: Raw dose-response data files were parsed to extract normalized Label-Free Quantification (LFQ) intensities. Relative intensities (percent of control) were calculated by normalizing drug-treated intensities to the corresponding DMSO control values within each batch. Median values were calculated across replicates for the 1µM dose.
2.  **Harmonization**: Relative binding intensities and apparent Kd values for LINCS and Klaeger were extracted from their respective published tidy datasets. All relative intensities were converted to a consistent 0-100 scale (% inhibition).
3.  **Apparent Kd Calculation**: For dose-response data where Kd was not explicitly provided, binding curves were integrated or retrieved from source metadata to provide the `mean_kd_apparent` values shown in this dashboard.
4.  **Integration**: Datasets were merged into a single master interaction table. Only kinases with sufficient coverage across the full integrated data (n > 100 measurements) were retained to ensure robustness.

---
*Created by Chinmaya Joisa for the Gomez Lab.*
