# gtdb-itol-tools
This project generates interactive phylogenetic tree visualizations that color-code bacterial taxa based on **carbon pathway ratios** (e.g., CO₂:C1 metabolic ratios). It processes GTDB bacterial taxonomy data and generates iTOL-compatible datasets with:

### Quick Start

Run the main pipeline:

```bash
bac120_taxonomy.tsv` — GTDB taxonomy mapping \
c_by_e_donor_phylum.csv` — Carbon pathway data by phylum \
bac120.tree_stripped` — Newick format phylogenetic tree \
```