# Splitpea

...
This documentation is old I will update!
...

This is the Python package implementation of **Splitpea**: *SPLicing InTeractions PErsonAlized*.

Original repository: [https://github.com/ylaboratory/splitpea](https://github.com/ylaboratory/splitpea)\
Repository for this package's code: [https://github.com/ylaboratory/splitpea_package](https://github.com/ylaboratory/splitpea_package)\
This pip package is an extended version with easier installation and additional functions.

Splitpea quantifies rewiring in protein-protein interaction (PPI) networks driven by alternative splicing events. It integrates differential exon usage (PSI values) with domain-domain interactions (DDIs) and PPIs to generate condition-specific networks.

## Installation

Install from PyPI:

```bash
# Note: Splitpea is not on the main PyPI yet.
# Eventually you will be able to install it with:
# pip install splitpea

# For now, install from the test repository:
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple splitpea
```

Some functions require `tabix` to be installed; to do so, simply run:

```bash
sudo apt-get install tabix
```

Alternatively, if you are using conda, you can install tabix via:

```bash
conda install -c bioconda tabix
```

Functions that depend on `tabix` will raise an error if `tabix` is not found.

---

## Usage

Splitpea can be run **from the command line** or **within Python**.

### Command-line

```bash
splitpea <in_file> <out_file_prefix> [options]
```

Examples:

- **Splitpea format**:

  ```bash
  splitpea psi.txt sample1
  ```

- **rMATS input**:

  ```bash
  splitpea psi_rmats.txt sample1 --input_format rmats
  ```

- **SUPPA2 input**:

  ```bash
  splitpea diffSplice.psivec diffSplice.dpsi sample1 --input_format suppa2
  ```

**Notes:**

- Currently, only **skipped exon (SE)** events are supported across all input types.
- rMATS support works for **JCEC and JC outputs**.
- For SUPPA2 input, provide both `.psivec` and `.dpsi` files (order does not matter as long as one ends with `.psivec` and the other ends with `.dpsi`).
- The default splitpea format (i.e. regular) requires a `.txt` file (tab-seperated) with the following header:
  ```
  ensembl.id  symbol  chr  strand  exon.start  exon.end  psi.gtex  psi.tcga  delta.psi  pval
  ```

> **Output:**\
> Splitpea automatically generates a `.dat` file (edge list) and a `.pickle` file (NetworkX graph) based on the provided output prefix. Additional optional outputs files also use the same output prefix.

---

### Python

```python
import splitpea

# Default run
net = splitpea.run("psi.txt", "sample1")

# rMATS input
net = splitpea.run("psi_rmats.txt", "sample1", input_format="rmats")

# SUPPA2 input
net = splitpea.run(["diffSplice.psivec", "diffSplice.dpsi"], "sample1", input_format="suppa2")
```

---

## Features

- Supports multiple input formats: Splitpea, rMATS, SUPPA2
- Species-specific reference handling:
  - Human (default) or mouse
  - Option to specify custom reference files, but default reference data is pre-built into the package 
- Adjustable PSI and significance thresholds
- Outputs:
  - `.edges.dat`: Edge list with weights and if an edge is chaos 
  - `.edges.pickle`: NetworkX graph object
- Optional outputs:
  - Gene degree statistics
  - Edge rewiring summaries
  - Auto-generated network plots
  - Gephi-compatible TSV files
  - Cytoscape-compatible gml files

---

## Requirements

- Python >= 3.8
- Packages:
  - `numpy`, `networkx`, `intervaltree`, `ipykernel`, `matplotlib`, `plotnine`, `scikit-learn`, `adjustText`, `importlib_resources`

---

## Full Parameters

### Python 

```python
splitpea.run(
    in_file, out_file_prefix: str, skip: int = 1, dpsi_cut: float = 0.05, sigscore_cut: float = 0.05,
    include_nas: bool = True, verbose: bool = False, input_format: str = "regular",
    ppif: str = None, ddif: str = None, entrezpfamf: str = None, pfamcoordsf: str = None, tbf: str = None,
    species: str = "human", index: int = None, edge_stats_file: str = None,
    gene_degree_stats: bool = False, plot_net: bool = False, gephi_tsv: bool = False, 
    cytoscape_gml: bool = False, map_path: str = None
)
```

### Command Line 

To view all command-line options:
```bash
splitpea -h
```
which lists: `in_file`, `out_file_prefix`, `--input_format`, `--skip`, `--dpsi_cut`, `--sigscore_cut`, `--include_nas`, `--verbose`, `--ppif`, `--ddif`, `--entrezpfamf`, `--pfamcoordsf`, `--tbf`, `--species`, `--index`, `--edge_stats_file`, `--gene_degree_stats`, `--plot_net`, `--gephi_tsv`, `--cytoscape_gml`, `--map_path` and explanations for each paramter.  


---

### Important Notes

- **Species Selection:**

  - Default: `species="human"`
  - To use mouse references: `species="mouse"`

- **Reference Files:**

  - Defaults are bundled for human and mouse.
  - Custom references can be supplied manually: `ppif` (PPI file), `ddif` (DDI file), and mappings ( including `entrezpfamf`, `pfamcoordsf`, and `tbf`) (note: `tbf` is often the same file as `pfamcoordsf`). `map_path` is a path to a .txt file containing gene ID mappings between `symbol`, `entrez`, `ensembl`, and `uniprot` identifiers (where those are the header/first line of the tab-seperated text file). There are default mappings built into the package, but you can define your own/update the mapping if you so wish.

- **Input Handling:**

  - `skip=1` by default to skip the header.
  - `index`: coordinate system adjustment:
    - 0-based for Splitpea format and rMATS input
    - 1-based for SUPPA2 input

- **Additional Outputs:**

  - `edge_stats_file`:  Appends rewired edge statistics per sample.  If file doesn't exist, it creates a new one with header:

    ```
    sample	num_gain_edges	num_loss_edges	num_chaos_edges
    ```

    (`sample` is simply the output prefix.)

  - `gene_degree_stats=True`:  Saves node degree stats including:

    - Node degree
    - Normalized degree = (degree in rewired network) / (degree in background PPI)
      - background network represent PPI network unaffected by splicing changes

  - `plot_net=True`:  Quick network plot with matplotlib (may crash for very large networks; it is recommended to use dedicated network visualization software, such as Gephi, for better handling and visualization of  the rewired networks).

  - `gephi_tsv=True`:  Saves a Gephi-compatible TSV file.

  - `cytoscape_gml=True`: Saves a Cytoscape-compatible gml file

---

## Extra functionality: 

(running the pipeline as shown in the [original repo example](https://github.com/ylaboratory/splitpea/tree/master))

In addition to Splitpea itself, this package includes several extra utility functions, derived from the original GitHub repository, that can be run (only via Python code). These functions enable preprocessing steps such as combining exon usage files, calculating differential splicing (delta PSI) in "Splitpea format", building consensus networks across patients, and constructing background PPI networks.

Note: These preprocessing functions are mainly designed for use with rMATS files, particularly PSI-based alternative splicing (skipped exon) matrices.

Example of expected input format:

```
AC	GeneName	chr	strand	exonStart	exonEnd	upstreamEE	downstreamES	sample1	sample2	sample3	...
```

Example pipeline: 

```python
# 1. Combine spliced exons across a folder
splitpea.combine_spliced_exon("/path/to/psi_files")

# This function averages exon PSI values across multiple samples in a folder.
# For each spliced exon input file, it computes the mean PSI value across samples,
# and outputs a summarized file with the mean PSI.
# Inputs:
# - `in_dir`: directory containing differential splicing PSI files (txt format).
# Output:
# - For each input file, generates a `{sample}_combined_mean.txt` summarizing PSI per exon.

# 2. Calculate delta PSI between normal and tumor
splitpea.calculate_delta_psi("normal_combined.txt", "normal.txt", "tumor.txt", "/path/to/psi_files")

# This function calculates the delta PSI (difference in exon inclusion levels) and associated p-values
# between a background condition (e.g., normal tissues) and a target condition (e.g., tumors).
# It compares average PSI values, performs statistical testing, and outputs files formatted
# for direct input into Splitpea for network construction.
# Inputs:
# - `sum_bg_file`: summarized background PSI values.
# - `bg_file`: raw sample-wise background PSI data.
# - `target_file`: raw sample-wise target PSI data.
# - `outdir`: output directory for the resulting delta PSI + p-value files.
# Output:
# - Files named {sample}-psi.txt with exon-level delta PSI and p-values, ready for network building.

# 3. Run Splitpea (see below)

# 4. Create consensus network across multiple patients
splitpea.get_consensus_network("/path/to/splitpea_output")

# This function loads individual patient-specific Splitpea networks from a directory,
# separates edges into positive (gained) and negative (lost) rewiring edges,
# and constructs two consensus graphs: one for consensus positive edges and one for consensus negative edges.
# It counts how many times each edge appears across the networks and sums edge weights accordingly.
# Outputs:
# - `consensus_network_neg.pickle`: consensus network of lost (negative) interactions.
# - `consensus_network_pos.pickle`: consensus network of gained (positive) interactions.

# 5. Construct the "background" PPI network
splitpea.get_background(ppif, ddif, entrezpfamf)

# This function constructs the "background" PPI network needed for downstream Splitpea analysis.
# It requires the following reference files:
# - `ppif`: protein-protein interaction file
# - `ddif`: domain-domain interaction file
# - `entrezpfamf`: gene-to-Pfam domain mapping file
# The resulting background network represents the baseline set of protein interactions assuming no alternative splicing-induced changes.
```

### Running Splitpea on Multiple Samples Example (Python)

```python
import os
import splitpea

data_dir = "/content/psis"
output_dir = "/content/output"
os.makedirs(output_dir, exist_ok=True)

for file in os.listdir(data_dir):
    if file.endswith(".txt"):
        input_path = os.path.join(data_dir, file)
        output_prefix = os.path.join(output_dir, file.replace("-psi.txt", ""))
        splitpea.run(input_path, output_prefix)
```

### Command Line Version

```bash
#!/bin/bash

data_dir="/content/psis"
output_dir="/content/output"
mkdir -p "$output_dir"

for filepath in "$data_dir"/*.txt; do
  filename=$(basename "$filepath" -psi.txt)
  splitpea "$filepath" "$output_dir/$filename"
done
```

---

## Citation

If you use Splitpea, please cite:

```bibtex
@inproceedings{dannenfelser2023splitpea,
  title={Splitpea: quantifying protein interaction network rewiring changes due to alternative splicing in cancer},
  author={Dannenfelser, Ruth and Yao, Vicky},
  booktitle={Pacific Symposium on Biocomputing 2024},
  pages={579--593},
  year={2023},
  organization={World Scientific}
}
```

---



