# Splitpea

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

## Requirements

- Python >= 3.8
- Packages:
  - `numpy`, `networkx`, `intervaltree`, `ipykernel`, `matplotlib`, `plotnine`, `scikit-learn`, `adjustText`, `importlib_resources`

---

## Features

- Adjustable PSI and significance thresholds
- Bundled reference files, with options for customization
- Species-specific references: human (default) or mouse
- Main Outputs:
  - `.edges.dat`: Edge list with weights and if an edge is chaos 
  - `.edges.pickle`: NetworkX graph object
- Additional features and outputs:
  - Edge rewiring summaries
  - Gene level statistics
  - Auto-generated network plots
  - Gephi-compatible TSV files
  - Cytoscape-compatible gml files

---

## Two modes

1. **Sample-specific mode** (`differential_format=sample_specific`)  
   - Input: a **differential exon table** for a single sample.  
   - Tab-delimited header:
     ```
     ensembl.id  symbol  chr  strand  exon.start  exon.end  psi.gtex  psi.tcga  delta.psi  pval
     ```
   - Can be generated via `preprocess_pooled` (see below) from two splicing matrices: one for background samples and another for samples to compare individually against it.

2. **Condition-specific mode** (`differential_format=suppa2` or `rmats`)  
   - **SUPPA2**: provide both `.psivec` and `.dpsi` (order agnostic)  
   - **rMATS**: use `SE.MATS.JC.txt` or `SE.MATS.JCEC.txt`

> Currently, only **skipped exon (SE)** events are supported.

---

## Quick start

Splitpea can be run **from the command line** or **within Python**.

Main subcommands:
- **`run`** — Build a rewired network
- **`plot`** — Visualize a saved network
- **`stats`** — Compute edge/gene statistics
- **`preprocess_pooled`** — Generate the differential exon table for sample-specific mode


### Command line

**Examples:**
```bash
# Sample-specific mode
splitpea run sample1-psi.txt out/sample1

# Condition-specific mode (SUPPA2)
splitpea run diffSplice.psivec diffSplice.dpsi out/condA \
  --differential_format suppa2

# Condition-specific mode (rMATS)
splitpea run SE.MATS.JCEC.txt out/condB \
  --differential_format rmats
```
> **Note:** The second argument (`out/sample1`, `out/condA`, etc.) is the **output file prefix**.

### Python API

**Examples:**
```python
import splitpea

# Sample-specific mode
G = splitpea.run("sample1-psi.txt", "out/sample1")

# Condition-specific mode (SUPPA2)
G = splitpea.run(["diffSplice.psivec", "diffSplice.dpsi"], "out/condA",
                 differential_format="suppa2")

# Condition-specific mode (rMATS)
G = splitpea.run("SE.MATS.JCEC.txt", "out/condB",
                 differential_format="rmats")
```
---

## `run` — Parameters and defaults

**Required:**
- **`in_file`** —  
  - `sample_specific`: one file path  
  - `suppa2`: two file paths (`.psivec` and `.dpsi`)  
  - `rmats`: one JC or JCEC file path
- **`out_file_prefix`** — Prefix for all output files (directory + base name)

**Options:**
- `--differential_format {sample_specific,suppa2,rmats}` (default: `sample_specific`)
- `--skip` *(int)* — Number of lines to skip in input file (default: `1`)
- `--dpsi_cut` *(float)* — Delta PSI cutoff (default: `0.05`)
- `--sigscore_cut` *(float)* — Significance score cutoff (default: `0.05`)
- `--include_nas` *(bool)* — Include NAs in significance testing (default: `True`)
- `--index {0,1}` — Coordinate indexing (default: auto-set per format)
- `--species {human,mouse}` — Which bundled references to use (default: `human`)
- `--map_path` — Custom mapping file between IDs and symbols (tab-delimited with headers: `symbol  entrez  ensembl  uniprot`)
- `--edge_stats_file` — Append gain/loss/chaos counts to a stats file

**Reference files**
Splitpea ships bundled reference datasets for human and mouse. These load automatically based on `--species` unless you override paths.

- `ppif` — protein–protein interactions
- `ddif` — domain–domain interactions
- `entrezpfamf` - Entrez–Pfam mapping*
- `pfamcoordsf` - Pfam genome coordinates
- `tbf` - Tabix index ( usually the .tbi for pfamcoordsf)
- `map_path` - Gene ID mapping: symbol/entrez/ensembl/uniprot

Overriding references:
- Pass custom paths to the corresponding parameters. 

---



### Outputs

For each run (`out_file_prefix`), Splitpea produces:

- `<prefix>.edges.dat` — edge list with `node1  node2  weight  chaos`
- `<prefix>.edges.pickle` — NetworkX graph

---

## Other Subcommands

### `plot`
Load a saved network (`.edges.pickle`) and export it in various formats or as a plot.

**Example (Python):**
```python
splitpea.plot(
    pickle_path="out/sample1.edges.pickle",
    pdf_path="plot/sample1_plot.pdf",         # optional
    gephi_path="plot/sample1_gephi.csv",      # optional
    cytoscape_path="plot/sample1_cyto.gml",   # optional
    with_labels=True,
    symbol=True,
    map_path=None,        # use bundled mapping if None
    species="human",
    self_edges=False,
    lcc=True
)
```

**Parameters:**
- **`pickle_path`** *(required)* — Path to `.edges.pickle` file generated by `splitpea run`.
- **`--with_labels`** — Draw node labels in plots. Default: `False`.
- **`--pdf_path`** — Path to save a PDF of the plotted network via Matplotlib. Omit to skip PDF.
- **`--gephi_path`** — Path to save a Gephi-compatible TSV file. Omit to skip.
- **`--cytoscape_path`** — Path to save a Cytoscape-compatible `.gml` file. Omit to skip.
- **`--symbol`** — If `True` (default), replace Entrez IDs with gene symbols.
- **`--map_path`** — Path to custom gene ID mapping file. Defaults to bundled mapping if not provided.
- **`--species`** — Species for mapping defaults (`human` or `mouse`).
- **`--self_edges`** — If `True`, keep self-loop edges in output. Default: `False`.
- **`--lcc`** — If `True` (default), plot only the largest connected component.
- **`--max_nodes`** — Maximum nodes allowed in matplotlib plotting; larger graphs skip matplotlib. Default: `2000`.
- **`--max_edges`** — Maximum edges allowed in matplotlib plotting; larger graphs skip matplotlib. Default: `10000`.

**Outputs (depending on args):**
- PDF plot, Gephi TSV, and/or Cytoscape GML.

--- 

### `stats`
Summarize edge counts and write per-gene degree stats.

Example (Python):
```python 
splitpea.stats(
    dat_file="out/sample1.edges.dat",
    rewire_net="out/sample1.edges.pickle",
    out_file_prefix="stats/sample1",
    species="human",
    map_path=None,      # use bundled if None
    ppif=None, ddif=None, entrezpfamf=None  # override to use custom references
)
```

**Parameters:**
- **`dat_file`** *(required)* — Path to `.edges.dat` file generated by `splitpea.run`.
- **`rewire_net`** *(required)* — Path to `.edges.pickle` file for the same network.
- **`out_file_prefix`** *(required)* — Prefix for output files (e.g., `stats/sample1`).
- **`ppif`** — Path to PPI reference file. Defaults to bundled species reference.
- **`ddif`** — Path to DDI reference file. Defaults to bundled species reference.
- **`entrezpfamf`** — Path to Entrez–Pfam mapping file. Defaults to bundled species reference.
- **`map_path`** — Path to gene ID mapping file. Defaults to bundled mapping if not provided.
- **`species`** — Species for selecting default references (`human` or `mouse`).

**Outputs:**
- Console printout of:
  - **Gain** — number of edges gained in rewired network
  - **Loss** — number of edges lost in rewired network
  - **Chaos** — number of chaos edges
- `<out_file_prefix>_gene_stats.csv` containing:
  - Entrez ID
  - Gene symbol
  - Node degree in rewired network
  - Normalized degree relative to background PPI network
  - Counts of gain/loss/chaos edges incident to the gene

### `preprocess_pooled`

End-to-end helper that builds sample-specific Splitpea inputs by comparing each target sample to a*pooled normal background.  
It either downloads a normal splicing matrix from IRIS (GTEx) for a chosen tissue or uses a normal matrix you provide.

Pipeline: 
1. Load a target splicing matrix (`compare_file`) with exon rows and sample columns (PSI values).
2. Obtain a normal/background splicing matrix:
   - Either download GTEx `<Tissue>` from IRIS (if `tissue` is given), or
   - Use your local file (`normal_path`).
3. Compute mean PSI per exon for the background (pooled normal).
4. For each target sample, compute delta PSI vs. background and a p-value.
5. Write per-sample Splitpea-ready files (`{sample}-psi.txt`) into `out_psi_dir`.

Example (Python):

```python
# Option A: download and use IRIS/GTEx matrix
splitpea.preprocess_pooled(
    compare_file="compare.txt",
    out_psi_dir="out_psi",
    tissue="Brain",
    tissue_download_root="/data/iris_cache"   # creates GTEx_<Tissue>/splicing_matrix/ here
)

# Option B: use a local normal matrix
splitpea.preprocess_pooled(
    compare_file="compare.txt",
    out_psi_dir="out_psi",
    normal_path="/path/to/GTEx_Brain/splicing_matrix.tsv"
)
```

**Inputs:**
- **`compare_file`** *(required)* — Target (case) splicing matrix to compare against the pooled normal.  
  - Tab-delimited; rows = exons, columns = samples; values = PSI \[0–1\].  
  - Expected header example:
    ```
    AC  GeneName  chr  strand  exonStart  exonEnd  upstreamEE  downstreamES  sample1  sample2  ...
    ```
- **Normal/background matrix** — Provide **one** of:
  - **`tissue`** — IRIS/GTEx tissue name to auto-download the background matrix (e.g., `Brain`, `AdiposeTissue`, …).  
    Optional: **`tissue_download_root`** to control where the `GTEx_<Tissue>/splicing_matrix/` folder is created.
  - **`normal_path`** — Path to a pre-downloaded normal splicing matrix (same exon schema; PSI values).

> You must supply **either** `tissue` **or** `normal_path`. If both are given, `normal_path` is used.

**Parameters:**
- **`compare_file`** *(str, required)* — Path to the target splicing matrix to be compared.
- **`out_psi_dir`** *(str, required)* — Output directory for per-sample Splitpea files (`{sample}-psi.txt`). Will be created if missing.
- `tissue` *(str)* — IRIS/GTEx tissue name to auto-download the normal matrix.  
    **Must be one of:**
    `AdiposeTissue`, `AdrenalGland`, `Bladder`, `Blood`, `BloodVessel`, `Brain`, `Breast`, `CervixUteri`, `Colon`, `Esophagus`, `FallopianTube`, `Heart`, `Kidney`, `Liver`, `Lung`, `Muscle`, `Nerve`, `Ovary`, `Pancreas`, `Pituitary`, `Prostate`, `SalivaryGland`, `Skin`, `SmallIntestine`, `Spleen`, `Stomach`, `Testis`, `Thyroid`, `Uterus`, `Vagina`. 
- **`tissue_download_root`** *(str, optional)* — Root directory for the IRIS download cache (creates `GTEx_<Tissue>/splicing_matrix/` under this path).
- **`normal_path`** *(str, optional)* — Path to an existing normal/background splicing matrix file; bypasses download.

**Outputs (written to `out_psi_dir`):**
- **`{sample}-psi.txt`** files (one per target sample), each a **Splitpea (sample-specific) format** table:

---

## Tutorials

- **Condition-specific mode (SUPPA2 / rMATS)** — _notebook/colab link (add later)_
- **Sample-specific mode (pooled normal)** — _notebook/colab link (add later)_

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



