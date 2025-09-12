# Splitpea

This is the Python package implementation of **Splitpea**: *SPLicing InTeractions PErsonAlized*.

Original repository: [https://github.com/ylaboratory/splitpea](https://github.com/ylaboratory/splitpea)\
Repository for this package's code: [https://github.com/ylaboratory/splitpea_package](https://github.com/ylaboratory/splitpea_package)\
This pip package is an extended version with easier installation and additional functions.

Splitpea quantifies rewiring in protein-protein interaction (PPI) networks driven by alternative splicing events. It integrates differential exon usage (PSI values) with domain-domain interactions (DDIs) and PPIs to generate condition or sample specific networks.

## Installation

Install from PyPI:

```bash
pip install splitpea
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


## Requirements

- Python >= 3.8
- Packages:
  - `numpy`, `networkx`, `intervaltree`, `ipykernel`, `matplotlib`, `plotnine`, `scikit-learn`, `adjustText`, `importlib_resources`, `requests`


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
     ensembl.id  symbol  chr  strand  exon.start  exon.end  psi.background  psi.sample  delta.psi  pval
     ```
   - Can be generated via `preprocess_pooled` (see below) from two splicing matrices: one for background samples and another for samples to compare individually against it. Supports SUPPA2 (.psi) and rMATS formats (SE.MATS.JC.txt or SE.MATS.JCEC.txt).

2. **Condition-specific mode** (`differential_format=suppa2` or `rmats`)  
   - **SUPPA2**: provide both `.psivec` and `.dpsi` (order agnostic)
   - **rMATS**: use `SE.MATS.JC.txt` or `SE.MATS.JCEC.txt`

> Currently, only **skipped exon (SE)** events are supported.

---

## Examples/Tutorials: 

- **Sample-specific mode:** https://colab.research.google.com/drive/1ktDhV5QSEm0em5R_Da9JhN1PGqMzrIsY?usp=sharing
- **Sample-specific mode (rMATS):** https://colab.research.google.com/drive/18bn9BRt8XPj-oUWceocKJkIVSeabWwOG?usp=sharing
- **Condition-specific mode (SUPPA2 / rMATS):** https://colab.research.google.com/drive/1tGua5zQYXhKEknvtRcIJ0cboFUrYKnWc?usp=sharing

---

## Quick start

Splitpea can be run **from the command line** or **within Python**.

Main subcommands:
- **`run`** — Build a rewired network
- **`plot`** — Visualize a saved network
- **`stats`** — Compute edge/gene statistics
- **`preprocess_pooled`** — Generate the differential exon table for sample-specific mode
- **`get_consensus_network`**- Generate summary consensus networks from a directory of rewired networks
- **`analyze_consensus_threshold`** — Threshold consensus neg/pos graphs across user-defined cutoffs and generate summary statistics and visualizations.

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

## `run` 

**Required:**
- **`in_file`** —  
  - `sample_specific`: one file path  
  - `suppa2`: two file paths (`.psivec` and `.dpsi`, typically outputed from diffSplice)  
  - `rmats`: one JC or JCEC file path

**Options:**
- **`out_file_prefix`** — Prefix for all output files (directory + base name); if no prefix is given, the output prefix defaults to `out_rewired_network` 
- `--differential_format {sample_specific,suppa2,rmats}` (default: `sample_specific`)
- `--skip` *(int)* — Number of lines to skip in input file (default: `1`)
- `--dpsi_cut` *(float)* — Delta PSI cutoff (default: `0.05`)
- `--sigscore_cut` *(float)* — Significance score cutoff (default: `0.05`)
- `--include_nas` *(bool)* — Include NAs in significance testing (default: `True`)
- `--index {0,1}` — Coordinate indexing (default: auto-set per format)
- `--map_path` — Custom mapping file between IDs and symbols (tab-delimited with headers: `symbol  entrez  ensembl  uniprot`)
- `--edge_stats_file` — Append gain/loss/chaos counts to a stats file

**Reference files**
Splitpea ships bundled reference datasets for human and mouse. These load automatically based on `--species` unless you override paths.

- `--species {human,mouse}` — Which bundled references to use (default: `human`)
- `ppif` — protein–protein interactions
- `ddif` — domain–domain interactions
- `entrezpfamf` - Entrez–Pfam mapping*
- `pfamcoordsf` - Pfam genome coordinates
- `tbf` - Tabix index ( usually the .tbi for pfamcoordsf)
- `map_path` - Gene ID mapping: symbol/entrez/ensembl/uniprot

> To override default bundled references, pass custom paths to the corresponding parameters. 

---



### Outputs

For each run (`out_file_prefix`), Splitpea produces:

- `<prefix>.edges.dat` — edge list with `node1  node2  weight  chaos`
- `<prefix>.edges.pickle` — NetworkX graph
-
> Python API: `splitpea.run(...)` returns the rewired `networkx.Graph`.
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
- **`--threshold`** — If set, only include edges with that pass threshold. Only change parameter when inputing a consensus network. Default: `0`.

**Outputs (depending on args):**
- PDF plot, Gephi TSV, and/or Cytoscape GML.

--- 

### `stats`
Summarize edge counts and write per-gene statistics.

Example (Python):
```python 
splitpea.stats(
    rewire_net="out/sample1.edges.pickle",
    out_file_prefix="stats/sample1",
    species="human",
    map_path=None,      # use bundled if None
    ppif=None, ddif=None, entrezpfamf=None  # override to use custom references
)
```

**Parameters:**
- **`rewire_net`** *(required)* — Path to `.edges.pickle` file for the same network.
- **`dat_file`** — Path to `.edges.dat` file generated by `splitpea.run`. If omitted, it will be inferred from `rewired_net`, but providing it directly will speed up execution.
- **`out_file_prefix`** *(required)* — Prefix for output files (e.g., `stats/sample1`).
- **`ppif`** — Path to PPI reference file. Defaults to bundled species reference.
- **`ddif`** — Path to DDI reference file. Defaults to bundled species reference.
- **`entrezpfamf`** — Path to Entrez–Pfam mapping file. Defaults to bundled species reference.
- **`map_path`** — Path to gene ID mapping file. Defaults to bundled mapping if not provided.
- **`species`** — Species for selecting default references (`human` or `mouse`).

**Outputs:**
- Console printout of:
  - Gain — number of edges gained in rewired network
  - Loss — number of edges lost in rewired network
  - Chaos — number of chaos edges
- `<out_file_prefix>_gene_stats.csv` containing:
  - Entrez ID
  - Gene symbol
  - Node degree in rewired network
  - Normalized degree relative to background PPI network
  - Counts of gain/loss/chaos edges incident to the gene

--- 

### `preprocess_pooled`

Helper function that builds sample-specific Splitpea inputs by comparing each target sample to a pooled normal background.  
It either downloads a normal splicing matrix from IRIS (GTEx) for a chosen tissue or a user can provide their own normal matrix.
You can provide either a .txt file from rMATS output files (SE.MATS.JC.txt or SE.MATS.JCEC.txt) for a single sample or a folder containing multiple samples and rMATS output files (SE.MATS.JC.txt or SE.MATS.JCEC.txt). The function will automatically process these inputs. For rMATS, it may be helpful to rename each file in the folder to match your sample names.

What the function does: 
1. Load a target splicing matrix (`compare_path`) with exon rows and sample columns (PSI values).
2. Obtain a normal/background splicing matrix:
   - Either download GTEx `<Tissue>` from IRIS (if `background` is given), or
   - Use your local file (`background_path`).
3. Compute mean PSI per exon for the background.
4. For each target sample, compute delta PSI vs. background and a p-value.
5. Write per-sample Splitpea-ready files (`{sample}-psi.txt`) into `out_psi_dir`.

Example (Python):

```python
# Option A: download and use IRIS/GTEx matrix
splitpea.preprocess_pooled(
    compare_path="compare.txt",
    background="Brain",
    background_download_root="/data/iris_cache",   # creates GTEx_<Tissue>/splicing_matrix/ here
    out_psi_dir="out_psi"
)

# Option B: use a local normal matrix
splitpea.preprocess_pooled(
    compare_path="compare.txt",
    background_path="/path/to/GTEx_Brain/splicing_matrix.txt",
    out_psi_dir="out_psi"
)
```

**Inputs:**
- **`compare_path`** *(required)* — Target (case) splicing matrix to compare against the pooled normal.  
  - Can be:  
    - a `.txt` file from rMATS (`SE.MATS.JC.txt` or `SE.MATS.JCEC.txt` format)
    - a directory of rMATS files (`SE.MATS.JC.txt` or `SE.MATS.JCEC.txt` format), with each file renamed to match its sample, or  
    - a `.txt` file in tab-delimited matrix format (rows = exons, columns = samples, values = PSI [0–1]).  
      - Expected header example:  
        ```
        GeneID  geneSymbol  chr  strand  exonStart  exonEnd  upstreamEE  downstreamES  sample1  sample2  ...
        ```
- **Normal/background matrix** — Provide **one** of:
  - **`background_`** — IRIS/GTEx tissue name to auto-download the background matrix (e.g., `Brain`, `AdiposeTissue`, …).  
    Optional: **`tissue_download_root`** to control where the `GTEx_<Tissue>/splicing_matrix/` folder is created.
  - **`background_path`** — Path to a pre-downloaded normal splicing data. Same accepted formats/data files as `compare_path`.

> You must supply **either** `background_` **or** `background__path`. 

**Parameters:**
- **`compare_path`** *(str, required)* — Path to the target splicing data to be compared.
- **`out_psi_dir`** *(str)* — Output directory for per-sample Splitpea files (`{sample}-psi.txt`). Will be created if missing. If none is given, defaults to create a out_psi folder in the current working directory.
- **`background`** *(str)* — IRIS/GTEx tissue name to auto-download the normal data.  
    **Must be one of:**
    `AdiposeTissue`, `AdrenalGland`, `Bladder`, `Blood`, `BloodVessel`, `Brain`, `Breast`, `CervixUteri`, `Colon`, `Esophagus`, `FallopianTube`, `Heart`, `Kidney`, `Liver`, `Lung`, `Muscle`, `Nerve`, `Ovary`, `Pancreas`, `Pituitary`, `Prostate`, `SalivaryGland`, `Skin`, `SmallIntestine`, `Spleen`, `Stomach`, `Testis`, `Thyroid`, `Uterus`, `Vagina`. 
- **`background_download_root`** *(str)* — Root directory for the IRIS download cache (creates `GTEx_<Tissue>/splicing_matrix/` under this path). If not given, defaults to current working directory.
- **`background_path`** *(str)* — Path to an existing normal/background splicing data file; bypasses download.
- **`map_path`** *(str)* — Gene ID mapping: symbol/entrez/ensembl/uniprot. If none is given, uses bundled mapping.
- **`species`** *(str)* — Species identifier (e.g. `human`, `mouse`); determines which default map path to use.   
- **`single_rMATS_compare`** *(bool, optional)* — If True, treat `compare_path` as a single rMATS file instead of a directory.  
- **`single_rMATS_background`** *(bool, optional)* — If True, treat `background_path` as a single rMATS file instead of a directory.  
- **`inclevel`** *(int, optional)* — Which inclusion-level field to use when parsing rMATS (1 or 2). Defaults to 1.  

**Outputs (written to `out_psi_dir`):**
- **`{sample}-psi.txt`** files (one per target sample), each a Splitpea (sample-specific) format table.

---

### `get_consensus_network`

Build consensus loss and gain summary networks from a directory of sample networks (`*.pickle` from `splitpea.run`). Chaos edges are excluded; edge weights and counts are accumulated across samples.

Example (Python):

```python
cons_neg, cons_pos = splitpea.get_consensus_network(
    net_dir="output"
)
```

Parameters: 
- **`net_dir`** *(str, required)* — Directory containing sample .pickle graphs.

**Outputs:**
- `(cons_neg, cons_pos)` — two NetworkX graphs:
  - **`cons_neg`**: edges with `weight <= 0`, attributes `weight` (sum) and `num_neg` (count).
  - **`cons_pos`**: edges with `weight > 0`, attributes `weight` (sum) and `num_pos` (count).
  - Both graphs have `graph['num_graphs']` = number of sample graphs combined.
- **Files written:**
  - `consensus_network_neg.pickle`
  - `consensus_network_pos.pickle` 

---

### `analyze_consensus_threshold`

Builds threshold consensus negative and/or positive rewired networks (from `splitpea.get_consensus_network`) across user-defined cutoffs and outputs basic statistics on them. Optionally save per-threshold graphs, and plot the # nodes and proportion of nodes retained across thresholds.

Example (Python):

```python
consensus_threshold_sizes = splitpea.analyze_consensus_threshold(
    neg_path="consensus_network_neg.pickle",
    pos_path="consensus_network_pos.pickle",
)
```

**Parameters:**
- **`neg_path`** *(str, optional)* — Path to consensus negative graph (`*_consensus_neg.pickle`).
- **`pos_path`** *(str, optional)* — Path to consensus positive graph (`*_consensus_pos.pickle`).
- **`thresholds`** *(list[float], optional)* — List of thresholds to evaluate (default: `[0.1, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0]`).
- **`label`** *(str, optional)* — Label/group name used in output and plots (default: `"sample"`).
- **`pickles_dir`** *(str, optional)* — Directory to save per-threshold graphs (default: `./threshold_networks/<label>/`).
- **`save_pickles`** *(bool, default `True`)* — Save per-threshold graphs as `.pickle`.
- **`write_txt`** *(bool, default `True`)* — Write TSV of sizes for each threshold.
- **`txt_path`** *(str, optional)* — Path for txt (default: `<pickles_dir>/consensus_threshold.lcc_sizes.txt`).
- **`save_pdf_prefix`** *(str, optional)* — If set, saves plots as `<prefix>_num_nodes.pdf` and `<prefix>_prop_nodes.pdf`.
- **`title_prefix`** *(str, optional)* — Optional text to prepend to plot titles.

**Outputs:**
- **`consensus_threshold_sizes`** *(pandas.DataFrame)* — Table of sizes per threshold with columns:  
  `label, direction, threshold, num_nodes, num_edges, prop_nodes`
- **Files written (optional):**
  - `threshold_networks/<label>/<direction>/<label>_consensus_<neg|pos>.thres_<t>.pickle` — thresholded graphs.
  - `threshold_networks/<label>/consensus_threshold.lcc_sizes.txt` — TSV summary of sizes.
  - `<prefix>_num_nodes.pdf` and `<prefix>_prop_nodes.pdf` — plots of #nodes and proportion-of-nodes vs threshold.


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



