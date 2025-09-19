# Splitpea

This is the Python package implementation of **Splitpea**: *SPLicing InTeractions PErsonAlized*, a method for calculating and analyzing protein-protein interaction changes due to alternative splicing. It integrates differential exon usage (PSI values) with domain-domain interactions (DDIs) and protein-protein interactions (PPIs) to generate rewired PPI networks compairing two conditions.

## Requirements and Installation

- Python >= 3.8
- Packages:
  - `numpy`, `networkx`, `intervaltree`, `ipykernel`, `matplotlib`, `plotnine`, `scikit-learn`, `adjustText`, `importlib_resources`, `requests`

This package also requires `tabix` which can be installed system wide:

```bash
sudo apt-get install tabix
```

or in a conda environment:

```bash
conda install -c bioconda tabix
```

The Splitpea package itself is installable directly from PyPI:

```bash
pip install splitpea
```

## Usage

Splitpea takes as input skipped exon data (SE) at the sample or differential expression level from SUPPA2 or rMATS and maps potential changes onto PPI networks resulting in a NetworkX graph object. The package provides two forms of operation for building rewired networks: _condition-specific mode_, which takes differential exon usage results directly from SUPPA2 or rMATS, and _sample-specific mode_, which takes rMATS individual sample files and compares them with a multi-sample reference background using the algorithm outlined in our [2023 Splitpea paper](https://pubmed.ncbi.nlm.nih.gov/38160308/). We outline these different modes in detail in the example Jupyter Notebooks found in the `examples` directory.

When running from SUPPA2, Splitpea requires both the `.psivec` and the `.dpsi` files as input. Starting from rMATS, Splitpea accepts either the `SE.MATS.JC.txt` or `SE.MATS.JCEC.txt` file. In sample-specific mode, users can choose from a single rMATS input file of SE events (example ending in `single`) or  a matrix of SE events where each sample in the matrix is individually compared against the background to create a sample-specific rewired network (example ending in `multi`).

Splitpea ships bundled reference datasets (protein-protein interaction, domain-domain interaction, Entrez-Pfam mappings, and gene symbol conversion) for human and mouse. These load automatically based on `--species` unless you override paths when using the `run` command. The package also includes normal tissue splicing backgrounds for tissues in GTEx as assembled in the [IRIS data set](https://www.pnas.org/doi/10.1073/pnas.2221116120). The full list of supported tissues are as follows, and will be downloaded if used in sample-specific mode with the `preprocess_pooled` command:

> `AdiposeTissue`, `AdrenalGland`, `Bladder`, `Blood`, `BloodVessel`, `Brain`, `Breast`, `CervixUteri`, `Colon`, `Esophagus`, `FallopianTube`, `Heart`, `Kidney`, `Liver`, `Lung`, `Muscle`, `Nerve`, `Ovary`, `Pancreas`, `Pituitary`, `Prostate`, `SalivaryGland`, `Skin`, `SmallIntestine`, `Spleen`, `Stomach`, `Testis`, `Thyroid`, `Uterus`, `Vagina`. 


### Examples / tutorials

We include three sample workflows to help users get started.

1. [Condition-specific mode](./examples/condition_mode.ipynb): takes the output of rMATS and SUPPA2 to produce a single rewired network. 

2. [Sample-specific mode - single sample](./examples/sample_mode_single.ipynb): which runs Splitpea on rMATS data for a single sample versus a background collection of healthy normal tissue samples to produce a single rewired network.

3. [Sample-specific mode - multiple samples](./examples/sample_mode_multi.ipynb): which runs Splitpea on group of rMATS samples versus a background collection of healthy normal tissue samples to produce one network for each rMATS sample. 

### Quick start

Splitpea can be run from the command line or within Python.

Main subcommands:

- **`run`** — Build a rewired network
- **`plot`** — Visualize a saved network
- **`stats`** — Compute edge/gene statistics
- **`preprocess_pooled`** — Generate the differential exon table for sample-specific mode
- **`get_consensus_network`** - Generate summary consensus networks across multiple rewired networks
- **`analyze_consensus_threshold`** — Analyze network properties of a consensus network

### Command line examples

```bash
# Condition-specific mode (SUPPA2)
splitpea run diffSplice.psivec diffSplice.dpsi suppa_example_output --differential_format suppa2

# Condition-specific mode (rMATS)
splitpea run SE.MATS.JCEC.txt rmats_example_output --differential_format rmats
```
> **Note:** The second argument is the output file prefix.

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




