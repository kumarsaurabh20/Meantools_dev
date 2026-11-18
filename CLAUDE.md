# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**MEANtools** (Metabolite ANticipation tools) is a Python bioinformatics research pipeline for predicting metabolic pathways by integrating multi-omics data (metabolomics + transcriptomics) with biochemical reaction databases (RetroRules, MetaNetX). The primary workflow is documented in `meantools.ipynb`.

## Running the Pipeline

There is no build system. Scripts are run directly with Python. Key invocation patterns:

```bash
# Query mass databases for metabolite annotation
python queryMassNPDB_mod.py -add data/ESI-MS-adducts.csv -ms metabolome.csv ...

# Correlation analysis (multi-omics integration)
python corrMultiomics_mod.py -ft feature_table.csv -qm expression.csv ...

# Cluster merging
python merge_clusters.py -ft metabolome.csv -qm transcriptome.csv ...

# Mass transition integration
python pathMassTransitions_mod.py -c clusters.csv -t transitions.csv ...

# Pathway prediction (core engine)
python heraldPathways_mod.py -c clusters.csv -r rules.csv ...

# Visualization and pathway filtering
python paveWays.py -sp structures.csv -r reactions.csv ...
```

The full end-to-end workflow with example data and parameter explanations is in `meantools.ipynb`.

## Architecture

### Four-Stage Pipeline

1. **Database Preparation** — `queryMassNPDB*.py`, `fill_db.py`, `pathMassTransitions*.py`
   - Parses RetroRules/MetaNetX databases, generates mass transition tables, validates SMILES

2. **Omics Data Preparation** — `corrMultiomics*.py`, `mutual_ranks.py`
   - Loads metabolomics/transcriptomics CSVs, queries mass DBs, computes Pearson/Spearman/log-Pearson correlations and mutual ranks, stores results in SQLite

3. **Cluster Analysis** — `merge_clusters.py`
   - Runs ClusterOne (`cluster_one-1.0.jar`) on correlation graphs, then merges overlapping functional clusters using fingerprinting, overlap, or co-expression strategies

4. **Pathway Prediction & Visualization** — `heraldPathways*.py`, `paveWays.py`
   - Applies SMARTS reaction rules to expand metabolite structures iteratively
   - Filters candidates by gene expression/correlation
   - Scores reaction sites using the GNN-SOM ensemble (`gnn_som/`, `models/model0.pt`–`model9.pt`)
   - Generates SVG/PNG molecular structure images and pathway networks via ReactionDecoder (`ReactionDecoder.jar`)

### Key Modules

| File | Role |
|------|------|
| `gizmos.py` | Core utilities: SQL queries, data conversion, statistics, heatmap generation |
| `heraldPathways.py` | Main prediction engine: reaction testing, virtual molecule generation |
| `paveWays.py` | Visualization: SVG generation, pathway networks, enzyme annotation |
| `corrMultiomics.py` | Correlation matrix computation and mutual rank calculations |
| `merge_clusters.py` | Functional cluster merging (overlap/fingerprint/co-expression) |
| `pathMassTransitions.py` | Links metabolites to reaction candidates via mass transitions |
| `gnn_som/__init__.py` | GNN-SOM model loading (`createGnnSom`, `loadGnnSomState`) |
| `gnn_som/MolFromKcf.py` | KEGG KCF format parser for molecule graph construction |
| `reaction_likelihoods.py` | GNN + ReactionDecoder integration for reaction scoring |

### Data Flow

```
mz values → queryMassNPDB → metabolite SMILES (SQLite)
Transcriptomics + Metabolomics → corrMultiomics → correlation matrices
Correlation matrices → merge_clusters → functional clusters (CSV)
Clusters + mass transitions → pathMassTransitions → reaction candidates
Reaction candidates + gene annotations → heraldPathways → predicted pathways
Predicted pathways → paveWays → SVG/PNG images + CSV/Excel results
```

### External Dependencies

- **JAR tools**: `ReactionDecoder.jar` (reaction visualization/scoring), `cluster_one-1.0.jar` (clustering)
- **KEGG HTTP API**: `gnn_som/MolFromKcf.py` calls the KEGG `/mol2kcf/` endpoint to fetch KCF-format molecules
- **PyTorch models**: 10 pre-trained GNN models in `models/` for site-of-metabolism prediction; configuration in `models/config_ori.json` and `models/config_mod.json`

### Script Variants

Many scripts have `*_mod.py` and `*_mod2.py` variants representing experimental iterations. The `_mod` variants are generally more current than the base versions.

## Dependencies

No `requirements.txt` or `setup.py` exists. Core libraries (inferred from imports):

- **Cheminformatics**: `rdkit`
- **ML/GNN**: `torch`, `torch_geometric`
- **Data**: `pandas`, `numpy`, `scipy`, `scikit-learn`
- **Graph**: `networkx`, `markov_clustering`
- **Visualization**: `matplotlib`, `seaborn`, `svgutils`, `svgelements`
- **DB**: `sqlalchemy`, `sqlite3` (stdlib)
- **Utilities**: `tqdm`, `requests`

## Input/Output Formats

**Inputs**: Metabolomics CSVs (`metabolite_id,mz`), transcriptomics expression matrices, gene annotation CSVs (`gene,pfam1;pfam2`), SMARTS reaction rules  
**Outputs**: SQLite databases, CSV/Excel tables, SVG/PNG structure images, pathway network files
