# Meantools — Development Version

This repository contains the development version of Meantools. It is under active development, with new functionality and improvements being added. Features and interfaces may change as they are tested and refined.

Current development includes:

1. **Faster correlation analysis**, with support for running on HPC environments using SLURM.
2. **New HPC scripts** that let users choose RAM-based or disk-based mutual-rank operations to suit available memory.
3. **Improvements to structural database searches.**
4. **Improvements to reaction prioritization scoring.**
5. **Improvements to reaction visualization.**

More functionality is coming. For HPC setup and usage, see the [HPC README](HPC_test/README.md).

## Dependencies

For the full development script collection, install the Python dependencies in a separate Python 3.10 or 3.11 environment:

```bash
python -m pip install -r requirements.txt
python -m pip check
```

For the standalone HPC correlation script and its tests, the smaller dependency set is sufficient:

```bash
python -m pip install -r requirements-hpc.txt
```

The files cover imports in the Python scripts, local modules, and notebook. They are development manifests, not fully validated environment lockfiles. The full set limits NumPy and pandas to versions below 2 because some older scripts still use removed APIs; the standalone HPC set has no such upper limits. Other legacy APIs may need updates. For example, visualization scripts import RDKit's `DrawingOptions`; that import was verified locally with RDKit 2024.3.5, but the full chemistry and visualization workflows have not been validated against every current dependency release.

For GPU reaction scoring, select a PyTorch build appropriate for the cluster before installing the full requirements. PyTorch Geometric extensions may also depend on the PyTorch/CUDA combination; follow the [official installation instructions](https://pytorch-geometric.readthedocs.io/en/stable/install/installation.html). The local model loader handles PyTorch Geometric 2.x checkpoints, so the requirements constrain that package to 2.x. Model files under `models/` are also needed for reaction-likelihood scoring.

External requirements are installed separately:

- **Java**, with `java` available on `PATH`, for `cluster_one-1.0.jar` and `ReactionDecoder.jar`.
- **SLURM** and the site's environment modules for HPC submission; local execution does not require SLURM.
- **SQLite support in Python**; disk-based HPC mutual ranking requires SQLite 3.25 or newer. `sqlite3` is part of Python, not a pip dependency.
- **Reference databases and input data**, supplied separately. Some annotation and reaction-scoring workflows also access remote services.

The package names `scikit-learn`, `umap-learn`, and `svg.py` provide the imports `sklearn`, `umap`, and `svg`, respectively. Project modules such as `gizmos`, `graphics`, `mutual_ranks`, and `gnn_som` are included in this repository and should not be installed from pip.
