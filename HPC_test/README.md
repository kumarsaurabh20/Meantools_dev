# Running correlations on Athena

The optimized `corrMultiomics_HPC.py` is standalone: copy it and the submission script to your HPC installation. It needs Python >=3.8, NumPy, SciPy, pandas, threadpoolctl, and SQLite (>=3.25 for disk MR). ClusterOne additionally needs Java and the existing JAR. It does not import `gizmos` or `mutual_ranks`, so you do not need to copy those modules for this HPC runner.

Install the one additional dependency in the activated environment if necessary:

```bash
python -m pip install threadpoolctl
python -c 'import sqlite3; print(sqlite3.sqlite_version)'
java -version
```

## What the supplied log establishes

- The batch script started at 20:22:35; correlation progress began around 20:23:55. The log has no submission or scheduling timestamps, so it cannot establish queue wait time.
- The first reported correlation result took 92 seconds. A chunk of eight feature rows was processed before results became visible.
- Correlations completed in 3:40:31 (7,007 feature rows). The excerpt continues until 09:10 the next morning, over nine hours after correlation progress reached 100%, without another stage appearing.
- The old monitor measured only the Python parent, excluding its worker processes and Java. Its low CPU readings do not prove that the allocation was idle; its ~2 GB RSS cannot establish the job's peak memory requirement. Thread CPU percentages were also divided by all CPUs on the node, obscuring usage within the allocation.
- The old code created two process pools, forked while a monitoring thread was active, recycled workers, and ran pandas plus SciPy separately for every pair. The observed pause is consistent with a worker shutdown/lifecycle problem, but the log alone cannot prove the exact cause.

## Implementation

- Float64 row normalization plus block matrix products replace per-pair pandas/SciPy calls. P-values use the Pearson beta distribution or Spearman t approximation used by SciPy. Constant rows are excluded because their correlations are undefined. `pearsonlog` removes zero rows, rejects negative values, and transforms each input once.
- A bounded thread executor shares the normalized arrays, uses one BLAS thread per worker, and holds at most one block per worker. There are no process pools, managers, worker copies, or fork/monitor interactions. `-t` defaults to `SLURM_CPUS_PER_TASK` and cannot exceed it.
- One writer streams retained correlation output into SQLite. Both MR backends visit retained edges instead of constructing all node combinations or repeatedly scanning the correlation DataFrame for every node. Disk mode stores intermediate ranks on scratch; RAM mode retains compact NumPy edge arrays and sorts them in memory (see below).
- MR retains the previous signed descending correlation ranking (even though the correlation cutoff uses absolute magnitude). It includes both endpoint roles in ranks, retains the last rank for duplicated neighbors, and emits one undirected edge per distinct node pair. For transcript-vs-transcript data, identical labels in the two matrices represent the same node, as before. Self correlations participate in ranks but are not emitted as MR edges, as before. Ties now have deterministic neighbor-ID ordering; old pandas tie ordering was unspecified, so exact tied MR values may change. Spearman can now also run MR when explicitly requested.
- Numeric node IDs are used only for ClusterOne input, avoiding errors for names containing spaces. The `<prefix>_nodes` table stores the mapping. Cluster output keeps `Members` labels and adds `Member_IDs`; use the latter with the lookup when labels contain whitespace. Edge and correlation tables keep original labels. Numeric ID ordering follows feature rows then previously unseen transcript rows.
- Intermediates have job-unique names on scratch. SQL cache is capped at 256 MiB; sorting and edge/rank storage can still require several to tens of GB of scratch for dense networks. The submission sets `SQLITE_TMPDIR` there too. Check node-local free space for larger datasets. No GPU or multi-node distribution is implemented.
- Java is checked before computing, its heap is capped, failed subprocesses stop with the actual diagnostic, and empty outputs produce useful errors. Input CSV errors, common samples, finite numeric values, and IDs are validated before computation.
- Each stage logs wall time; correlation blocks log cumulative progress. `/usr/bin/time -v` and SLURM accounting provide measurements without a noisy background monitor.

The correlation schema remains `metabolite, gene, correlation, P`; MR/weight and cluster table suffixes retain the original names for `-mdr` runs. For backward compatibility, omitting `-mdr` still uses the original default decay list `[5,10,25]`. The `metabolite` column is the feature-table ID even for transcript-vs-transcript analysis.

Filtering now applies to shared samples, rather than filtering on extra samples that would subsequently be discarded. Nonfinite inputs are rejected instead of propagating undefined values. The unused `-p` and `-a` options now report that they are unsupported. Input readers accept comma/semicolon text CSVs; spreadsheets must be exported first.

## Choose RAM or disk for mutual ranks

The same ranking rules, tie handling, weights, and output schemas apply in both modes:

| Option | Intermediate storage | Tradeoff |
| --- | --- | --- |
| `--mr-storage disk` (default) | SQLite edge/rank tables on scratch | Lower RAM requirement; sorting/joins depend on disk speed |
| `--mr-storage ram` | NumPy edge/rank arrays in memory | Avoids intermediate database I/O; uses more RAM |

RAM mode stores retained edges during correlations, sorts directed connections in NumPy, assigns ranks within each node, retains the largest rank for repeated neighbors, and aligns forward/reverse ranks. It does **not** simply put SQLite on a RAM disk. Ranking is primarily single-threaded in both modes; `-t` controls correlation workers. Final correlation/MR/weight tables still write to the requested SQLite database in both modes. ClusterOne still uses small network/output files on scratch.

For a job allocated 32 GB, an initial RAM working budget of 20 GiB leaves room for other Python data and buffers:

```bash
# Add to a direct Python invocation:
--mr-storage ram --mr-memory-gb 20

# Or choose with the provided submission script:
sbatch --export=ALL,MR_STORAGE=ram,MR_MEMORY_GB=20 HPC_test/corrMultiomics.sh

# Disk mode:
sbatch --export=ALL,MR_STORAGE=disk HPC_test/corrMultiomics.sh
```

For RAM mode, the script estimates ranking working memory as **256 bytes per retained correlation edge + 128 MiB**. For example, 10 million retained edges estimate about 2.51 GiB and 49 million about 11.81 GiB. This is a conservative planning estimate, not a measured peak or hard RSS cap; other process memory, library allocations, and allocator behavior are additional. The retained edge count includes self edges and both orientations, since those affect ranks. An optional `--mr-memory-gb` budget is checked while collecting edges and again before sorting, so there is no need to allocate a worst-case full network upfront. The normalized correlation matrices are released before ranking.

If no explicit budget is given in RAM mode, `SLURM_MEM_PER_NODE`, when available as a positive numeric MiB value, supplies a budget of 60% of the allocation. This is based on **allocated** RAM, not free RAM across the node. Outside that case, no automatic budget is assumed; set one explicitly. Choose a budget smaller than the job memory allocation. A budget is not a SLURM memory request: changing it does not increase `#SBATCH --mem`.

If the estimate exceeds the budget, the run stops with an actionable error rather than silently changing modes. A Python `MemoryError` also suggests disk mode, but a scheduler/cgroup OOM kill cannot be caught by Python. Partial output tables can remain after a failed run; rerun using a fresh `-tn`. ClusterOne's Java heap must also fit within the allocation, although it starts after NumPy ranking completes and the rank arrays have been released.

The choice applies from the beginning of a **new run**; it cannot switch the already-running HPC job. Disk mode remains the default so existing submissions keep their prior memory behavior.

## Resource starting point

Use **one node, one task, 8 physical cores, 32 GB RAM, 4 hours, `hmq`** for the first full optimized run of these inputs. This is a conservative starting request, not a measured optimum or a runtime guarantee. The submission uses `hmq` as requested. This partition choice does not itself require increasing CPU or memory requests. Account access, actual queue wait, and partition wall-time limits depend on local policy.

The standard nodes have 40 physical cores / 80 hardware threads. `--hint=nomultithread` asks SLURM for one thread per physical core. The original request reserved 10 CPUs while Python used 8 workers, with another setting allowing 3 OpenMP threads. The new request and worker count agree and avoid nested BLAS threading.

More cores will help only the matrix block phase. Input reading, SQLite writing/ranking, and much of ClusterOne will not scale with Python `-t`. A small, well-used allocation can therefore finish sooner overall than a large request that queues longer.

Start with `--block-size 128`. Approximate live RAM includes the two normalized input matrices, raw input/preparation copies while loading, SQL buffers, and up to `workers * block_size * transcript_rows` correlation entries plus index arrays, selected coefficients, p-values, and DataFrame buffers. These latter arrays can cost several times the raw 8-byte correlation block. The script logs matrix and raw-block sizes; these are not a complete peak estimate. If block memory becomes large, lower `--block-size` before increasing the memory allocation. The Java heap defaults to 8 GB and must leave space for the Python process, Java native memory, and accounting overhead.

## Submit and measure

Review the paths in `corrMultiomics.sh`, then:

```bash
sbatch corrMultiomics.sh
```

The original submission uses two transcript matrices. This is preserved intentionally. If this should be metabolite-vs-transcript, change `-ft` to the metabolite CSV before submitting.

The table prefix includes the job ID so reruns do not duplicate existing results. Direct Python invocations default to failing if any destination table exists. Explicit `--if-exists replace` replaces destination tables; `append` adds duplicates and requires an identical node mapping. Output commits incrementally; a failed job can leave partial tables. Use a new prefix for a fresh run.

After completion, inspect all accounting rows (including `.batch` and the `srun` step):

```bash
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem,ExitCode
scontrol show job JOBID
```

`TotalCPU / (Elapsed * AllocCPUS)` estimates allocation-wide CPU efficiency when taken from the appropriate aggregate/step record; do not add parent and step records that account for the same work. MaxRSS interpretation depends on site accounting: check step/task records, and use cgroup peak memory if your site exposes it. `/usr/bin/time` reports process resource usage; it is not a simultaneous sum of Python and Java RSS.

To compare scaling on the same input and parameters, run **4, 8, and 16 cores** with 32 GB initially (each job writes its own prefix):

```bash
sbatch --cpus-per-task=4 corrMultiomics.sh
sbatch --cpus-per-task=8 corrMultiomics.sh
sbatch --cpus-per-task=16 corrMultiomics.sh
```

Choose the smallest CPU count near the fastest total runtime, accounting for queue time. If measured full-run peak memory is comfortably below 12 GB, try 16 GB on a subsequent run; otherwise request measured peak plus roughly 25–50% headroom. Raise Java heap only if ClusterOne reports heap exhaustion and the SLURM allocation has sufficient room. Adjust wall time from measured full-run time with headroom. No definitive best CPU/memory pair can be derived from the parent-only log.

## Validation

```bash
python -m unittest discover -s HPC_test/tests -v
bash -n HPC_test/corrMultiomics.sh
```

Seven tests cover Pearson/log-Pearson/Spearman agreement with SciPy, input formats, MR with overlapping node roles, end-to-end streamed SQLite output and rerun rejection, ClusterOne errors/name restoration, exact RAM/disk output equivalence on ties/duplicates/self edges, RAM budgets, and the RAM pipeline. Passed on local Python 3.13 and on the Python 3.8 meantools interpreter (with threadpoolctl supplied temporarily for that check; install it on HPC).

A local synthetic kernel benchmark (64 x 256 pairs, 250 samples) measured 0.902 s for a SciPy-per-pair loop versus 0.012 s for the block implementation with 4 workers. Coefficients agreed within 1e-12. This excludes I/O, ranks, and ClusterOne and is not an Athena full-run speedup prediction. Full SLURM execution and real Java clustering have not been tested here.

References: [SLURM sbatch resource semantics](https://slurm.schedmd.com/sbatch.html), [SLURM CPU allocation and binding](https://slurm.schedmd.com/cpu_management.html), [SciPy Pearson p-value definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html).
