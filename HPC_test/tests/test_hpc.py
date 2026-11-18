"""Numerical regression and pipeline tests; run: python -m unittest discover -s HPC_test/tests -v."""
import importlib.util
import sqlite3
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import subprocess
import numpy as np
import pandas as pd
from scipy import stats

spec = importlib.util.spec_from_file_location('hpc', Path(__file__).resolve().parents[1] / 'corrMultiomics_HPC.py')
hpc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hpc)


def options(**kwargs):
    defaults = dict(method='pearson', mad_filter=False, remove_zeros=False,
                    threads=2, block_size=3, correlation_cutoff=0.2, edge_weight_cutoff=0.01,
                    multi_decay_rates=[10, 25], decay_rate=5, clusterone=False, mutual_rank=True,
                    plot=False, annotation='', scratch_dir=None, mr_storage='disk', mr_memory_gb=None, sqlite_table_name='test', if_exists='fail')
    defaults.update(kwargs)
    return SimpleNamespace(**defaults)


class NumericalTests(unittest.TestCase):
    def test_methods_match_scipy(self):
        rng = np.random.RandomState(22)
        raw_x, raw_y = rng.rand(7, 13) + 0.1, rng.rand(11, 13) + 0.1
        # Include ties for Spearman.
        for method in ['pearson', 'pearsonlog', 'spearman']:
            a = np.round(raw_x, 1) if method == 'spearman' else raw_x
            b = np.round(raw_y, 1) if method == 'spearman' else raw_y
            opt = options(method=method, correlation_cutoff=0)
            x, _ = hpc.prepare(pd.DataFrame(a), opt)
            y, _ = hpc.prepare(pd.DataFrame(b), opt)
            for _, _, i, j, r, p in hpc.blocks(x, y, opt):
                for k in range(len(i)):
                    left, right = a[i[k]], b[j[k]]
                    if method == 'pearsonlog':
                        left, right = np.log10(left), np.log10(right)
                    expected = stats.spearmanr(left, right) if method == 'spearman' else stats.pearsonr(left, right)
                    np.testing.assert_allclose([r[k], p[k]], expected, atol=1e-12)

    def test_input_formats(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'matrix.csv'
            for sep in [',', ';']:
                path.write_text(f'id{sep}s1{sep}s2\na{sep}1{sep}2\n')
                self.assertEqual(hpc.read_matrix(path).shape, (1, 2))
            path.write_bytes(b'id,s1\na,\x82')
            with self.assertRaisesRegex(ValueError, 'UTF-8'):
                hpc.read_matrix(path)

    def test_sql_ranks_match_reference_with_overlapping_roles(self):
        # Repeated orientations emulate transcript-vs-transcript inputs.
        edges = [(0, 1, .9), (1, 0, .9), (0, 2, .8), (1, 2, -.7), (2, 2, 1.)]
        with tempfile.TemporaryDirectory() as tmp, sqlite3.connect(':memory:') as scratch, sqlite3.connect(':memory:') as out:
            scratch.execute('CREATE TABLE edges (a INTEGER,b INTEGER,correlation REAL)')
            scratch.executemany('INSERT INTO edges VALUES (?,?,?)', edges)
            hpc.mutual_rank(scratch, out, options(), np.array(['a', 'b', 'c']), tmp)
            result = pd.read_sql_query('SELECT * FROM test_MR_edges', out)
            rank = {}
            for node in range(3):
                neighbors = [(b, r) for a, b, r in edges if a == node] + [(a, r) for a, b, r in edges if b == node]
                rank[node] = {neighbor: i for i, (neighbor, _) in enumerate(sorted(neighbors, key=lambda v: (-v[1], v[0])), 1)}
            for row in result.itertuples():
                a, b = 'abc'.index(row.metabolite), 'abc'.index(row.gene)
                self.assertAlmostEqual(row.MR, np.sqrt(rank[a][b] * rank[b][a]))
            self.assertEqual(len(result), 3)

    def test_pipeline_and_rerun_protection(self):
        rng = np.random.RandomState(4)
        with tempfile.TemporaryDirectory() as tmp:
            x, y = Path(tmp)/'x.csv', Path(tmp)/'y.csv'
            pd.DataFrame(rng.randn(6, 20), index=['a','b','c','d','e','f']).to_csv(x)
            pd.DataFrame(rng.randn(9, 20), index=[f'g{i}' for i in range(9)]).to_csv(y)
            opt = options(feature_table=x, quantitation_matrix=y, sqlite_db_name=str(Path(tmp)/'out.sqlite'))
            hpc.main(opt)
            with sqlite3.connect(opt.sqlite_db_name) as conn:
                df = pd.read_sql_query('SELECT * FROM test', conn)
                self.assertGreater(len(df), 0)
                self.assertTrue((df.correlation.abs() >= .2).all())
                count = conn.execute('SELECT count(*) FROM test_MR_edges').fetchone()[0]
                self.assertEqual(count, len(df))
            with self.assertRaisesRegex(ValueError, 'already exist'):
                hpc.main(opt)

    def test_ram_disk_equivalence(self):
        rng = np.random.RandomState(52)
        # Include ties, different duplicate scores, both orientations, self edges,
        # isolated nodes, and negative correlations.
        edges = [(int(a), int(b), float(r)) for a, b, r in
                 zip(rng.randint(0, 20, 200), rng.randint(0, 20, 200),
                     rng.choice([-.9, -.2, .2, .8, 1.], 200))]
        for data in [edges, [(0, 0, 1.)], [(1, 4, -.4)], []]:
            with self.subTest(edges=len(data)), tempfile.TemporaryDirectory() as tmp, \
                    sqlite3.connect(':memory:') as scratch, sqlite3.connect(':memory:') as disk, \
                    sqlite3.connect(':memory:') as ram:
                scratch.execute('CREATE TABLE edges (a INTEGER,b INTEGER,correlation REAL)')
                scratch.executemany('INSERT INTO edges VALUES (?,?,?)', data)
                nodes = np.array([f'node {i}' for i in range(21)])
                hpc.mutual_rank(scratch, disk, options(), nodes, tmp)
                # Deliberately split edges across chunks, like correlation blocks.
                chunks = [np.array(data[i:i+13], dtype=hpc.EDGE_DTYPE) for i in range(0, len(data), 13)]
                hpc.ram_mutual_rank(chunks, ram, options(mr_storage='ram'), nodes, tmp)
                self.assertEqual(chunks, [])
                tables = [r[0] for r in disk.execute("SELECT name FROM sqlite_master WHERE type='table'")]
                self.assertEqual(tables, [r[0] for r in ram.execute("SELECT name FROM sqlite_master WHERE type='table'")])
                for table in tables:
                    expected = pd.read_sql_query('SELECT * FROM ' + table, disk).sort_values(['metabolite','gene']).reset_index(drop=True)
                    actual = pd.read_sql_query('SELECT * FROM ' + table, ram).sort_values(['metabolite','gene']).reset_index(drop=True)
                    pd.testing.assert_frame_equal(expected, actual)

    def test_ram_budget_and_ram_pipeline(self):
        with self.assertRaisesRegex(ValueError, 'exceeding'):
            hpc.check_ram_budget(1000000, .2)
        rng = np.random.RandomState(16)
        with tempfile.TemporaryDirectory() as tmp:
            x, y = Path(tmp)/'x.csv', Path(tmp)/'y.csv'
            pd.DataFrame(rng.randn(8, 20)).to_csv(x)
            pd.DataFrame(rng.randn(10, 20)).to_csv(y)
            opt = options(feature_table=x, quantitation_matrix=y, sqlite_db_name=str(Path(tmp)/'out.sqlite'),
                          mr_storage='ram', mr_memory_gb=1.)
            hpc.main(opt)
            with sqlite3.connect(opt.sqlite_db_name) as conn:
                self.assertGreater(conn.execute('SELECT count(*) FROM test_MR_edges').fetchone()[0], 0)
            # RAM storage must not create a scratch SQLite edge database.
            self.assertFalse(list(Path(tmp).glob('**/ranks.sqlite')))

    def test_cluster_failures_and_name_mapping(self):
        with tempfile.TemporaryDirectory() as tmp, sqlite3.connect(':memory:') as conn:
            out = Path(tmp)/'out.csv'
            opt = options(clusterone_jar='jar', java_heap='1g')
            with patch.object(hpc.subprocess, 'run', return_value=subprocess.CompletedProcess([], 1, stderr='bad network')):
                with self.assertRaisesRegex(RuntimeError, 'bad network'):
                    hpc.cluster(opt, 'weights', out, np.array(['Orotic acid', 'Ankrd49']), conn, 'clusters')
            def fake_run(cmd, stdout, **kwargs):
                stdout.write('Cluster,Members\n1,0 1\n')
                return subprocess.CompletedProcess(cmd, 0, stderr='')
            with patch.object(hpc.subprocess, 'run', side_effect=fake_run):
                hpc.cluster(opt, 'weights', out, np.array(['Orotic acid', 'Ankrd49']), conn, 'clusters')
            self.assertEqual(conn.execute('SELECT Members, Member_IDs FROM clusters').fetchone(), ('Orotic acid Ankrd49', '0 1'))


if __name__ == '__main__':
    unittest.main()
