# rtk2

`rtk2` is a C++ command-line toolkit for building, transforming, summarizing,
normalizing, and rarefying large gene-abundance matrices. It also contains
specialized modes for functional aggregation, gene-cluster decluttering, and
abundance-aware discovery of mosaic gene clusters.

The current source reports version 2.28.

## Build

The build requires a C++11 compiler, `make`, zlib, and pthread support.

```text
make -j
./rtk2 -h
./rtk2 -v
```

`make clean` removes the executable and object files.

## Matrix convention

Most matrix modes expect a delimited text table with:

- a header row containing sample names;
- a row identifier, normally a gene or functional identifier, in column 1;
- one non-negative abundance or count per sample.

Tab-delimited input is the normal convention. Several high-volume paths accept
gzip-compressed input and produce indexed, parallel gzip output. Exact secondary
inputs and output naming depend on the selected mode; run `rtk2 -h` for the
built-in option summary.

## Modes

| Mode | Purpose |
| --- | --- |
| `memory`, `swap`, `rarefaction`, `rare_inmat` | Rarefy count matrices and calculate diversity measures. `memory` keeps the matrix in memory; `swap` reduces peak memory with temporary files. |
| `geneMat` | Build a gene-abundance matrix from CD-HIT-style clusters and mapping/coverage inputs. |
| `normalize` | Normalize each sample column, using sparse in-memory processing by default or `-normalizeReadPass` for two-pass streaming. |
| `mergeMat` | Outer-join matrices by row ID and append their sample columns. |
| `sumMat` | Aggregate rows through a supplied functional or taxonomic hierarchy. |
| `lineExtr`, `rowSubset` | Extract a requested set of matrix rows. |
| `submatrices` | Create multiple row-selected matrices from a reference definition. |
| `colSums` | Report sample-column totals in original and sorted order. |
| `rowOcc` | Report row occurrence statistics. |
| `module` | Estimate module abundance and completeness from a KO matrix and module database. |
| `decluter`, `decluter2` | Consolidate related gene clusters using abundance co-exclusion evidence. |
| `mosaic` | Report sequence-related representatives whose complementary abundance supports a split locus. See [Mosaic mode](docs/MOSAIC_MODE.md). |

## Common examples

Normalize a matrix:

```text
rtk2 normalize -i Matrix.mat.gz -o Matrix.normalized -gz -t 8
```

Outer-join matrices:

```text
rtk2 mergeMat -i first.tsv.gz,second.tsv.gz -o merged.tsv.gz
```

Rarefy a table to 10,000 counts per sample:

```text
rtk2 memory -i counts.tsv -o rarefied/run. -d 10000 -r 10 -w 1 -t 8
```

Find candidate mosaic loci:

```text
rtk2 mosaic \
  -i Matrix.mat.gz \
  -reference SB.clusters.core \
  -paf catalog.self.paf.gz \
  -groups gene_to_nog.tsv \
  -o mosaic/run
```

The mosaic mode is report-only: it does not rewrite the abundance matrix. Its
detailed input contracts, evidence model, thresholds, outputs, and MATAFILER4
integration are documented in [docs/MOSAIC_MODE.md](docs/MOSAIC_MODE.md).

## Development and tests

Build with warnings enabled using the supplied `Makefile`:

```text
make -j2
```

Run the mosaic regression:

```text
tests/test_mosaic.sh
```

`testCases.txt` contains historical, environment-specific command examples. It is
retained as test/development material rather than user documentation.
