`postal` is a command line program useful for creating pipelines for the
initial post-analysis of seqFish data.   It consists of the following modules:

1.  `mkad`  Make an anndata object from seqFish output.
2.  `qc`  Run QC on anndata object.
3. `filter`   Filter cells and genes.
4. `latent`  Calcuate latent representation of cell data.
5. `cluster`  Cluster cells.

## Install postal

1.  Install [poetry](https://python-poetry.org/docs/)

```
curl -sSL https://install.python-poetry.org | python3 -

poetry --version
```

2.  Install postal

```
poetry install
```

3.  Test

```
pytest
```

## Run postal

The input to postal is a configuration yaml file.  The yaml files consists of a list
modules, one module for each subcommand you want to run.  For example,
to run the `mkad` subcommand to make an anndata file, the yaml file needs to have the
following entry:

```yaml
- module: mkad
  arguments:
    counts_file:  your_counts.csv
    cell_data_file:  your_cell_data.csv
    transcript_data_file:  your_transcript_data.csv
    decode_file: your_decode.hdf5
    outs:  your_output_directory
	anndata_file:  your_output_anndata.h5ad
```

List the available subcommands:

```
postal list
```

Run a particular subcommand:

```
postal mkad config.yaml
```


