# NEEAT benchmark

This subfolder contains data and a Snakemake workflow to run the NEEAT filter with various parameter settings.

## Install the software dependencies

The software required to run the NEEAT benchmark are specified in the `pixi.toml` file. To make use of this file, first install [pixi](https://pixi.sh/latest/). Then run the following command to start a bash shell with the software environment:

```bash
pixi shell
```

### Alternative installation

You can also create a conda environment with the required software:

```bash
conda env create -f environment.yml
```

Then activate the environment with 
    
```bash 
conda activate neeat-benchmark
```

## Prepare input data

The input data was created from a dataset of ASVs generated from samples collected in Sweden, clustered with swarm (`d=15` setting).

1. Extract the tar archive at `data/benchmakr-data.tar.gz` to the `data` folder.
```bash
tar -xvzf data/benchmark-data.tar.gz -C data
```

Also edit the `config.yml` file so that `taxonomy_file` points to the cluster taxonomy file for the correct data.

## Run the workflow

1. If running locally on a laptop:

```bash
snakemake --cores 4 -k -p --rerun-triggers mtime --rerun-incomplete 
```

2. If running on a cluster with the SLURM workload manager, edit the `myprofile/config.yaml` file to add your SLURM account:

```yaml
default-resources:
  slurm_account: "your_slurm_account"
```

Then run the workflow with the following command:

```bash
snakemake --profile myprofile 
```

## Output

- **raw_results/<order>/** contains lists of discarded OTUs/clusters for each order/parameter combination, _e.g_ `raw_results/Diptera/evo_run1_discarded_otus.tsv` contain discarded OTUs from run 1 of the EVO filter for the Diptera order.
- **results/<order>/** contains evaluation results for each order/parameter combination, _e.g_ `results/Diptera/evo_res_1.tsv` contain evaluation results from run 1 of the EVO filter for the Diptera order.
- **results/<method>_res.tsv** are aggregated results for each method (`evo`, `echo`, `abundance`, `ee`, `eee`, `eeea`) including the combined **Hexapoda** results.