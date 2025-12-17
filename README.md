# nextstrain-hpai-north-america

This repository contains files for constructing nextstrain builds with additional metadata designations relevant to the current high path avian influenza in North America found [here](https://nextstrain.org/groups/moncla-lab/h5nx/north-america/ha). For a walk through of these data and some of these additional metadata analyses, check out our [nexstrain narrative](https://nextstrain.org/groups/moncla-lab/narratives/h5nx/north-america-2021-present)!

## Installation

Make sure you have installed [Bioconda](https://bioconda.github.io/) and configured it correctly as described in the link.

Install dependencies:

```
conda env create -f environment.yml
```

Clone with submodules:

```
git clone https://github.com/moncla-lab/nextstrain-hpai
cd nextstrain-hpai-north-america
git submodule update --init
```

## Usage

### Activate the environment

```bash
conda activate nextstrain-hpai
```

### Run pipeline:

To set the number of jobs to run, use the `-j $NUMBER_OF_JOBS` flag. To set the number of cores to run in parallel, use `--cores` with 4, 8, or 10 depending on your computer. 

```bash
snakemake --cores $NUMBER_OF_CORES all
```

### View results:
```bash
auspice view --datasetDir auspice
```

### Upload results

```
nextstrain remote upload nextstrain.org/groups/moncla-lab/h5nx-northamerica ./auspice/*
```

## Customizing Reference Sequences

The pipeline automatically downloads reference sequences from NCBI based on the configuration in `config/references.tsv`.

### Default References

By default, the pipeline uses A/Goose/Guangdong/1/96 (H5N1) RefSeq sequences (NC_* accessions), which are curated by NCBI and recommended for robust phylodynamics analysis.

### Using Different References

To use different reference sequences:

1. Edit `config/references.tsv` and replace the `ncbi_accession` values with your desired NCBI accession numbers
2. Re-run the pipeline - it will automatically download the new references

Example `config/references.tsv`:
```
segment	ncbi_accession
pb2	YOUR_ACCESSION_1
pb1	YOUR_ACCESSION_2
pa	YOUR_ACCESSION_3
ha	YOUR_ACCESSION_4
np	YOUR_ACCESSION_5
na	YOUR_ACCESSION_6
mp	YOUR_ACCESSION_7
ns	YOUR_ACCESSION_8
```

The pipeline will download these sequences and save them as GenBank files in the `results/` directory following the pattern: `results/{segment}/reference_sequence.gb`

## Submodule Management

This repository utilizes [Git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to reuse code that is used in several projects.

### Updating submodules when data changes:
```
git submodule update --remote
```
Then add, commit, and push like usual.

### For maintainers:
- **h5-data-updates**: A submodule that contains shared Genoflu analysis functions. Changes here affect both CEIRR and North America pipelines.
- **metadata\_mod\_scripts**: Contains flyway and species classification data specific to North America analysis.

## Metadata Enhancement

Metadata files have been modified to add migratory flyway assignments, as well as higher resolution species and order data, using the script and csv files available within the metadata\_mod\_scripts folder. 
"species.csv" provide order assigments to sequences based on parsing out the "animal" from the strain name (meta_mod scripts provide code to do this). Further groupings are done by creating lists of orders and defining these via lists in this script and a new column called "species-grouped". 

2 versions of the metadata modifying scripts are available here for if you are looping through with a separate metadata file for each gene (looped) or a single metadata file (1file). Wildcards and input files will need to be adjusted accordingly in the Snakefile.

Within this script, domestic status is updated via matching by sequence ID to a file "NA-H5Nx-2021-2023-seqmerge.tsv" with high resolution of domestic status data. This file however cannot be shared due to strain information derived from GISAID.

## Domesticity Classification Strategy

We implement a two-tiered strategy to classify avian samples as domestic or wild:

 1. Per-Sample GISAID Annotations (Highest Priority)
	- Uses existing domestic_status field from Fauna/GISAID metadata
	- Only accepts clean "wild"/"domestic" values (skips "?" markers)
	- Covers < 10% of samples, highest accuracy when present

2. Species-Level Classifications (Fallback)
  - Uses hand-curated species lookup table (config/species_lookup.tsv)
  - Maps animal names (e.g., "goose" → "wild", "chicken" → "domestic")
  - Covers more samples, strategically targetting high-frequency species
