# nextstrain-hpai-north-america

This repository contains files for constructing nextstrain builds with additional metadata designations relevant to the current high path avian influenza in North America found [here](https://nextstrain.org/groups/moncla-lab/h5nx/north-america/ha). For a walk through of these data and some of these additional metadata analyses, check out our [nexstrain narrative](https://nextstrain.org/groups/moncla-lab/narratives/h5nx/north-america-2021-present)!

Metadata files have been modified to add migratory flyway assignments, as well as higher resolution species and order data, using the script and csv files available within the metadata_mod_scripts folder. 
"species.csv" provide order assigments to sequences based on parsing out the "animal" from the strain name (meta_mod scripts provide code to do this). Further groupings are done by creating lists of orders and defining these via lists in this script and a new column called "species-grouped". 

2 versions of the metadata modifying scripts are available here for if you are looping through with a separate metadata file for each gene (looped) or a single metadata file (1file). Wildcards and input files will need to be adjusted accordingly in the Snakefile.

Within this script, domestic status is updated via matching by sequence ID to a file "NA-H5Nx-2021-2023-seqmerge.tsv" with high resolution of domestic status data. This file however cannot be shared due to strain information derived from GISAID.
