from nextstrain_ceirr.ceirr import *
from na_hpai import *

"""This file specifies the entire avian-flu pipeline that will be run, with
specific parameters for subsampling, tree building, and visualization. In this
build, you will generate 1 tree: an H5N1 tree for the HA genes. In this simple
build, clade annotation has been removed. This template should provide a
reasonable starting place to customize your own build. Simply edit and add
components to this Snakefile."""


"""Here, define your wildcards. To include more subtypes or gene segments, simply
add those to these lists, separated by commas"""
# it will take too long to run all 8 gene segments, but I've included all the files to do so.
# the other segments are pb1, pa, ha, np, mp, and ns
SEGMENTS = ["pb1","pb2","na","pa","ha","np","mp","ns"]
NUMBER_OF_GENOTYPES = 15
GENOTYPES_TO_INCLUDE = ['D1.1', 'D1.2']

"""This rule tells Snakemake that at the end of the pipeline, you should have
generated JSON files in the auspice folder for each subtype and segment."""
rule all:
    input:
        auspice_json = expand("auspice/h5nx_{segment}.json", segment=SEGMENTS)

"""Specify all input files here. For this build, you'll start with input sequences
from the example_data folder, which contain metadata information in the
sequence header. Specify here files denoting specific strains to include or drop,
references sequences, and files for auspice visualization (colors)"""
rule files:
    params:
        input_sequences = "data/h5nx/{segment}/sequences.fasta",
        input_metadata = "data/h5nx/metadata-with-clade.tsv",
        dropped_strains = "config/exclude_strains.txt",
        include_strains = "config/include_strains.txt",
        reference = "config/reference_sequence_{segment}_A_goose_CR_2021.gb", #H3N8 from 1997
        auspice_config = "config/auspice_config.json",
        colors = "config/colors.tsv",
        description = "config/description.md"


files = rules.files.params

"""This rule unzips the h5nx.zip from the h5-data-updates submodule.
"""
rule unzip_h5_data:
    input:
        "h5-data-updates/h5nx.zip"
    output:
        metadata=files.input_metadata,
        sequences=expand(files.input_sequences, segment=SEGMENTS)
    shell:
        'unzip -o h5-data-updates/h5nx.zip -d data/'

"""GenoFlu requires all segments be present in order to proceed with an annotation.
This organizes the data into a subset of complete genomes accordingly.
"""
rule genoflu_dataflow:
    input:
        rules.unzip_h5_data.output.sequences
    output:
        expand("data/genoflu/{segment}.fasta", segment=SEGMENTS)
    run:
        genoflu_dataflow()

"""Compute GenoFlu annotations by running GenoFlu-multi on subset of complete genomes.
"""
rule genoflu_run:
    input:
        rules.genoflu_dataflow.output
    output:
        'data/genoflu/results/results.tsv'
    shell:
        '''
            # this avoids a quirk of the GenoFlu package... avoids UnboundLocalError related to excel_stats
            rm -rf data/genoflu/temp/
            python GenoFLU-multi/bin/genoflu-multi.py -n 12 -f data/genoflu
        '''

"""This post-processes the GenoFlu data to include only lineages that either dominate
or are specifically included to avoid visual clutter.
"""
rule genoflu_postprocess:
    input:
        metadata=files.input_metadata,
        genoflu=rules.genoflu_run.output[0]
    output:
        metadata='results/metadata-with-genoflu.tsv',
        counts='results/genoflu/results/counts.tsv'
    run:
        genoflu_postprocess(
            input.metadata, input.genoflu, output.metadata, output.counts,
            NUMBER_OF_GENOTYPES, GENOTYPES_TO_INCLUDE
        )

rule metadata_annotation:
    input:
        rules.genoflu_postprocess.output.metadata
    output:
        'results/metadata.tsv',
    run:
        metadata_annotation(input[0], output[0])

"""In this section of the Snakefile, rules are specified for each step of the pipeline.
Each rule has inputs, outputs, parameters, and the specific text for the commands in
bash. Rules reference each other, so altering one rule may require changing another
if they depend on each other for inputs and outputs. Notes are included for
specific rules."""


"""The parse rule is used to separate out sequences and metadata into 2 distinct
files. This rule assumes an input fasta file that contains metadata information
in the header. By specifying the order of those fields in the `fasta_fields` line,
`augur parse` will separate those fields into labeled columns in the output metadata
file."""


"""The minimum length required for sequences. Sequences shorter than these will be
subsampled out of the build. Here, we're requiring all segments to be basically
complete. To include partial genomes, shorten these to your desired length"""

def min_length(w):
    len_dict = {"pb2": 2100, "pb1": 2100, "pa": 2000, "ha":1600, "np":1400, "na":1270, "mp":900, "ns":800}
    length = len_dict[w.segment]
    return(length)


"""This rule specifies how to subsample data for the build, which is highly
customizable based on your desired tree."""
rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
          - samples with missing region and country metadata
          - excluding strains prior to {params.min_date}
        """
    input:
        sequences = files.input_sequences,
        metadata = rules.metadata_annotation.output[0],
        include = files.include_strains,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered_{segment}.fasta"
    params:
        group_by = "month host location",
        sequences_per_group = 25,
        min_date = 2021,
        min_length = min_length,  # instead of specifying one parameter value, we can use a function to specify minimum lengths that are unique to each segment
        exclude_where = "host=laboratoryderived host=ferret host=unknown host=other host=host country=? region=? region!='north america' h5_label_clade!=2.3.4.4b"

    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --exclude-where {params.exclude_where} \
            --min-length {params.min_length} \
            --non-nucleotide
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{segment}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --nthreads 1
        """


rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{segment}.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads 1
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.metadata_annotation.output[0]
    output:
        tree = "results/tree_{segment}.nwk",
        node_data = "results/branch-lengths_{segment}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt-muts_{segment}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}\
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa-muts_{segment}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.metadata_annotation.output[0]
    output:
        node_data = "results/traits_{segment}.json",
    params:
        columns = "host region country division flyway Domestic_Status",
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

"""This makes a segment specific config for GenoFlu segment lineages.
"""
rule auspice_config:
    input:
        files.auspice_config
    output:
        "config/{segment}/auspice_config.json"
    run:
        auspice_segment_config(input[0], output[0], wildcards.segment)

"""This rule exports the results of the pipeline into JSON format, which is required
for visualization in auspice. To make changes to the categories of metadata
that are colored, or how the data is visualized, alter the auspice_config files"""
rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.metadata_annotation.output[0],
        node_data = [rules.refine.output.node_data,rules.traits.output.node_data,rules.ancestral.output.node_data,rules.translate.output.node_data],
        auspice_config = rules.auspice_config.output[0],
        colors = files.colors,
        description = files.description

    output:
        auspice_json = "auspice/h5nx_{segment}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data}\
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --include-root-sequence \
            --colors {input.colors} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
