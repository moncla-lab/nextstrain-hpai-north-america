from na_hpai import *

wildcard_constraints:
    region="[^/]+",
    subset="[^/]+"

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

"""This rule tells Snakemake that at the end of the pipeline, you should have
generated JSON files in the auspice folder for each subtype and segment."""
rule all:
    input:
        auspice_json = expand(
            "auspice/{region}_{segment}.json",
            region=['north-america-only', 'global-context'],
            segment=SEGMENTS
        ),
        tip_frequencies = expand(
            "auspice/{region}_{segment}_tip-frequencies.json",
            region=['north-america-only', 'global-context'],
            segment=SEGMENTS
        )

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
        reference = "results/{segment}/reference_sequence.gb",
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

rule download_references:
    message: "Downloading RefSeq reference sequences from NCBI"
    input:
        references_config = "config/references.tsv"
    output:
        reference_files = expand("results/{segment}/reference_sequence.gb", segment=SEGMENTS)
    run:
        from na_hpai import download_refseq_references
        download_refseq_references(
            references_tsv=input.references_config,
            output_dir="results"
        )

rule metadata_annotation:
    input:
        files.input_metadata
    output:
        'results/metadata.tsv',
    run:
        metadata_annotation(input[0], output[0])

rule generate_genotype_colors:
    input:
        rules.metadata_annotation.output[0]
    output:
        'results/genotype_colors.tsv'
    run:
        generate_genotype_colors(input[0], output[0])

rule combine_colors:
    input:
        static = files.colors,
        genotype = rules.generate_genotype_colors.output[0]
    output:
        'results/colors.tsv'
    shell:
        'cat {input.static} {input.genotype} > {output}'

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


each_exclude = "host=laboratoryderived host=ferret host=unknown host=other host=host country=? region=?"
asia_exclude = "region='north america' region='africa' region='antarctica' region='south america' region='europe'"

exclude_where = {
    'north-america-only': {
        'na-cattle-resolved': None,  # Uses --query instead
        'na-cattle-unresolved': None,  # Uses --query instead
        'na-noncattle': "region!='north america' host='Cattle'",
    },
    'global-context': {
        'h5n5': None,  # Uses --query instead
        'na-cattle-resolved': "subtype='h5n5'",  # Exclude H5N5 (handled by h5n5 bucket)
        'na-cattle-unresolved': "subtype='h5n5'",
        'na-noncattle': "region!='north america' host='Cattle' subtype='h5n5'",
        'sa': "region!='south america' subtype='h5n5'",
        'europe': "region!='europe' subtype='h5n5'",
        'asia': asia_exclude + " subtype='h5n5'"
    }
}

# Define which subsets to include in each build
SUBSETS_BY_REGION = {
    'north-america-only': ['na-cattle-resolved', 'na-cattle-unresolved', 'na-noncattle'],
    'global-context': ['h5n5', 'na-cattle-resolved', 'na-cattle-unresolved', 'na-noncattle', 'sa', 'europe', 'asia']
}


def exclude_by_region(wildcards):
    subset_exclude = exclude_where[wildcards.region][wildcards.subset]
    if subset_exclude is None:
        return each_exclude
    return each_exclude + ' ' + subset_exclude


def sequences_per_group(wildcards):
    spg_dict = {
        'h5n5': 500,                  # High sampling for rare H5N5 subtype
        'na-cattle-resolved': 100,    # High sampling for well-resolved cattle
        'na-cattle-unresolved': 1000, # Grab most poorly-resolved cattle
        'na-noncattle': 25,           # Normal sampling for other North American sequences
        'sa': 10,
        'europe': 1,
        'asia': 5
    }
    return spg_dict[wildcards.subset]


def group_by_strategy(wildcards):
    """Use different grouping for unresolved cattle (no month/location)"""
    if wildcards.subset == 'na-cattle-unresolved':
        return "host subtype"  # One big bucket for poorly-resolved cattle
    elif wildcards.subset == 'h5n5':
        return "month host location"  # No subtype grouping needed (already filtered to h5n5)
    else:
        return "month host location subtype"  # Normal grouping for others


def query_param(wildcards):
    """Return full --query argument string. All subsets get clade filter, some get additional filters."""
    clade_filter = "h5_label_clade.isin(['2.3.4.4b', 'Am-nonGsGD'])"
    subset_queries = {
        'h5n5': "subtype == 'h5n5'",
        'na-cattle-resolved': "host == 'Cattle' and region == 'North America' and location != 'Usa' and not date.str.contains('XX')",
        'na-cattle-unresolved': "host == 'Cattle' and region == 'North America' and (location == 'Usa' or date.str.contains('XX'))",
    }
    subset_query = subset_queries.get(wildcards.subset)
    if subset_query:
        return f'--query "{subset_query} and {clade_filter}"'
    else:
        return f'--query "{clade_filter}"'


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
        sequences = "results/{region}/{subset}/filtered_{segment}.fasta"
    params:
        group_by = group_by_strategy,  # Dynamic grouping strategy
        sequences_per_group = sequences_per_group,
        min_date = 2021,
        min_length = min_length,  # instead of specifying one parameter value, we can use a function to specify minimum lengths that are unique to each segment
        exclude_where = exclude_by_region,
        query = query_param  # Returns full --query arg for cattle subsets, empty string for others

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
            --non-nucleotide \
            {params.query} || true
        """

rule concatenate:
    message:
        """
        Concatenating {wildcards.region} to full FASTA
        """
    input:
        lambda wildcards: expand(
            "results/{{region}}/{subset}/filtered_{{segment}}.fasta",
            subset=SUBSETS_BY_REGION[wildcards.region]
        )
    output:
        sequences="results/{region}/filtered_{segment}.fasta"
    shell:
        "cat {input} | seqkit rmdup -n > {output.sequences}"

rule metadata_post_filter:
    """Extract metadata for strains that survived filtering and concatenation."""
    input:
        sequences = rules.concatenate.output.sequences,
        metadata = rules.metadata_annotation.output[0]
    output:
        metadata = "results/{region}/metadata_post_filter_{segment}.tsv"
    run:
        extract_metadata_post_filter(input.sequences, input.metadata, output.metadata)

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.concatenate.output.sequences,
        reference = files.reference
    output:
        alignment = "results/{region}/aligned_{segment}.fasta"
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
        tree = "results/{region}/tree-raw_{segment}.nwk"
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
        tree = "results/{region}/tree_{segment}.nwk",
        node_data = "results/{region}/branch-lengths_{segment}.json"
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

rule metadata_post_refine:
    """Extract metadata for strains that survived clock filtering."""
    input:
        tree = rules.refine.output.tree,
        metadata = rules.metadata_annotation.output[0]
    output:
        metadata = "results/{region}/metadata_post_refine_{segment}.tsv"
    run:
        extract_metadata_post_refine(input.tree, input.metadata, output.metadata)

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/{region}/nt-muts_{segment}.json"
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
        node_data = "results/{region}/aa-muts_{segment}.json"
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
        node_data = "results/{region}/traits_{segment}.json",
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
        "results/{segment}/auspice_config.json"
    run:
        auspice_segment_config(input[0], output[0], wildcards.segment)

"""Calculate tip frequencies for visualization"""
rule tip_frequencies:
    message: "Estimating tip frequencies for {wildcards.segment}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.metadata_annotation.output[0]
    output:
        "auspice/{region}_{segment}_tip-frequencies.json"
    params:
        min_date = "2022",
        pivot_interval = 1
    shell:
        """
        augur frequencies \
            --method kde \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --min-date {params.min_date} \
            --pivot-interval {params.pivot_interval} \
            --output {output}
        """

"""This is a custom rule developed for the avian influenza builds and is not part
of the Nextstrain architecture. It uses custom python scripts to determine the
sequence of amino acids at the HA cleavage site, and annotate those sequences
for whether they contain a furin cleavage site."""
rule cleavage_site:
    message: "determining sequences that harbor furin cleavage sites"
    input:
        alignment = "results/{region}/aligned_ha.fasta"
    output:
        cleavage_site_annotations = "results/{region}/cleavage-site_ha.json",
        cleavage_site_sequences = "results/{region}/cleavage-site-sequences_ha.json"
    shell:
        """
        python scripts/annotate-ha-cleavage-site.py \
            --alignment {input.alignment} \
            --furin_site_motif {output.cleavage_site_annotations} \
            --cleavage_site_sequence {output.cleavage_site_sequences}
        """


"""This function allows us to annotate HA sequences with cleavage site information,
without trying to apply it to the other segments"""
def node_data_by_wildcards(w):
    """for ha, include cleavage site data during export; for other segments, do not"""
    if w.segment == "ha":
        node_data = [rules.refine.output.node_data,rules.traits.output.node_data,rules.ancestral.output.node_data,rules.translate.output.node_data,rules.cleavage_site.output.cleavage_site_annotations,rules.cleavage_site.output.cleavage_site_sequences]
    else:
        node_data = [rules.refine.output.node_data,rules.traits.output.node_data,rules.ancestral.output.node_data,rules.translate.output.node_data]
    return(node_data)


"""This rule exports the results of the pipeline into JSON format, which is required
for visualization in auspice. To make changes to the categories of metadata
that are colored, or how the data is visualized, alter the auspice_config files"""
rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.metadata_annotation.output[0],
        metadata_post_filter = rules.metadata_post_filter.output.metadata,
        metadata_post_refine = rules.metadata_post_refine.output.metadata,
        node_data = node_data_by_wildcards,
        auspice_config = rules.auspice_config.output[0],
        colors = rules.combine_colors.output,
        description = files.description

    output:
        auspice_json = "auspice/{region}_{segment}.json"
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

onstart:
    shell("date '+Start: %Y-%m-%d %H:%M:%S' > time.txt")

onsuccess:
    shell("date '+End: %Y-%m-%d %H:%M:%S' >> time.txt")
