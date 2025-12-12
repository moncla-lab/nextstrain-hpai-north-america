import json

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
import matplotlib.colors as mcolors

# Import species annotation functions from standalone script
from scripts.annotate_species_metadata import annotate_metadata as _annotate_metadata

SEGMENTS = ["pb1", "pb2", "na", "pa", "ha", "np", "mp", "ns"]


def metadata_annotation(input_metadata_tsv, output_metadata_tsv):
    """
    Annotate metadata with species information, flyways, and domesticity.

    This is a wrapper around the standalone annotate_species_metadata.py script.
    The standalone script can be used by other builds as a submodule.

    Pipeline-specific processing (GenoFLU) is added after species annotation.
    """
    # Use standalone script for species annotation
    _annotate_metadata(
        input_path=input_metadata_tsv,
        output_path=output_metadata_tsv,
        config_dir='config',
        flyway_path='metadata_mod_scripts/flyway_regions.csv'
    )

    # Add pipeline-specific GenoFLU postprocessing
    metadata = pd.read_csv(output_metadata_tsv, sep="\t")
    metadata = genoflu_postprocess(metadata)
    metadata.to_csv(output_metadata_tsv, index=False, sep="\t")


def auspice_segment_config(input_config_path, output_config_path, segment):
    with open(input_config_path) as input_file:
        config = json.load(input_file)

    segment_key = segment.upper()
    genoflu_segment_config = {
        "key": f"genoflu_{segment_key}",
        "title": f"GenoFlu {segment_key} lineage",
        "type": "categorical",
    }
    config["colorings"] += [genoflu_segment_config]

    # Add FCS colorings for HA segment only
    if segment == "ha":
        fcs_annotation_config = {
            "key": "furin_cleavage_motif",
            "title": "Furin cleavage site",
            "type": "categorical",
        }
        fcs_sequence_config = {
            "key": "cleavage_site_sequence",
            "title": "Cleavage site sequence",
            "type": "categorical",
        }
        config["colorings"] += [fcs_annotation_config, fcs_sequence_config]

    with open(output_config_path, "w") as output_file:
        json.dump(config, output_file, indent=2)


def genoflu_refine_genotype(row):
    """Refine genotype: keep A*/B*/C*/D* patterns, collapse Minor*, map unassigned by region"""
    genotype = row["genoflu"]

    # Collapse all Minor* to just "Minor"
    if isinstance(genotype, str) and genotype.startswith("Minor"):
        return "Minor"

    # Check if it's empty or unassigned
    if (
        pd.isna(genotype)
        or genotype == ""
        or "Not assigned" in str(genotype)
        or "Unseen constellation" in str(genotype)
    ):
        # Check region to determine Americas vs not
        region = row.get("region", "")
        if "America" in str(region):
            return "Unassigned-Americas"
        else:
            return "Unassigned-Not Americas"

    # Keep A*, B*, C*, D* patterns as-is
    return genotype


def genoflu_postprocess(metadata_df):
    """Add genotype_ml column by refining raw genoflu annotations"""
    # Apply genotype refinement
    metadata_df["genotype_ml"] = metadata_df.apply(genoflu_refine_genotype, axis=1)

    # Print genotype counts for debugging
    counts = metadata_df["genotype_ml"].value_counts()
    print("Refined genotype counts:", counts.to_string())

    return metadata_df


def generate_genotype_colors(metadata_tsv, output_tsv):
    """Generate colors for all genotypes in genotype_ml column"""
    metadata = pd.read_csv(metadata_tsv, sep="\t")
    unique_genotypes = metadata["genotype_ml"].dropna().unique()

    # Categorize genotypes
    major = sorted([g for g in unique_genotypes if g[0] in "ABCD" and g != "Minor"])
    has_minor = "Minor" in unique_genotypes
    unassigned_americas = "Unassigned-Americas" in unique_genotypes
    unassigned_not_americas = "Unassigned-Not Americas" in unique_genotypes

    colors = []

    # Major genotypes: continuous gradient
    if major:
        color_hex_list = [
            "#4042C7",
            "#4274CE",
            "#5199B7",
            "#69B091",
            "#88BB6C",
            "#ADBD51",
            "#CEB541",
            "#E39B39",
            "#E56C2F",
            "#DC2F24",
        ]
        custom_cmap = mcolors.LinearSegmentedColormap.from_list(
            "custom_gradient", color_hex_list
        )
        normalized_values = np.linspace(0, 1, len(major))

        for i, genotype in enumerate(major):
            rgba_color = custom_cmap(normalized_values[i])
            hex_color = mcolors.to_hex(rgba_color)
            colors.append(f"genotype_ml\t{genotype}\t{hex_color}")

    # Minor: very light grey (almost white)
    if has_minor:
        colors.append(f"genotype_ml\tMinor\t#F0F0F0")

    # Unassigned categories
    if unassigned_americas:
        colors.append("genotype_ml\tUnassigned-Americas\t#404040")
    if unassigned_not_americas:
        colors.append("genotype_ml\tUnassigned-Not Americas\t#808080")

    # Write to file
    with open(output_tsv, "w") as f:
        f.write("\n".join(colors) + "\n")


def extract_metadata_post_filter(sequences_fasta, metadata_tsv, output_tsv):
    """Extract metadata for strains that survived filtering and concatenation."""
    strains = set()
    with open(sequences_fasta) as f:
        for line in f:
            if line.startswith(">"):
                strains.add(line[1:].strip())
    meta = pd.read_csv(metadata_tsv, sep="\t", low_memory=False)
    filtered = meta[meta["strain"].isin(strains)]
    filtered.to_csv(output_tsv, sep="\t", index=False)


def extract_metadata_post_refine(tree_nwk, metadata_tsv, output_tsv):
    """Extract metadata for strains that survived clock filtering."""
    tree = Phylo.read(tree_nwk, "newick")
    tips = {tip.name for tip in tree.get_terminals()}
    meta = pd.read_csv(metadata_tsv, sep="\t", low_memory=False)
    filtered = meta[meta["strain"].isin(tips)]
    filtered.to_csv(output_tsv, sep="\t", index=False)
