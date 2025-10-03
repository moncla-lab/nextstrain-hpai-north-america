import json

import pandas as pd
from Bio import SeqIO

SEGMENTS = ["pb1", "pb2", "na", "pa", "ha", "np", "mp", "ns"]


def extract_animal(item):
    parts = item.split("/")
    if len(parts) > 1:
        return parts[1].lower()
    return None


def add_orders(df):
    carnivore_list = [
        "skunk",
        "redfox",
        "fox",
        "bobcat",
        "harborseal",
        "raccoon",
        "blackbear",
        "stripedskunk",
        "cat",
        "vulpesvulpes",
        "coyote",
        "greyseal",
        "wildmink",
    ]
    marsup_list = ["virginiaopossum"]
    artiodactyl_list = ["bottlenosedolphin", "dolphin", "dairycattle", "goat"]
    anseriformes_list = ["lesserscaup"]
    passeriformes_list = ["greattailedgrackle", "americanraven", "commongrackle"]
    pelican_list = ["brownpelican", "snowyegret"]
    accipitriformes_list = ["osprey", "turkeyvulture", "coopershawk"]

    df.loc[df["Animal"].isin(carnivore_list), "order"] = "carnivora"
    df.loc[df["Animal"].isin(marsup_list), "order"] = "didelphimorphia"
    df.loc[df["Animal"].isin(artiodactyl_list), "order"] = "artiodactyl"
    df.loc[df["Animal"].isin(anseriformes_list), "order"] = "anseriformes"
    df.loc[df["Animal"].isin(passeriformes_list), "order"] = "passeriformes"
    df.loc[df["Animal"].isin(pelican_list), "order"] = "pelecaniformes"
    df.loc[df["Animal"].isin(accipitriformes_list), "order"] = "accipitriformes"


def add_species_group(df):
    df["species_group"] = "unknown"

    wild_terrest_list = [
        "skunk",
        "redfox",
        "fox",
        "bobcat",
        "raccoon",
        "blackbear",
        "stripedskunk",
        "cat",
        "domesticcat",
        "feline",
        "vulpesvulpes",
        "coyote",
        "wildmink",
    ]
    rum_list = ["dairycattle", "goat"]
    marine_list = ["harborseal", "greyseal", "bottlenosedolphin", "dolphin"]
    humans = ["Human"]

    anser_list = ["anseriformes"]
    gall_list = ["galliformes"]
    raptor_list = ["accipitriformes", "falconiformes", "strigiformes"]
    waterbird_list = [
        "charadriiformes",
        "pelecaniformes",
        "suliformes",
        "podicipediformes",
    ]
    passer_list = ["passeriformes"]
    other_avian_list = ["casuariiformes", "rheiformes", "avian"]

    df.loc[df["Animal"].isin(wild_terrest_list), "species_group"] = (
        "Mammal- Terrestrial"
    )
    df.loc[df["Animal"].isin(marine_list), "species_group"] = "Mammal- Marine"
    df.loc[df["Animal"].isin(rum_list), "species_group"] = "Ruminant"

    df.loc[df["host"].isin(humans), "species_group"] = "Human"

    df.loc[df["order"].isin(anser_list), "species_group"] = "Anseriformes"
    df.loc[df["order"].isin(gall_list), "species_group"] = "Galliformes"
    df.loc[df["order"].isin(raptor_list), "species_group"] = "Raptor"
    df.loc[df["order"].isin(waterbird_list), "species_group"] = "Other- Waterbird"
    df.loc[df["order"].isin(passer_list), "species_group"] = "Passerine"
    df.loc[df["order"].isin(other_avian_list), "species_group"] = "Other- Avian"


def metadata_annotation(input_metadata_tsv, output_metadata_tsv):
    metadata = pd.read_csv(input_metadata_tsv, sep="\t")

    metadata["Animal"] = metadata.strain.apply(extract_animal)
    species = pd.read_csv("metadata_mod_scripts/species.csv").rename(
        columns={"annotated": "Animal"}
    )
    metadata = pd.merge(metadata, species, how="left", on=["Animal"])

    flyways = pd.read_csv("metadata_mod_scripts/flyway_regions.csv")
    metadata = pd.merge(metadata, flyways, how="left", on=["location"])

    add_orders(metadata)

    add_species_group(metadata)

    # Add refined genotype column
    metadata = genoflu_postprocess(metadata)

    metadata.to_csv(output_metadata_tsv, index=False, sep="\t")


def auspice_segment_config(input_config_path, output_config_path, segment):
    with open(input_config_path) as input_file:
        config = json.load(input_file)

    genoflu_segment_config = {
        "key": f"genoflu_{segment}_lineage",
        "title": f"GenoFlu {segment.upper()} lineage",
        "type": "categorical",
    }
    config["colorings"] += [genoflu_segment_config]
    with open(output_config_path, "w") as output_file:
        json.dump(config, output_file, indent=2)


def genoflu_refine_genotype(row):
    """Refine genotype: keep A*/B*/C*/D*/Minor* patterns, map unassigned by region"""
    genotype = row["genoflu"]

    # Check if it's empty or unassigned
    if pd.isna(genotype) or genotype == "" or "Not assigned" in genotype or "Unseen constellation" in genotype:
        # Check region to determine Americas vs not
        region = row.get("region", "")
        if "America" in region:
            return "Unassigned-Americas"
        else:
            return "Unassigned-Not Americas"

    # Keep A*, B*, C*, D*, Minor* patterns as-is
    return genotype


def genoflu_postprocess(metadata_df):
    """Add genotype_ml column by refining raw genoflu annotations"""
    # Apply genotype refinement
    metadata_df["genotype_ml"] = metadata_df.apply(genoflu_refine_genotype, axis=1)

    # Print genotype counts for debugging
    counts = metadata_df["genotype_ml"].value_counts()
    print("Refined genotype counts:", counts.to_string())

    return metadata_df
