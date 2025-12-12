#!/usr/bin/env python3
"""
Standalone species metadata annotation script.

Annotates viral metadata with species taxonomy, flyway regions, and domesticity status.
Can be used as a submodule in other Nextstrain builds or run directly from this repository.

Usage:
    python scripts/annotate_species_metadata.py \\
        --input metadata.tsv \\
        --output metadata_annotated.tsv \\
        --config-dir config \\
        --flyway-regions metadata_mod_scripts/flyway_regions.csv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Cache for config files to avoid repeated disk reads
_CONFIG_CACHE = {}


def load_config(filename):
    """
    Load TSV config file with caching.

    Args:
        filename: Path to TSV file

    Returns:
        DataFrame with config data
    """
    if filename not in _CONFIG_CACHE:
        _CONFIG_CACHE[filename] = pd.read_csv(filename, sep='\t')
    return _CONFIG_CACHE[filename]


def extract_animal(strain):
    """
    Extract animal name from strain identifier.

    Args:
        strain: Strain name (e.g., "A/chicken/Oregon/123/2020")

    Returns:
        Animal name (lowercase) or None
    """
    parts = strain.split("/")
    if len(parts) > 1:
        return parts[1].lower()
    return None


def get_gisaid_domestic_status(value):
    """
    Extract valid domestic_status from GISAID metadata field.

    GISAID provides domestic_status for ~8% of samples.
    Returns the value if it's 'wild' or 'domestic', otherwise None.

    Args:
        value: domestic_status value (or pandas Series/scalar)

    Returns:
        'domestic', 'wild', or None
    """
    # Handle both Series and scalar values
    if isinstance(value, pd.Series):
        domestic_status = value.get('domestic_status', None)
    else:
        domestic_status = value

    if pd.notna(domestic_status):
        status = str(domestic_status).lower().strip()
        if status in ['wild', 'domestic']:
            return status

    return None


def add_orders(df, order_overrides_path):
    """
    Fill missing orders using override lookup table.

    Args:
        df: DataFrame with 'Animal' and 'order' columns
        order_overrides_path: Path to order_overrides.tsv
    """
    order_overrides = load_config(order_overrides_path)

    # Create mapping dict for faster lookup
    animal_to_order = dict(zip(order_overrides['animal_name'], order_overrides['order']))

    # Apply overrides only where order is missing
    for animal_name, order in animal_to_order.items():
        mask = (df['Animal'] == animal_name) & (df['order'].isna())
        df.loc[mask, 'order'] = order


def add_species_group(df, species_groups_path, order_groups_path):
    """
    Assign species groups based on animal name and taxonomic order.

    Priority: Animal name lookup > Order lookup

    Args:
        df: DataFrame with 'Animal', 'order', and 'host' columns
        species_groups_path: Path to species_groups.tsv
        order_groups_path: Path to order_groups.tsv
    """
    df["species_group"] = "unknown"

    # Load config files
    species_groups = load_config(species_groups_path)
    order_groups = load_config(order_groups_path)

    # Apply order-based groups first (lower priority)
    order_to_group = dict(zip(order_groups['order'], order_groups['species_group']))
    for order, group in order_to_group.items():
        df.loc[df['order'] == order, 'species_group'] = group

    # Apply animal name-based groups second (higher priority, will override)
    animal_to_group = dict(zip(species_groups['animal_name'], species_groups['species_group']))
    for animal_name, group in animal_to_group.items():
        df.loc[df['Animal'] == animal_name, 'species_group'] = group

    # Special case: Human via host column (highest priority)
    if 'host' in df.columns:
        df.loc[df['host'] == 'Human', 'species_group'] = 'Human'


def annotate_metadata(input_path, output_path, config_dir, flyway_path=None):
    """
    Annotate metadata with species information, flyways, and domesticity.

    Args:
        input_path: Path to input metadata TSV
        output_path: Path to output metadata TSV
        config_dir: Directory containing config TSV files
        flyway_path: Path to flyway_regions.csv (optional)
    """
    config_dir = Path(config_dir)

    # Load metadata
    metadata = pd.read_csv(input_path, sep="\t")

    # Extract animal name from strain
    metadata["Animal"] = metadata.strain.apply(extract_animal)

    # Load unified species lookup table
    species_lookup_path = config_dir / 'species_lookup.tsv'
    species_lookup = load_config(str(species_lookup_path)).rename(
        columns={'strain_host': 'Animal'}
    )

    # Merge species data (corrected_host, latin_name, category, order, domesticstatus)
    metadata = pd.merge(metadata, species_lookup, how="left", on=["Animal"])

    # Merge flyway data (optional)
    if flyway_path:
        flyways = pd.read_csv(flyway_path)
        metadata = pd.merge(metadata, flyways, how="left", on=["location"])

    # Fill missing orders using overrides
    order_overrides_path = config_dir / 'order_overrides.tsv'
    add_orders(metadata, str(order_overrides_path))

    # Assign species groups
    species_groups_path = config_dir / 'species_groups.tsv'
    order_groups_path = config_dir / 'order_groups.tsv'
    add_species_group(metadata, str(species_groups_path), str(order_groups_path))

    # Merge domesticity from two sources with intelligent priority:
    # 1. GISAID domestic_status (if wild/domestic - specific per-sample annotation)
    # 2. Clara's species_lookup domesticstatus (if not "unknown" - expert species-level)
    # 3. Clara's species_lookup domesticstatus (even if "unknown" - ambiguous species)
    # This ensures GISAID's specific annotations override Clara's generic "unknown" assignments

    metadata['gisaid_domestic'] = metadata['domestic_status'].apply(get_gisaid_domestic_status)

    # Use GISAID if it exists and is definitive (wild/domestic)
    # Otherwise use Clara if she has a definitive answer (not "unknown")
    # Otherwise fall back to whatever we have
    metadata['domesticstatus'] = metadata.apply(
        lambda row: row['gisaid_domestic'] if pd.notna(row['gisaid_domestic'])
                    else (row['domesticstatus'] if pd.notna(row['domesticstatus']) and row['domesticstatus'] != 'unknown'
                          else row['domesticstatus']),
        axis=1
    )

    metadata = metadata.drop('gisaid_domestic', axis=1)

    # Write output
    metadata.to_csv(output_path, index=False, sep="\t")


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Annotate viral metadata with species taxonomy and domesticity",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input metadata TSV file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output metadata TSV file'
    )
    parser.add_argument(
        '--config-dir',
        default='config',
        help='Directory containing config TSV files (default: config)'
    )
    parser.add_argument(
        '--flyway-regions',
        help='Path to flyway_regions.csv (optional)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.config_dir).is_dir():
        print(f"Error: Config directory not found: {args.config_dir}", file=sys.stderr)
        sys.exit(1)

    if args.flyway_regions and not Path(args.flyway_regions).exists():
        print(f"Error: Flyway regions file not found: {args.flyway_regions}", file=sys.stderr)
        sys.exit(1)

    # Run annotation
    annotate_metadata(
        input_path=args.input,
        output_path=args.output,
        config_dir=args.config_dir,
        flyway_path=args.flyway_regions
    )

    print(f"Annotated metadata written to {args.output}")


if __name__ == '__main__':
    main()
