#!/usr/bin/env python3
"""
Merge species metadata from multiple sources into config/species_lookup.tsv

Based on MERGE_STRATEGY.md:
Priority: Clara NAm > Clara Birds > Jordan > Clara Global > Anna (skip)
"""

import pandas as pd
import sys

def load_clara_nam():
    """Load Clara's North America data (highest priority)"""
    df = pd.read_csv('clara/NAm_host_list_order-domesticstat.csv', encoding='latin-1')
    df = df.rename(columns={'latinname': 'latin_name', 'host_category': 'category'})
    df['source'] = 'clara_nam'
    return df[['strain_host', 'corrected_host', 'latin_name', 'category', 'order', 'domesticstatus', 'source']]

def load_clara_birds():
    """Load Clara's birdspecies (Lambo's work)"""
    df = pd.read_csv('clara/birdspecies.csv')
    df = df.rename(columns={
        'annotated': 'strain_host',
        'correction': 'corrected_host',
        'broad': 'category',
        'Host': 'latin_name'
    })
    df['source'] = 'clara_birds'
    # Note: No domesticstatus column in this file
    return df[['strain_host', 'corrected_host', 'latin_name', 'category', 'order', 'source']]

def load_jordan():
    """Load Jordan's updates (new species additions)"""
    df = pd.read_csv('jordan/updates-jo.csv')
    df = df.rename(columns={'annotation': 'strain_host'})
    df['source'] = 'jordan'
    return df[['strain_host', 'corrected_host', 'latin_name', 'category', 'order', 'domesticstatus', 'source']]

def load_clara_global():
    """Load Clara's global data (gap filler)"""
    df = pd.read_csv('clara/global_host_list_order-domesticstat.csv', encoding='latin-1')
    # Column already named latin_name in this file
    df['source'] = 'clara_global'
    return df[['strain_host', 'corrected_host', 'latin_name', 'category', 'order', 'domesticstatus', 'source']]

def merge_with_priority(sources_list):
    """
    Merge multiple sources with priority order.

    For each strain_host (case-insensitive), use first non-null value
    from highest priority source for each field.

    Args:
        sources_list: List of (priority, name, dataframe) tuples, sorted by priority

    Returns:
        Merged dataframe
    """

    all_data = []
    for priority, name, df in sources_list:
        df_copy = df.copy()
        df_copy['priority'] = priority
        df_copy['source_name'] = name
        all_data.append(df_copy)

    # Combine all sources
    combined = pd.concat(all_data, ignore_index=True)

    # Normalize strain_host for deduplication
    combined['strain_host_key'] = combined['strain_host'].str.lower().str.strip()

    print(f"Total rows before deduplication: {len(combined)}")
    print(f"Unique strain_host keys: {combined['strain_host_key'].nunique()}")

    # Sort by priority (higher = better)
    combined = combined.sort_values('priority', ascending=False)

    # Group by strain_host_key and take first non-null for each field
    def first_valid(series):
        """Return first non-null value, or None if all null"""
        valid = series.dropna()
        return valid.iloc[0] if len(valid) > 0 else None

    merged = combined.groupby('strain_host_key', as_index=False).agg({
        'strain_host': 'first',  # Use the form from highest priority source
        'corrected_host': first_valid,
        'latin_name': first_valid,
        'category': first_valid,
        'order': first_valid,
        'domesticstatus': first_valid,
        'source_name': lambda x: ','.join(x.unique())  # Track all contributors
    })

    # Rename source_name back to source
    merged = merged.rename(columns={'source_name': 'source'})

    # Drop temporary key
    merged = merged.drop('strain_host_key', axis=1)

    # Sort by strain_host for readability
    merged = merged.sort_values('strain_host')

    return merged

def print_statistics(df):
    """Print summary statistics about the merged data"""

    print(f"\n=== MERGE RESULTS ===")
    print(f"Total species: {len(df)}")

    print(f"\n=== FIELD COVERAGE ===")
    for col in ['corrected_host', 'latin_name', 'category', 'order', 'domesticstatus']:
        count = df[col].notna().sum()
        pct = (count / len(df)) * 100
        print(f"  {col}: {pct:.1f}% ({count}/{len(df)})")

    print(f"\n=== SOURCE DISTRIBUTION (Top 20) ===")
    source_counts = df['source'].value_counts()
    for source, count in source_counts.head(20).items():
        print(f"  {source}: {count}")

    print(f"\n=== DOMESTICSTATUS BREAKDOWN ===")
    status_counts = df['domesticstatus'].value_counts(dropna=False)
    for status, count in status_counts.items():
        print(f"  {status}: {count}")

    print(f"\n=== SPECIES WITH INCOMPLETE DATA ===")
    incomplete = df[
        df['corrected_host'].isna() |
        df['category'].isna() |
        df['order'].isna()
    ]
    print(f"  Missing critical fields: {len(incomplete)}")
    if len(incomplete) > 0:
        print("\nFirst 10 incomplete species:")
        print(incomplete[['strain_host', 'corrected_host', 'category', 'order', 'source']].head(10).to_string(index=False))

def main():
    print("Loading data sources...")

    clara_nam = load_clara_nam()
    print(f"  Clara NAm: {len(clara_nam)} rows")

    clara_birds = load_clara_birds()
    print(f"  Clara Birds: {len(clara_birds)} rows")

    jordan = load_jordan()
    print(f"  Jordan: {len(jordan)} rows")

    clara_global = load_clara_global()
    print(f"  Clara Global: {len(clara_global)} rows")

    # Define priority order (higher number = higher priority)
    sources = [
        (4, 'clara_nam', clara_nam),
        (3, 'clara_birds', clara_birds),
        (2, 'jordan', jordan),
        (1, 'clara_global', clara_global),
    ]

    print("\nMerging with priority: clara_nam > clara_birds > jordan > clara_global")
    merged = merge_with_priority(sources)

    # Reorder columns for output
    column_order = ['strain_host', 'corrected_host', 'latin_name', 'category', 'order', 'domesticstatus', 'source']
    merged = merged[column_order]

    # Print statistics
    print_statistics(merged)

    # Write output
    output_path = 'config/species_lookup.tsv'
    merged.to_csv(output_path, sep='\t', index=False)
    print(f"\n✓ Wrote {len(merged)} species to {output_path}")

if __name__ == '__main__':
    main()
