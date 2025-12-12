#!/usr/bin/env python3
"""
Analyze all species data sources to understand overlap and gaps before merging.
"""

import pandas as pd

def load_all_sources():
    """Load all data sources"""

    # Anna's data (legacy)
    anna = pd.read_csv('metadata_mod_scripts/species.csv')
    anna = anna.rename(columns={
        'annotated': 'strain_host',
        'correction': 'corrected_host',
        'broad': 'category'
    })
    anna['strain_host_key'] = anna['strain_host'].str.lower().str.strip()
    print(f"Anna (species.csv): {len(anna)} rows")

    # Clara's North America data (authoritative for NAm)
    clara_nam = pd.read_csv('clara/NAm_host_list_order-domesticstat.csv', encoding='latin-1')
    clara_nam['strain_host_key'] = clara_nam['strain_host'].str.lower().str.strip()
    clara_nam = clara_nam.rename(columns={'latinname': 'latin_name', 'host_category': 'category'})
    print(f"Clara NAm: {len(clara_nam)} rows")

    # Clara's global data (fill gaps)
    clara_global = pd.read_csv('clara/global_host_list_order-domesticstat.csv', encoding='latin-1')
    clara_global['strain_host_key'] = clara_global['strain_host'].str.lower().str.strip()
    clara_global = clara_global.rename(columns={'latin_name': 'latin_name', 'host_category': 'category'})
    print(f"Clara Global: {len(clara_global)} rows")

    # Clara's birdspecies (Lambo's work, incorporated by Jordan per Slack)
    clara_birds = pd.read_csv('clara/birdspecies.csv')
    clara_birds = clara_birds.rename(columns={
        'annotated': 'strain_host',
        'correction': 'corrected_host',
        'broad': 'category',
        'Host': 'latin_name'
    })
    clara_birds['strain_host_key'] = clara_birds['strain_host'].str.lower().str.strip()
    print(f"Clara birdspecies: {len(clara_birds)} rows")

    # Jordan's updates (348 additions to Clara's data)
    jordan = pd.read_csv('jordan/updates-jo.csv')
    jordan = jordan.rename(columns={'annotation': 'strain_host'})
    jordan['strain_host_key'] = jordan['strain_host'].str.lower().str.strip()
    print(f"Jordan updates: {len(jordan)} rows")

    return {
        'anna': anna,
        'clara_nam': clara_nam,
        'clara_global': clara_global,
        'clara_birds': clara_birds,
        'jordan': jordan
    }

def analyze_overlaps(sources):
    """Analyze overlaps between sources"""

    anna_keys = set(sources['anna']['strain_host_key'])
    clara_nam_keys = set(sources['clara_nam']['strain_host_key'])
    clara_global_keys = set(sources['clara_global']['strain_host_key'])
    clara_birds_keys = set(sources['clara_birds']['strain_host_key'])
    jordan_keys = set(sources['jordan']['strain_host_key'])

    print("\n=== OVERLAP ANALYSIS ===")
    print(f"Total unique species across all sources: {len(anna_keys | clara_nam_keys | clara_global_keys | clara_birds_keys | jordan_keys)}")

    print(f"\nAnna ONLY (not in any Clara/Jordan): {len(anna_keys - clara_nam_keys - clara_global_keys - clara_birds_keys - jordan_keys)}")
    print(f"Clara NAm ONLY: {len(clara_nam_keys - anna_keys - clara_global_keys - clara_birds_keys - jordan_keys)}")
    print(f"Clara Global ONLY: {len(clara_global_keys - anna_keys - clara_nam_keys - clara_birds_keys - jordan_keys)}")
    print(f"Clara Birds ONLY: {len(clara_birds_keys - anna_keys - clara_nam_keys - clara_global_keys - jordan_keys)}")
    print(f"Jordan ONLY: {len(jordan_keys - anna_keys - clara_nam_keys - clara_global_keys - clara_birds_keys)}")

    # Species in Anna but missing from Clara/Jordan
    anna_exclusive = anna_keys - clara_nam_keys - clara_global_keys - clara_birds_keys - jordan_keys
    if anna_exclusive:
        print(f"\n=== SPECIES IN ANNA BUT NOT IN CLARA/JORDAN ({len(anna_exclusive)}) ===")
        anna_only_df = sources['anna'][sources['anna']['strain_host_key'].isin(anna_exclusive)]
        print(anna_only_df[['strain_host', 'corrected_host', 'category', 'order']].head(20).to_string(index=False))

        # Save full list for review
        anna_only_df[['strain_host', 'corrected_host', 'category', 'order']].to_csv(
            'anna_exclusive_species.csv', index=False
        )
        print(f"\nFull list saved to: anna_exclusive_species.csv")

def analyze_field_coverage(sources):
    """Analyze which sources have which fields populated"""

    print("\n=== FIELD COVERAGE ===")
    for name, df in sources.items():
        print(f"\n{name}:")
        for col in ['corrected_host', 'latin_name', 'category', 'order', 'domesticstatus']:
            if col in df.columns:
                pct = (df[col].notna().sum() / len(df)) * 100
                print(f"  {col}: {pct:.1f}% ({df[col].notna().sum()}/{len(df)})")
            else:
                print(f"  {col}: NOT PRESENT")

def check_conflicts(sources):
    """Check for conflicts between sources on the same species"""

    print("\n=== CHECKING FOR CONFLICTS ===")

    # Compare Clara NAm vs Jordan on shared species
    clara_nam = sources['clara_nam']
    jordan = sources['jordan']

    shared_keys = set(clara_nam['strain_host_key']) & set(jordan['strain_host_key'])
    print(f"\nSpecies in both Clara NAm and Jordan: {len(shared_keys)}")

    conflicts = []
    for key in shared_keys:
        c_row = clara_nam[clara_nam['strain_host_key'] == key].iloc[0]
        j_row = jordan[jordan['strain_host_key'] == key].iloc[0]

        # Check for disagreements
        conflict_fields = []
        for field in ['corrected_host', 'category', 'order', 'domesticstatus']:
            if field in c_row and field in j_row:
                c_val = c_row[field] if pd.notna(c_row[field]) else None
                j_val = j_row[field] if pd.notna(j_row[field]) else None

                if c_val and j_val and str(c_val).lower() != str(j_val).lower():
                    conflict_fields.append(f"{field}:C={c_val}/J={j_val}")

        if conflict_fields:
            conflicts.append({
                'strain_host': c_row['strain_host'],
                'conflicts': ', '.join(conflict_fields)
            })

    if conflicts:
        print(f"\nFound {len(conflicts)} conflicts between Clara NAm and Jordan:")
        for c in conflicts[:10]:
            print(f"  {c['strain_host']}: {c['conflicts']}")
        if len(conflicts) > 10:
            print(f"  ... and {len(conflicts) - 10} more")

def main():
    sources = load_all_sources()
    analyze_overlaps(sources)
    analyze_field_coverage(sources)
    check_conflicts(sources)

if __name__ == '__main__':
    main()
