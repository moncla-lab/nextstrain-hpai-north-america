import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import sys


def extract_and_categorize_genotypes(metadata_file, number_of_genotypes=9, included_genotypes=[]):
    print(f"Reading metadata from {metadata_file}...")
    metadata_df = pd.read_csv(metadata_file, sep='\t')

    if 'Genotype' not in metadata_df.columns:
        print("Error: 'Genotype' column not found in metadata")
        sys.exit(1)

    # Count genotypes
    genotype_counts = metadata_df['Genotype'].value_counts()
    print(f"\nTotal unique genotypes found: {len(genotype_counts)}")
    print("\nTop 15 genotypes by count:")
    print(genotype_counts.head(15).to_string())

    # Define special values to exclude
    special_values = {
        "Not Americas",
        "Missing segments",
        "Unassigned",
        "Unassigned-US",
        "Not dominant genotype"
    }

    # Get top N genotypes (excluding special values)
    valid_genotypes = [g for g in genotype_counts.index if g not in special_values]
    top_genotypes = valid_genotypes[:number_of_genotypes]

    # Combine with included genotypes
    desired_genotypes = set(top_genotypes) | set(included_genotypes)

    # Categorize all genotypes
    all_genotypes = set(metadata_df['Genotype'].unique())
    dominant = []
    nondominant = []

    for genotype in all_genotypes:
        if genotype in special_values:
            continue
        elif genotype in desired_genotypes:
            dominant.append(genotype)
        else:
            nondominant.append(genotype)

    # Sort alphabetically
    dominant.sort()
    nondominant.sort()

    return dominant, nondominant, genotype_counts

def generate_colors(dominant, nondominant, output_file='proof_of_concept_colors.tsv'):
    colors = []

    print(f"\n=== Color Assignment ===")
    print(f"Dominant genotypes ({len(dominant)}): {dominant[:5]}{'...' if len(dominant) > 5 else ''}")
    print(f"Non-dominant genotypes ({len(nondominant)}): {nondominant[:5]}{'...' if len(nondominant) > 5 else ''}")

    # Rainbow colors for dominant genotypes
    if dominant:
        rainbow_cmap = cm.get_cmap('hsv')
        # Generate values from 0 to 5/6 (300°) to avoid purple-red overlap
        hues = np.linspace(0, 5/6, len(dominant))

        print("\nRainbow colors for dominant genotypes:")
        for i, genotype in enumerate(dominant):
            rgb_color = rainbow_cmap(hues[i])
            hex_color = mcolors.to_hex(rgb_color)
            colors.append(['genotype', genotype, hex_color])
            if i < 5:  # Show first 5 as examples
                print(f"  {genotype:20s} -> {hex_color}")

    # Grayscale colors for non-dominant genotypes
    if nondominant:
        gray_cmap = cm.get_cmap('gray')
        # Gray range from 0.8 (light) to 0.3 (dark)
        gray_start, gray_end = 0.8, 0.3
        gray_values = np.linspace(gray_start, gray_end, len(nondominant))

        print("\nGrayscale colors for non-dominant genotypes:")
        for i, genotype in enumerate(nondominant):
            rgb_color = gray_cmap(gray_values[i])
            hex_color = mcolors.to_hex(rgb_color)
            colors.append(['genotype', genotype, hex_color])
            if i < 5:  # Show first 5 as examples
                print(f"  {genotype:20s} -> {hex_color}")

    # Write to TSV file
    with open(output_file, 'w') as f:
        for row in colors:
            f.write('\t'.join(row) + '\n')

    print(f"\nColors written to {output_file}")
    print(f"Total entries: {len(colors)}")

    return colors

def main():
    """
    Main execution - proof of concept for dynamic color generation.
    """
    print("=" * 60)
    print("PROOF OF CONCEPT: Dynamic Genotype Color Generation")
    print("=" * 60)

    # Configuration (matching the Snakefile defaults)
    NUMBER_OF_GENOTYPES = 9
    GENOTYPES_TO_INCLUDE = []  # Could add specific genotypes here if needed

    # Check if metadata file exists
    metadata_file = 'results/metadata-with-genoflu.tsv'

    # Extract and categorize genotypes
    dominant, nondominant, counts = extract_and_categorize_genotypes(
        metadata_file,
        NUMBER_OF_GENOTYPES,
        GENOTYPES_TO_INCLUDE
    )

    # Generate colors
    colors = generate_colors(dominant, nondominant)

    print("\n" + "=" * 60)
    print("SUCCESS: Proof of concept complete!")
    print("=" * 60)

    # Show some statistics
    print(f"\nStatistics:")
    print(f"  - Dominant genotypes: {len(dominant)}")
    print(f"  - Non-dominant genotypes: {len(nondominant)}")
    print(f"  - Total colored genotypes: {len(dominant) + len(nondominant)}")

    # Validation: Check that we're covering most of the data
    total_samples_with_genotypes = counts[~counts.index.isin({
        "Not Americas", "Missing segments", "Unassigned", "Unassigned-US"
    })].sum()
    samples_in_dominant = counts[counts.index.isin(dominant)].sum()
    coverage = (samples_in_dominant / total_samples_with_genotypes) * 100

    print(f"\nCoverage analysis:")
    print(f"  - Total samples with valid genotypes: {total_samples_with_genotypes}")
    print(f"  - Samples with dominant genotypes: {samples_in_dominant}")
    print(f"  - Coverage by dominant genotypes: {coverage:.1f}%")

if __name__ == "__main__":
    main()