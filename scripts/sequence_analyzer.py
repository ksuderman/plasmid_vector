"""
Basic sequence analysis for C. acnes promoters
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from utils import analyze_promoter_dataset, load_fasta_sequences


def plot_feature_distributions(df: pd.DataFrame, output_dir: str = "../results"):
    """Plot distributions of promoter features."""
    Path(output_dir).mkdir(exist_ok=True)

    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")

    # Feature distribution plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # GC content
    axes[0, 0].hist(df['gc_content'], bins=20, alpha=0.7)
    axes[0, 0].set_title('GC Content Distribution')
    axes[0, 0].set_xlabel('GC Content')
    axes[0, 0].set_ylabel('Frequency')

    # -10 box scores
    axes[0, 1].hist(df['minus_10_score'], bins=20, alpha=0.7)
    axes[0, 1].set_title('-10 Box Consensus Scores')
    axes[0, 1].set_xlabel('Consensus Score')
    axes[0, 1].set_ylabel('Frequency')

    # -35 box scores
    axes[0, 2].hist(df['minus_35_score'], bins=20, alpha=0.7)
    axes[0, 2].set_title('-35 Box Consensus Scores')
    axes[0, 2].set_xlabel('Consensus Score')
    axes[0, 2].set_ylabel('Frequency')

    # Spacer lengths
    valid_spacers = df[df['spacer_length'] > 0]['spacer_length']
    axes[1, 0].hist(valid_spacers, bins=range(10, 25), alpha=0.7)
    axes[1, 0].set_title('Spacer Length Distribution')
    axes[1, 0].set_xlabel('Spacer Length (bp)')
    axes[1, 0].set_ylabel('Frequency')

    # Sequence lengths
    axes[1, 1].hist(df['length'], bins=20, alpha=0.7)
    axes[1, 1].set_title('Sequence Length Distribution')
    axes[1, 1].set_xlabel('Length (bp)')
    axes[1, 1].set_ylabel('Frequency')

    # Correlation heatmap
    numeric_cols = ['gc_content', 'minus_10_score', 'minus_35_score', 'spacer_length', 'length']
    correlation_data = df[numeric_cols].corr()
    sns.heatmap(correlation_data, annot=True, cmap='coolwarm', center=0, ax=axes[1, 2])
    axes[1, 2].set_title('Feature Correlations')

    plt.tight_layout()
    plt.savefig(f"{output_dir}/promoter_feature_analysis.png", dpi=300, bbox_inches='tight')
    plt.show()


def analyze_consensus_sequences(df: pd.DataFrame):
    """Analyze consensus sequences for -10 and -35 boxes."""
    print("\\n=== Consensus Sequence Analysis ===")

    # -10 box sequences
    minus_10_seqs = df[df['minus_10_sequence'] != '']['minus_10_sequence'].tolist()
    if minus_10_seqs:
        print(f"\\nFound {len(minus_10_seqs)} -10 boxes:")
        print("Top -10 sequences:")
        for seq, score in zip(df[df['minus_10_sequence'] != '']['minus_10_sequence'][:10],
                             df[df['minus_10_sequence'] != '']['minus_10_score'][:10]):
            print(f"  {seq} (score: {score:.3f})")

    # -35 box sequences
    minus_35_seqs = df[df['minus_35_sequence'] != '']['minus_35_sequence'].tolist()
    if minus_35_seqs:
        print(f"\\nFound {len(minus_35_seqs)} -35 boxes:")
        print("Top -35 sequences:")
        for seq, score in zip(df[df['minus_35_sequence'] != '']['minus_35_sequence'][:10],
                             df[df['minus_35_sequence'] != '']['minus_35_score'][:10]):
            print(f"  {seq} (score: {score:.3f})")


def generate_summary_report(df: pd.DataFrame, output_file: str = "../results/analysis_summary.txt"):
    """Generate a summary report of the analysis."""
    Path(output_file).parent.mkdir(exist_ok=True)

    with open(output_file, 'w') as f:
        f.write("C. acnes Promoter Analysis Summary\\n")
        f.write("=" * 40 + "\\n\\n")

        f.write(f"Total sequences analyzed: {len(df)}\\n")
        f.write(f"Average sequence length: {df['length'].mean():.1f} bp\\n")
        f.write(f"Average GC content: {df['gc_content'].mean():.3f}\\n\\n")

        # -10 box statistics
        valid_minus_10 = df[df['minus_10_score'] > 0]
        f.write(f"-10 boxes found: {len(valid_minus_10)} ({len(valid_minus_10)/len(df)*100:.1f}%)\\n")
        if len(valid_minus_10) > 0:
            f.write(f"Average -10 consensus score: {valid_minus_10['minus_10_score'].mean():.3f}\\n")
            f.write(f"Best -10 score: {valid_minus_10['minus_10_score'].max():.3f}\\n")

        # -35 box statistics
        valid_minus_35 = df[df['minus_35_score'] > 0]
        f.write(f"\\n-35 boxes found: {len(valid_minus_35)} ({len(valid_minus_35)/len(df)*100:.1f}%)\\n")
        if len(valid_minus_35) > 0:
            f.write(f"Average -35 consensus score: {valid_minus_35['minus_35_score'].mean():.3f}\\n")
            f.write(f"Best -35 score: {valid_minus_35['minus_35_score'].max():.3f}\\n")

        # Spacer statistics
        valid_spacers = df[df['spacer_length'] > 0]
        if len(valid_spacers) > 0:
            f.write(f"\\nSpacer length statistics:\\n")
            f.write(f"Average spacer length: {valid_spacers['spacer_length'].mean():.1f} bp\\n")
            f.write(f"Most common spacer length: {valid_spacers['spacer_length'].mode().iloc[0]} bp\\n")

    print(f"Summary report saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Analyze promoter sequences")
    parser.add_argument("fasta_file", help="Input FASTA file with promoter sequences")
    parser.add_argument("--output-dir", default="../results", help="Output directory for results")
    parser.add_argument("--no-plots", action="store_true", help="Skip generating plots")

    args = parser.parse_args()

    print(f"Analyzing promoter sequences from: {args.fasta_file}")

    # Analyze sequences
    df = analyze_promoter_dataset(args.fasta_file)

    # Print basic statistics
    print(f"\\nAnalyzed {len(df)} sequences")
    print(f"Average GC content: {df['gc_content'].mean():.3f}")
    print(f"Sequences with -10 boxes: {(df['minus_10_score'] > 0).sum()}")
    print(f"Sequences with -35 boxes: {(df['minus_35_score'] > 0).sum()}")

    # Generate detailed analysis
    analyze_consensus_sequences(df)

    # Save detailed results
    output_file = f"{args.output_dir}/detailed_analysis.csv"
    Path(args.output_dir).mkdir(exist_ok=True)
    df.to_csv(output_file, index=False)
    print(f"\\nDetailed results saved to: {output_file}")

    # Generate plots
    if not args.no_plots:
        plot_feature_distributions(df, args.output_dir)

    # Generate summary report
    generate_summary_report(df, f"{args.output_dir}/analysis_summary.txt")


if __name__ == "__main__":
    main()