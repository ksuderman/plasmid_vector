"""
Utility functions for promoter sequence analysis
"""

import re
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


def load_fasta_sequences(fasta_file: str) -> Dict[str, str]:
    """
    Load sequences from a FASTA file.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        Dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


def find_minus_10_box(sequence: str, position: int = None) -> List[Tuple[int, str, float]]:
    """
    Find potential -10 boxes (Pribnow box) in the sequence.

    Args:
        sequence: DNA sequence
        position: Expected position of -10 box (if None, searches whole sequence)

    Returns:
        List of (position, sequence, score) tuples
    """
    consensus = "TATAAT"
    minus_10_candidates = []

    if position:
        # Search around expected position ±5 bp
        start = max(0, position - 5)
        end = min(len(sequence), position + 6)
        search_seq = sequence[start:end]
        search_start = start
    else:
        search_seq = sequence
        search_start = 0

    for i in range(len(search_seq) - 5):
        hexamer = search_seq[i:i+6]
        score = calculate_consensus_score(hexamer, consensus)
        if score >= 0.5:  # Threshold for potential -10 box
            minus_10_candidates.append((search_start + i, hexamer, score))

    return sorted(minus_10_candidates, key=lambda x: x[2], reverse=True)


def find_minus_35_box(sequence: str, minus_10_pos: int = None) -> List[Tuple[int, str, float]]:
    """
    Find potential -35 boxes in the sequence.

    Args:
        sequence: DNA sequence
        minus_10_pos: Position of -10 box (to constrain search to ~17bp upstream)

    Returns:
        List of (position, sequence, score) tuples
    """
    consensus = "TTGACA"
    minus_35_candidates = []

    if minus_10_pos:
        # Search for -35 box 15-20 bp upstream of -10 box
        start = max(0, minus_10_pos - 25)
        end = max(0, minus_10_pos - 10)
        search_seq = sequence[start:end]
        search_start = start
    else:
        search_seq = sequence
        search_start = 0

    for i in range(len(search_seq) - 5):
        hexamer = search_seq[i:i+6]
        score = calculate_consensus_score(hexamer, consensus)
        if score >= 0.5:  # Threshold for potential -35 box
            minus_35_candidates.append((search_start + i, hexamer, score))

    return sorted(minus_35_candidates, key=lambda x: x[2], reverse=True)


def calculate_consensus_score(sequence: str, consensus: str) -> float:
    """
    Calculate similarity score between a sequence and consensus.

    Args:
        sequence: Input sequence
        consensus: Consensus sequence

    Returns:
        Score between 0 and 1 (1 = perfect match)
    """
    if len(sequence) != len(consensus):
        return 0.0

    matches = sum(1 for i, (s, c) in enumerate(zip(sequence, consensus)) if s == c)
    return matches / len(consensus)


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence) if len(sequence) > 0 else 0.0


def extract_promoter_features(sequence: str) -> Dict[str, float]:
    """
    Extract various features from a promoter sequence.

    Args:
        sequence: Promoter sequence

    Returns:
        Dictionary of features
    """
    features = {}

    # Basic composition features
    features['length'] = len(sequence)
    features['gc_content'] = calculate_gc_content(sequence)
    features['at_content'] = 1 - features['gc_content']

    # Find best -10 and -35 boxes
    minus_10_results = find_minus_10_box(sequence)
    minus_35_results = find_minus_35_box(sequence)

    if minus_10_results:
        best_minus_10 = minus_10_results[0]
        features['minus_10_score'] = best_minus_10[2]
        features['minus_10_position'] = best_minus_10[0]
        features['minus_10_sequence'] = best_minus_10[1]

        # Look for -35 box relative to best -10
        minus_35_relative = find_minus_35_box(sequence, best_minus_10[0])
        if minus_35_relative:
            best_minus_35 = minus_35_relative[0]
            features['minus_35_score'] = best_minus_35[2]
            features['minus_35_position'] = best_minus_35[0]
            features['minus_35_sequence'] = best_minus_35[1]

            # Calculate spacer length
            features['spacer_length'] = best_minus_10[0] - best_minus_35[0] - 6
        else:
            features['minus_35_score'] = 0.0
            features['minus_35_position'] = -1
            features['minus_35_sequence'] = ''
            features['spacer_length'] = -1
    else:
        # No -10 box found
        features['minus_10_score'] = 0.0
        features['minus_10_position'] = -1
        features['minus_10_sequence'] = ''
        features['minus_35_score'] = 0.0
        features['minus_35_position'] = -1
        features['minus_35_sequence'] = ''
        features['spacer_length'] = -1

    return features


def analyze_promoter_dataset(fasta_file: str) -> pd.DataFrame:
    """
    Analyze a dataset of promoter sequences.

    Args:
        fasta_file: Path to FASTA file with promoter sequences

    Returns:
        DataFrame with promoter features
    """
    sequences = load_fasta_sequences(fasta_file)

    results = []
    for seq_id, sequence in sequences.items():
        features = extract_promoter_features(sequence)
        features['sequence_id'] = seq_id
        features['sequence'] = sequence
        results.append(features)

    return pd.DataFrame(results)


if __name__ == "__main__":
    # Example usage
    test_sequence = "TTGACAGATCGATCGATCGATCGATCGATATAAT"
    features = extract_promoter_features(test_sequence)
    print("Test sequence features:")
    for key, value in features.items():
        print(f"  {key}: {value}")