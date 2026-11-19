"""Protein Analyzer: calculate properties of proteins from their sequences."""

from protein_analyzer.fasta import parse_fasta, read_fasta
from protein_analyzer.sequence import ProteinSequence

__all__ = ["ProteinSequence", "parse_fasta", "read_fasta"]
