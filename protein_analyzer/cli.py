"""Command-line interface for the protein analyzer.

Examples:
    python -m protein_analyzer MQIFVKTLTGK
    python -m protein_analyzer --fasta examples/ubiquitin.fasta --ph 7.4
"""

from __future__ import annotations

import argparse
import sys

from protein_analyzer.fasta import read_fasta
from protein_analyzer.sequence import ProteinSequence


def build_parser() -> argparse.ArgumentParser:
    """Return the parser that reads command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="protein_analyzer",
        description="Calculate properties of a protein from its sequence.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "sequence", nargs="?",
        help="amino acid sequence in one-letter codes (e.g. MQIFVK)",
    )
    source.add_argument(
        "--fasta", metavar="FILE", help="path to a FASTA file",
    )
    parser.add_argument(
        "--ph", type=float, default=7.0,
        help="pH for the net-charge calculation (default: 7.0)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the program. Return 0 on success and 1 on an error."""
    args = build_parser().parse_args(argv)

    if not 0 <= args.ph <= 14:
        print("Error: pH must be between 0 and 14.", file=sys.stderr)
        return 1

    try:
        if args.fasta:
            proteins = read_fasta(args.fasta)
        else:
            proteins = [ProteinSequence(args.sequence)]
    except (ValueError, OSError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    print("\n\n".join(protein.summary(args.ph) for protein in proteins))
    return 0
