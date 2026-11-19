"""Reading protein sequences from FASTA files.

A FASTA file looks like this:

    >sp|P0CG48|UBC_HUMAN Ubiquitin
    MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ
    KESTLHLVLRLRGG

A line starting with '>' is a header (the name). The lines after it are the
sequence, until the next header.
"""

from __future__ import annotations

from protein_analyzer.sequence import ProteinSequence


def parse_fasta(text: str) -> list[ProteinSequence]:
    """Return every protein found in FASTA-formatted text.

    Raises:
        ValueError: If sequence lines appear before any header, or a header
            has no sequence after it.

    >>> proteins = parse_fasta(">pep1\\nGG\\n>pep2\\nAC\\nDE\\n")
    >>> [(p.name, p.residues) for p in proteins]
    [('pep1', 'GG'), ('pep2', 'ACDE')]
    """
    proteins: list[ProteinSequence] = []
    name: str | None = None
    chunks: list[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                proteins.append(_build(name, chunks))
            name = line[1:].strip() or "unnamed"
            chunks = []
        elif name is None:
            raise ValueError("Found sequence before any '>' header line.")
        else:
            chunks.append(line)

    if name is not None:
        proteins.append(_build(name, chunks))
    return proteins


def read_fasta(path: str) -> list[ProteinSequence]:
    """Return every protein in the FASTA file at <path>."""
    with open(path, encoding="utf-8") as file:
        return parse_fasta(file.read())


def _build(name: str, chunks: list[str]) -> ProteinSequence:
    """Join sequence lines into one ProteinSequence."""
    if not chunks:
        raise ValueError(f"Header '{name}' has no sequence.")
    return ProteinSequence("".join(chunks), name=name)
