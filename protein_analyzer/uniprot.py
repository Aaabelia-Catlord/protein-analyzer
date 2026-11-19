"""Downloading protein sequences from UniProt (https://www.uniprot.org).

Uses UniProt's REST API. For example, this URL returns human insulin
in FASTA format:
    https://rest.uniprot.org/uniprotkb/P01308.fasta
"""

from __future__ import annotations

import re
import urllib.error
import urllib.request

from protein_analyzer.fasta import parse_fasta
from protein_analyzer.sequence import ProteinSequence

BASE_URL = "https://rest.uniprot.org/uniprotkb"
TIMEOUT_SECONDS = 15

# The official format of a UniProt accession (ID), from UniProt's docs.
# Examples: P01308, Q9Y261, A0A023GPI8
ACCESSION_PATTERN = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


class UniProtError(Exception):
    """Raised when a protein can't be downloaded from UniProt."""


def fetch_by_accession(accession: str) -> ProteinSequence:
    """Download the protein with this UniProt accession (e.g. 'P01308').

    Raises:
        UniProtError: If the accession is badly formatted, doesn't exist,
            UniProt can't be reached, or the sequence can't be analyzed.
    """
    accession = accession.strip().upper()
    if not ACCESSION_PATTERN.match(accession):
        raise UniProtError(
            f"'{accession}' is not a valid UniProt accession "
            "(example: P01308)."
        )

    text = _download(f"{BASE_URL}/{accession}.fasta")
    if not text.strip():
        raise UniProtError(f"No protein found for accession '{accession}'.")

    try:
        protein = parse_fasta(text)[0]
    except ValueError as error:
        # e.g. the sequence has selenocysteine (U), which we don't support.
        raise UniProtError(f"Can't analyze {accession}: {error}") from error

    protein.name = short_name(protein.name)
    return protein


def short_name(header: str) -> str:
    """Turn a long UniProt FASTA header into a readable name.

    >>> short_name("sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS")
    'Insulin (P01308, Homo sapiens)'
    """
    ids, _, rest = header.partition(" ")
    parts = ids.split("|")
    accession = parts[1] if len(parts) >= 2 else ids

    description = rest.split(" OS=")[0].strip() or accession
    organism = re.search(r"OS=(.+?)(?: [A-Z]{2}=|$)", rest)
    if organism:
        return f"{description} ({accession}, {organism.group(1)})"
    return f"{description} ({accession})"


def _download(url: str) -> str:
    """Return the text at <url>, turning network problems into UniProtError."""
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT_SECONDS) as response:
            return response.read().decode("utf-8")
    except urllib.error.HTTPError as error:
        if error.code in (400, 404):
            raise UniProtError("No protein found with that accession.") from error
        raise UniProtError(
            f"UniProt returned an error (HTTP {error.code}). Try again later."
        ) from error
    except (urllib.error.URLError, TimeoutError) as error:
        raise UniProtError(
            "Could not connect to UniProt. Check your internet connection."
        ) from error
