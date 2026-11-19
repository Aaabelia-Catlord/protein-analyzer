"""Downloading protein sequences from UniProt (https://www.uniprot.org).

Uses UniProt's REST API. For example, this URL returns human insulin
in FASTA format:
    https://rest.uniprot.org/uniprotkb/P01308.fasta
"""

from __future__ import annotations

import re
import urllib.error
import urllib.parse
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


# Common names for organisms, mapped to their NCBI taxonomy IDs.
# Users can also type any taxonomy ID directly (e.g. 9913 for cow).
ORGANISMS: dict[str, int] = {
    "human": 9606,
    "mouse": 10090,
    "rat": 10116,
    "zebrafish": 7955,
    "fly": 7227,      # Drosophila melanogaster
    "worm": 6239,     # Caenorhabditis elegans
    "yeast": 559292,  # Saccharomyces cerevisiae (strain S288c)
    "ecoli": 83333,   # Escherichia coli (strain K12)
}

MAX_GENE_RESULTS = 5


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


def fetch_by_gene(gene: str, organism: str = "human"
                  ) -> list[ProteinSequence]:
    """Search UniProt for proteins made by <gene> in <organism>.

    Only reviewed (Swiss-Prot) entries are returned, because they are
    checked by curators. Returns at most MAX_GENE_RESULTS proteins, best
    match first.

    Raises:
        UniProtError: If the gene name or organism is invalid, nothing
            matches, or UniProt can't be reached.
    """
    gene = gene.strip()
    # Gene names only use letters, numbers, '-', '_' and '.'.
    if not re.fullmatch(r"[A-Za-z0-9._-]+", gene):
        raise UniProtError(f"'{gene}' is not a valid gene name.")

    taxonomy_id = organism_to_taxonomy_id(organism)
    query = (f"gene_exact:{gene} AND organism_id:{taxonomy_id} "
             "AND reviewed:true")
    params = urllib.parse.urlencode(
        {"query": query, "format": "fasta", "size": MAX_GENE_RESULTS}
    )
    text = _download(f"{BASE_URL}/search?{params}")

    if not text.strip():
        raise UniProtError(
            f"No reviewed UniProt entry found for gene '{gene}' "
            f"in organism '{organism}'."
        )

    proteins = []
    for protein in parse_fasta_safely(text):
        protein.name = short_name(protein.name)
        proteins.append(protein)
    if not proteins:
        raise UniProtError(f"Can't analyze any protein for gene '{gene}'.")
    return proteins


def organism_to_taxonomy_id(organism: str) -> int:
    """Return the NCBI taxonomy ID for a common name or a numeric ID.

    >>> organism_to_taxonomy_id("Mouse")
    10090
    >>> organism_to_taxonomy_id("9913")
    9913
    """
    organism = organism.strip().lower()
    if organism.isdigit():
        return int(organism)
    if organism in ORGANISMS:
        return ORGANISMS[organism]
    raise UniProtError(
        f"Unknown organism '{organism}'. Use one of: "
        f"{', '.join(ORGANISMS)}, or a taxonomy ID number."
    )


def parse_fasta_safely(text: str) -> list[ProteinSequence]:
    """Parse FASTA text, skipping records this tool can't analyze.

    For example, proteins with selenocysteine (U) are skipped instead of
    stopping the whole search.
    """
    proteins = []
    record = ""
    for line in text.splitlines(keepends=True):
        if line.startswith(">") and record:
            proteins.extend(_try_parse(record))
            record = ""
        record += line
    if record:
        proteins.extend(_try_parse(record))
    return proteins


def _try_parse(record: str) -> list[ProteinSequence]:
    """Return the protein in one FASTA record, or [] if it can't be parsed."""
    try:
        return parse_fasta(record)
    except ValueError:
        return []


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
            raise UniProtError(
                "No protein found. Check the ID or gene name.") from error
        raise UniProtError(
            f"UniProt returned an error (HTTP {error.code}). Try again later."
        ) from error
    except (urllib.error.URLError, TimeoutError) as error:
        raise UniProtError(
            "Could not connect to UniProt. Check your internet connection."
        ) from error
