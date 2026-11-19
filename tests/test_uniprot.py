"""Tests for downloading from UniProt.

These tests never use the internet. Instead, they replace the real
download function with a fake one ("mocking"), so the tests are fast and
still pass when you're offline or UniProt is down.
"""

import io
import urllib.error

import pytest

from protein_analyzer import uniprot
from protein_analyzer.cli import main
from protein_analyzer.uniprot import UniProtError, fetch_by_accession

# A real response from https://rest.uniprot.org/uniprotkb/P01308.fasta
INSULIN_FASTA = (
    ">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1\n"
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN\n"
)


@pytest.fixture
def fake_download(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Replace the real download with one that returns insulin.

    Returns a list that records every URL that was "downloaded".
    """
    requested_urls: list[str] = []

    def fake(url: str) -> str:
        requested_urls.append(url)
        return INSULIN_FASTA

    monkeypatch.setattr(uniprot, "_download", fake)
    return requested_urls


# ---------- fetch_by_accession ----------

def test_fetch_uses_correct_url(fake_download: list[str]) -> None:
    fetch_by_accession("p01308")  # lowercase should work too
    assert fake_download == ["https://rest.uniprot.org/uniprotkb/P01308.fasta"]


def test_fetch_returns_protein(fake_download: list[str]) -> None:
    protein = fetch_by_accession("P01308")
    assert len(protein) == 110
    assert protein.residues.startswith("MALWMRLL")
    assert protein.name == "Insulin (P01308, Homo sapiens)"


def test_invalid_accession_does_not_download(fake_download: list[str]) -> None:
    with pytest.raises(UniProtError, match="not a valid"):
        fetch_by_accession("insulin")
    assert fake_download == []


def test_unsupported_residues(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(uniprot, "_download", lambda url: ">sp|P1|X\nMUA\n")
    with pytest.raises(UniProtError, match="Can't analyze"):
        fetch_by_accession("P12345")


# ---------- Network errors ----------

def _raise(error: Exception):
    """Return a fake urlopen that always raises <error>."""
    def fake_urlopen(*args, **kwargs):
        raise error
    return fake_urlopen


def test_not_found(monkeypatch: pytest.MonkeyPatch) -> None:
    error = urllib.error.HTTPError("url", 404, "Not Found", {}, io.BytesIO())
    monkeypatch.setattr(uniprot.urllib.request, "urlopen", _raise(error))
    with pytest.raises(UniProtError, match="No protein found"):
        fetch_by_accession("P99999")


def test_server_error(monkeypatch: pytest.MonkeyPatch) -> None:
    error = urllib.error.HTTPError("url", 500, "Error", {}, io.BytesIO())
    monkeypatch.setattr(uniprot.urllib.request, "urlopen", _raise(error))
    with pytest.raises(UniProtError, match="HTTP 500"):
        fetch_by_accession("P01308")


def test_no_internet(monkeypatch: pytest.MonkeyPatch) -> None:
    error = urllib.error.URLError("no connection")
    monkeypatch.setattr(uniprot.urllib.request, "urlopen", _raise(error))
    with pytest.raises(UniProtError, match="internet"):
        fetch_by_accession("P01308")


# ---------- Command line ----------

def test_cli_uniprot(fake_download: list[str],
                     capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--uniprot", "P01308"]) == 0
    assert "Insulin (P01308, Homo sapiens)" in capsys.readouterr().out


def test_cli_uniprot_error(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--uniprot", "not-an-id"]) == 1
    assert "not a valid" in capsys.readouterr().err
