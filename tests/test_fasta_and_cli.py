"""Tests for FASTA reading and the command-line interface."""

from pathlib import Path

import pytest

from protein_analyzer import parse_fasta, read_fasta
from protein_analyzer.cli import main

EXAMPLE = Path(__file__).parent.parent / "examples" / "ubiquitin.fasta"


def test_parse_multiple_records() -> None:
    proteins = parse_fasta(">one\nGG\n\n>two\nAC\nDE\n")
    assert [p.name for p in proteins] == ["one", "two"]
    assert [p.residues for p in proteins] == ["GG", "ACDE"]


def test_sequence_before_header_raises() -> None:
    with pytest.raises(ValueError):
        parse_fasta("GGG\n>name\nAAA")


def test_header_without_sequence_raises() -> None:
    with pytest.raises(ValueError):
        parse_fasta(">empty\n>next\nAAA")


def test_read_example_file() -> None:
    proteins = read_fasta(str(EXAMPLE))
    assert len(proteins) == 1
    assert len(proteins[0]) == 76


def test_cli_with_sequence(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["GG"]) == 0
    assert "Length: 2 residues" in capsys.readouterr().out


def test_cli_with_fasta(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--fasta", str(EXAMPLE), "--ph", "7.4"]) == 0
    output = capsys.readouterr().out
    assert "Human ubiquitin" in output
    assert "Net charge at pH 7.4" in output


def test_cli_bad_sequence(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["HELLOWORLD"]) == 1
    assert "Invalid" in capsys.readouterr().err


def test_cli_bad_ph(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["GG", "--ph", "20"]) == 1
    assert "pH" in capsys.readouterr().err
