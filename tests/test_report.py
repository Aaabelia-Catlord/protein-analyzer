"""Tests for the Markdown report."""

from pathlib import Path

import pytest

from protein_analyzer import ProteinSequence
from protein_analyzer.cli import main
from protein_analyzer.report import build_report, format_sequence, write_report


def make_protein() -> ProteinSequence:
    protein = ProteinSequence("MKWVTFISLLLLFSSAYS", name="Test protein")
    protein.source = "https://www.uniprot.org/uniprotkb/P02768"
    return protein


def test_report_contains_main_properties() -> None:
    text = build_report([make_protein()], "curve.png", ph=7.4)
    assert text.startswith("# Protein report: Test protein")
    assert "| Length | 18 residues |" in text
    assert "Molecular weight" in text
    assert "Isoelectric point (pI)" in text
    assert "Net charge at pH 7.4" in text
    assert "GRAVY" in text


def test_report_links_image_and_source() -> None:
    text = build_report([make_protein()], "curve.png")
    assert "![Titration curve](curve.png)" in text
    assert "Source: https://www.uniprot.org/uniprotkb/P02768" in text


def test_report_composition_uses_full_names() -> None:
    text = build_report([make_protein()], "curve.png")
    assert "| Leucine (L) | 4 |" in text


def test_report_without_source() -> None:
    text = build_report([ProteinSequence("GG")], "curve.png")
    assert "Source:" not in text


def test_report_with_several_proteins() -> None:
    proteins = [ProteinSequence("GG", name="one"),
                ProteinSequence("KK", name="two")]
    text = build_report(proteins, "curve.png")
    assert text.startswith("# Protein report: 2 proteins")
    assert "## one" in text
    assert "## two" in text


def test_format_sequence_line_numbers() -> None:
    lines = format_sequence(ProteinSequence("A" * 130)).splitlines()
    assert len(lines) == 3
    assert lines[1].startswith("  61 ")
    assert lines[2] == " 121 AAAAAAAAAA"


def test_write_report_creates_both_files(tmp_path: Path) -> None:
    report = tmp_path / "insulin.md"
    image = write_report([make_protein()], str(report))
    assert report.exists()
    assert image == tmp_path / "insulin_titration.png"
    assert image.exists()
    assert "insulin_titration.png" in report.read_text(encoding="utf-8")


def test_cli_report(tmp_path: Path) -> None:
    report = tmp_path / "report.md"
    assert main(["KKKDDD", "--report", str(report)]) == 0
    assert report.exists()


def test_cli_report_bad_folder(tmp_path: Path,
                               capsys: pytest.CaptureFixture[str]) -> None:
    report = tmp_path / "missing_folder" / "report.md"
    assert main(["KKKDDD", "--report", str(report)]) == 1
    assert "couldn't save report" in capsys.readouterr().err
