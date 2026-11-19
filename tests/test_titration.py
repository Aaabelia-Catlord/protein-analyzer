"""Tests for the titration curve feature."""

from pathlib import Path

from protein_analyzer import ProteinSequence
from protein_analyzer.cli import main
from protein_analyzer.plot import plot_titration_curves, titration_points


def test_titration_points_cover_ph_0_to_14() -> None:
    phs, charges = titration_points(ProteinSequence("KDE"))
    assert len(phs) == len(charges) == 141
    assert phs[0] == 0
    assert phs[-1] == 14


def test_charge_goes_from_positive_to_negative() -> None:
    _, charges = titration_points(ProteinSequence("KDE"))
    assert charges[0] > 0
    assert charges[-1] < 0


def test_plot_creates_image(tmp_path: Path) -> None:
    image = tmp_path / "curve.png"
    plot_titration_curves([ProteinSequence("KDE", name="test")], str(image))
    assert image.exists()


def test_cli_plot_option(tmp_path: Path) -> None:
    image = tmp_path / "curve.png"
    assert main(["KKKDDD", "--plot", str(image)]) == 0
    assert image.exists()
