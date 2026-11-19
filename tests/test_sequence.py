"""Tests for ProteinSequence. Run with: python -m pytest"""

import pytest

from protein_analyzer import ProteinSequence

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
    "KESTLHLVLRLRGG"
)


# ---------- Creating a sequence ----------

def test_cleans_whitespace_and_case() -> None:
    protein = ProteinSequence(" mq if\nvk ")
    assert protein.residues == "MQIFVK"
    assert len(protein) == 6


def test_empty_sequence_raises() -> None:
    with pytest.raises(ValueError):
        ProteinSequence("   ")


def test_invalid_letter_raises() -> None:
    with pytest.raises(ValueError, match="B, X"):
        ProteinSequence("ACXB")


# ---------- Composition ----------

def test_composition_counts() -> None:
    assert ProteinSequence("GAGG").composition() == {"A": 1, "G": 3}


def test_percent_composition_adds_to_100() -> None:
    percents = ProteinSequence(UBIQUITIN).percent_composition()
    assert sum(percents.values()) == pytest.approx(100)


# ---------- Molecular weight ----------

def test_dipeptide_weight() -> None:
    # Two glycine residues + one water.
    assert ProteinSequence("GG").molecular_weight() == pytest.approx(
        132.12, abs=0.01)


def test_ubiquitin_weight_matches_known_value() -> None:
    # Human ubiquitin is about 8564.8 Da (ExPASy ProtParam).
    assert ProteinSequence(UBIQUITIN).molecular_weight() == pytest.approx(
        8564.8, abs=1)


# ---------- Charge and pI ----------

def test_charge_is_positive_at_low_ph_and_negative_at_high_ph() -> None:
    protein = ProteinSequence(UBIQUITIN)
    assert protein.net_charge(1.0) > 0
    assert protein.net_charge(13.0) < 0


def test_charge_is_zero_at_pi() -> None:
    protein = ProteinSequence(UBIQUITIN)
    assert protein.net_charge(protein.isoelectric_point()) == pytest.approx(
        0, abs=0.01)


def test_pi_with_no_charged_side_chains() -> None:
    # Only the two ends are charged, so pI is halfway between their pKas.
    assert ProteinSequence("GGG").isoelectric_point() == pytest.approx(
        (8.6 + 3.6) / 2, abs=0.01)


def test_basic_peptide_has_higher_pi_than_acidic_peptide() -> None:
    basic = ProteinSequence("KKKRRR")
    acidic = ProteinSequence("DDDEEE")
    assert basic.isoelectric_point() > 10
    assert acidic.isoelectric_point() < 4


# ---------- Extinction coefficient ----------

def test_extinction_coefficient() -> None:
    protein = ProteinSequence("WYCC")
    assert protein.extinction_coefficient() == 5500 + 1490 + 125
    assert protein.extinction_coefficient(all_cystines=False) == 5500 + 1490


def test_odd_cysteine_count_rounds_down_pairs() -> None:
    assert ProteinSequence("CCC").extinction_coefficient() == 125


# ---------- GRAVY ----------

def test_gravy() -> None:
    assert ProteinSequence("AAA").gravy() == pytest.approx(1.8)
    assert ProteinSequence("IR").gravy() == pytest.approx((4.5 - 4.5) / 2)
