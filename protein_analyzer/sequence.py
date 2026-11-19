"""The ProteinSequence class: calculates properties of a protein from its
amino acid sequence."""

from __future__ import annotations

from collections import Counter

from protein_analyzer.constants import (
    C_TERMINUS_PKA,
    CYSTINE_EXTINCTION,
    HYDROPATHY,
    N_TERMINUS_PKA,
    NEGATIVE_SIDE_CHAIN_PKA,
    POSITIVE_SIDE_CHAIN_PKA,
    RESIDUE_MASSES,
    TRP_EXTINCTION,
    TYR_EXTINCTION,
    WATER_MASS,
)


class ProteinSequence:
    """A protein written as one-letter amino acid codes.

    Attributes:
        name: A label for the protein (for example, a FASTA header).
        residues: The sequence in uppercase one-letter codes.
        source: A link to where the sequence came from (e.g. its UniProt
            page), or None if it was typed in or read from a file.

    >>> protein = ProteinSequence("gg")
    >>> protein.residues
    'GG'
    >>> round(protein.molecular_weight(), 2)
    132.12
    """

    name: str
    residues: str
    source: str | None

    def __init__(self, residues: str, name: str = "unnamed") -> None:
        """Create a protein from a sequence string.

        Whitespace is removed and letters are made uppercase.

        Raises:
            ValueError: If the sequence is empty or has a letter that is not
                one of the 20 standard amino acids.
        """
        cleaned = "".join(residues.split()).upper()
        if not cleaned:
            raise ValueError("Sequence is empty.")

        invalid = sorted(set(cleaned) - set(RESIDUE_MASSES))
        if invalid:
            raise ValueError(
                f"Invalid amino acid code(s): {', '.join(invalid)}"
            )

        self.residues = cleaned
        self.name = name
        self.source = None
        self._counts = Counter(cleaned)

    def __len__(self) -> int:
        """Return the number of residues."""
        return len(self.residues)

    def __repr__(self) -> str:
        return f"ProteinSequence(name={self.name!r}, length={len(self)})"

    def composition(self) -> dict[str, int]:
        """Return how many times each amino acid appears, in A-Z order.

        Amino acids that do not appear are left out.
        """
        return {aa: self._counts[aa] for aa in sorted(self._counts)}

    def percent_composition(self) -> dict[str, float]:
        """Return the percentage of the sequence made up by each amino acid."""
        total = len(self)
        return {aa: 100 * n / total for aa, n in self.composition().items()}

    def molecular_weight(self) -> float:
        """Return the average molecular weight in daltons (Da).

        Sum of the residue masses plus one water for the two free ends.
        """
        residue_total = sum(
            RESIDUE_MASSES[aa] * n for aa, n in self._counts.items()
        )
        return residue_total + WATER_MASS

    def net_charge(self, ph: float) -> float:
        """Return the net charge of the protein at the given pH.

        Uses the Henderson-Hasselbalch equation to find the fraction of each
        ionizable group that is charged:
            positive group charged fraction = 1 / (1 + 10^(pH - pKa))
            negative group charged fraction = 1 / (1 + 10^(pKa - pH))
        """
        positive = 1 / (1 + 10 ** (ph - N_TERMINUS_PKA))
        for aa, pka in POSITIVE_SIDE_CHAIN_PKA.items():
            positive += self._counts[aa] / (1 + 10 ** (ph - pka))

        negative = 1 / (1 + 10 ** (C_TERMINUS_PKA - ph))
        for aa, pka in NEGATIVE_SIDE_CHAIN_PKA.items():
            negative += self._counts[aa] / (1 + 10 ** (pka - ph))

        return positive - negative

    def isoelectric_point(self, tolerance: float = 0.001) -> float:
        """Return the isoelectric point (pI): the pH where net charge is 0.

        Net charge always goes down as pH goes up, so we can use binary
        search (bisection) between pH 0 and pH 14.
        """
        low, high = 0.0, 14.0
        while high - low > tolerance:
            middle = (low + high) / 2
            if self.net_charge(middle) > 0:
                low = middle   # still positive, so the pI is higher
            else:
                high = middle  # negative, so the pI is lower
        return (low + high) / 2

    def extinction_coefficient(self, all_cystines: bool = True) -> int:
        """Return the molar extinction coefficient at 280 nm (M^-1 cm^-1).

        Use with the Beer-Lambert law (A = e * c * l) to find protein
        concentration from absorbance.

        Args:
            all_cystines: If True, assume every pair of cysteines forms a
                disulfide bond (cystine). If False, assume none do.
        """
        total = (self._counts["W"] * TRP_EXTINCTION
                 + self._counts["Y"] * TYR_EXTINCTION)
        if all_cystines:
            total += (self._counts["C"] // 2) * CYSTINE_EXTINCTION
        return total

    def gravy(self) -> float:
        """Return the GRAVY score (grand average of hydropathy).

        Positive means hydrophobic overall; negative means hydrophilic.
        """
        total = sum(HYDROPATHY[aa] * n for aa, n in self._counts.items())
        return total / len(self)

    def summary(self, ph: float = 7.0) -> str:
        """Return a readable report of this protein's properties."""
        lines = [
            f"Protein: {self.name}",
            f"Length: {len(self)} residues",
            f"Molecular weight: {self.molecular_weight():.2f} Da",
            f"Isoelectric point (pI): {self.isoelectric_point():.2f}",
            f"Net charge at pH {ph}: {self.net_charge(ph):+.2f}",
            "Extinction coefficient (280 nm): "
            f"{self.extinction_coefficient()} M^-1 cm^-1 "
            "(all Cys as cystines), "
            f"{self.extinction_coefficient(all_cystines=False)} "
            "(all Cys reduced)",
            f"GRAVY: {self.gravy():.3f}",
            "Composition:",
        ]
        for aa, percent in self.percent_composition().items():
            lines.append(f"  {aa}: {self._counts[aa]:>4}  ({percent:.1f}%)")
        return "\n".join(lines)
