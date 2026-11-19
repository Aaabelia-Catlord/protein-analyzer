"""Reference values used to analyze protein sequences.

Sources:
- Average residue masses: ExPASy ProtParam (average isotopic masses, in Da).
- pKa values: EMBOSS set (commonly used by sequence-analysis tools).
- Hydropathy: Kyte & Doolittle (1982), J. Mol. Biol. 157:105-132.
- Extinction coefficients at 280 nm: Pace et al. (1995), Protein Sci. 4:2411-2423.
"""

# Mass of each amino acid *residue* (the amino acid minus one water,
# because a water is lost every time a peptide bond forms).
RESIDUE_MASSES: dict[str, float] = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}

# One water is added back for the free N-terminus (H) and C-terminus (OH).
WATER_MASS: float = 18.01524

# pKa values for ionizable groups.
N_TERMINUS_PKA: float = 8.6
C_TERMINUS_PKA: float = 3.6

# Side chains that are positively charged when protonated.
POSITIVE_SIDE_CHAIN_PKA: dict[str, float] = {"K": 10.8, "R": 12.5, "H": 6.5}

# Side chains that are negatively charged when deprotonated.
NEGATIVE_SIDE_CHAIN_PKA: dict[str, float] = {
    "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1,
}

# Kyte-Doolittle hydropathy (positive = hydrophobic, negative = hydrophilic).
HYDROPATHY: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Molar extinction coefficients at 280 nm (M^-1 cm^-1).
TRP_EXTINCTION: int = 5500
TYR_EXTINCTION: int = 1490
CYSTINE_EXTINCTION: int = 125  # per disulfide bond (pair of cysteines)

# Full names of the amino acids, for reports.
AMINO_ACID_NAMES: dict[str, str] = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartate",
    "C": "Cysteine", "E": "Glutamate", "Q": "Glutamine", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
}
