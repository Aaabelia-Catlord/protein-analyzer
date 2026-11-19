# Protein Analyzer

A Python command-line tool that calculates the physical and chemical properties of a protein from its amino acid sequence.

I built this to connect what I'm learning in biochemistry (protein structure, acid–base chemistry, the Beer–Lambert law) with what I'm learning in computer science (object-oriented design, algorithms and testing).

## Features

Give it a sequence (typed in, or from a FASTA file) and it reports:

| Property | How it's calculated |
|---|---|
| **Molecular weight** | Sum of average residue masses + one water for the free termini |
| **Net charge at any pH** | Henderson–Hasselbalch equation applied to every ionizable group |
| **Isoelectric point (pI)** | Binary search (bisection) for the pH where net charge = 0 |
| **Extinction coefficient (280 nm)** | Trp, Tyr and cystine contributions (Pace et al., 1995), for use with the Beer–Lambert law |
| **GRAVY score** | Average Kyte–Doolittle hydropathy |
| **Amino acid composition** | Count and percentage of each residue |

## Example

```
$ python -m protein_analyzer --fasta examples/ubiquitin.fasta --ph 7.4
Protein: Human ubiquitin
Length: 76 residues
Molecular weight: 8564.84 Da
Isoelectric point (pI): 7.54
Net charge at pH 7.4: +0.05
Extinction coefficient (280 nm): 1490 M^-1 cm^-1 (all Cys as cystines), 1490 (all Cys reduced)
GRAVY: -0.489
Composition:
  A:    2  (2.6%)
  ...
```

You can also pass a sequence directly:

```
python -m protein_analyzer MQIFVKTLTGKTITLEVEPSDTIENVK
```

## Getting started

Requires Python 3.10+.

```
git clone https://github.com/<your-username>/protein-analyzer.git
cd protein-analyzer
pip install -r requirements.txt   # only needed to run the tests
python -m protein_analyzer --fasta examples/ubiquitin.fasta
```

Run the tests:

```
python -m pytest --doctest-modules
```

## Project structure

```
protein_analyzer/
    constants.py   # residue masses, pKa values, hydropathy scale
    sequence.py    # ProteinSequence class with all calculations
    fasta.py       # FASTA file parsing
    cli.py         # command-line interface (argparse)
tests/             # pytest unit tests
examples/          # sample FASTA file
```

## Validation

I checked the results for human ubiquitin against [ExPASy ProtParam](https://web.expasy.org/protparam/):

| Property | This tool | ExPASy |
|---|---|---|
| Molecular weight | 8564.84 Da | 8564.84 Da |
| GRAVY | −0.489 | −0.489 |
| pI | 7.54 | 6.56 |

The pI is different because pI depends heavily on which pKa values you use. This tool uses the EMBOSS pKa set. ExPASy uses a different method, with pKa values that change depending on neighbouring residues. Published pKa tables don't agree with each other, so calculated pI values are only estimates.

## Possible future improvements

- Let the user choose a pKa set (EMBOSS, Lehninger, ExPASy-style)
- Plot net charge against pH (titration curve) with matplotlib
- Simulate a trypsin digest and predict peptide masses for mass spectrometry
- Support non-standard residue codes (e.g. X, B, Z)

## References

- Kyte J, Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. *J Mol Biol* 157:105–132.
- Pace CN et al. (1995). How to measure and predict the molar absorption coefficient of a protein. *Protein Sci* 4:2411–2423.
- Gasteiger E et al. (2005). Protein identification and analysis tools on the ExPASy server. *The Proteomics Protocols Handbook*, 571–607.
