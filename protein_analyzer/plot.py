"""Plotting titration curves (net charge vs. pH) for proteins."""

from __future__ import annotations

from matplotlib.figure import Figure

from protein_analyzer.sequence import ProteinSequence


def titration_points(protein: ProteinSequence, step: float = 0.1
                     ) -> tuple[list[float], list[float]]:
    """Return (pH values, net charges) from pH 0 to pH 14.

    >>> phs, charges = titration_points(ProteinSequence("GG"), step=7.0)
    >>> phs
    [0.0, 7.0, 14.0]
    """
    count = round(14 / step)
    phs = [i * step for i in range(count + 1)]
    charges = [protein.net_charge(ph) for ph in phs]
    return phs, charges


def plot_titration_curves(proteins: list[ProteinSequence],
                          output_path: str) -> None:
    """Draw the titration curve of each protein and save it as an image.

    Each protein's pI (where the curve crosses zero) is marked with a dot.
    """
    fig = Figure(figsize=(7, 4.5))
    ax = fig.subplots()

    for protein in proteins:
        phs, charges = titration_points(protein)
        pi = protein.isoelectric_point()
        line, = ax.plot(phs, charges, linewidth=2,
                        label=f"{protein.name} (pI = {pi:.2f})")
        ax.plot(pi, 0, "o", color=line.get_color())

    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlim(0, 14)
    ax.set_xlabel("pH")
    ax.set_ylabel("Net charge")
    ax.set_title("Titration curve")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
