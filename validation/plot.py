"""Render the recorded Pb–Pb geometry and TGlauberMC validation figures."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


NATIVE = "#146c83"
REFERENCE = "#b45726"
NEUTRAL = "#526980"


def geometry(events, report, output):
    b = events["b_fm"]
    collisions = events["ncoll"]
    central_edge = report["bins"][1]["b_high_fm"]["native"]
    b_edges = np.arange(0, 20.0001, .4)
    minimum_bias, _ = np.histogram(b, b_edges, density=True)
    hard, _ = np.histogram(b, b_edges, weights=collisions, density=True)

    fig, (impact, overlap) = plt.subplots(1, 2, figsize=(11.7, 4.7), layout="constrained")
    fig.suptitle("Pb–Pb at 5.02 TeV · minimum-bias geometry", fontsize=16, weight="bold")
    impact.axvspan(0, central_edge, color=NATIVE, alpha=.10, label="0–10% class")
    impact.stairs(minimum_bias, b_edges, color=NATIVE, lw=2.1, label="Minimum bias")
    impact.stairs(hard, b_edges, color=REFERENCE, lw=2.1, label=r"$N_{\mathrm{coll}}$-weighted")
    impact.set(xlim=(0, 20), ylim=(0, None), xlabel="Impact parameter $b$ [fm]",
               ylabel="Probability density [fm$^{-1}$]", title="Impact parameter and hard-event weighting")
    impact.legend(frameon=False, loc="upper right")
    impact.grid(axis="y", alpha=.18)

    density = overlap.hexbin(events["npart"], collisions, gridsize=48, mincnt=1, bins="log",
                             cmap="viridis", linewidths=0)
    overlap.set(xlabel=r"Participant nucleons $N_{\mathrm{part}}$",
                ylabel=r"Binary collisions $N_{\mathrm{coll}}$", title="Participants and binary collisions")
    fig.colorbar(density, ax=overlap, label="Events per hexagon", pad=.025)
    fig.savefig(output / "pbpb-geometry.png", dpi=150)
    plt.close(fig)


def validation(report, output):
    bins = report["bins"]
    centers = np.array([np.mean(row["centrality_percent"]) for row in bins])
    fig, panels = plt.subplots(1, 3, figsize=(13.4, 4.6), layout="constrained")
    fig.suptitle("Glauber C++ vs TGlauberMC v3.2 · Pb–Pb at 5.02 TeV", fontsize=16, weight="bold")
    fields = ("npart_mean", "taa_mean_mb_inverse")
    labels = (r"$\langle N_{\mathrm{part}}\rangle$", r"$\langle T_{\mathrm{AA}}\rangle$ [mb$^{-1}$]")
    for ax, field, label in zip(panels[:2], fields, labels):
        for key, name, color, marker in (
            ("native", "Glauber C++", NATIVE, "o"),
            ("reference", "TGlauberMC", REFERENCE, "s"),
        ):
            ax.errorbar(centers, [row[field][key] for row in bins],
                        yerr=[row[field][key + "_se"] for row in bins],
                        fmt=marker, lw=1, ms=4, fillstyle="none" if key == "reference" else "full",
                        capsize=2, color=color, label=name)
        ax.set(ylabel=label, title=label)
        ax.legend(frameon=False, fontsize=9)
    panels[1].set_yscale("log")

    pull = panels[2]
    pull.axhspan(-4, 4, color=NEUTRAL, alpha=.09, label="Acceptance band")
    pull.axhline(0, color=NEUTRAL, lw=.8)
    boundaries = [(row["centrality_percent"][1], row["b_high_fm"]["pull"])
                  for row in bins if row["centrality_percent"][1] < 100]
    boundary_x, boundary_y = zip(*boundaries)
    pull.plot(boundary_x, boundary_y, ".", color=NEUTRAL, label=r"$b$ upper edge")
    pull.plot(centers, [row[fields[0]]["pull"] for row in bins], "o",
              ms=3.3, color=NATIVE, label=r"$\langle N_{\mathrm{part}}\rangle$")
    pull.plot(centers, [row[fields[1]]["pull"] for row in bins], "s",
              ms=3.3, color=REFERENCE, label=r"$\langle T_{\mathrm{AA}}\rangle$")
    pull.set(ylim=(-4.5, 4.5), ylabel="Difference / combined MC standard error",
             title="Monte Carlo pulls")
    pull.legend(frameon=False, fontsize=8, loc="lower left")
    for ax in panels:
        ax.set(xlim=(0, 100), xlabel="Centrality [%]")
        ax.grid(axis="y", alpha=.17)
    fig.savefig(output / "pbpb-validation.png", dpi=150)
    plt.close(fig)


def main():
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=folder.parent / "docs/figures")
    args = parser.parse_args()

    path = folder / "native-run/events.csv"
    report = json.loads((folder / "pbpb-5020.json").read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == report["native"]["events_sha256"]
    events = np.genfromtxt(path, names=True, delimiter=",")
    assert len(events) == report["native"]["inelastic_events"]
    assert report["passed"] and len(report["bins"]) == 20

    plt.rcParams.update({"font.size": 10, "axes.titleweight": "bold", "figure.facecolor": "white"})
    args.output.mkdir(parents=True, exist_ok=True)
    geometry(events, report, args.output)
    validation(report, args.output)


if __name__ == "__main__":
    main()
