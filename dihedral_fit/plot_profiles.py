#!/usr/bin/env python3
"""Plot one or more torsional energy profiles on the same axes.

Every profile is a two-column file, angle(deg) and energy, as written by
gaussian_dihedral.py, gromacs_dihedral.py or any other source. Each curve is
shifted so that its own minimum sits at zero.

Examples
--------
  # the QM profile alone
  plot_profiles.py qm_scan.dat --labels "QM, M062X/6-311++G(d,p)" \\
      -o dihedral_m062x.png

  # QM against the MM profile with the torsion removed
  plot_profiles.py qm_scan.dat mm_scan.dat \\
      --labels "QM" "MM, torsion removed" -o dihedral_qm_mm.png

  # and after the fitted terms are put back in
  plot_profiles.py qm_scan.dat mm_scan.dat mm_corrected.dat \\
      --labels "QM" "MM, torsion removed" "MM, fitted torsion" \\
      --rmsd-against 1 -o dihedral_qm_mm_corr.png
"""

import argparse
import sys

import numpy as np

KJ_PER_KCAL = 4.184
COLOURS = ["black", "#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
STYLES = ["-", "--", "-.", ":", "-", "--"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("files", nargs="+", help="two-column profiles, angle(deg) and energy")
    p.add_argument("--labels", nargs="+", default=None,
                   help="legend entry of each file (default: the file names)")
    p.add_argument("-o", "--output", default="profile.png",
                   help="image to write (default: %(default)s)")
    p.add_argument("--dpi", type=int, default=600,
                   help="resolution of the saved image (default: %(default)s)")
    p.add_argument("--units", choices=["kcal", "kJ"], default="kcal",
                   help="units of the energy column of the input files "
                        "(default: %(default)s)")
    p.add_argument("--plot-units", choices=["kcal", "kJ"], default=None,
                   help="units of the vertical axis (default: the input units)")
    p.add_argument("--title", default=None, help="title of the figure")
    p.add_argument("--xlabel", default="dihedral angle (degrees)",
                   help="label of the horizontal axis (default: %(default)s)")
    p.add_argument("--rmsd-against", type=int, default=None, metavar="I",
                   help="add the RMSD of every other curve against file number I "
                        "(1-based) to the legend")
    p.add_argument("--figsize", nargs=2, type=float, default=[7.0, 4.4],
                   metavar=("W", "H"), help="figure size in inches (default: 7.0 4.4)")
    p.add_argument("--ymax", type=float, default=None, help="upper limit of the energy axis")
    p.add_argument("--legend-loc", default="outside",
                   help="where the legend goes: a matplotlib location, in which case "
                        "the energy axis is extended so that it clears the curves, or "
                        "'outside', the default, to put it above the axes (default: %(default)s)")
    p.add_argument("--no-markers", action="store_true", help="draw lines only")
    p.add_argument("--no-show", action="store_true",
                   help="only write the file, do not open a window")
    return p.parse_args(argv)


def read_profile(path, to_kj):
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise SystemExit(f"{path}: expected two columns, angle(deg) and energy")
    angles = data[:, 0].astype(float)
    energies = data[:, 1].astype(float) * to_kj
    order = np.argsort(angles)
    angles, energies = angles[order], energies[order]
    return angles, energies - energies.min()


def rmsd(reference, other):
    ra, re_ = reference
    oa, oe = other
    interpolated = np.interp(np.mod(ra, 360.0), np.mod(oa, 360.0), oe, period=360.0)
    return float(np.sqrt(np.mean((interpolated - re_) ** 2)))


def make_room_for(fig, ax, legend, data_max, fixed_top):
    """Raise the top of the axes until a legend drawn inside clears the curves.

    The legend is laid out in points, so its height as a fraction of the axes
    does not change when the limits do; one pass is enough.
    """
    if fixed_top is not None or data_max <= 0:
        return
    fig.canvas.draw()
    try:
        box = legend.get_window_extent().transformed(ax.transAxes.inverted())
    except Exception:
        return
    bottom, top = ax.get_ylim()
    free = 1.0 - box.height - 0.05           # fraction of the axes left for the data
    if free <= 0.2:
        return                               # legend too tall to fit inside; leave it
    needed = bottom + (data_max - bottom) / free
    if needed > top:
        ax.set_ylim(bottom, needed)


def main(argv=None):
    args = parse_args(argv)
    labels = args.labels or args.files
    if len(labels) != len(args.files):
        raise SystemExit("--labels needs one entry per file")

    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    to_kj = KJ_PER_KCAL if args.units == "kcal" else 1.0
    plot_units = args.plot_units or args.units
    out_factor = 1.0 if plot_units == "kJ" else 1.0 / KJ_PER_KCAL
    unit_label = "kJ mol$^{-1}$" if plot_units == "kJ" else "kcal mol$^{-1}$"

    profiles = [read_profile(f, to_kj) for f in args.files]

    if args.rmsd_against:
        i = args.rmsd_against - 1
        if not 0 <= i < len(profiles):
            raise SystemExit("--rmsd-against is out of range")
        for j, profile in enumerate(profiles):
            if j == i:
                continue
            value = rmsd(profiles[i], profile) * out_factor
            labels[j] = f"{labels[j]} (RMSD {value:.2f})"

    fig, ax = plt.subplots(figsize=tuple(args.figsize), facecolor="white")
    ax.set_facecolor("white")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.set_axisbelow(True)

    for n, ((angles, energies), label) in enumerate(zip(profiles, labels)):
        ax.plot(angles, energies * out_factor,
                marker="" if args.no_markers else MARKERS[n % len(MARKERS)],
                linestyle=STYLES[n % len(STYLES)],
                color=COLOURS[n % len(COLOURS)],
                markersize=3.6, markerfacecolor="white", markeredgewidth=0.9,
                linewidth=1.4, label=label)

    ax.set_xlim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    if args.ymax is not None:
        ax.set_ylim(top=args.ymax)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(f"relative energy ({unit_label})")

    # a single curve needs no legend: its label reads better as a title
    title = args.title
    if title is None and len(profiles) == 1 and args.labels:
        title = labels[0]
    if title:
        ax.set_title(title, fontsize=11)

    if len(profiles) > 1:
        if args.legend_loc == "outside":
            legend = ax.legend(frameon=False, fontsize=9, loc="lower left",
                               bbox_to_anchor=(0.0, 1.01), borderaxespad=0.0)
        else:
            legend = ax.legend(frameon=False, fontsize=9, loc=args.legend_loc)
            make_room_for(fig, ax, legend, max(e.max() for _a, e in profiles) * out_factor,
                          args.ymax)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, facecolor="white", bbox_inches="tight")
    print(f"{args.output} ({args.dpi} dpi)", file=sys.stderr)
    if not args.no_show:
        try:
            plt.show()
        except Exception as err:          # no display available
            print(f"note: cannot open a window ({err})", file=sys.stderr)
    plt.close(fig)


if __name__ == "__main__":
    main()
