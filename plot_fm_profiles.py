#!/usr/bin/env python3
"""Plot radial Fishbone-Moncrief profiles from a reproducible VTK selection."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


DEFAULT_PATTERN = "FM_KS_a0.9d0_data/weno5_hlle/*.vtk"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay radial profiles from Fishbone-Moncrief VTK frames."
    )
    parser.add_argument(
        "pattern",
        nargs="?",
        default=DEFAULT_PATTERN,
        help=f"Glob pattern for VTK frames (default: {DEFAULT_PATTERN})",
    )
    parser.add_argument("--field", default="Density", help="VTK point-data field to plot.")
    parser.add_argument("--nx", type=int, default=400, help="Radial cells per frame.")
    parser.add_argument("--frames", type=int, default=36, help="Maximum number of frames.")
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Plot the field without normalization to its initial maximum.",
    )
    parser.add_argument("--xlim", type=float, nargs=2, metavar=("MIN", "MAX"))
    parser.add_argument("--ylim", type=float, nargs=2, metavar=("MIN", "MAX"))
    parser.add_argument("--output", type=Path, help="Optional image path to save.")
    parser.add_argument("--no-show", action="store_true", help="Do not open an interactive window.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    files = sorted(Path().glob(args.pattern))
    if not files:
        raise SystemExit(f"No VTK files match: {args.pattern}")

    if args.frames < 1 or args.nx < 1:
        raise SystemExit("--frames and --nx must both be positive.")

    files = files[: args.frames]
    fig, axis = plt.subplots(figsize=(12, 7))
    scale = None

    for index, path in enumerate(files):
        mesh = pv.read(path)
        if args.field not in mesh.point_data:
            raise SystemExit(f"{path}: missing point-data field {args.field!r}")

        values = mesh.point_data[args.field][: args.nx]
        coordinates = mesh.points[: args.nx]
        if len(values) != args.nx or len(coordinates) != args.nx:
            raise SystemExit(f"{path}: fewer than {args.nx} radial points.")

        radius = np.linalg.norm(coordinates, axis=1)
        if scale is None:
            scale = float(np.max(np.abs(values)))
            if not np.isfinite(scale) or scale <= 0.0:
                raise SystemExit(f"{path}: invalid initial normalization scale {scale}.")

        profile = values if args.raw else values / scale
        label = "Temporal evolution" if index == 0 else None
        axis.plot(radius, profile, color="teal", alpha=0.5, label=label)

    axis.set_xlabel("Radius (r)", fontsize=14)
    ylabel = args.field if args.raw else f"{args.field} (normalized)"
    axis.set_ylabel(ylabel, fontsize=14)
    axis.set_title("Equatorial Fishbone-Moncrief torus evolution", fontsize=16)
    axis.grid(True, linestyle="--", alpha=0.5)
    if args.xlim:
        axis.set_xlim(args.xlim)
    if args.ylim:
        axis.set_ylim(args.ylim)
    if axis.get_legend_handles_labels()[0]:
        axis.legend(loc="upper right")
    fig.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=200)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
