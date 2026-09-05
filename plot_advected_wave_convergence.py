#!/usr/bin/env python3
"""Grafica normas y ordenes observados para la onda advectada con HLLC."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


METHOD_ORDER = ["godunov", "mc", "weno3", "mp5", "weno5"]
METHOD_LABEL = {
    "godunov": "Godunov",
    "mc": "TVD-MC",
    "weno3": "WENO3",
    "mp5": "MP5",
    "weno5": "WENO5",
}
NORM_NAMES = ("L1", "L2", "Linf")


def read_results(data_directory: Path) -> dict[str, np.ndarray]:
    results: dict[str, np.ndarray] = {}
    for method in METHOD_ORDER:
        path = data_directory / f"convergence_{method}_hllc.dat"
        if not path.is_file():
            raise FileNotFoundError(f"Missing convergence data: {path}")
        values = np.loadtxt(path, comments="#", ndmin=2)
        if values.shape[1] != 4 or values.shape[0] < 2:
            raise ValueError(f"Unexpected convergence table shape in {path}: {values.shape}")
        values = values[np.argsort(values[:, 0])]
        if np.any(values[:, 1:] <= 0.0) or not np.all(np.isfinite(values)):
            raise ValueError(f"Non-positive or non-finite norm in {path}")
        results[method] = values
    return results


def observed_orders(values: np.ndarray) -> np.ndarray:
    resolutions = values[:, 0]
    refinement = resolutions[1:] / resolutions[:-1]
    return np.log(values[:-1, 1:] / values[1:, 1:]) / np.log(refinement[:, None])


def plot_norms(results: dict[str, np.ndarray], output_path: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for norm_index, (axis, norm_name) in enumerate(zip(axes, NORM_NAMES), start=1):
        for method in METHOD_ORDER:
            values = results[method]
            axis.loglog(
                values[:, 0],
                values[:, norm_index],
                marker="o",
                linewidth=1.8,
                label=METHOD_LABEL[method],
            )
        axis.set_title(rf"Norma $L_{{{norm_name[1:]}}}$" if norm_name != "Linf" else r"Norma $L_\infty$")
        axis.set_xlabel("Numero de celdas N")
        axis.set_ylabel("Error en densidad")
        axis.grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=9)
    figure.suptitle("Onda sinusoidal advectada: convergencia con HLLC")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_orders(results: dict[str, np.ndarray], output_path: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for norm_index, (axis, norm_name) in enumerate(zip(axes, NORM_NAMES)):
        for method in METHOD_ORDER:
            values = results[method]
            orders = observed_orders(values)
            axis.semilogx(
                values[1:, 0],
                orders[:, norm_index],
                marker="o",
                linewidth=1.8,
                label=METHOD_LABEL[method],
            )
        for reference_order in (1, 2, 3, 5):
            axis.axhline(reference_order, color="0.65", linewidth=0.8, linestyle="--")
        axis.set_title(rf"Orden observado en $L_{{{norm_name[1:]}}}$" if norm_name != "Linf" else r"Orden observado en $L_\infty$")
        axis.set_xlabel("Resolucion fina N")
        axis.set_ylabel(r"$p=\log(E_h/E_{h/2})/\log(2)$")
        axis.set_ylim(0.0, 6.0)
        axis.grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=9)
    figure.suptitle("Orden de convergencia observado con SSP-RK3 + HLLC")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_order_table(results: dict[str, np.ndarray], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("method", "N_coarse", "N_fine", "order_L1", "order_L2", "order_Linf"))
        for method in METHOD_ORDER:
            values = results[method]
            orders = observed_orders(values)
            for index, order_values in enumerate(orders):
                writer.writerow(
                    (
                        method,
                        int(values[index, 0]),
                        int(values[index + 1, 0]),
                        *(f"{value:.10e}" for value in order_values),
                    )
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_directory", type=Path)
    parser.add_argument("--output-directory", type=Path)
    arguments = parser.parse_args()

    data_directory = arguments.data_directory.resolve()
    output_directory = (arguments.output_directory or data_directory.parent).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    results = read_results(data_directory)

    plot_norms(results, output_directory / "advected_wave_hllc_norms.png")
    plot_orders(results, output_directory / "advected_wave_hllc_orders.png")
    write_order_table(results, output_directory / "advected_wave_hllc_observed_orders.csv")

    for method in METHOD_ORDER:
        values = results[method]
        last_orders = observed_orders(values)[-1]
        print(
            f"{METHOD_LABEL[method]:8s} N={int(values[-1, 0]):4d} "
            f"L1={values[-1, 1]:.6e} L2={values[-1, 2]:.6e} "
            f"Linf={values[-1, 3]:.6e} "
            f"p=({last_orders[0]:.3f}, {last_orders[1]:.3f}, {last_orders[2]:.3f})"
        )


if __name__ == "__main__":
    main()
