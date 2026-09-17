#!/usr/bin/env python3
"""Compare the binary scalar fields of two GRHD legacy VTK files."""

from __future__ import annotations

import argparse
import array
import math
import pathlib
import sys


def read_line(stream) -> str:
    raw = stream.readline()
    if not raw:
        raise ValueError("fin inesperado del archivo VTK")
    return raw.decode("ascii").strip()


def read_big_endian_doubles(stream, count: int) -> array.array:
    raw = stream.read(8 * count)
    if len(raw) != 8 * count:
        raise ValueError("bloque binario VTK incompleto")
    values = array.array("d")
    values.frombytes(raw)
    if sys.byteorder == "little":
        values.byteswap()
    if stream.read(1) != b"\n":
        raise ValueError("separador VTK invalido despues de un bloque binario")
    return values


def read_vtk_fields(path: pathlib.Path) -> tuple[int, dict[str, array.array]]:
    with path.open("rb") as stream:
        if not read_line(stream).startswith("# vtk DataFile"):
            raise ValueError(f"{path}: encabezado VTK invalido")
        read_line(stream)
        if read_line(stream) != "BINARY":
            raise ValueError(f"{path}: se esperaba un VTK binario")
        if read_line(stream) != "DATASET STRUCTURED_GRID":
            raise ValueError(f"{path}: este comparador espera STRUCTURED_GRID")

        dimensions = read_line(stream).split()
        if len(dimensions) != 4 or dimensions[0] != "DIMENSIONS":
            raise ValueError(f"{path}: DIMENSIONS invalido")

        points = read_line(stream).split()
        if len(points) != 3 or points[0] != "POINTS" or points[2] != "double":
            raise ValueError(f"{path}: POINTS invalido")
        read_big_endian_doubles(stream, 3 * int(points[1]))

        cell_data = read_line(stream).split()
        if len(cell_data) != 2 or cell_data[0] != "CELL_DATA":
            raise ValueError(f"{path}: CELL_DATA invalido")
        cell_count = int(cell_data[1])
        fields: dict[str, array.array] = {}

        while True:
            raw = stream.readline()
            if not raw:
                break
            line = raw.decode("ascii").strip()
            if not line:
                continue
            scalar = line.split()
            if len(scalar) != 4 or scalar[0] != "SCALARS" or scalar[2:] != ["double", "1"]:
                raise ValueError(f"{path}: declaracion SCALARS invalida: {line}")
            if read_line(stream) != "LOOKUP_TABLE default":
                raise ValueError(f"{path}: LOOKUP_TABLE invalido")
            fields[scalar[1]] = read_big_endian_doubles(stream, cell_count)

    return cell_count, fields


def compare(reference: pathlib.Path, candidate: pathlib.Path, relative_tolerance: float,
            absolute_tolerance: float) -> int:
    reference_count, reference_fields = read_vtk_fields(reference)
    candidate_count, candidate_fields = read_vtk_fields(candidate)
    if reference_count != candidate_count:
        raise ValueError("los archivos tienen numeros de celdas distintos")
    if reference_fields.keys() != candidate_fields.keys():
        raise ValueError("los archivos no contienen los mismos campos escalares")

    print(f"Celdas comparadas: {reference_count}")
    print(f"{'campo':<24} {'L1 abs':>13} {'L2 abs':>13} {'Linf abs':>13} {'L2 rel':>13} estado")
    failed = False
    for name in reference_fields:
        left = reference_fields[name]
        right = candidate_fields[name]
        absolute_sum = 0.0
        squared_sum = 0.0
        reference_squared_sum = 0.0
        maximum = 0.0
        for reference_value, candidate_value in zip(left, right):
            if not math.isfinite(reference_value) or not math.isfinite(candidate_value):
                raise ValueError(f"el campo {name} contiene valores no finitos")
            difference = candidate_value - reference_value
            absolute_difference = abs(difference)
            absolute_sum += absolute_difference
            squared_sum += difference * difference
            reference_squared_sum += reference_value * reference_value
            maximum = max(maximum, absolute_difference)

        l1 = absolute_sum / reference_count
        l2 = math.sqrt(squared_sum / reference_count)
        reference_l2 = math.sqrt(reference_squared_sum / reference_count)
        relative_l2 = l2 / reference_l2 if reference_l2 > 0.0 else l2
        passed = maximum <= absolute_tolerance or relative_l2 <= relative_tolerance
        failed = failed or not passed
        print(f"{name:<24} {l1:13.5e} {l2:13.5e} {maximum:13.5e} "
              f"{relative_l2:13.5e} {'OK' if passed else 'REVISAR'}")

    print(f"Tolerancias: L2_rel <= {relative_tolerance:.3e} o Linf_abs <= {absolute_tolerance:.3e}")
    return 2 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=pathlib.Path, help="VTK de referencia CPU")
    parser.add_argument("candidate", type=pathlib.Path, help="VTK candidato GPU")
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-12)
    arguments = parser.parse_args()
    try:
        return compare(arguments.reference, arguments.candidate,
                       arguments.relative_tolerance, arguments.absolute_tolerance)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
