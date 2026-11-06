#!/usr/bin/env python3
"""Generate a tensor-product graded Hex8 mesh for CB-Geo MPM.

The generator is intentionally standalone and uses only the Python standard
library. It writes the ASCII files consumed by mpm::IOMeshAscii together with
an entity-set file and a reusable MPM input fragment.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence


EPS = 1.0e-12


class ConfigError(ValueError):
    """Raised when the generator configuration is invalid."""


@dataclass(frozen=True)
class AxisResult:
    """One-dimensional graded grid data."""

    coordinates: tuple[float, ...]
    widths: tuple[float, ...]
    core_start: int
    core_end: int
    requested_fine: float
    requested_coarse: float

    @property
    def cell_count(self) -> int:
        return len(self.widths)

    @property
    def actual_fine(self) -> float:
        return self.widths[self.core_start]

    @property
    def max_adjacent_ratio(self) -> float:
        if len(self.widths) < 2:
            return 1.0
        return max(
            self.widths[index + 1] / self.widths[index]
            for index in range(len(self.widths) - 1)
        )


@dataclass(frozen=True)
class OutputFiles:
    mesh: Path
    particles: Path
    particles_cells: Path
    particles_volumes: Path
    entity_sets: Path
    mpm_fragment: Path
    report: Path
    axis_spacing: Path
    center_ids: Path


def _require_mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{name} must be a JSON object")
    return value


def _number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{name} must be a number")
    result = float(value)
    if not math.isfinite(result):
        raise ConfigError(f"{name} must be finite")
    return result


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigError(f"{name} must be a positive integer")
    return value


def _triplet(value: Any, name: str, *, positive: bool = False) -> tuple[float, float, float]:
    if not isinstance(value, list) or len(value) != 3:
        raise ConfigError(f"{name} must be an array of three numbers")
    result = tuple(_number(item, f"{name}[{index}]") for index, item in enumerate(value))
    if positive and any(item <= 0.0 for item in result):
        raise ConfigError(f"{name} entries must be positive")
    return result  # type: ignore[return-value]


def _int_triplet(value: Any, name: str) -> tuple[int, int, int]:
    if not isinstance(value, list) or len(value) != 3:
        raise ConfigError(f"{name} must be an array of three positive integers")
    result = tuple(
        _positive_int(item, f"{name}[{index}]") for index, item in enumerate(value)
    )
    return result  # type: ignore[return-value]


def _as_list(value: Any, name: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ConfigError(f"{name} must be an array")
    return value


def _smoothstep(value: float) -> float:
    return value * value * (3.0 - 2.0 * value)


def _transition_sizes(
    length: float,
    fine: float,
    coarse: float,
    mode: str,
    max_adjacent_ratio: float,
) -> list[float]:
    """Return sizes ordered from the fine core outward to the coarse boundary.

    For the linear mode this is the requested discrete sequence

        N = round(2 L / (fine + coarse))
        h_i = fine + (coarse - fine) i / (N - 1)

    followed by a uniform scaling that makes sum(h_i) exactly equal to L.
    """

    if length <= EPS:
        return []

    if mode not in {"linear", "smoothstep"}:
        raise ConfigError(
            "spacing.transition must be either 'linear' or 'smoothstep'"
        )
    if coarse + EPS < fine:
        raise ConfigError("spacing.coarse must be greater than or equal to spacing.fine")

    count = max(1, int(round(2.0 * length / (fine + coarse))))
    if count == 1:
        sizes = [length]
    else:
        weights = [
            index / (count - 1) for index in range(count)
        ]
        if mode == "smoothstep":
            weights = [_smoothstep(weight) for weight in weights]
        sizes = [fine + (coarse - fine) * weight for weight in weights]
        scale = length / sum(sizes)
        sizes = [size * scale for size in sizes]

    ratio = (
        1.0
        if len(sizes) < 2
        else max(sizes[index + 1] / sizes[index] for index in range(len(sizes) - 1))
    )
    if ratio > max_adjacent_ratio + 1.0e-6:
        raise ConfigError(
            "The requested fine/coarse size range is too large for the available "
            "transition length while satisfying spacing.max_adjacent_ratio "
            f"({ratio:.4g} > {max_adjacent_ratio:.4g}). Increase the distance "
            "between the core and boundary, reduce spacing.coarse, or relax "
            "spacing.max_adjacent_ratio."
        )
    return sizes


def _build_axis(
    low: float,
    high: float,
    core_low: float,
    core_high: float,
    fine: float,
    coarse: float,
    transition_mode: str,
    max_adjacent_ratio: float,
) -> AxisResult:
    if low >= high:
        raise ConfigError("domain minimum must be smaller than domain maximum")
    if core_low >= core_high:
        raise ConfigError("core minimum must be smaller than core maximum")
    if core_low < low - EPS or core_high > high + EPS:
        raise ConfigError("core bounds must lie inside the domain")

    core_length = core_high - core_low
    core_intervals = max(1, int(round(core_length / fine)))
    actual_fine = core_length / core_intervals
    core_nodes = [
        core_low + actual_fine * index for index in range(core_intervals + 1)
    ]
    core_nodes[-1] = core_high

    left_sizes = _transition_sizes(
        core_low - low,
        actual_fine,
        coarse,
        transition_mode,
        max_adjacent_ratio,
    )
    right_sizes = _transition_sizes(
        high - core_high,
        actual_fine,
        coarse,
        transition_mode,
        max_adjacent_ratio,
    )

    left_nodes: list[float] = []
    if left_sizes:
        distance = 0.0
        distances = [0.0]
        for size in left_sizes:
            distance += size
            distances.append(distance)
        left_nodes = [core_low - value for value in reversed(distances)]
        left_nodes[0] = low
        left_nodes[-1] = core_low

    right_nodes: list[float] = []
    if right_sizes:
        distance = 0.0
        right_nodes.append(core_high)
        for size in right_sizes:
            distance += size
            right_nodes.append(core_high + distance)
        right_nodes[-1] = high
        right_nodes = right_nodes[1:]

    coordinates = tuple(
        (left_nodes[:-1] if left_nodes else []) + core_nodes + right_nodes
    )
    widths = tuple(
        coordinates[index + 1] - coordinates[index]
        for index in range(len(coordinates) - 1)
    )

    if any(width <= 0.0 for width in widths):
        raise ConfigError("axis generation produced a non-positive cell size")
    if abs(coordinates[0] - low) > 1.0e-9 * max(1.0, abs(high - low)):
        raise ConfigError("axis generation did not reach the lower domain bound")
    if abs(coordinates[-1] - high) > 1.0e-9 * max(1.0, abs(high - low)):
        raise ConfigError("axis generation did not reach the upper domain bound")

    core_start = max(0, len(left_nodes) - 1)
    core_end = core_start + core_intervals
    return AxisResult(
        coordinates=coordinates,
        widths=widths,
        core_start=core_start,
        core_end=core_end,
        requested_fine=fine,
        requested_coarse=coarse,
    )


def _point_in_box(
    point: Sequence[float], minimum: Sequence[float], maximum: Sequence[float]
) -> bool:
    return all(
        minimum[axis] - EPS <= point[axis] <= maximum[axis] + EPS
        for axis in range(3)
    )


def _point_in_shape(point: Sequence[float], shape: dict[str, Any], name: str) -> bool:
    shape_type = str(shape.get("type", "")).lower()
    if shape_type == "box":
        minimum = _triplet(shape.get("min"), f"{name}.min")
        maximum = _triplet(shape.get("max"), f"{name}.max")
        return _point_in_box(point, minimum, maximum)

    if shape_type == "sphere":
        center = _triplet(shape.get("center"), f"{name}.center")
        radius = _number(shape.get("radius"), f"{name}.radius")
        if radius <= 0.0:
            raise ConfigError(f"{name}.radius must be positive")
        distance_squared = sum((point[axis] - center[axis]) ** 2 for axis in range(3))
        return distance_squared <= radius * radius + EPS

    if shape_type == "cylinder":
        center = _triplet(shape.get("center"), f"{name}.center")
        radius = _number(shape.get("radius"), f"{name}.radius")
        axis_name = str(shape.get("axis", "z")).lower()
        if radius <= 0.0:
            raise ConfigError(f"{name}.radius must be positive")
        if axis_name not in {"x", "y", "z"}:
            raise ConfigError(f"{name}.axis must be x, y, or z")

        axis_index = {"x": 0, "y": 1, "z": 2}[axis_name]
        if "height" in shape:
            height = _number(shape.get("height"), f"{name}.height")
            if height <= 0.0:
                raise ConfigError(f"{name}.height must be positive")
            axis_minimum = center[axis_index] - 0.5 * height
            axis_maximum = center[axis_index] + 0.5 * height
        else:
            axis_minimum = _number(shape.get("min"), f"{name}.min")
            axis_maximum = _number(shape.get("max"), f"{name}.max")
        if axis_minimum > axis_maximum:
            raise ConfigError(f"{name} axis bounds are reversed")
        if not axis_minimum - EPS <= point[axis_index] <= axis_maximum + EPS:
            return False

        radial_squared = 0.0
        for index in range(3):
            if index != axis_index:
                radial_squared += (point[index] - center[index]) ** 2
        return radial_squared <= radius * radius + EPS

    raise ConfigError(f"{name}.type must be box, sphere, or cylinder")


def _removed_by_shapes(
    point: Sequence[float], remove_regions: Sequence[dict[str, Any]]
) -> bool:
    return any(
        _point_in_shape(point, shape, f"remove_regions[{index}]")
        for index, shape in enumerate(remove_regions)
    )


def _format(value: float) -> str:
    return f"{value:.17g}"


def _write_mesh(
    path: Path, x_axis: AxisResult, y_axis: AxisResult, z_axis: AxisResult
) -> tuple[int, int]:
    nx = x_axis.cell_count
    ny = y_axis.cell_count
    nz = z_axis.cell_count
    nnodes = (nx + 1) * (ny + 1) * (nz + 1)
    ncells = nx * ny * nz

    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("! Generated graded orthorhombic Hex8 mesh\n")
        stream.write("! elementShape hexahedron\n")
        stream.write("! elementNumPoints 8\n")
        stream.write(f"{nnodes}\t{ncells}\n")

        for z_value in z_axis.coordinates:
            for y_value in y_axis.coordinates:
                for x_value in x_axis.coordinates:
                    stream.write(
                        f"{_format(x_value)}\t{_format(y_value)}\t{_format(z_value)}\n"
                    )

        nx1 = nx + 1
        ny1 = ny + 1
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    n000 = k * nx1 * ny1 + j * nx1 + i
                    n100 = n000 + 1
                    n110 = n000 + nx1 + 1
                    n010 = n000 + nx1
                    n001 = n000 + nx1 * ny1
                    n101 = n001 + 1
                    n111 = n001 + nx1 + 1
                    n011 = n001 + nx1
                    stream.write(
                        f"{n000}\t{n100}\t{n110}\t{n010}\t"
                        f"{n001}\t{n101}\t{n111}\t{n011}\n"
                    )
    return nnodes, ncells


def _copy_body_with_count(source: Path, destination: Path, count: int) -> None:
    with destination.open("w", encoding="utf-8", newline="\n") as output:
        output.write(f"{count}\n")
        with source.open("r", encoding="utf-8") as body:
            shutil.copyfileobj(body, output)


def _generate_particles(
    output_dir: Path,
    x_axis: AxisResult,
    y_axis: AxisResult,
    z_axis: AxisResult,
    per_cell: tuple[int, int, int],
    particle_domain_min: tuple[float, float, float],
    particle_domain_max: tuple[float, float, float],
    remove_regions: Sequence[dict[str, Any]],
    files: OutputFiles,
    include_all_particle_set: bool,
    include_core_particle_set: bool,
) -> tuple[
    int,
    int,
    int,
    int,
    list[int] | None,
    list[int] | None,
    tuple[int, tuple[float, float, float], float] | None,
]:
    nx = x_axis.cell_count
    ny = y_axis.cell_count
    nz = z_axis.cell_count
    nx1 = nx + 1
    ny1 = ny + 1
    px, py, pz = per_cell
    particles_per_cell = px * py * pz

    body_particles = output_dir / ".particles.body.tmp"
    body_cells = output_dir / ".particles-cells.tmp"
    body_volumes = output_dir / ".particles-volumes.tmp"

    particle_count = 0
    removed_count = 0
    outside_particle_domain_count = 0
    active_cell_count = 0
    all_particle_ids: list[int] | None = [] if include_all_particle_set else None
    core_particle_ids: list[int] | None = [] if include_core_particle_set else None
    core_center = (
        0.5 * (x_axis.coordinates[x_axis.core_start] + x_axis.coordinates[x_axis.core_end]),
        0.5 * (y_axis.coordinates[y_axis.core_start] + y_axis.coordinates[y_axis.core_end]),
        0.5 * (z_axis.coordinates[z_axis.core_start] + z_axis.coordinates[z_axis.core_end]),
    )
    closest_particle: tuple[int, tuple[float, float, float], float] | None = None
    closest_particle_distance_squared = math.inf

    with (
        body_particles.open("w", encoding="utf-8", newline="\n") as particle_stream,
        body_cells.open("w", encoding="utf-8", newline="\n") as cells_stream,
        body_volumes.open("w", encoding="utf-8", newline="\n") as volume_stream,
    ):
        cell_id = 0
        for k in range(nz):
            z0 = z_axis.coordinates[k]
            dz = z_axis.widths[k]
            in_core_z = z_axis.core_start <= k < z_axis.core_end
            for j in range(ny):
                y0 = y_axis.coordinates[j]
                dy = y_axis.widths[j]
                in_core_y = y_axis.core_start <= j < y_axis.core_end
                for i in range(nx):
                    x0 = x_axis.coordinates[i]
                    dx = x_axis.widths[i]
                    in_core_x = x_axis.core_start <= i < x_axis.core_end
                    cell_volume = dx * dy * dz
                    particle_volume = cell_volume / particles_per_cell

                    cell_is_core = in_core_x and in_core_y and in_core_z
                    cell_has_particle = False
                    for pzz in range(pz):
                        z_value = z0 + (pzz + 0.5) * dz / pz
                        for pyy in range(py):
                            y_value = y0 + (pyy + 0.5) * dy / py
                            for pxx in range(px):
                                x_value = x0 + (pxx + 0.5) * dx / px
                                point = (x_value, y_value, z_value)
                                if not _point_in_box(
                                    point, particle_domain_min, particle_domain_max
                                ):
                                    outside_particle_domain_count += 1
                                    continue
                                if _removed_by_shapes(point, remove_regions):
                                    removed_count += 1
                                    continue

                                particle_stream.write(
                                    f"{_format(x_value)}\t{_format(y_value)}\t"
                                    f"{_format(z_value)}\n"
                                )
                                cells_stream.write(f"{particle_count}\t{cell_id}\n")
                                volume_stream.write(
                                    f"{particle_count}\t{_format(particle_volume)}\n"
                                )
                                if all_particle_ids is not None:
                                    all_particle_ids.append(particle_count)
                                if core_particle_ids is not None and cell_is_core:
                                    core_particle_ids.append(particle_count)
                                if cell_is_core:
                                    distance_squared = sum(
                                        (point[index] - core_center[index]) ** 2
                                        for index in range(3)
                                    )
                                    if distance_squared < closest_particle_distance_squared:
                                        closest_particle = (
                                            particle_count,
                                            point,
                                            math.sqrt(distance_squared),
                                        )
                                        closest_particle_distance_squared = distance_squared
                                particle_count += 1
                                cell_has_particle = True
                    if cell_has_particle:
                        active_cell_count += 1
                    cell_id += 1

    _copy_body_with_count(body_particles, files.particles, particle_count)
    body_cells.replace(files.particles_cells)
    body_volumes.replace(files.particles_volumes)
    body_particles.unlink(missing_ok=True)

    return (
        particle_count,
        removed_count,
        outside_particle_domain_count,
        active_cell_count,
        all_particle_ids,
        core_particle_ids,
        closest_particle,
    )


def _node_id(i: int, j: int, k: int, nx: int, ny: int) -> int:
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i


def _closest_node(
    x_axis: AxisResult,
    y_axis: AxisResult,
    z_axis: AxisResult,
    center: tuple[float, float, float],
) -> tuple[int, tuple[float, float, float]]:
    index_x = min(
        range(len(x_axis.coordinates)),
        key=lambda index: abs(x_axis.coordinates[index] - center[0]),
    )
    index_y = min(
        range(len(y_axis.coordinates)),
        key=lambda index: abs(y_axis.coordinates[index] - center[1]),
    )
    index_z = min(
        range(len(z_axis.coordinates)),
        key=lambda index: abs(z_axis.coordinates[index] - center[2]),
    )
    return (
        _node_id(index_x, index_y, index_z, x_axis.cell_count, y_axis.cell_count),
        (
            x_axis.coordinates[index_x],
            y_axis.coordinates[index_y],
            z_axis.coordinates[index_z],
        ),
    )


def _boundary_node_sets(
    nx: int, ny: int, nz: int, set_ids: dict[str, int]
) -> list[dict[str, Any]]:
    x_min: list[int] = []
    x_max: list[int] = []
    y_min: list[int] = []
    y_max: list[int] = []
    z_min: list[int] = []
    z_max: list[int] = []

    for k in range(nz + 1):
        for j in range(ny + 1):
            x_min.append(_node_id(0, j, k, nx, ny))
            x_max.append(_node_id(nx, j, k, nx, ny))
    for k in range(nz + 1):
        for i in range(nx + 1):
            y_min.append(_node_id(i, 0, k, nx, ny))
            y_max.append(_node_id(i, ny, k, nx, ny))
    for j in range(ny + 1):
        for i in range(nx + 1):
            z_min.append(_node_id(i, j, 0, nx, ny))
            z_max.append(_node_id(i, j, nz, nx, ny))

    return [
        {"id": set_ids["x_min"], "set": sorted(x_min)},
        {"id": set_ids["x_max"], "set": sorted(x_max)},
        {"id": set_ids["y_min"], "set": sorted(y_min)},
        {"id": set_ids["y_max"], "set": sorted(y_max)},
        {"id": set_ids["z_min"], "set": sorted(z_min)},
        {"id": set_ids["z_max"], "set": sorted(z_max)},
    ]


def _write_entity_sets(
    path: Path,
    axis_x: AxisResult,
    axis_y: AxisResult,
    axis_z: AxisResult,
    entity_config: dict[str, Any],
    all_particle_ids: list[int] | None,
    core_particle_ids: list[int] | None,
    center_node_id: int,
    center_particle_id: int | None,
) -> None:
    set_ids = {
        "x_min": int(entity_config.get("x_min_id", 1)),
        "x_max": int(entity_config.get("x_max_id", 2)),
        "y_min": int(entity_config.get("y_min_id", 3)),
        "y_max": int(entity_config.get("y_max_id", 4)),
        "z_min": int(entity_config.get("z_min_id", 5)),
        "z_max": int(entity_config.get("z_max_id", 6)),
    }
    if len(set(set_ids.values())) != 6:
        raise ConfigError("node-set IDs must be unique")

    entity_sets: dict[str, Any] = {
        "node_sets": _boundary_node_sets(
            axis_x.cell_count, axis_y.cell_count, axis_z.cell_count, set_ids
        )
    }

    particle_sets: list[dict[str, Any]] = []
    if all_particle_ids is not None:
        particle_sets.append(
            {"id": int(entity_config.get("all_particles_id", 0)), "set": all_particle_ids}
        )
    if core_particle_ids is not None:
        particle_sets.append(
            {"id": int(entity_config.get("core_particles_id", 10)), "set": core_particle_ids}
        )
    if center_particle_id is not None:
        particle_sets.append(
            {
                "id": int(entity_config.get("center_particle_set_id", 100)),
                "set": [center_particle_id],
            }
        )
    if particle_sets:
        entity_sets["particle_sets"] = particle_sets

    entity_sets["node_sets"].append(
        {
            "id": int(entity_config.get("center_node_set_id", 100)),
            "set": [center_node_id],
        }
    )

    if bool(entity_config.get("include_core_cell_set", False)):
        core_cells: list[int] = []
        nx = axis_x.cell_count
        ny = axis_y.cell_count
        for k in range(axis_z.core_start, axis_z.core_end):
            for j in range(axis_y.core_start, axis_y.core_end):
                for i in range(axis_x.core_start, axis_x.core_end):
                    core_cells.append(k * nx * ny + j * nx + i)
        entity_sets["cell_sets"] = [
            {"id": int(entity_config.get("core_cells_id", 10)), "set": core_cells}
        ]

    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(entity_sets, stream, indent=2)
        stream.write("\n")


def _write_axis_spacing(path: Path, axes: Sequence[AxisResult]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("axis,node_index,coordinate,cell_size\n")
        for axis_name, axis in zip(("x", "y", "z"), axes):
            for index, coordinate in enumerate(axis.coordinates):
                size = axis.widths[index] if index < len(axis.widths) else ""
                stream.write(
                    f"{axis_name},{index},{_format(coordinate)},{size}\n"
                )


def _write_center_ids(
    path: Path,
    center: tuple[float, float, float],
    node_id: int,
    node_coordinates: tuple[float, float, float],
    particle: tuple[int, tuple[float, float, float], float] | None,
    node_set_id: int,
    particle_set_id: int,
) -> None:
    node_distance = math.sqrt(
        sum((node_coordinates[index] - center[index]) ** 2 for index in range(3))
    )
    payload: dict[str, Any] = {
        "core_center": list(center),
        "closest_node": {
            "id": node_id,
            "coordinates": list(node_coordinates),
            "distance": node_distance,
            "node_set_id": node_set_id,
        },
        "closest_particle": None,
    }
    if particle is not None:
        particle_id, particle_coordinates, particle_distance = particle
        payload["closest_particle"] = {
            "id": particle_id,
            "coordinates": list(particle_coordinates),
            "distance": particle_distance,
            "particle_set_id": particle_set_id,
        }
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, indent=2)
        stream.write("\n")


def _resolve_output_files(
    output_dir: Path, output_config: dict[str, Any]
) -> OutputFiles:
    return OutputFiles(
        mesh=output_dir / str(output_config.get("mesh_file", "mesh-3d.txt")),
        particles=output_dir / str(output_config.get("particles_file", "particles-3d.txt")),
        particles_cells=output_dir
        / str(output_config.get("particles_cells_file", "particles-cells.txt")),
        particles_volumes=output_dir
        / str(output_config.get("particles_volumes_file", "particles-volumes.txt")),
        entity_sets=output_dir
        / str(output_config.get("entity_sets_file", "entity_sets.json")),
        mpm_fragment=output_dir
        / str(output_config.get("mpm_fragment_file", "mpm-input-fragment.json")),
        report=output_dir / str(output_config.get("report_file", "mesh-generation-report.json")),
        axis_spacing=output_dir
        / str(output_config.get("axis_spacing_file", "axis-spacing.csv")),
        center_ids=output_dir
        / str(output_config.get("center_ids_file", "core-center-ids.json")),
    )


def _all_output_paths(files: OutputFiles) -> list[Path]:
    return [
        files.mesh,
        files.particles,
        files.particles_cells,
        files.particles_volumes,
        files.entity_sets,
        files.mpm_fragment,
        files.report,
        files.axis_spacing,
        files.center_ids,
    ]


def _write_mpm_fragment(
    path: Path,
    files: OutputFiles,
    particle_config: dict[str, Any],
) -> None:
    fragment = {
        "mesh": {
            "mesh": files.mesh.name,
            "entity_sets": files.entity_sets.name,
            "io_type": "Ascii3D",
            "check_duplicates": False,
            "isoparametric": False,
            "node_type": "N3D",
            "cell_type": "ED3H8",
            "particle_cells": files.particles_cells.name,
            "particles_volumes": files.particles_volumes.name,
        },
        "particles": [
            {
                "generator": {
                    "type": "file",
                    "material_id": int(particle_config.get("material_id", 0)),
                    "pset_id": int(particle_config.get("pset_id", 0)),
                    "io_type": "Ascii3D",
                    "particle_type": str(particle_config.get("particle_type", "P3D")),
                    "check_duplicates": False,
                    "location": files.particles.name,
                }
            }
        ],
    }
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(fragment, stream, indent=2)
        stream.write("\n")


def generate(
    config: dict[str, Any],
    *,
    config_dir: Path,
    output_override: Path | None = None,
    force: bool = False,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Generate all output files and return the report dictionary."""

    domain = _require_mapping(config.get("domain"), "domain")
    core = _require_mapping(config.get("core"), "core")
    spacing = _require_mapping(config.get("spacing"), "spacing")

    domain_min = _triplet(domain.get("min"), "domain.min")
    domain_max = _triplet(domain.get("max"), "domain.max")
    core_min = _triplet(core.get("min"), "core.min")
    core_max = _triplet(core.get("max"), "core.max")
    particle_domain_config = config.get("particle_domain")
    if particle_domain_config is None:
        particle_domain_min = domain_min
        particle_domain_max = domain_max
    else:
        particle_domain = _require_mapping(
            particle_domain_config, "particle_domain"
        )
        particle_domain_min = _triplet(
            particle_domain.get("min"), "particle_domain.min"
        )
        particle_domain_max = _triplet(
            particle_domain.get("max"), "particle_domain.max"
        )
        for axis in range(3):
            if (
                particle_domain_min[axis] < domain_min[axis] - EPS
                or particle_domain_max[axis] > domain_max[axis] + EPS
                or particle_domain_min[axis] >= particle_domain_max[axis]
            ):
                raise ConfigError(
                    "particle_domain must be a non-empty box inside domain"
                )
    fine = _triplet(spacing.get("fine"), "spacing.fine", positive=True)
    coarse = _triplet(spacing.get("coarse"), "spacing.coarse", positive=True)
    transition_mode = str(spacing.get("transition", "linear")).lower()
    max_adjacent_ratio = _number(
        spacing.get("max_adjacent_ratio", 1.20),
        "spacing.max_adjacent_ratio",
    )
    if max_adjacent_ratio <= 1.0:
        raise ConfigError("spacing.max_adjacent_ratio must be greater than 1")

    particle_config = _require_mapping(config.get("particles", {}), "particles")
    per_cell = _int_triplet(
        particle_config.get("per_cell", [2, 2, 2]),
        "particles.per_cell",
    )
    remove_regions = [
        _require_mapping(shape, f"remove_regions[{index}]")
        for index, shape in enumerate(_as_list(config.get("remove_regions"), "remove_regions"))
    ]

    x_axis = _build_axis(
        domain_min[0],
        domain_max[0],
        core_min[0],
        core_max[0],
        fine[0],
        coarse[0],
        transition_mode,
        max_adjacent_ratio,
    )
    y_axis = _build_axis(
        domain_min[1],
        domain_max[1],
        core_min[1],
        core_max[1],
        fine[1],
        coarse[1],
        transition_mode,
        max_adjacent_ratio,
    )
    z_axis = _build_axis(
        domain_min[2],
        domain_max[2],
        core_min[2],
        core_max[2],
        fine[2],
        coarse[2],
        transition_mode,
        max_adjacent_ratio,
    )

    nnodes = (x_axis.cell_count + 1) * (y_axis.cell_count + 1) * (
        z_axis.cell_count + 1
    )
    ncells = x_axis.cell_count * y_axis.cell_count * z_axis.cell_count
    max_particle_slots = ncells * per_cell[0] * per_cell[1] * per_cell[2]

    output_config = _require_mapping(config.get("output", {}), "output")
    if output_override is not None:
        output_dir = output_override
    else:
        output_dir = Path(str(output_config.get("directory", "generated_graded_hex")))
        if not output_dir.is_absolute():
            output_dir = config_dir / output_dir
    output_dir = output_dir.resolve()
    files = _resolve_output_files(output_dir, output_config)

    min_size = min(min(x_axis.widths), min(y_axis.widths), min(z_axis.widths))
    time_step_config = _require_mapping(config.get("time_step", {}), "time_step")
    suggested_dt: float | None = None
    if "wave_speed" in time_step_config:
        wave_speed = _number(time_step_config["wave_speed"], "time_step.wave_speed")
        max_velocity = _number(
            time_step_config.get("max_velocity", 0.0),
            "time_step.max_velocity",
        )
        cfl = _number(time_step_config.get("cfl", 0.5), "time_step.cfl")
        if wave_speed + max_velocity <= 0.0:
            raise ConfigError("time_step.wave_speed + max_velocity must be positive")
        suggested_dt = cfl * min_size / (wave_speed + max_velocity)

    report: dict[str, Any] = {
        "config": {
            "domain": {"min": domain_min, "max": domain_max},
            "core": {"min": core_min, "max": core_max},
            "particle_domain": {
                "min": particle_domain_min,
                "max": particle_domain_max,
            },
            "spacing": {
                "requested_fine": fine,
                "requested_coarse": coarse,
                "transition": transition_mode,
                "max_adjacent_ratio": max_adjacent_ratio,
            },
            "particles_per_cell": per_cell,
        },
        "mesh": {
            "nnodes": nnodes,
            "ncells": ncells,
            "cell_type": "ED3H8",
            "isoparametric": False,
            "min_cell_size": min_size,
            "axes": {
                "x": {
                    "cell_count": x_axis.cell_count,
                    "actual_fine_size": x_axis.actual_fine,
                    "min_size": min(x_axis.widths),
                    "max_size": max(x_axis.widths),
                    "max_adjacent_ratio": x_axis.max_adjacent_ratio,
                },
                "y": {
                    "cell_count": y_axis.cell_count,
                    "actual_fine_size": y_axis.actual_fine,
                    "min_size": min(y_axis.widths),
                    "max_size": max(y_axis.widths),
                    "max_adjacent_ratio": y_axis.max_adjacent_ratio,
                },
                "z": {
                    "cell_count": z_axis.cell_count,
                    "actual_fine_size": z_axis.actual_fine,
                    "min_size": min(z_axis.widths),
                    "max_size": max(z_axis.widths),
                    "max_adjacent_ratio": z_axis.max_adjacent_ratio,
                },
            },
        },
        "particles": {
            "max_slots": max_particle_slots,
        },
        "time_step": {
            "suggested_dt": suggested_dt,
        },
        "output_directory": str(output_dir),
    }

    if dry_run:
        return report

    output_dir.mkdir(parents=True, exist_ok=True)
    if not force:
        existing = [path for path in _all_output_paths(files) if path.exists()]
        if existing:
            names = ", ".join(path.name for path in existing)
            raise ConfigError(
                f"output files already exist ({names}); pass --force to overwrite"
            )

    for temporary in output_dir.glob(".*.tmp"):
        temporary.unlink(missing_ok=True)

    _, actual_ncells = _write_mesh(files.mesh, x_axis, y_axis, z_axis)
    if actual_ncells != ncells:
        raise RuntimeError("internal error: mesh cell count mismatch")

    entity_config = _require_mapping(config.get("entity_sets", {}), "entity_sets")
    include_all_particle_set = bool(entity_config.get("include_all_particle_set", False))
    include_core_particle_set = bool(entity_config.get("include_core_particle_set", False))
    (
        particle_count,
        removed_count,
        outside_particle_domain_count,
        active_cell_count,
        all_particle_ids,
        core_particle_ids,
        closest_particle,
    ) = _generate_particles(
        output_dir,
        x_axis,
        y_axis,
        z_axis,
        per_cell,
        particle_domain_min,
        particle_domain_max,
        remove_regions,
        files,
        include_all_particle_set,
        include_core_particle_set,
    )

    core_center = (
        0.5
        * (
            x_axis.coordinates[x_axis.core_start]
            + x_axis.coordinates[x_axis.core_end]
        ),
        0.5
        * (
            y_axis.coordinates[y_axis.core_start]
            + y_axis.coordinates[y_axis.core_end]
        ),
        0.5
        * (
            z_axis.coordinates[z_axis.core_start]
            + z_axis.coordinates[z_axis.core_end]
        ),
    )
    center_node_id, center_node_coordinates = _closest_node(
        x_axis, y_axis, z_axis, core_center
    )
    center_particle_id = closest_particle[0] if closest_particle is not None else None
    center_node_set_id = int(entity_config.get("center_node_set_id", 100))
    center_particle_set_id = int(entity_config.get("center_particle_set_id", 100))

    _write_entity_sets(
        files.entity_sets,
        x_axis,
        y_axis,
        z_axis,
        entity_config,
        all_particle_ids,
        core_particle_ids,
        center_node_id,
        center_particle_id,
    )
    _write_center_ids(
        files.center_ids,
        core_center,
        center_node_id,
        center_node_coordinates,
        closest_particle,
        center_node_set_id,
        center_particle_set_id,
    )
    _write_axis_spacing(files.axis_spacing, (x_axis, y_axis, z_axis))
    _write_mpm_fragment(files.mpm_fragment, files, particle_config)

    report["particles"].update(
        {
            "written": particle_count,
            "removed_by_boolean": removed_count,
            "removed_outside_particle_domain": outside_particle_domain_count,
            "active_cells": active_cell_count,
            "empty_cells": ncells - active_cell_count,
        }
    )
    report["center_selection"] = {
        "core_center": list(core_center),
        "closest_node": {
            "id": center_node_id,
            "coordinates": list(center_node_coordinates),
            "node_set_id": center_node_set_id,
        },
        "closest_particle": (
            {
                "id": closest_particle[0],
                "coordinates": list(closest_particle[1]),
                "distance": closest_particle[2],
                "particle_set_id": center_particle_set_id,
            }
            if closest_particle is not None
            else None
        ),
    }
    report["files"] = {
        "mesh": files.mesh.name,
        "particles": files.particles.name,
        "particles_cells": files.particles_cells.name,
        "particles_volumes": files.particles_volumes.name,
        "entity_sets": files.entity_sets.name,
        "mpm_fragment": files.mpm_fragment.name,
        "axis_spacing": files.axis_spacing.name,
        "center_ids": files.center_ids.name,
    }

    with files.report.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    return report


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a graded orthorhombic Hex8 mesh and particle inputs for "
            "CB-Geo MPM."
        )
    )
    parser.add_argument("config", type=Path, help="JSON generator configuration")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="override output.directory from the configuration",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite generated files if they already exist",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate the configuration and print the estimated mesh counts",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    try:
        with config_path.open("r", encoding="utf-8") as stream:
            config = json.load(stream)
        report = generate(
            config,
            config_dir=config_path.parent,
            output_override=args.output_dir,
            force=args.force,
            dry_run=args.dry_run,
        )
    except (OSError, json.JSONDecodeError, ConfigError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if args.dry_run:
        print(json.dumps(report, indent=2))
    else:
        print(
            "Generated "
            f"{report['mesh']['nnodes']} nodes, "
            f"{report['mesh']['ncells']} Hex8 cells, "
            f"{report['particles']['written']} particles "
            f"in {report['output_directory']}"
        )
        print(f"Suggested dt: {report['time_step']['suggested_dt']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
