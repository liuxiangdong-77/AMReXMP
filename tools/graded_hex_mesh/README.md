# Graded Hex8 Mesh Generator

This tool generates an initialization-only, tensor-product, graded orthorhombic
Hex8 mesh for CB-Geo MPM. It also generates particles, particle-to-cell mapping,
particle volumes, boundary node sets, and a reusable MPM input fragment.

The generator uses only the Python standard library.

See `USAGE.zh-CN.md` for the complete Chinese workflow and configuration guide.

## Quick Start

```bash
python3 generate_graded_hex_mesh.py graded_hex.example.json --dry-run
python3 generate_graded_hex_mesh.py graded_hex.example.json
python3 validate_graded_hex_output.py generated_graded_hex
```

Use `--force` to overwrite previously generated files.

## Configuration

The domain and core are axis-aligned boxes:

```json
{
  "domain": {"min": [0, 0, 0], "max": [30000, 30000, 35000]},
  "core": {"min": [0, 0, 25000], "max": [1000, 1000, 28000]}
}
```

The core uses `spacing.fine`. Outside the core, spacing grows to
`spacing.coarse` independently along x, y, and z.

The optional `particle_domain` defines where particles are generated. The Hex8
mesh always spans the full `domain`, so cells outside `particle_domain` remain
present but empty:

```json
"particle_domain": {
  "min": [1000, 1000, 1000],
  "max": [29000, 29000, 34000]
}
```

Empty cells let particles move into a surrounding buffer region without leaving
the computational mesh. The mesh must still cover the complete expected motion;
empty cells do not prevent errors beyond the outer mesh boundary.

`spacing.transition` supports:

- `linear`: uses `N = round(2L/(hf+hc))` and
  `h_i = hf + (hc-hf) * i/(N-1)`, followed by a uniform scaling that preserves
  the transition length.
- `smoothstep`: applies a smoothstep weight to the same discrete cell sequence.

There is no separate transition subregion: the full distance from the core
boundary to the domain boundary participates in the grading.

`spacing.max_adjacent_ratio` validates the largest ratio between adjacent cell
sizes. If the specified transition distance cannot satisfy both the requested
endpoint sizes and this limit, generation stops with an actionable error.

The small example uses a relaxed ratio because its transition region is short.
For wave-propagation studies, start with a ratio near 1.1 to 1.3 and provide
enough transition length for the requested fine-to-coarse size range.

## Particle Generation

`particles.per_cell` defines a structured lattice inside every cell. The default
`[2, 2, 2]` produces approximately two particles per cell in each direction.

The particle volume is:

```text
cell_volume / total_particle_slots_per_cell
```

If a Boolean cutter removes a particle, the remaining particles keep this
volume, so removed material reduces total particle volume.

## Boolean Removal

`remove_regions` supports:

- `box`: `min` and `max`
- `sphere`: `center` and `radius`
- `cylinder`: `center`, `radius`, `axis`, and either `height` or axis `min`/`max`

The generated background mesh remains a conforming tensor-product Hex8 grid.
Boolean geometry is represented by removing particles, not by cutting mesh
cells.

Boolean removal is disabled when `remove_regions` is omitted or is an empty
array. The default example has it disabled.

## Core-Center Selection

The generator always writes `core-center-ids.json` with:

- the exact core-center coordinates;
- the nearest node ID and coordinates;
- the nearest retained particle ID and coordinates.

The nearest node is also written as a one-node entity set, and the nearest
particle is written as a one-particle entity set. Their IDs can be configured:

```json
"entity_sets": {
  "center_node_set_id": 100,
  "center_particle_set_id": 100
}
```

This makes the center point directly usable for nodal or particle loading.

## Generated Files

- `mesh-3d.txt`: CB-Geo MPM ASCII Hex8 mesh.
- `particles-3d.txt`: particle coordinates.
- `particles-cells.txt`: particle ID and cell ID.
- `particles-volumes.txt`: particle ID and volume.
- `entity_sets.json`: boundary node sets and optional particle/cell sets.
- `mpm-input-fragment.json`: mesh and particle configuration for `mpm.json`.
- `mesh-generation-report.json`: counts, spacing, ratios, and suggested dt.
- `axis-spacing.csv`: coordinates and sizes along each axis.
- `core-center-ids.json`: nearest core-center node and particle IDs.

## MPM Notes

Use the mesh and particle blocks from `mpm-input-fragment.json` in the MPM input.
The generated mesh is intended for:

```json
"cell_type": "ED3H8",
"node_type": "N3D",
"io_type": "Ascii3D",
"isoparametric": false
```

MPM does not automatically derive a stable time step from the smallest cell.
Use `analysis.dt` no larger than the value in
`mesh-generation-report.json`.

To write the normalized graded-mesh particle radius to VTK, include
`radius_factor` in `post_processing.vtk`:

```json
"post_processing": {
  "vtk": ["radius_factor"]
}
```

MPM computes it as the local cell equivalent size divided by the largest cell
equivalent size in the model. The coarsest particles have value 1 and finer
particles have values below 1.

A combined file `particle_data<step>.vtp` is also written. It contains the
selected stress/velocity/displacement fields plus `radius_factor` and the
three-component `particle_spacing` vector.

For anisotropic cells, use ParaView `3D Glyphs` with a `Box` source, select
`Scale by Vector Components`, and use `particle_spacing` as the scale array.
This fills each orthorhombic cell uniformly without the single-scalar
limitation of Point Gaussian.

The tensor-product grid cannot create an isolated fine hexahedral box while
keeping all cells axis-aligned and conforming. Refinement propagates as
transition bands along x, y, and z. This is intentional and avoids hanging
nodes.
