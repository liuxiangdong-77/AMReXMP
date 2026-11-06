#!/usr/bin/env python3

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("generate_graded_hex_mesh.py")
SPEC = importlib.util.spec_from_file_location("graded_hex_generator", MODULE_PATH)
GENERATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = GENERATOR
SPEC.loader.exec_module(GENERATOR)


class GradedHexGeneratorTest(unittest.TestCase):
    def test_axis_is_monotonic_and_smooth(self) -> None:
        axis = GENERATOR._build_axis(
            0.0,
            8.0,
            2.0,
            6.0,
            0.25,
            1.0,
            "smoothstep",
            3.0,
        )
        self.assertEqual(axis.coordinates[0], 0.0)
        self.assertEqual(axis.coordinates[-1], 8.0)
        self.assertTrue(all(width > 0.0 for width in axis.widths))
        self.assertLessEqual(axis.max_adjacent_ratio, 3.0 + 1.0e-6)

    def test_generated_mesh_and_particle_counts(self) -> None:
        config = {
            "domain": {"min": [0.0, 0.0, 0.0], "max": [2.0, 2.0, 2.0]},
            "core": {"min": [0.5, 0.5, 0.5], "max": [1.5, 1.5, 1.5]},
            "particle_domain": {
                "min": [0.5, 0.5, 0.5],
                "max": [1.5, 1.5, 1.5],
            },
            "spacing": {
                "fine": [0.25, 0.25, 0.25],
                "coarse": [1.0, 1.0, 1.0],
                "transition": "linear",
                "max_adjacent_ratio": 3.0,
            },
            "particles": {"per_cell": [2, 2, 2]},
            "remove_regions": [
                {
                    "type": "sphere",
                    "center": [1.0, 1.0, 1.0],
                    "radius": 0.2,
                }
            ],
            "output": {"directory": "generated"},
        }

        with tempfile.TemporaryDirectory() as temporary:
            report = GENERATOR.generate(
                config,
                config_dir=Path(temporary),
                force=True,
            )
            output_dir = Path(report["output_directory"])
            mesh_path = output_dir / "mesh-3d.txt"
            particles_path = output_dir / "particles-3d.txt"

            with mesh_path.open("r", encoding="utf-8") as stream:
                lines = [
                    line.strip()
                    for line in stream
                    if line.strip() and not line.lstrip().startswith("!")
                ]
            nnodes, ncells = (int(value) for value in lines[0].split())
            self.assertEqual(len(lines), 1 + nnodes + ncells)

            with particles_path.open("r", encoding="utf-8") as stream:
                nparticles = int(stream.readline().strip())
                coordinates = [
                    tuple(float(value) for value in line.split())
                    for line in stream
                    if line.strip()
                ]
            self.assertEqual(nparticles, len(coordinates))
            self.assertGreater(nparticles, 0)
            self.assertLess(nparticles, report["particles"]["max_slots"])
            self.assertGreater(report["particles"]["empty_cells"], 0)
            self.assertTrue(
                all(
                    (x - 1.0) ** 2 + (y - 1.0) ** 2 + (z - 1.0) ** 2
                    > 0.2**2 + 1.0e-10
                    for x, y, z in coordinates
                )
            )

            with (output_dir / "entity_sets.json").open(
                "r", encoding="utf-8"
            ) as stream:
                entity_sets = json.load(stream)
            self.assertEqual(len(entity_sets["node_sets"]), 7)
            self.assertFalse((output_dir / "cell-levels.txt").exists())

            with (output_dir / "core-center-ids.json").open(
                "r", encoding="utf-8"
            ) as stream:
                center_ids = json.load(stream)
            self.assertIn("closest_node", center_ids)
            self.assertIn("closest_particle", center_ids)
            self.assertIsNotNone(center_ids["closest_particle"])

            validator_path = MODULE_PATH.with_name("validate_graded_hex_output.py")
            validator_spec = importlib.util.spec_from_file_location(
                "graded_hex_validator", validator_path
            )
            validator = importlib.util.module_from_spec(validator_spec)
            assert validator_spec.loader is not None
            sys.modules[validator_spec.name] = validator
            validator_spec.loader.exec_module(validator)
            summary = validator.validate(output_dir, 1.0e-8)
            self.assertEqual(summary["nodes"], report["mesh"]["nnodes"])
            self.assertEqual(summary["cells"], report["mesh"]["ncells"])
            self.assertEqual(summary["particles"], nparticles)


if __name__ == "__main__":
    unittest.main()
