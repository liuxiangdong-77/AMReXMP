"""Configure the AMR standard case in ParaView.

Run with ParaView's Python shell or pvpython:

    pvpython setup_paraview_amr_standard.py --data-dir DATA_DIRECTORY
"""

import argparse
import glob
import os
import sys

from paraview.simple import (
    Clip,
    CreateView,
    GetActiveViewOrCreate,
    Hide,
    Render,
    SaveState,
    Show,
    XMLPolyDataReader,
    XMLUnstructuredGridReader,
)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        default=os.environ.get("MPM_PARAVIEW_DATA_DIR", ""),
        help="Directory containing grid_data*.vtu and particle_data*.vtp",
    )
    parser.add_argument(
        "--state",
        default="",
        help="Optional output ParaView state file (.pvsm)",
    )
    parser.add_argument("--plane-x", type=float, default=5.0)
    parser.add_argument("--particle-radius", type=float, default=0.12)
    return parser.parse_args()


def configure(data_dir, state="", plane_x=5.0, particle_radius=0.12):
    data_dir = os.path.abspath(data_dir)
    if not os.path.isdir(data_dir):
        raise RuntimeError("Data directory does not exist: " + data_dir)

    particle_files = sorted(
        glob.glob(os.path.join(data_dir, "particle_data*.vtp"))
    )
    grid_files = sorted(glob.glob(os.path.join(data_dir, "grid_data*.vtu")))
    if not particle_files:
        raise RuntimeError("No particle_data*.vtp files were found")
    if not grid_files:
        raise RuntimeError("No grid_data*.vtu files were found")

    view = GetActiveViewOrCreate("RenderView")
    if view is None:
        view = CreateView("RenderView")
    view.ViewSize = [1280, 800]
    view.OrientationAxesVisibility = 1

    grid = XMLUnstructuredGridReader(FileName=grid_files[0])
    grid_display = Show(grid, view)
    grid_display.Representation = "Wireframe"
    grid_display.ColorArrayName = ["CELLS", "amr_level"]
    grid_display.LineWidth = 1.0

    particles = XMLPolyDataReader(FileName=particle_files)
    particles_display = Show(particles, view)
    particles_display.Representation = "Point Gaussian"
    particles_display.GaussianRadius = particle_radius
    particles_display.ScaleArray = ["POINTS", "radius_factor"]
    particles_display.ScaleFactor = 1.0

    particle_clip = Clip(Input=particles)
    particle_clip.ClipType = "Plane"
    particle_clip.ClipType.Origin = [plane_x, 5.0, 5.0]
    particle_clip.ClipType.Normal = [1.0, 0.0, 0.0]
    particle_clip.Invert = 0
    Hide(particles, view)
    clipped_display = Show(particle_clip, view)
    clipped_display.Representation = "Point Gaussian"
    clipped_display.GaussianRadius = particle_radius
    clipped_display.ScaleArray = ["POINTS", "radius_factor"]
    clipped_display.ScaleFactor = 1.0

    grid_clip = Clip(Input=grid)
    grid_clip.ClipType = "Plane"
    grid_clip.ClipType.Origin = [plane_x, 5.0, 5.0]
    grid_clip.ClipType.Normal = [1.0, 0.0, 0.0]
    grid_clip.Invert = 0
    Hide(grid, view)
    clipped_grid_display = Show(grid_clip, view)
    clipped_grid_display.Representation = "Wireframe"
    clipped_grid_display.ColorArrayName = ["CELLS", "amr_level"]
    clipped_grid_display.LineWidth = 1.0

    view.CameraPosition = [-22.0, 5.0, 5.0]
    view.CameraFocalPoint = [5.0, 5.0, 5.0]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelProjection = 1
    Render(view)

    if state:
        SaveState(state, view)


if __name__ == "__main__" and "--data-dir" in sys.argv:
    arguments = parse_arguments()
    configure(
        arguments.data_dir,
        arguments.state,
        arguments.plane_x,
        arguments.particle_radius,
    )
