import csv
import gsd.fl
import gsd.hoomd
import os
import signac
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import argparse
import itertools
import operator
from scipy.signal import argrelextrema
from scipy.ndimage import gaussian_filter
from rhaco.simulate import AVOGADRO, BOLTZMANN, KCAL_TO_J, AMU_TO_KG, ANG_TO_M


"""
This module plots the residence time distributions for each job in the workspace.
Residency is defined as the time that a particular reactant molecule is located within
a 1nm distance of the surface crystal (configurable by the --tolerance argument).
"""


def get_job_frame(job):
    gsd_file_name = os.path.join(job.ws, "output_traj.gsd")
    try:
        gsd_file = gsd.fl.GSDFile(gsd_file_name, "rb")
    except OSError:
        raise SystemError(" ".join([gsd_file_name, "not found. Skipping..."]))
    trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
    initial_frame = trajectory[0]
    return initial_frame


def get_surface_atoms(job, morphology):
    # In the first frame, we can use the definitions of the reactant boxes to find
    # all of the non-reactant atoms, giving us the crystal atoms
    # In rhaco/generate.py the top box is defined as:
    # box_top = mb.Box(mins=[-(args.crystal_x * args.dimensions[0]) / 2.0,
    #                        -(args.crystal_y * args.dimensions[1]) / 2.0,
    #                        args.crystal_separation / 20.0
    #                        + (args.crystal_z * args.dimensions[2])],
    #                  maxs=[(args.crystal_x * args.dimensions[0]) / 2.0,
    #                        (args.crystal_y * args.dimensions[1]) / 2.0,
    #                        args.z_reactor_size / 2.0])
    # The surface we want is everything between z = 0 and the bottom of the box_top
    z_min = 0.0
    z_max = float(job.sp["crystal_separation"]) / 2.0 + (
        float(job.sp["crystal_z"]) * 10.0 * float(job.sp["dimensions"].split("x")[2])
    )
    z_vals = morphology.particles.position[:,2]
    return morphology.particles.position[(z_vals > z_min) * (z_vals < z_max)]


def create_mesh(morphology, lattice_spacing, z_position):
    box_dims = morphology.configuration.box[:3]
    # Work out all of the x, y, and z values
    x_vals = np.arange(-box_dims[0]/2.0, box_dims[0]/2.0, lattice_spacing)
    y_vals = np.arange(-box_dims[1]/2.0, box_dims[1]/2.0, lattice_spacing)
    z_vals = np.array([z_position])
    # Create the mesh coordinates using vstack and meshgrid
    mesh_coords = np.vstack(np.meshgrid(x_vals, y_vals, z_vals)).reshape(3, -1).T
    # Now centre the mesh over the middle
    centroid = np.array(
        [
            sum(mesh_coords[:,axis]) / len(mesh_coords[:,axis]) for axis in range(2)
        ]
        + [0.0]
    )
    mesh_coords -= centroid
    centroid = np.array(
        [
            sum(mesh_coords[:,axis]) / len(mesh_coords[:,axis]) for axis in range(2)
        ]
        + [0.0]
    )
    return mesh_coords


def flatten(input_list):
    return [item for sublist in input_list for item in sublist]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-j",
        "--job",
        type=str,
        required=False,
        default=None,
        help=(
            "If present, only consider the job in the current directory's workspace"
        ),
    )
    parser.add_argument(
        "-i",
        "--interactive",
        required=False,
        action="store_true",
        help=(
            "Open the heatmap in interactive mode"
        ),
    )
    parser.add_argument(
        "-z",
        "--z_val",
        type=float,
        required=False,
        default=20.0,
        help=(
            "The z value to use when interactive mode is not activated"
        ),
    )
    parser.add_argument(
        "-m",
        "--mesh",
        type=float,
        required=False,
        default=10.0,
        help=(
            "The spacing between mesh probes (grid resolution) in angstroems"
        ),
    )
    args, directory_list = parser.parse_known_args()
    project = signac.get_project("../")
    schema = project.detect_schema()
    plt.figure()
    # Find the crystal extents and append them to the job documents
    for job in project:
        if (args.job is not None) and (job.get_id() != args.job):
            continue
        job_frame = get_job_frame(job)
        crystal_posns = get_surface_atoms(job, job_frame)
        mesh_posns = create_mesh(job_frame, args.mesh, args.z_val)
    plt.close()
