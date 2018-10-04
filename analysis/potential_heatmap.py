import gsd.fl
import gsd.hoomd
import os
import signac
import matplotlib.colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import argparse
import freud.locality
import freud.box
import rhaco.simulate
from matplotlib.widgets import Slider

import time as T

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
    surface_posns = morphology.particles.position[(z_vals > z_min) * (z_vals < z_max)]
    surface_type_ID = morphology.particles.typeid[(z_vals > z_min) * (z_vals < z_max)]
    type_lookup = {ID: type_name for ID, type_name in enumerate(morphology.particles.types)}
    surface_type = np.vectorize(type_lookup.get)(surface_type_ID)
    return surface_posns, surface_type


def create_mesh(morphology, z_range, args):
    lattice_spacing = args.mesh_spacing
    mesh_type = args.mesh_type
    box_dims = morphology.configuration.box[:3]
    # Work out all of the x, y, and z values
    x_vals = np.arange(-box_dims[0]/2.0, box_dims[0]/2.0, lattice_spacing)
    y_vals = np.arange(-box_dims[1]/2.0, box_dims[1]/2.0, lattice_spacing)
    z_vals = z_range
    mesh_shape = (len(x_vals), len(y_vals), len(z_vals))
    print("Creating mesh of shape", mesh_shape, "...")
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
    mesh_types = [mesh_type] * len(mesh_coords)
    return mesh_coords, mesh_types, mesh_shape


def create_freud_nlist(job_frame, crystal_mesh_posns, mesh_size, cut_off):
    # Firstly create the simulation box so that Freud can do the periodic stuff
    simulation_box = freud.box.Box(*job_frame.configuration.box)
    # Get the set of the mesh_IDs for later
    mesh_IDs = set(np.arange(mesh_size))
    # Calculate the cell list based on the input cut_off (same as used in the sim)
    cell_list = freud.locality.LinkCell(simulation_box, cut_off)
    # Compute the neighbourlist for all particles in the crystal and mesh
    cell_list.compute(simulation_box, crystal_mesh_posns)
    neighbour_list = cell_list.nlist
    # Create a neighbourlist dictionary of crystal_IDs for each probe atom in the mesh
    nlist = {}
    for probe_ID in range(mesh_size):
        # Set manipulation
        neighbour_IDs = set(neighbour_list.index_j[neighbour_list.index_i == probe_ID])
        nlist[probe_ID] = list(neighbour_IDs - mesh_IDs)

        # # Original
        # crystal_IDs = neighbour_list.index_j[neighbour_list.index_i == probe_ID]
        # nlist[probe_ID] = [int(ID) for ID in crystal_IDs if ID >= mesh_size]

        # # Array masking
        # neighbour_IDs = neighbour_list.index_j[neighbour_list.index_i == probe_ID]
        # nlist[probe_ID] = neighbour_IDs[neighbour_IDs >= mesh_size]
    return nlist


def calculate_potentials(
    job, nlist, n_probes, crystal_mesh_posns, crystal_mesh_types, u_max
):
    # Vectorised version
    ff_coeffs = rhaco.simulate.get_coeffs(os.path.join(job.ws, "output.hoomdxml"))
    pair_coeffs = {coeff[0]: [coeff[1], coeff[2]] for coeff in ff_coeffs["pair_coeffs"]}
    if len(ff_coeffs["external_forcefields"]) > 0:
        raise SystemError("EXTERNAL FORCEFIELD DETECTED (E.G. EAM), CANNOT INTERPRET")
    probe_atoms = {}
    for probe_index in range(n_probes):
        probe_atoms[probe_index] = [crystal_mesh_posns[probe_index]]
    for probe_ID, neighbours in nlist.items():
        type_1 = crystal_mesh_types[probe_ID]
        pos_1 = crystal_mesh_posns[probe_ID]
        types_2 = crystal_mesh_types[neighbours]
        posns_2 = crystal_mesh_posns[neighbours]
        potential = np.sum(LJ_pair_potential_vec(pair_coeffs, type_1, types_2, pos_1, posns_2))
        if u_max is not None:
            potential = np.clip(potential, None, u_max)
        probe_atoms[probe_ID].append(np.sum(potential))

    # # For loop version
    # ff_coeffs = rhaco.simulate.get_coeffs(os.path.join(job.ws, "output.hoomdxml"))
    # pair_coeffs = {coeff[0]: [coeff[1], coeff[2]] for coeff in ff_coeffs["pair_coeffs"]}
    # if len(ff_coeffs["external_forcefields"]) > 0:
    #     raise SystemError("EXTERNAL FORCEFIELD DETECTED (E.G. EAM), CANNOT INTERPRET")
    # probe_atoms = {}
    # for probe_index in range(n_probes):
    #     probe_atoms[probe_index] = [crystal_mesh_posns[probe_index]]
    # # NOTE If this is slow, then consider vectorizing it by passing in all neighbours
    # # at once
    # for probe_ID, neighbours in nlist.items():
    #     potential = 0.0
    #     potential_values = []
    #     type_1 = crystal_mesh_types[probe_ID]
    #     pos_1 = crystal_mesh_posns[probe_ID]
    #     for neighbour_ID in neighbours:
    #         type_2 = crystal_mesh_types[neighbour_ID]
    #         pos_2 = crystal_mesh_posns[neighbour_ID]
    #         potential_increment = LJ_pair_potential(pair_coeffs, type_1, type_2, pos_1, pos_2)
    #         potential_values.append(potential_increment)
    #         potential += potential_increment
    #     # If the potential is more than about 5 then there's no chance of the particle
    #     # being there, so set this as a cap
    #     if (u_max is not None) and (potential > u_max):
    #         potential = u_max
    #     probe_atoms[probe_ID].append(potential)
    return probe_atoms


def LJ_pair_potential(coeffs, type1, type2, type1_posn, type2_posn):
    # Using the geometric mixing rules
    epsilon = np.sqrt(coeffs[type1][0] * coeffs[type2][0])
    sigma = np.sqrt(coeffs[type1][1] * coeffs[type2][1])
    r = np.sqrt(np.sum((np.array(type1_posn) - np.array(type2_posn)) ** 2))
    # Lennard-Jones 12-6 formula
    potential = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
    return potential


def LJ_pair_potential_vec(coeffs, type1, type2, type1_posn, type2_posn):
    # Create a vectorized dictionary lookup
    v_lookup = np.vectorize(lambda dic, val, ix: dic[val][ix])
    # Obtain the epsilon, sigma and separation matrices
    epsilon_atom2 = v_lookup(coeffs, type2, 0)
    sigma_atom2 = v_lookup(coeffs, type2, 1)
    r = np.sqrt(np.sum((np.array(type1_posn) - np.array(type2_posn)) ** 2, axis=1))
    # Using the geometric mixing rules
    epsilon = np.sqrt(coeffs[type1][0] * epsilon_atom2)
    sigma = np.sqrt(coeffs[type1][1] * sigma_atom2)
    # Lennard-Jones 12-6 formula
    potential = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
    return potential


def create_potential_array(potential_dict, mesh_shape, args):
    potential_array = np.zeros(mesh_shape)
    lattice_spacing = args.mesh_spacing
    # Obtain the mesh_positions
    mesh_details = np.array([np.append(val[0], [val[1]]) for val in potential_dict.values()])
    # We're not going to need the z values so we can drop those here.
    mesh_posns = mesh_details[:,[0,1]]
    mesh_potentials = mesh_details[:,3]
    # Since we are no longer operating in coordinate space, shift the positions
    # so that the bottom left (-large_x, -large_y) is now the origin
    offset = np.array([np.min(mesh_posns[:,axis]) for axis in range(2)])
    mesh_posns -= offset
    # Divide by the lattice spacing
    mesh_posns /= lattice_spacing
    # Drop the z_coordinate as we don't need it, and remap the matrix as integers
    mesh_posns = np.array(mesh_posns[:,[0,1]], dtype=int)
    for probe_ID, probe_posn in enumerate(mesh_posns):
        potential_array[probe_posn[0], probe_posn[1]] = mesh_potentials[probe_ID]
    # As our origin is currently the bottom left (-large_x, -large_y), but array
    # plotting subroutines have the origin in the top left, we need to flip the array
    potential_array = np.flip(potential_array, 0)
    return potential_array


def get_z_range(job, z_step):
    crystal_offset = float(job.sp["crystal_separation"]) / 2.0
    crystal_z = float(job.sp["crystal_z"]) * 10.0 # convert from nm to ang
    dim_z = int(job.sp["dimensions"].split("x")[2])
    # Calculated as the job.document["crystal_top_layer"]
    z_min = crystal_offset + ((dim_z - 2) * crystal_z / 2.0)
    # Calculated as the job.document["crystal_top_edge"]
    top_edge = crystal_offset + (dim_z * crystal_z / 2.0)
    # FF pair cutoff is 10.0 (NOTE: Hardcoded into rhaco, if soft-coded then updat3
    # this)
    z_max = top_edge + 10.0
    return np.arange(z_min, z_max, z_step)


def plot_heatmap(input_array, z_range, args):
    fig, ax = plt.subplots()
    # Reverse the colour map so that the `hot spots' show places where particles are
    # more likely to reside (i.e. lower potential energy)
    colour_map = plt.get_cmap("inferno_r")
    c_norm = matplotlib.colors.Normalize(
        vmin=np.min(input_array),
        vmax=np.max(input_array),
    )
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=colour_map)
    heatmap = plt.imshow(input_array[0], cmap=colour_map, interpolation='nearest')
    # Create the colour bar
    scalar_map.set_array(input_array.flatten())
    cbar = plt.colorbar(scalar_map, aspect=20)
    ax_zval = plt.axes([0.2, 0.03, 0.4, 0.03], axisbg='black')
    z_slider = create_slider(ax_zval, z_range, input_array, heatmap, args)
    plt.show()


class create_slider(object):
    def __init__(self, axis, z_range, input_array, heatmap, args):
        z_min = np.min(z_range)
        z_max = np.max(z_range)
        self.z_range = z_range
        self.input_array = input_array
        self.heatmap = heatmap
        self.slider = Slider(axis, r'$z$', z_min, z_max, valinit=z_min)
        self.slider.on_changed(self.update)

    def update(self, val):
        slice_index, discrete_val = find_nearest(self.z_range, val)
        new_slice = self.input_array[slice_index]
        self.heatmap.set_data(new_slice)


def find_nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index, array[index]



def update(zval):
    global input_array
    z_val = s_zval.val
    z_slice = input_array[zval]
    img_plt.set_data(z_slice)




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
    # parser.add_argument(
    #     "-i",
    #     "--interactive",
    #     required=False,
    #     action="store_true",
    #     help=(
    #         "Open the heatmap in interactive mode"
    #     ),
    # )
    parser.add_argument(
        "-z",
        "--z_step",
        type=float,
        required=False,
        default=0.5,
        help=(
            "The z step to use (in ang)"
        ),
    )
    parser.add_argument(
        "-ms",
        "--mesh_spacing",
        type=float,
        required=False,
        default=10.0,
        help=(
            "The spacing between mesh probes (grid resolution) in angstroems"
        ),
    )
    parser.add_argument(
        "-mt",
        "--mesh_type",
        type=str,
        required=False,
        default="C",
        help=(
            "The atom type to use for the probe atoms in the mesh. Required for"
            " calculating the potential energy at each probe site."
        ),
    )
    parser.add_argument(
        "-r",
        "--r_cut",
        type=float,
        required=False,
        default=10.0,
        help=(
            "The pair interaction cut-off in angstroems"
        ),
    )
    parser.add_argument(
        "-u",
        "--u_max",
        type=float,
        required=False,
        default=None,
        help=(
            "The maximum potential to consider, in energy units. If the probe"
            " potential is more than this value, it will be truncated to this value."
            " The force is expected to be sufficiently high at and beyond this value"
            " that a particle could not reside there for any length of time."
        ),
    )
    args, directory_list = parser.parse_known_args()
    project = signac.get_project("../")
    schema = project.detect_schema()
    potential_array_3d = []
    # Find the crystal extents and append them to the job documents
    for job in project:
        if (args.job is not None) and (job.get_id() != args.job):
            continue
        z_range = get_z_range(job, args.z_step)
        job_frame = get_job_frame(job)
        crystal_posns, crystal_types = get_surface_atoms(job, job_frame)
        mesh_posns, mesh_types, mesh_shape = create_mesh(job_frame, z_range, args)
        crystal_mesh_posns = np.vstack([mesh_posns, crystal_posns])
        crystal_mesh_types = np.hstack([mesh_types, crystal_types])
        print("Calculating the LinkedCell neighbourlist...")
        nlist = create_freud_nlist(
            job_frame,
            crystal_mesh_posns,
            np.prod(mesh_shape),
            args.r_cut
        )
        for z_index, z_val in enumerate(z_range):
            print("\rCalculating potentials for z_slice {:d} of {:d}".format(z_index, len(z_range)), end="")
            potential_dict = calculate_potentials(
                job, nlist, np.prod(mesh_shape), crystal_mesh_posns, crystal_mesh_types,
                args.u_max,
            )
            potential_array = create_potential_array(potential_dict, mesh_shape, args)
            potential_array_3d.append(potential_array)
        exit()
        plot_heatmap(np.array(potential_array_3d), z_range, args)
