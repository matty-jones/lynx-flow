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


def find_crystal_extents_z(project, type_names, args):
    """
    This function finds the z coordinates of the outermost layer of crystal,
    i.e. the one that is exposed to the reactant and will actually undergo the
    catalysed reaction.
    This is important when calculating the RDFs to make sure that the features
    do not get washed out as more crystal layers are included.
    """
    for job in project:
        if args.job is not None:
            if job.get_id() != args.job:
                continue
        # Skip if this is a parent job
        if "job_type" in job.sp:
            if job.sp.job_type == "parent":
                continue
        # Skip the calculation if the crystal extents already exist
        try:
            job.document["crystal_min_z"]
            job.document["crystal_max_z"]
            print("Crystal dimensions already calculated, skipping...")
            continue
        except KeyError:
            pass
        print("\nConsidering job", job.ws)
        print("".join(["Calculating Z-distribution for ", type_names, "..."]))
        # print(job.statepoint.dimensions)
        # if job.statepoint.dimensions.split('x')[2] != '1':
        #     continue
        # save_dir = os.path.join(job.ws, "RDFs")
        # if not os.path.exists(save_dir):
        #     os.makedirs(save_dir)
        gsd_file_name = os.path.join(job.ws, "output_traj.gsd")
        try:
            gsd_file = gsd.fl.GSDFile(gsd_file_name, "rb")
        except OSError:
            print(gsd_file_name, "not found. Skipping...")
            continue
        trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
        atom_types = [trajectory[0].particles.types.index(type_name) for type_name in type_names]
        final_frame = trajectory[-1]
        atom_posns = final_frame.particles.position[
            np.where(final_frame.particles.typeid in atom_types)
        ]
        n_troughs = 0
        target_troughs = (int(job.statepoint.dimensions.split("x")[2]) - 1) * 2
        n_bins = 20
        bins = None
        central_bins = None
        n = None
        while n_troughs < target_troughs:
            print("Currently found", n_troughs, "of", target_troughs, "troughs...")
            n, bins = np.histogram(
                atom_posns[:, 2],
                bins=np.linspace(
                    -float(job.statepoint.crystal_separation),
                    float(job.statepoint.crystal_separation),
                    n_bins,
                ),
            )
            # n, bins, patches = plt.hist(atom_posns[:,2], bins = np.linspace(-job.statepoint.crystal_separation, job.statepoint.crystal_separation, n_bins))
            central_bins = (bins[1:] + bins[:-1]) / 2.0
            troughs = argrelextrema(n, np.less)[0]
            n_troughs = len(troughs)
            smoothed_n = gaussian_filter(n, 1.0)
            # plt.figure()
            # plt.title(" ".join(["Z-separation of", type_name]))
            # plt.plot(central_bins, n)
            # plt.plot(central_bins, smoothed_n, c='r')
            # for trough in troughs:
            #     plt.axvline(central_bins[trough], c='k')
            # plt.xlabel('Z-separation (Ang)')
            # plt.ylabel('Frequency (Arb. U.)')
            # plt.show()
            # #plt.savefig(os.path.join(save_dir, av_rdf_title + '.pdf'))
            # plt.close()

            print("Increasing n_bins from", n_bins, "to", n_bins * 2)
            n_bins *= 2
        if n_troughs == 0:
            # Only one layer, so include everything
            job.document["crystal_min_z"] = 0.0
            job.document["crystal_max_z"] = 0.0
        else:
            print("Found all", n_troughs, "troughs:")
            print(troughs)
            # smoothed_n = gaussian_filter(n, 1.0)
            # plt.title(" ".join(["Z-separation of", type_name]))
            # plt.plot(central_bins, n)
            # plt.plot(central_bins, smoothed_n, c='r')
            # for trough in troughs:
            #     plt.axvline(central_bins[trough], c='k')
            # plt.xlabel('Z-separation (Ang)')
            # plt.ylabel('Frequency (Arb. U.)')
            # plt.show()
            # #plt.savefig(os.path.join(save_dir, av_rdf_title + '.pdf'))
            # plt.close()
            trough_positions = [central_bins[x] for x in troughs]
            print("Trough positions =", [central_bins[x] for x in troughs])
            job.document["crystal_min_z"] = min(trough_positions)
            job.document["crystal_max_z"] = max(trough_positions)


def calc_COM(list_of_positions, list_of_atom_types=None, list_of_atom_masses=None):
    """
    This function calculates the centre of mass of a collection of sites/atoms
    (listOfPositions) with corresponding type (listOfAtomTypes) or mass (listOfMasses)
    if list_of_atom_masses is not specified, then list_of_atom_types must be.
    """
    if len(list_of_positions) == 1:
        return np.array(list_of_positions[0])
    mass_weighted = np.array([0.0, 0.0, 0.0])
    if list_of_atom_masses is None:
        list_of_atom_masses = []
        for atom_type in list_of_atom_types:
            # Masses obtained from nist.gov, for the atoms we are likely to
            # simulate the most.
            # Add in new atoms here if your molecule requires it!
            if atom_type.lower()[:2] == "si":
                list_of_atom_masses.append(27.976926)
            elif atom_type.lower()[0] == "c":
                list_of_atom_masses.append(12.000000)
            elif atom_type.lower()[0] == "h":
                list_of_atom_masses.append(1.0078250)
            elif atom_type.lower()[0] == "s":
                list_of_atom_masses.append(31.972071)
            elif atom_type.lower()[0] == "o":
                list_of_atom_masses.append(15.994914)
            elif atom_type.lower()[0] == "n":
                list_of_atom_masses.append(14.003074)
            elif atom_type.lower()[:2] == "mo":
                list_of_atom_masses.append(95.960000)
            elif atom_type.lower()[:2] == "nb":
                list_of_atom_masses.append(92.906380)
            elif atom_type.lower()[:2] == "te":
                list_of_atom_masses.append(127.60000)
            elif atom_type.lower()[:2] == "v":
                list_of_atom_masses.append(50.941500)
            elif atom_type.lower()[:2] == "ni":
                list_of_atom_masses.append(140.91120)
            elif atom_type.lower()[:2] == "ga":
                list_of_atom_masses.append(69.723000)
            elif atom_type.lower()[:2] == "mn":
                list_of_atom_masses.append(54.938045)
            elif atom_type.lower()[:2] == "cu":
                list_of_atom_masses.append(63.546000)
            elif atom_type.lower()[:2] == "ag":
                list_of_atom_masses.append(107.86820)
            elif atom_type.lower()[:2] == "au":
                list_of_atom_masses.append(196.96657)
            elif (atom_type.lower()[0] == "d") or (atom_type.lower()[0] == "a"):
                list_of_atom_masses.append(1.0)
            else:
                raise SystemError(
                    "Unknown atomic mass " + str(atom_type) + ". Please hardcode"
                    " into calc_COM."
                )
    total_mass = np.sum(list_of_atom_masses)
    for atom_ID, position in enumerate(list_of_positions):
        for axis in range(3):
            mass_weighted[axis] += position[axis] * list_of_atom_masses[atom_ID]
    if total_mass == 0.0:
        print("List of posns", list_of_positions)
        print("List of masses", list_of_atom_masses)
        print("Mass weighted", mass_weighted)
        print("Total mass", total_mass)
        raise SystemError("Error in mass calculation.")
    return mass_weighted / float(total_mass)


def obtain_bonded_list(bond_list):
    # Create a lookup table `neighbour list' for all connected atoms called
    # {bondedAtoms}
    bonded_atoms = {}
    for bond in bond_list:
        if bond[0] not in bonded_atoms:
            bonded_atoms[bond[0]] = [bond[1]]
        else:
            bonded_atoms[bond[0]].append(bond[1])
        if bond[1] not in bonded_atoms:
            bonded_atoms[bond[1]] = [bond[0]]
        else:
            bonded_atoms[bond[1]].append(bond[0])
    return bonded_atoms


def split_molecules(frame, atom_type_ID):
    # Split the full morphology into individual molecules
    # Create a lookup table `neighbour list' for all connected atoms called
    # {bondedAtoms}
    bonded_atoms = obtain_bonded_list(frame.bonds.group)
    molecule_list = [i for i in range(len(frame.particles.typeid))]
    # Recursively add all atoms in the neighbour list to this molecule
    for mol_ID in range(len(molecule_list)):
        molecule_list = update_molecule(mol_ID, molecule_list, bonded_atoms)
    # Here we have a list of len(atoms) where each index gives the molID
    mol_ID_dict = {}
    for atom_ID, mol_ID in enumerate(molecule_list):
        if frame.particles.typeid[atom_ID] != atom_type_ID:
            continue
        if mol_ID not in mol_ID_dict:
            mol_ID_dict[mol_ID] = [atom_ID]
        else:
            mol_ID_dict[mol_ID].append(atom_ID)
    mol_containing_type_dict = {}
    # Now only return molecules that contain the atom we care about
    for mol_ID, AAIDs in mol_ID_dict.items():
        if np.any(np.array([frame.particles.typeid[AAID] for AAID in AAIDs]) == atom_type_ID):
            mol_containing_type_dict[mol_ID] = AAIDs
    return mol_containing_type_dict


def update_molecule(atom_ID, molecule_list, bonded_atoms):
    # Recursively add all neighbours of atom number atomID to this molecule
    try:
        for bonded_atom in bonded_atoms[atom_ID]:
            # If the moleculeID of the bonded atom is larger than that of the
            # current one, update the bonded atom's ID to the current one's to
            # put it in this molecule, then iterate through all of the bonded
            # atom's neighbours
            if molecule_list[bonded_atom] > molecule_list[atom_ID]:
                molecule_list[bonded_atom] = molecule_list[atom_ID]
                molecule_list = update_molecule(
                    bonded_atom, molecule_list, bonded_atoms
                )
            # If the moleculeID of the current atom is larger than that of the
            # bonded one, update the current atom's ID to the bonded one's to
            # put it in this molecule, then iterate through all of the current
            # atom's neighbours
            elif molecule_list[bonded_atom] < molecule_list[atom_ID]:
                molecule_list[atom_ID] = molecule_list[bonded_atom]
                molecule_list = update_molecule(atom_ID, molecule_list, bonded_atoms)
            # Else: both the current and the bonded atom are already known to
            # be in this molecule, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single molecule)
        pass
    return molecule_list


def get_type_positions(AAID_list, frame, crystal_min_z=None, crystal_max_z=None):
    type_positions = []
    for AAID_sublist in AAID_list:
        sublist_positions = []
        sublist_types = []
        for AAID in AAID_sublist:
            position = frame.particles.position[AAID]
            if (crystal_min_z is not None) and (crystal_max_z is not None):
                if (position[2] < crystal_max_z) and (position[2] > crystal_min_z):
                    continue
            sublist_positions.append(position)
            type_ID = frame.particles.typeid[AAID]
            atom_type = frame.particles.types[type_ID]
            sublist_types.append(atom_type)
        type_positions.append(calc_COM(sublist_positions, list_of_atom_types=sublist_types))
    return np.array(type_positions)


def plot_residence_time_per_job(project, args):
    for job in project:
        if args.job is not None:
            if job.get_id() != args.job:
                continue
        # Skip if this is a parent job
        if "job_type" in job.sp:
            if job.sp.job_type == "parent":
                continue
        print("\nConsidering job", job.ws)
        print(
            "".join(["Calculating residence times for carbon-containing molecules..."])
        )
        gsd_file_name = os.path.join(job.ws, "output_traj.gsd")
        try:
            gsd_file = gsd.fl.GSDFile(gsd_file_name, "rb")
        except OSError:
            print(gsd_file_name, "not found. Skipping...")
            continue
        # Skip the calculation if the residence times have already been calculated
        try:
            job.document["mean_residence_time"]
            job.document["mean_residence_time_error"]
            print("Mean residence time already calculated, skipping...")
            continue
        except KeyError:
            pass
        trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
        type1_ID = trajectory[0].particles.types.index(args.atom_type)
        molID_to_AAIDs = split_molecules(trajectory[0], type1_ID)
        type1_AAIDs = list(molID_to_AAIDs.values())
        z_max = job.document["crystal_max_z"] + args.tolerance
        z_min = job.document["crystal_min_z"] - args.tolerance

        residence_dict = {}
        for frame_no, frame in enumerate(trajectory):
            print(
                "".join(
                    [
                        "\rCalculating molecule residence for frame ",
                        str(frame_no + 1),
                        " of ",
                        str(len(trajectory)),
                    ]
                ),
                end=" ",
            )
            type1_pos = get_type_positions(type1_AAIDs, frame)
            # Brute force approach, check all mols every frame
            for mol_ID, position in enumerate(type1_pos):
                if (position[2] > z_min) and (position[2] < z_max):
                    # Molecule is residing
                    if mol_ID in residence_dict:
                        residence_dict[mol_ID].append(frame_no)
                    else:
                        residence_dict[mol_ID] = [frame_no]
        print("")
        # Turn the residence dictionary into a histogram of residence frames
        residence_frames = []
        for frames in residence_dict.values():
            for _, group in itertools.groupby(enumerate(frames), lambda ix: ix[0] - ix[1]):
                # Slice the list into consecutive chunks
                slice_list = list(map(operator.itemgetter(1), group))
                # Append the residence_frames with the length of each chunk
                residence_frames.append(len(slice_list))
        # Convert frames to times
        residence_times = calculate_residence_times(residence_frames, job)
        mean_residence_time = np.mean(residence_times)
        residence_time_error = np.std(residence_times) / np.sqrt(len(residence_times))
        mean_string = r"{:.2E} $\pm$ {:.1E}".format(mean_residence_time, residence_time_error)
        print("Mean residence time ==", mean_string)
        # Now plot the residence histogram
        plt.clf()
        plt.hist(residence_times / 1E-9, bins=20)
        plt.xlabel("Residence (ns)")
        plt.ylabel("Frequency (Molecules)")
        plt.title(mean_string)
        plt.savefig(os.path.join(job.ws, "residence_time.pdf"))
        job.document["mean_residence_time"] = mean_residence_time
        job.document["mean_residence_time_error"] = residence_time_error


def calculate_residence_times(residence_frames, job):
    # Firstly, how many timesteps does one frame equal?
    # Rhaco dumps either 500 evenly-distributed frames across the whole simulation
    # or one frame per timestep - whichever is longer
    frame_period_timesteps = max([int(int(job.statepoint["run_time"]) / 500), 1])
    # Now multiply this by the timestep to get the period in hoomd time units
    frame_period_dimless = frame_period_timesteps * job.statepoint["timestep"]
    # Now calculate the dimensionless time in SI to get residence time
    mass_factor = 1.0 * AMU_TO_KG
    distance_factor = job.statepoint.get("distance_scale_unit", 1.0) * ANG_TO_M
    energy_factor = job.statepoint.get("energy_scale_unit", 1.0) * KCAL_TO_J / AVOGADRO
    dimless_time_factor = np.sqrt((mass_factor * (distance_factor)**2) / energy_factor)
    # Multiply the factors to get the frame_period in SI
    frame_period_SI = frame_period_dimless * dimless_time_factor
    # Multiply this by residence_frames to get the residence_times
    return frame_period_SI * np.array(residence_frames)


def plot_residence_time_vs_temp(project):
    temperatures = {}
    residence_times = {}
    for job in project:
        if "job_type" in job.sp:
            if job.sp.job_type == "parent":
                continue
        dimensions = job.sp["dimensions"]
        temperature = job.sp["temperature"]
        z_size = job.sp["z_reactor_size"]
        mean_res_time = job.document["mean_residence_time"]
        mean_res_time_err = job.document["mean_residence_time_error"]
        identify_string = "".join(["dims", str(dimensions), "z_size", str(z_size)])
        if identify_string not in temperatures:
            temperatures[identify_string] = []
            residence_times[identify_string] = []
        temperatures[identify_string].append(temperature)
        residence_times[identify_string].append([mean_res_time, mean_res_time_err])
    for key, temp_vals in temperatures.items():
        res_time_vals = residence_times[key]
        temp_vals, res_time_vals_combined = zip(*sorted(zip(temp_vals, res_time_vals)))
        mean_time, error = zip(*res_time_vals_combined)
        plt.clf()
        plt.errorbar(temp_vals, np.array(mean_time)/1E-12, yerr=np.array(error)/1E-12)
        plt.xlabel("Temperature (K)")
        plt.ylabel("Mean residence time (ps)")
        fig_name = "".join(["../outputs/", key, "_res_time.pdf"])
        plt.savefig(fig_name)
        print("Figure saved as", fig_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-j",
        "--job",
        type=str,
        required=False,
        default=None,
        help=(
            "If present, only consider the job in the current directory's workspace."
        ),
    )
    parser.add_argument(
        "-a",
        "--atom_type",
        type=str,
        required=False,
        default="C",
        help=(
            "Find the residence times of molecules containing this type of atom."
        ),
    )
    parser.add_argument(
        "-t",
        "--tolerance",
        type=float,
        required=False,
        default=10.0,
        help=(
            "If a molecule containing args.atom_type has a centre-of-mass location"
            " within args.tolerance of the surface, it is considered as residing on"
            " the surface."
        ),
    )
    args, directory_list = parser.parse_known_args()
    project = signac.get_project("../")
    plt.figure()
    surface_atom_types = []
    # schema = project.detect_schema()
    # print(schema["reactant_composition"])
    # reactant_atom_types = next(iter(list(schema["reactant_composition"].values())[0]))
    # print(reactant_atom_types, type(reactant_atom_types))
    # exit()
    for stoic_dict_str, _ in project.groupby('stoichiometry'):
        surface_atom_types += list(eval(stoic_dict_str).keys())
    surface_atom_types = list(set(surface_atom_types))
    # Plot distribution of z-values for atoms
    find_crystal_extents_z(project, surface_atom_types, args)
    plot_residence_time_per_job(project, args)
    if args.job is None:
        plot_residence_time_vs_temp(project)
