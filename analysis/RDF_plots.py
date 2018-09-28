import csv
import freud
import gsd.fl
import gsd.hoomd
import os
import signac
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import argparse
from scipy.signal import argrelextrema
from scipy.ndimage import gaussian_filter


"""
This module plots several RDFs for each job in the workspace.
10 RDFs are plotted over time as the simulation progresses (data taken from the
trajectory GSD), as well as an aggregated RDF describing the average over the
entire simulation.
Additionally, csv files are written for every RDF to permit subsequent analysis
All 22 files (2 * 11) are written for each surface atom in the system (Mo, Nb,
Te, V).
"""


def find_crystal_extents_z(project, args):
    """
    This function finds the z coordinates of the outermost layer of crystal,
    i.e. the one that is exposed to the reactant and will actually undergo the
    catalysed reaction.
    This is important when calculating the RDFs to make sure that the features
    do not get washed out as more crystal layers are included.
    It also calculates the extreme edge of the crystal at the top and bottom, which is
    important when determining the residence times of the reactants.
    """
    for job in project:
        if args.job is not None:
            if job.get_id() != args.job:
                continue
        # Skip if this is a parent job
        if "job_type" in job.sp:
            if job.sp.job_type == "parent":
                continue
        if not args.overwrite:
            # Skip the calculation if the crystal extents already exist
            try:
                job.document["crystal_top_edge"]
                job.document["crystal_bot_edge"]
                job.document["crystal_top_layer"]
                job.document["crystal_bot_layer"]
                print("Crystal dimensions already calculated, skipping...")
                continue
            except KeyError:
                pass
        print("\nConsidering job", job.ws)
        crystal_offset = float(job.sp["crystal_separation"]) / 2.0
        crystal_z = float(job.sp["crystal_z"]) * 10.0 # convert from nm to ang
        dim_z = int(job.sp["dimensions"].split("x")[2])
        job.document["crystal_top_edge"] = crystal_offset + (dim_z * crystal_z / 2.0)
        job.document["crystal_bot_edge"] = -(crystal_offset + (dim_z * crystal_z / 2.0))
        job.document["crystal_top_layer"] = crystal_offset + (
            (dim_z - 2) * crystal_z / 2.0
        )
        job.document["crystal_bot_layer"] = -(crystal_offset + (
            (dim_z - 2) * crystal_z / 2.0
        ))


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
        raise SystemError("WOBBEY")
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
        # Skip calculating the COM if none of the atoms in this molecule/sublist
        # have satisfied the crystal max and min conditions
        if len(sublist_positions) > 0:
            type_positions.append(
                calc_COM(sublist_positions, list_of_atom_types=sublist_types)
            )
    return np.array(type_positions)


def plot_rdf(project, type1_name, type2_name, args, r_max=20, stride=50, type1_by_mol=False, type2_by_mol=False):
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
            "".join(["Calculating RDFs between ", type1_name, "-", type2_name, "..."])
        )
        save_dir = os.path.join(job.ws, "RDFs")
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        gsd_file_name = os.path.join(job.ws, "output_traj.gsd")
        try:
            gsd_file = gsd.fl.GSDFile(gsd_file_name, "rb")
        except OSError:
            print(gsd_file_name, "not found. Skipping...")
            continue
        trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
        sim_box = trajectory[0].configuration.box[:3]
        av_rdf = freud.density.RDF(rmax=r_max, dr=0.1)
        av_rdf.resetRDF()
        type1_ID = trajectory[0].particles.types.index(type1_name)
        type2_ID = trajectory[0].particles.types.index(type2_name)

        print("Splitting Molecules...")
        # AAID_to_molID, molID_to_AAIDs = split_molecules(trajectory[0])
        # Split_molecules is super fast (because the bond list is already a lookup table).
        # instead of running it once and then iterating through every AAID to find the molecules
        # we want, change split_molecules to also accept an atom type.
        # Then just have it return [[mol1_atom1, mol1_atom2,...], [mol2_atom1, mol2_atom2,...], ...]
        # which is type1_AAIDs.



        if type1_by_mol is True:
            print("Determining molecules containing type_1...")
            molID_to_AAIDs = split_molecules(trajectory[0], type1_ID)
            type1_AAIDs = list(molID_to_AAIDs.values())
            # OLD CODE, SUPER SLOW
            # # Return AAIDs for all atoms in the same molecule, if that molecule includes
            # # a type1_ID atom
            # type1_AAIDs = []
            # mols_considered = []
            # for AAID, molID in enumerate(AAID_to_molID):
            #     if trajectory[0].particles.typeid[AAID] == type1_ID:
            #         if molID in mols_considered:
            #             continue
            #         mols_considered.append(molID)
            #         AAIDs_in_mol = molID_to_AAIDs[molID]
            #         type1_AAIDs.append(AAIDs_in_mol)
        else:
            type1_AAIDs = [[AAID] for AAID in np.where(trajectory[0].particles.typeid == type1_ID)[0]]

        if type2_by_mol is True:
            print("Determining molecules containing type_2...")
            molID_to_AAIDs = split_molecules(trajectory[0], type2_ID)
            type2_AAIDs = list(molID_to_AAIDs.values())
            # Return AAIDs for all atoms in the same molecule, if that molecule includes
            # a type1_ID atom
            # OLD CODE, SUPER SLOW
            # type2_AAIDs = []
            # mols_considered = []
            # for AAID, molID in enumerate(AAID_to_molID):
            #     if trajectory[0].particles.typeid[AAID] == type2_ID:
            #         if molID in mols_considered:
            #             continue
            #         mols_considered.append(molID)
            #         AAIDs_in_mol = molID_to_AAIDs[molID]
            #         type2_AAIDs.append(AAIDs_in_mol)
        else:
            type2_AAIDs = [[AAID] for AAID in np.where(trajectory[0].particles.typeid == type2_ID)[0]]

        for frame_no, frame in enumerate(trajectory):
            print(
                "".join(
                    [
                        "\rCalculating RDF for frame ",
                        str(frame_no + 1),
                        " of ",
                        str(len(trajectory)),
                    ]
                ),
                end=" ",
            )
            frame_rdf = freud.density.RDF(rmax=r_max, dr=0.1)
            frame_rdf.resetRDF()
            # In case box size changes
            box = frame.configuration.box
            freud_box = freud.box.Box(Lx=box[0], Ly=box[1], Lz=box[2])

            type1_pos = get_type_positions(type1_AAIDs, frame)
            type2_pos = get_type_positions(
                type2_AAIDs,
                frame,
                crystal_min_z=job.document["crystal_bot_layer"],
                crystal_max_z=job.document["crystal_top_layer"],
            )
            av_rdf.accumulate(freud_box, type1_pos, type2_pos)
            if frame_no % stride == 0:
                print(
                    "".join(
                        [
                            "\rSaving RDF image for frame ",
                            str(frame_no + 1),
                            " of ",
                            str(len(trajectory)),
                        ]
                    ),
                    end=" ",
                )
                frame_rdf.compute(freud_box, type1_pos, type2_pos)
                frame_rdf_title = "".join(
                    ["RDF_", type1_name, "-", type2_name, "_{:03d}".format(frame_no)]
                )
                with open(
                    os.path.join(save_dir, frame_rdf_title + ".csv"), "w+"
                ) as csv_file:
                    writer = csv.writer(csv_file)
                    writer.writerow(["r", "g(r)"])
                    for r, g_r in zip(frame_rdf.getR(), frame_rdf.getRDF()):
                        writer.writerow([r, g_r])
                plt.figure()
                plt.title(frame_rdf_title)
                plt.plot(frame_rdf.getR(), frame_rdf.getRDF())
                plt.xlabel("r (Ang)")
                plt.ylabel("RDF (Arb. U.)")
                plt.savefig(os.path.join(save_dir, frame_rdf_title + ".pdf"))
                plt.close()
        print("\rCalculating RDF averaged over all frames", end=" ")
        av_rdf_title = "RDF_" + type1_name + "-" + type2_name + "_Av"
        with open(os.path.join(save_dir, av_rdf_title + ".csv"), "w+") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["r", "g(r)"])
            for r, g_r in zip(av_rdf.getR(), av_rdf.getRDF()):
                writer.writerow([r, g_r])
        plt.figure()
        plt.title(av_rdf_title)
        plt.plot(av_rdf.getR(), av_rdf.getRDF())
        plt.xlabel("r (Ang)")
        plt.ylabel("RDF (Arb. U.)")
        plt.savefig(os.path.join(save_dir, av_rdf_title + ".pdf"))
        plt.close()
        print("")


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
        "-o",
        "--overwrite",
        required=False,
        action="store_true",
        help=(
            "Recalculate any job.document properties and update them, regardless of"
            " whether they have already been calculated or not."
        ),
    )
    args, directory_list = parser.parse_known_args()
    project = signac.get_project("../")
    # Find crystal extents for each job in the project
    find_crystal_extents_z(project, args)
    surface_atom_types = []
    for stoic_dict_str, _ in project.groupby('stoichiometry'):
        surface_atom_types += list(eval(stoic_dict_str).keys())
    surface_atom_types = list(set(surface_atom_types))
    for atom_type in surface_atom_types:
        # Plot RDF variation
        plot_rdf(project, "C", atom_type, args, type1_by_mol=True)
