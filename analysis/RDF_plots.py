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
All 22 files (2 * 11) are written for each metal atom in the system (Mo, Nb,
Te, V).
"""


def find_crystal_extents_z(project, type_name, args):
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
        print("\nConsidering job", job.ws)
        print("".join(["Calculating Z-distribution for ", type_name, "..."]))
        # print(job.statepoint.dimensions)
        # if job.statepoint.dimensions.split('x')[2] != '1':
        #     continue
        save_dir = os.path.join(job.ws, "RDFs")
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        gsd_file_name = os.path.join(job.ws, "output_traj.gsd")
        gsd_file = gsd.fl.GSDFile(gsd_file_name, "rb")
        trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
        atom_type = trajectory[0].particles.types.index(type_name)
        final_frame = trajectory[-1]
        atom_posns = final_frame.particles.position[
            np.where(final_frame.particles.typeid == atom_type)
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
                    -job.statepoint.crystal_separation,
                    job.statepoint.crystal_separation,
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


def plot_rdf(project, type1_name, type2_name, args, r_max=20, stride=50):
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
        gsd_file = gsd.fl.GSDFile(gsd_file_name, "rb")
        trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
        sim_box = trajectory[0].configuration.box[:3]
        av_rdf = freud.density.RDF(rmax=r_max, dr=0.1)
        av_rdf.resetRDF()
        type1 = trajectory[0].particles.types.index(type1_name)
        type2 = trajectory[0].particles.types.index(type2_name)
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
            type1_pos = frame.particles.position[
                np.where(frame.particles.typeid == type1)
            ]
            # old_type2_pos = frame.particles.position[
            #     np.where(frame.particles.typeid == type2)]
            # New Type 2 pos needs to be all particles of type type2,
            # but with min_z <= z-position < max_z
            type2_pos = frame.particles.position[
                np.where(
                    (frame.particles.typeid == type2)
                    & (
                        (
                            frame.particles.position[:, 2]
                            >= job.document["crystal_max_z"]
                        )
                        | (
                            frame.particles.position[:, 2]
                            <= job.document["crystal_min_z"]
                        )
                    )
                )
            ]
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
            "If present, only consider the job in the current" " directory's workspace."
        ),
    )
    args, directory_list = parser.parse_known_args()
    project = signac.get_project("../")
    for type2 in ["Mo", "V", "Nb", "Te"]:
        # Plot distribution of z-values for atoms
        find_crystal_extents_z(project, type2, args)
        # Plot RDF variation
        plot_rdf(project, "C", type2, args)
