import os
import signac
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np

"""
This module plots the RDF first-peak progression for each job in the workspace.
Evolutions are hardcoded in and show the evolution in the first-peak as a
function of temperature for each z_reactor_size at a given crystal dimension.
Files are output to the output_figures directory in rhaco-flow's root.
"""


def plot_first_peaks(project):
    colours = pl.cm.coolwarm(np.linspace(0, 1, 4))
    all_dimensions = []
    all_z_sizes = []
    all_surface_atom_types = []
    # Iterate through all jobs in the project to find the unique dimensions
    # and reactor sizes
    for job in project:
        statepoint = job.sp()
        all_dimensions.append(statepoint["dimensions"])
        all_z_sizes.append(statepoint["z_reactor_size"])
        all_surface_atom_types += list(eval(statepoint["stoichiometry"]).keys())
    all_dimensions = list(set(all_dimensions))
    all_z_sizes = list(set(all_z_sizes))
    all_surface_atom_types = list(set(all_surface_atom_types))

    for plot_index, plot_type in enumerate(["Position (A)", "Magnitude (Arb. U.)"]):
        print(" ".join(["Creating", plot_type[:3].lower(), "plots..."]))
        for dimension in all_dimensions:
            print("Creating plots for dimensions =", dimension)
            # New figure for every dimension
            plt.figure()
            for z_size in all_z_sizes:
                # New line on the plot for each z_size
                print("Creating plots for reactor size =", z_size)
                for index, atom_type in enumerate(all_surface_atom_types):
                    temperatures = []
                    first_peak_data = []
                    for job in project.find_jobs(
                        {"dimensions": dimension,
                         "z_reactor_size": z_size,
                         "job_type": "child",
                        }
                    ):
                        try:
                            first_peak = job.document[
                                "_".join(["RDF_first_peak", atom_type])
                            ]
                            first_peak_data.append(first_peak[plot_index])
                            temperatures.append(job.sp()["temperature"])
                        except KeyError:
                            continue
                    plot_T, plot_peak = zip(*sorted(zip(temperatures, first_peak_data)))
                    plt.plot(plot_T, plot_peak, color=colours[index], label=atom_type)
                plt.xlabel("T (K)")
                plt.ylabel(plot_type)
                plt.title("".join(["Dims = ", dimension, ", Z = ", str(z_size)]))
                plt.legend(prop={"size": 10})
                plot_file_name = "".join(
                    [
                        "first_peak_",
                        plot_type[:3].lower(),
                        "_",
                        dimension,
                        "_Z",
                        str(z_size),
                        ".pdf",
                    ]
                )
                plot_save_loc = os.path.join("..", "output_figures", plot_file_name)
                plt.savefig(plot_save_loc)
                print("Figure saved to", plot_save_loc)
                plt.close()


if __name__ == "__main__":
    project = signac.get_project("../")
    # Plot RDF first-peak variation
    plot_first_peaks(project)
