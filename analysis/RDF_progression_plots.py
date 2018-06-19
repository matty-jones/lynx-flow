import os
import signac
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np

"""
This module plots the RDF first-peak progression for each job in the workspace.
Evolutions are hardcoded in and show the evolution in the first-peak as a
function of temperature for each z_reactor_size at a given crystal dimension.
Files are output to the output_figures directory in lynx-flow's root.
"""


def plot_first_peaks(project):
    colours = pl.cm.seismic(np.linspace(0, 1, 4))
    for plot_index, plot_type in enumerate(["Position (A)", "Magnitude (Arb. U.)"]):
        print(" ".join(["Creating", plot_type[:3].lower(), "plots..."]))
        for dimension in ["10x10x1", "10x10x2", "10x10x3"]:
            print("Creating plots for dimensions =", dimension)
            plt.figure()
            for z_size in [10, 15, 20, 25]:
                print("Creating plots for reactor size =", z_size)
                for index, metal in enumerate(["Mo", "Nb", "Te", "V"]):
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
                                "_".join(["RDF_first_peak", metal])
                            ]
                            first_peak_data.append(first_peak[plot_index])
                            temperatures.append(job.sp()["temperature"])
                        except KeyError:
                            continue
                    plot_T, plot_peak = zip(*sorted(zip(temperatures, first_peak_data)))
                    plt.plot(plot_T, plot_peak, color=colours[index], label=metal)
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
