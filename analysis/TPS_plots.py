import os
import signac
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np

"""
This module plots the TPS progression for each job in the workspace.
Evoliutions are hardcoded in and show the evolution in the TPS as a function of
temperature for each z_reactor_size at a given crystal dimension.
Files are output to the output_figures directory in rhaco-flow's root.
"""


def extract_av_tps(file_name):
    with open(file_name, "r") as o_file:
        lines = o_file.readlines()
    # Reverse the lines to make it faster
    lines = lines[::-1]
    for line in lines:
        if "Average TPS" in line:
            tps = float(line.split()[-1])
            return tps
    return None


def plot_tpses(project):
    colours = pl.cm.plasma(np.linspace(0, 1, 4))
    for dimension in ["10x10x1", "10x10x2", "10x10x3"]:
        print("Creating plot for dimensions =", dimension)
        plt.figure()
        for index, z_reactor_size in enumerate([10, 15, 20, 25]):
            print("Plotting line for reactor size =", z_reactor_size)
            temperatures = []
            tpses = []
            for job in project.find_jobs(
                {"dimensions": dimension, "z_reactor_size": z_reactor_size}
            ):
                tps = extract_av_tps(os.path.join(job.ws, "stdout.o"))
                temperatures.append(job.sp()["temperature"])
                tpses.append(tps)
            temperatures, tpses = zip(*sorted(zip(temperatures, tpses)))
            plt.plot(
                temperatures,
                tpses,
                color=colours[index],
                label="Z = " + str(z_reactor_size),
            )
        plt.xlabel("T (K)")
        plt.ylabel("TPS (Arb. U.)")
        plt.title("Dims = " + dimension)
        plt.legend(prop={"size": 10})
        plt.savefig("../output_figures/tps_" + dimension + ".pdf")
        plt.close()


if __name__ == "__main__":
    project = signac.get_project("../")
    # Plot TPS variation
    plot_tpses(project)
