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


def calculate_xy_profile(job, z_lim):
    print("\nCalculating XY-Plane profile for", job.ws)
    save_dir = os.path.join(job.ws, "XY_profile")
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    gsd_file_name = os.path.join(job.ws, "output_traj.gsd")
    try:
        gsd_file = gsd.fl.GSDFile(gsd_file_name, "rb")
    except OSError:
        print(gsd_file_name, "not found. Skipping...")
        return
    trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
    final_frame = trajectory[-1]
    position_matrix = final_frame.particles.position
    reactant_positions = position_matrix[np.where(position_matrix[:,2] > z_lim)]
    min_x = np.min(reactant_positions[:,0])
    max_x = np.max(reactant_positions[:,0])
    x_bins = np.linspace(np.floor(min_x), np.ceil(max_x), 50)
    central_bins = (x_bins[1:] + x_bins[:-1]) / 2.0
    plt.clf()
    n, _, _ = plt.hist(reactant_positions[:,0], bins=x_bins)
    smoothed_n = gaussian_filter(n, 1.0)
    troughs = argrelextrema(smoothed_n, np.less)[0]
    if len(troughs) > 1:
        # Multiple turning points, find the deepest minimum
        trough_val = 9E99
        deepest_trough = 0
        for trough in troughs:
            if smoothed_n[trough] <= trough_val:
                deepest_trough = trough
                trough_val = smoothed_n[trough]
    else:
        deepest_trough = troughs
    plt.plot(central_bins, smoothed_n)
    try:
        plt.title("Neck = {:.3f}".format(float(smoothed_n[deepest_trough])))
    except TypeError:
        print("NECK NOT FOUND FOR THIS JOB")
        plt.title("NECK NOT FOUND")
    plt.xlabel("X position (Ang)")
    plt.ylabel("Particles in Slice (Arb. U.)")
    plt.savefig(os.path.join(save_dir, "xy_profile.png"))
    return smoothed_n[troughs][0]


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
        "-z",
        "--z_lim",
        type=str,
        required=False,
        default=14.9150,
        help=(
            "The z-value of the lowest non-surface atom in the simulation."
        ),
    )
    args, directory_list = parser.parse_known_args()
    project = signac.get_project("../")
    job_temperatures = []
    job_neck_widths = []
    plt.figure()
    for job in project:
        if args.job is not None:
            if job.get_id() != args.job:
                continue
        if (job.sp.job_type == "parent"):
            continue
        temperature = job.sp.temperature
        neck_width = calculate_xy_profile(job, args.z_lim)
        job_temperatures.append(temperature)
        job_neck_widths.append(neck_width)
    job_temperatures, job_neck_widths = zip(*sorted(zip(job_temperatures, job_neck_widths)))
    plt.clf()
    plt.plot(job_temperatures, job_neck_widths)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Neck width (Ang)")
    save_file = "../outputs/neck_evolution.png"
    plt.savefig(save_file)
    print("Sintering Neck Evolution plot saved as", save_file)
