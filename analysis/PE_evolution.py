import argparse
import os
import signac
import numpy as np
import matplotlib.pyplot as plt


def calculate_PE(job):
    print("\nCalculating PE for", job.ws)
    log_file = os.path.join(job.ws, "output.log")
    log_data = np.genfromtxt(log_file)
    PE_vals = log_data[-100:, 6]
    return np.average(PE_vals)


def plot_PE_evolution(temperature_data, PE_data, no_plot_melting_point):
    temp_centres = (np.array(temperature_data[1:]) + np.array(temperature_data[:-1])) / 2.0
    dE_dT = np.diff(PE_data)
    temp_centres_2 = (np.array(temp_centres[1:]) + np.array(temp_centres[:-1])) / 2.0
    dE_dT_2 = np.diff(dE_dT)
    # Find the maximum value of dE_dT_2
    max_grad_change_index = np.argmax(dE_dT_2) + 1
    min_grad_change_index = np.argmin(dE_dT_2) + 1
    max_grad_change_temp = temperature_data[max_grad_change_index]
    min_grad_change_temp = temperature_data[min_grad_change_index]
    # Linearly interpolate between the two to get the transition temperature
    transition_temp = (max_grad_change_temp + min_grad_change_temp) / 2.0
    # Now plot the graph
    plt.plot(temperature_data, PE_data)
    if not no_plot_melting_point:
        plt.axvline(min_grad_change_temp, c="r")
        plt.axvline(max_grad_change_temp, c="r")
        plt.axvline(transition_temp, c="k")
        plt.title("".join(["MP = {:.1f}".format(transition_temp)]))
    plt.xlabel("Temperature (K)")
    plt.ylabel("PE (Arb. U.)")
    save_dir = "../outputs"
    save_file = "PE.png"
    try:
        plt.savefig(os.path.join(save_dir, save_file))
    except FileNotFoundError:
        os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, save_file))
    plt.clf()
    print(
        "Potential Energy Evolution plot saved as",
        os.path.join(save_dir, save_file)
    )
    # plt.figure()
    # plt.plot(temp_centres, dE_dT)
    # plt.title("".join(["Size = ln_", str(size), "_dEdT"]))
    # plt.show()
    # plt.figure()
    # plt.plot(temp_centres_2, dE_dT_2)
    # plt.axhline(0.0, c="k")
    # plt.title("".join(["Size = ln_", str(size), "_d2EdT2"]))
    # plt.show()
    return transition_temp


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
        "-m",
        "--no_melting_point",
        action="store_false",
        required=False,
        help=(
            "If present, skip plotting the melting point on the output figure."
        ),
    )
    args, directory_list = parser.parse_known_args()
    project = signac.get_project("../")
    job_temperatures = []
    job_potentials = []
    plt.figure()
    for job in project:
        if args.job is not None:
            if job.get_id() != args.job:
                continue
        if (job.sp.job_type == "parent"):
            continue
        temperature = job.sp.temperature
        potential_energy = calculate_PE(job)
        job_temperatures.append(temperature)
        job_potentials.append(potential_energy)
    job_temperatures, job_potentials = zip(*sorted(zip(job_temperatures, job_potentials)))
    melting_temp = plot_PE_evolution(job_temperatures, job_potentials, args.no_melting_point)
    print("Melting temperature =", melting_temp)
