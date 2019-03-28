import os
import numpy as np
import matplotlib.pyplot as plt


def get_T_PE(morphologies):
    temperature_data = []
    PE_data = []
    for morphology in morphologies:
        print("Examining", morphology, "...")
        temperature = int(morphology.split("-T_")[1][:4])
        temperature_data.append(temperature)
        try:
            log_file = os.path.join(morphology, "analyze_log.log")
            log_data = np.genfromtxt(log_file)
        except OSError:
            log_file = os.path.join(morphology, "analyze.log")
            log_data = np.genfromtxt(log_file)
        PE_vals = log_data[-100:, 4]
        PE_data.append(np.average(PE_vals))
    # Parallel sort the temperature and PE data
    temperature_data, PE_data = list(zip(*sorted(zip(temperature_data, PE_data))))
    return temperature_data, PE_data


def plot_PE_evolution(temperature_data, PE_data, size):
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
    plt.axvline(min_grad_change_temp, c="r")
    plt.axvline(max_grad_change_temp, c="r")
    plt.axvline(transition_temp, c="k")
    plt.title("".join(["Size = ln_", str(size), " mp = {:.1f}".format(transition_temp)]))
    plt.xlabel("Temperature (K)")
    plt.ylabel("PE (Arb. U.)")
    plt.savefig("".join(["PE_evolution_size_", str(size), ".png"]))
    plt.clf()
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
    melting_points = []
    particle_sizes = []
    particle_lookup = {
        10: 4.7,
        20: 10.3,
        30: 15.5,
        40: 20.0,
        50: 25.4,
    }
    bulk_morphologies = []
    plt.figure()
    for size in [10, 20, 30, 40, 50]:
        print("Gathering data for particles with size ln =", size, "=", particle_lookup[size], "nm...")
        morphologies = []
        for dir_name in os.listdir():
            if ("_".join(["ln", str(size)]) in dir_name):
                if ("lm_1.0" not in dir_name):
                    morphologies.append(dir_name)
                elif dir_name not in bulk_morphologies:
                    bulk_morphologies.append(dir_name)
        temperature_data, PE_data = get_T_PE(morphologies)
        transition_temp = plot_PE_evolution(temperature_data, PE_data, size)
        melting_points.append(transition_temp)
        particle_sizes.append(particle_lookup[size])
    # if len(bulk_morphologies) > 0:
    #     print("Examining the bulk data...")
    #     bulk_temp, bulk_PE = get_T_PE(bulk_morphologies)
    #     bulk_transition = plot_PE_evolution(bulk_temp, bulk_PE, "bulk")
    #     plt.axhline(bulk_transition, c="r")
    plt.plot(particle_sizes, melting_points)
    plt.axhline(1234.8, c="k")
    plt.xlabel("d (nm)")
    plt.ylabel("Temperature (K)")
    plt.savefig("melting_point_curve_evolution.png")
    plt.close()
