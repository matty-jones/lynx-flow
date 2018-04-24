import signac
import os
import freud
import gsd.fl
import gsd.hoomd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np


def extract_av_tps(file_name):
    with open(file_name, 'r') as o_file:
        lines = o_file.readlines()
    # Reverse the lines to make it faster
    lines = lines[::-1]
    for line in lines:
        if 'Average TPS' in line:
            tps = float(line.split()[-1])
            return tps
    return None


def plot_tpses(project):
    colours = pl.cm.plasma(np.linspace(0, 1, 4))
    for dimension in ['10x10x1', '10x10x2', '10x10x3']:
        print("Creating plot for dimensions =", dimension)
        plt.figure()
        for index, z_reactor_size in enumerate([10, 15, 20, 25]):
            print("Plotting line for reactor size =", z_reactor_size)
            temperatures = []
            tpses = []
            for job in project.find_jobs({'dimensions': dimension,
                                         'z_reactor_size': z_reactor_size}):
                tps = extract_av_tps(os.path.join(job.ws,
                                                  'stdout.o'))
                temperatures.append(job.sp()['temperature'])
                tpses.append(tps)
            temperatures, tpses = zip(*sorted(zip(temperatures, tpses)))
            plt.plot(temperatures, tpses, color=colours[index], label='Z = ' + str(z_reactor_size))
        plt.xlabel('T (K)')
        plt.ylabel('TPS (Arb. U.)')
        plt.title('Dims = ' + dimension)
        plt.legend(prop={'size':10})
        plt.savefig('./output_figures/tps_' + dimension + '.pdf')
        plt.close()


def plot_rdf(project, type1_name, type2_name, r_max=20, stride=50):
    for job in project:
        print("Considering job", job.ws)
        print("Calculating RDFs between " + type1_name + "-" + type2_name + "...")
        save_dir = os.path.join(job.ws, 'RDFs')
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        gsd_file_name = os.path.join(job.ws, 'output_traj.gsd')
        gsd_file = gsd.fl.GSDFile(gsd_file_name, 'rb')
        trajectory = gsd.hoomd.HOOMDTrajectory(gsd_file)
        sim_box = trajectory[0].configuration.box[:3]
        av_rdf = freud.density.RDF(rmax=r_max, dr=0.1)
        av_rdf.resetRDF()
        type1 = trajectory[0].particles.types.index(type1_name)
        type2 = trajectory[0].particles.types.index(type2_name)
        for frame_no, frame in enumerate(trajectory):
            print("\rCalculating RDF for frame", frame_no + 1, "of", len(trajectory), end=' ')
            frame_rdf = freud.density.RDF(rmax=r_max, dr=0.1)
            frame_rdf.resetRDF()
            # In case box size changes
            box = frame.configuration.box
            freud_box = freud.box.Box(Lx = box[0], Ly = box[1], Lz = box[2])
            type1_pos = frame.particles.position[
                np.where(frame.particles.typeid == type1)]
            type2_pos = frame.particles.position[
                np.where(frame.particles.typeid == type2)]
            av_rdf.accumulate(freud_box, type1_pos, type2_pos)
            if frame_no%stride == 0:
                print("\rSaving RDF image for frame", frame_no + 1, "of", len(trajectory), end=' ')
                frame_rdf.compute(freud_box, type1_pos, type2_pos)
                frame_rdf_title = 'RDF_' + type1_name + '-' + type2_name + '_{:03d}'.format(frame_no)
                plt.figure()
                plt.title(frame_rdf_title)
                plt.plot(frame_rdf.getR(), frame_rdf.getRDF())
                plt.xlabel('r (Ang)')
                plt.ylabel('RDF (Arb. U.)')
                plt.savefig(os.path.join(save_dir, frame_rdf_title + '.pdf'))
                plt.close()
        print("\rCalculating RDF averaged over all frames", end=' ')
        av_rdf_title = 'RDF_' + type1_name + '-' + type2_name + '_Av'
        plt.figure()
        plt.title(av_rdf_title)
        plt.plot(av_rdf.getR(), av_rdf.getRDF())
        plt.xlabel('r (Ang)')
        plt.ylabel('RDF (Arb. U.)')
        plt.savefig(os.path.join(save_dir, av_rdf_title + '.pdf'))
        plt.close()
        print("")


if __name__ == "__main__":
    project = signac.get_project('./')
    plot_tpses(project)
    for type2 in ['Mo', 'V', 'Nb', 'Te']:
        plot_rdf(project, 'C', type2)
