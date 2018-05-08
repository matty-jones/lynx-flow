import argparse
import csv
import freud
import gsd.fl
import gsd.hoomd
import os
import signac
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import pandas as pd
import numpy as np
from scipy.signal import argrelextrema
from scipy.ndimage import gaussian_filter


def get_first_peak(project, metal=None, layers=None, z_size=None):
    job_list = project.find_jobs({
        'dimensions': 'x'.join(['10', '10', str(layers)]),
        'z_reactor_size': z_size})
    for job in job_list:
        csv_file_location = os.path.join(job.ws, 'RDFs',
                                         ''.join(['RDF_C-', metal,
                                                 '_Av.csv']))
        RDF_Data = pd.read_csv(csv_file_location)

        smoothed_RDF = gaussian_filter(RDF_Data['g(r)'], 2.0)

        plt.figure()
        plt.plot(RDF_Data['r'], RDF_Data['g(r)'], c='b')
        plt.plot(RDF_Data['r'], smoothed_RDF, c='r')
        plt.xlabel('r (Ang)')
        plt.ylabel('RDF (Arb. U.)')
        figure_file = csv_file_location.replace('.csv', '_smooth.pdf')
        plt.savefig(figure_file)
        plt.close()
        peaks = argrelextrema(smoothed_RDF, np.greater)
        try:
            first_peak = [RDF_Data['r'].values[peaks[0][0]],
                          smoothed_RDF[peaks[0][0]]]
            # Update job document with the first peak data
            job.document['RDF_first_peak'] = first_peak
        except IndexError:
            print("No maximum found")
            print("Check", figure_file, "for more details")
            pass


if __name__ == "__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-m", "--metal",
    #                    type=str,
    #                    required=True,
    #                    help='''The transition metal in the M1 catalyst to
    #                    consider when calculating the RDF''')
    #parser.add_argument("-l", "--layers",
    #                    type=int,
    #                    required=True,
    #                    help='''Consider only RDFs that contain this number of
    #                    layers of M1 catalyst in the system''')
    #parser.add_argument("-z", "--z_reactor_size",
    #                    type=float,
    #                    required=True,
    #                    help='''Consider only RDFs that for systems that were
    #                    simulated in a reactor of this size''')
    #args = parser.parse_args()
    project = signac.get_project('../')
    for type2 in ['Mo', 'V', 'Nb', 'Te']:
        for layer in [1, 2, 3]:
            for reactor_size in [10, 15, 20, 25]:
                # Plot RDF variation
                get_first_peak(project, metal=type2, layers=layer,
                               z_size=reactor_size)
