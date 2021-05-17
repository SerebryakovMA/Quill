import sys 
from os import path as osp

import numpy as np
import pandas as pd
import qplot
import qplot_energy

def read_data(file_path):
    if not osp.exists(file_path):
        print("File doesn't exist: ", file_path)
        sys.exit(-1)

    data = np.loadtxt(file_path, delimiter=",", skiprows=0, dtype=str)
    header = np.loadtxt(file_path, delimiter=",", max_rows=1, dtype=str)

    return data, header

def energy_plots(data, header):
    for i in data:
        qplot.energy(df="../"+i, clf=True, save2=i+".png")

def treck_length(data, header):
    for i in data:
        if (i[-9:-6] = "not"):
            qplot_energy.energy_level(df="../"+i, species="t_np")
        else:
            qplot_energy.energy_level(df="../"+i, species="t")
    plt.savefig("treck_length_duo.png")


def main():
    data, header = read_data("results_folders.csv")
    #energy_plots(data, header)
    treck_length(data, header)

if __name__ == "__main__":
    main()
