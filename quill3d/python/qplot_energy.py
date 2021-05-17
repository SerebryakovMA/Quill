#!/usr/bin/python

from __future__ import print_function
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import resread
import tcmap
from collections import OrderedDict
import expression_parser


def __get_data_folder(data_folder, kwargs_dict):
    if data_folder is not None:
        return data_folder
    elif 'df' in kwargs_dict:
        return kwargs_dict.pop('df')
    else:
        return None


def tex_format(space_item):
    tmp = space_item
    if space_item in ['ex','ey','ez','bx','by','bz']:
        tmp = space_item[0].upper() + '_' + space_item[1]
    elif space_item in ['ux','uy','uz']:
        tmp = 'p_' + space_item[1]
    elif space_item in ['vx','vy','vz']:
        tmp = space_item[0] + '_' + space_item[1]
    elif space_item in ['xi','chi','theta']:
        tmp = '\\' + space_item
    elif space_item == 'g':
        tmp = '\\gamma'
    elif space_item == 'phi':
        tmp = '\\varphi'
    return '$' + tmp + '$'


def energy_level(data_folder=None, save2=None, include_deleted=True, clf=False, species='epgiwt', **kwargs):
    'Plots time of energy decrease of electrons, ions, etc. to mathched vs time'
    'species - what should be plotted: e, p, g, i - particles, w - field energy, t - total'
    #Important: several types of ions are not supported

    def mid_value_ind(arr):
        abs_val_arr = np.abs([a - 0.5*arr[0] for a in arr])
        smallest_difference_index = abs_val_arr.argmin()
        return smallest_difference_index
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    tmp = resread.t_data('energy')
    if clf:
        plt.clf()
    #total_energy = np.sum(tmp[:,1:5], axis=1)
    if 't' in species:
        plt.plot(tmp[mid_value_ind(tmp[:,1])], tmp[mid_value_ind(tmp[:,1]),0], ".", color="g")
    if "t_np" in species:
        plt.plot(tmp[mid_value_ind(tmp[:,1])], tmp[mid_value_ind(tmp[:,1]),0], ".", color="b") 

    if (resread.catching or (resread.dump_photons and 'g' in species)) and include_deleted:
        # deleted energy
        tmp_del = resread.t_data('energy_deleted')
        if 'e' in species:
            sum_e, min_len = safe_sum(tmp_del[:,1], tmp[:,2])
            plt.plot(tmp_del[:min_len,0], sum_e, '--g') # electrons
        if 'p' in species:
            sum_p, min_len = safe_sum(tmp_del[:,2], tmp[:,3])
            plt.plot(tmp_del[:min_len,0], sum_p, '--r') # positrons
        if 'g' in species:
            sum_g, min_len = safe_sum(tmp_del[:,3], tmp[:,4])
            plt.plot(tmp_del[:min_len,0], sum_g, '--b') # hard photons
        
        total_del_energy, min_len = safe_sum(total_energy, np.sum(tmp_del[:,1:4], axis=1))
        if (resread.n_ion_populations>0):
            total_del_energy, min_len = safe_sum(total_del_energy, tmp_del[:,4])
            sum_i, min_len = safe_sum(tmp_del[:,4], tmp[:,5])
            if 'i' in species:
                plt.plot(tmp_del[:min_len,0], sum_i, '--m') # ions
        if 't' in species:
            plt.plot(tmp[:min_len,0], total_del_energy, ':k') # sum energy
    
    plt.ylabel('$ct/\lambda$')
    plt.xlabel('Energy, J')

    if save2 is not None:
        plt.savefig(save2)

def N(data_folder=None, particles='gep', save2=None, clf=False, **kwargs):
    'Plots number of particles over time. Ions are not currently supported'
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    a = resread.t_data('N')
    t = a[:,0]
    Ne = a[:,1]
    Np = a[:,2]
    Ng = a[:,3]

    lw = 1.5 # linewidth
    if clf:
        plt.clf()
    if 'e' in particles:
        plt.plot(t, Ne, 'g--', linewidth = lw, label = r'$N_e$')
    if 'p' in particles:
        plt.plot(t, Np, 'r-', linewidth = lw, label = r'$N_p$')
    if 'g' in particles:
        plt.plot(t, Ng, 'b:', linewidth = lw, label = r'$N_g$')

    plt.legend(loc = 'best', fontsize = 'medium')
    plt.xlabel(r'$ct/\lambda$')
    plt.ylabel(r'$N$')
    plt.yscale('log')

    if save2 is not None:
        plt.savefig(save2)


def onaxis(t, particles='we', colors='rgbcmyk', norm='true',
        data_folder=None, save2=None, plotargs={}, rrargs={}, clf=False, **kwargs):
    'Plot of particle density, fields etc. along the x-axis.\n\
    *norm* can be set to \'optimal\', \'true\' or to array of desired maximal values.\n\
    For rrargs see help(qplot.resread.onaxis),\n\
    for plotargs see help(plt.plot).\n\
    \n\
    Examples:\n\
    qplot.onaxis(0, \'e\'),\n\
    qplot.onaxis(0, \'ep\', norm = \'true\'),\n\
    qplot.onaxis(0, \'we\', norm = [0.5, 1]),\n\
    qplot.onaxis(15,[\'ey\', \'bz\', \'e\'], \'rgb\', plotargs = {linewidth: 1.2})\n\
    qplot.onaxis(4, \'g\', rrargs = {\'av\': \'y\'}),\n\
    qplot.onaxis(4, \'p\', rrargs = {\'sx\': 10, \'sz\': 3}).'
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    x = resread.onaxis('x', **rrargs)
    a = []
    for p in particles:
        if p == 'e':
            filename = 'rho'
        elif p == 'g':
            filename = 'rho_ph'
        elif p == 'p':
            filename = 'rho_p'
        elif p == 'i':
            filename = 'irho_' + str(resread.icmr[0]) + '_'
        else:
            filename = p
        resread.t = '%g' % t
        if p == 'e':
            a.append(-resread.onaxis(filename, **rrargs))
        else:
            a.append(resread.onaxis(filename, **rrargs))
    ma = 0
    for i, b in enumerate(a):
        mb = max(b)
        print('max value for', particles[i], '=', mb)
        if mb > ma:
            ma = mb
    if norm == 'optimal':
        for i, b in enumerate(a):
            mb = max(b)
            if mb < ma / 4:
                a[i] = a[i] / mb * ma / 4
    elif norm != 'true':
        for i, m in enumerate(norm):
            a[i] = norm[i] * a[i] / max(a[i])

    if clf:
        plt.clf()
    for i, b in enumerate(a):
        plt.plot(x, b, color = colors[i % len(colors)], **plotargs)
    plt.xlabel('$x$')
    if save2 is not None:
        plt.savefig(save2)


def mwcoordinate(data_folder=None, save2=None, clf=False, **kwargs):
    """
    Plots the time dependence of the coordinate of the moving window

    Parameters
    ----------
    data_folder or df : path to folder with data
    save2 : file to save the image to
    clf : clear figure before plotting

    """
    data_folder = __get_data_folder(data_folder, kwargs)
    if data_folder is not None:
        resread.data_folder = data_folder
    resread.read_parameters()
    a = resread.t_data('mwcoordinate')
    t = a[:, 0]
    x = a[:, 1]

    linewidth = 1.5
    if clf:
        plt.clf()
    plt.plot(t, x, '-', linewidth=linewidth)

    plt.xlabel(r'$ct/\lambda$')
    plt.ylabel(r'$x/\lambda$')

    if save2 is not None:
        plt.savefig(save2)


def ion():
    """
    Invokes 'plt.ion()'
    """
    plt.ion()


def ioff():
    """
    Invokes 'plt.ioff()'
    """
    plt.ioff()


def reset_style():
    mpl.rcParams.update(rc_backup) 

import __main__
if not plt.isinteractive() and not hasattr(__main__, '__file__'):
    print("For the interactive regime use IPython in the Pylab mode ('ipython --pylab' or '%pylab'), 'plt.ion()', "
          "or 'qplot.ion()'")

lw = 0.7 # linewidth for border
lwl = 1.0 # linewidth for lines in plots
font = {'family' : 'serif', 'serif' : 'cmr10', 'size' : 9}
rc_backup = mpl.rcParams.copy()
mpl.rc('font', **font)
mpl.rc('lines', linewidth=lwl)
mpl.rc('axes', linewidth=lw)
mpl.rc('axes', unicode_minus=False)
mpl.rc('figure', figsize=(3.5,2.5), dpi=200, autolayout=True)
mpl.rc('mathtext' ,fontset='cm')
mpl.rc('savefig',dpi=300)

tcmap.red()
tcmap.green()
tcmap.blue()
tcmap.orange()
tcmap.purple()
