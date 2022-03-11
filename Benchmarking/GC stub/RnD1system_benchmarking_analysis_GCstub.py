# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:30:01 2022

@author: AnnaWiniwarter
"""

import scipy as sp
import numpy as np
import ixdat

from scipy import ndimage
from pathlib import Path
from matplotlib import pyplot as plt
from ixdat.techniques.ms import MSMeasurement
from ixdat.techniques.ec import ECMeasurement

# functions for exponential fitting of decay
# for now copy paste from here: https://stackoverflow.com/questions/3938042/fitting-exponential-decay-with-no-initial-guessing       
def model_func(t, A, K, C):
    return A * np.exp(K * t) + C
    
def fit_exp_linear(t, y, C=0):
    y = y - C
    y = np.log(y)
    K, A_log = np.polyfit(t, y, 1)
    A = np.exp(A_log)
    return A, K

def exp_fit_signal(data, signal, tspan, color=None):
        """
        Function to do an exponential fit on selection of data
        Parameters
        ----------
        data : ECMSMeasurement object (probably ECMeasurement or MSMeasurement 
                                       will also work)
        signal : which part of the ECMSMeasurement object to fit, i.e. column 
        that will then be selected using ECMSMeasurement.grab(signal)
        tspan : ixdat syntax for providing a time span [tstart,tend], time span
        over which the fit will be done
        color : color to use for plotting "signal", for MS signals, standard 
        is selected automatically
        
        Returns
        -------
        t_half, (fig, ax1)
        """
        
        if color is None:
            if "M" in signal:
                color = ixdat.plotters.ms_plotter.STANDARD_COLORS[signal]
            else: 
                color = "k"
        if "M" in signal:
            ylabel = "MS signal / [A]"
        else:
            ylabel = "signal"
                
        t, sig = data.grab(signal)
        #smooth signal (moving average)
        # sig = ndimage.uniform_filter1d(sig,30)
        # t_mask = np.where(np.logical_and((t>tspan[0]), (t<tspan[1])))
        # t_to_fit = t[t_mask]
        # sig_to_fit = sig[t_mask]
        
        data_to_fit = data.cut(tspan=tspan)
        t_to_fit, sig_to_fit = data_to_fit.grab(signal)
              
        # by plotting, check that everything makes sense
        fig, ax1 = plt.subplots()
        ax1.plot(t, sig, linestyle="", marker="o", markersize="3", markerfacecolor="w", markeredgecolor=color ,label=signal+ " signal")
        ax1.plot(t_to_fit, sig_to_fit, label="selected")
        ax1.set_xlim(tspan[0]-10,tspan[1]+20)
        ax1.set_ylim(0, 1.25*np.max(sig[np.where(np.logical_and((t>tspan[0]-10), (t<tspan[1]+20)))]))
        ax1.set_xlabel("time / [s]")
        ax1.set_ylabel(ylabel)
        
        sig_fit_params = fit_exp_linear(t_to_fit, sig_to_fit, C=0)
        sig_fit = model_func(t_to_fit, sig_fit_params[0], sig_fit_params[1], C=0)
        
        ax1.plot(t_to_fit, sig_fit, ":", label="fit")
        ax1.legend()
        
        t_half_sig = np.log(2)/-sig_fit_params[1]
        ax1.annotate(f"t_half={t_half_sig:.2f} s", (0.5,0.5), xycoords="subfigure fraction")
        
        print(signal + " t_half at " + str(time) + "s = " + str(t_half_sig))
        
        return t_half_sig, (fig, ax1)
    
def find_decay_edge(data, signal, gradient_cutoff=None):
        """
        Function to find time where a mass signal starts decaying 
        Parameters
        ----------
        data : ECMSMeasurement object (probably ECMeasurement or MSMeasurement 
                                       will also work)
        signal : which part of the ECMSMeasurement object to fit, i.e. column 
        that will then be selected using ECMSMeasurement.grab(signal)
        
        gradient_cutoff: defines gradient under which mass signal is "decreasing drastically", 
        depends on the absolute value of the signal -> needs to be chosen individually (a 
        good first guess is -1*((order of magnitude of signal)-1))
        
        Returns
        -------
        t_list: list of times where decay starts
        """
        t, m = data.grab(signal)
        # now select the timespan where we want to fit the decay curve
        # the way below is probably not the smartest or most pythonic way, but
        # it seems to work for this type of data
        # define a cutoff for select the range where the mass signal is decreasing drastically 
        # choice of this number is a bit arbitrary, but it seems to work like this
        # needs to be tested with other signal intensities!!
        
        if gradient_cutoff is None:
            gradient_cutoff = -np.max(m)/25
            
        #plot the gradient vs time to determine the cutoff value.
        fig, ax = plt.subplots()
        ax.plot(t, np.gradient(m))
        # ax.plot(t, ndimage.uniform_filter1d(np.gradient(m),30))
        
        mask1 = np.where(np.gradient(m)<gradient_cutoff)
        # mask1 = np.where(ndimage.uniform_filter1d(np.gradient(m),30)<gradient_cutoff)
        # select the range where there is a step in the time because values are removed 
        # by the previous mask
        mask2 = np.where(np.gradient(mask1[0])>1)
        t_list = [t[mask1][0]] + t[mask1][mask2].tolist()
        t_list_clean = t_list[0::2] # necessary because the second mask find 2 times for each step
        
        return t_list_clean

# define a common file name for all generated data files and plots
exp_name = "RnD1system_benchmark_GC"
data_directory = Path.home() / r"C:\Users\AnnaWiniwarter\Github_repo_for_sharing\application-notes-data\Benchmarking\data"

if False: # first time/fresh import from zilien
    full_data = ixdat.Measurement.read(data_directory / r"2022-03-10 08_39_41 GC benchmarking2.tsv", reader="zilien")
    #2022-02-24 09_50_59 benchmarking GC electrode/2022-02-24 09_50_59 benchmarking GC electrode.tsv
    axes_a = full_data.plot_measurement(tspan=[0,10000])
    full_plot= axes_a[0].get_figure()
    full_plot.savefig("./" + exp_name + "full_experiment.png")
    full_data.export("./" + exp_name + ".csv")

if True: # import from ixdat file
    full_data = ixdat.Measurement.read("./" + exp_name + ".csv", reader="ixdat")

if True: # plot the CV part
    cvs = full_data.cut(tspan=[2050,3500])
    cvs.tstamp += 2050
    cvs = cvs.as_cv()
    axes_b = cvs.plot_measurement()
    cvs_vs_time = axes_b[0].get_figure()
    cvs_vs_time.savefig("./" + exp_name + "CVs_vs_time.png")
    # plot one of the CVs vs potential. 
    axes_c = cvs[3].ec_plotter.plot_vs_potential()
    cvs_ec_vs_pot = axes_c.get_figure()
    cvs_ec_vs_pot.savefig("./" + exp_name + "CV_vs_potential_EC.png")
    
    if True: # To only plot averaged (less noisy) EC data, import biologic file directly 
        cvs_ec_only = ECMeasurement.read(data_directory / "2022-03-10 08_39_41 GC benchmarking2_09_01_CVA_DUSB0_C01.mpt", reader="biologic")
        cvs_ec_only = cvs_ec_only.as_cv()
        axes_c = cvs_ec_only[6].plot_vs_potential()
        cvs_ec_vs_pot = axes_c.get_figure()
        cvs_ec_vs_pot.savefig("./" + exp_name + "CV_vs_potential_EC_biologic.png")
    
if True: # plot and fit the HER QC
    her = full_data.cut(tspan=[3350,6030])
    her.tstamp += 3355
    axes_her = her.plot_measurement(tspan=[0,3000])
    fig_her = axes_her[0].get_figure()
    fig_her.tight_layout()
    fig_her.savefig("./" + exp_name + "_HER_CP.png")
    signal =  "M2"
    # signal_tlist = "raw current / [mA]"
    t_list_clean = find_decay_edge(her, signal, gradient_cutoff= -2E-11)  #for M2 gradient_cutoff= -1E-11 for raw current = 0.001
    t_half_h2_list = []        
    for time in t_list_clean:
        t_half_h2, fig = exp_fit_signal(her, signal=signal, tspan=[time, time+5])
        fig[0].savefig("./" + exp_name + "_" + signal + f"decay_at_{time:.0f}s.png")
        t_half_h2_list.append(t_half_h2)
    np.savetxt("./" + exp_name + "_" + signal + "_decay_times.csv", t_half_h2_list,
            delimiter=", ", fmt='%s')
if True: # plot and fit the gas exchange QC
    gas_exchange = full_data.cut(tspan=[6100, 8250])
    gas_exchange.tstamp += 6100
    axes_gas_ex = gas_exchange.ms_plotter.plot_measurement()
    fig_gas_ex = axes_gas_ex.get_figure()
    fig_gas_ex.savefig("./" + exp_name + "_gas_exchange.png")
    times_ar = find_decay_edge(gas_exchange, "M40") #, gradient_cutoff=-1E-10)
    times_he = find_decay_edge(gas_exchange, "M4") #, gradient_cutoff=-1E-10)
    t_half_list = []
    for time in times_ar:
        t_half_ar, fig = exp_fit_signal(gas_exchange, signal="M40", tspan=[time, time+5])
        fig[0].savefig("./" + exp_name + f"_M40_decay_at_{time:.0f}s.png")
        t_half_list.append(("M40", t_half_ar))
    for time in times_he:
        t_half_he, fig = exp_fit_signal(gas_exchange, signal="M4", tspan=[time, time+5])
        fig[0].savefig("./" + exp_name + f"_M4_decay_at_{time:.0f}s.png")
        t_half_list.append(("M4", t_half_he))
    np.savetxt("./" + exp_name + "_gas_exchange_decay_times.csv", t_half_list,
           delimiter =", ", fmt='%s')
    
    