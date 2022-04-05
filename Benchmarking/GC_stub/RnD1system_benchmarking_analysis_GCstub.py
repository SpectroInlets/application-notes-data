# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:30:01 2022

@author: AnnaWiniwarter
"""
  
import numpy as np
import ixdat

from pathlib import Path
from matplotlib import pyplot as plt
from ixdat.techniques.ms import MSMeasurement
from ixdat.techniques.ec import ECMeasurement


# -----------------------EDIT SETTINGS HERE----------------------------------
#Settings for choosing which part of the script is executed:
DATA_SOURCE = "ixdat"
# DATA_SOURCE can be "raw" (for importing EC and MS data from Zilien .tsv file),
# or "ixdat" (for importing ixdat .csv files as generated when using setting "raw")
# in either case, for plotting the CV only, the biologic CVA file is imported
# automatically
WHICH_PART = "plot_CVs"
# WHICH_PART can be "plot_CVs", "plot+fit_HER", "plot+fit_gas_exchange"

SAVE_FIGURES = True
# SAVE_FIGURES set to True: figures will be saved automatically.
FIGURE_TYPE = ".png"
# FIGURE_TYPE can be any format matplotlib accepts for saving figures, eg. 
# ".png", ".svg", ".pdf" etc.

# select the directory of the raw data.
EXP_NAME = "RnD1system_benchmark_GC"
THIS_DIR = Path(__file__).parents[1].resolve()
DATA_DIRECTORY = THIS_DIR / "data"

# select the filenames of the raw data.
ZILIEN_FILENAME = "2022-03-10 08_39_41 GC benchmarking2.tsv"
BIOLOGIC_FILENAME = "2022-03-10 08_39_41 GC benchmarking2_09_01_CVA_DUSB0_C01.mpt"


# In addition to the above general settings, when using a diffent dataset, 
# the time of when the different measurements occur in the dataset needs to be 
# adjusted in by changing the "tspan" in each section accordingly. To this end
# it is recommended to use full_data.plot_measurement(tspan=[t1,t2]) after importing
# and varying t1 and t2 to find the tspans of interest. (see README.md)
# ----------------------- END OF EDIT SETTINGS ------------------------------


# ----------- functions for exponential fitting of decay --------------------
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
    data_to_fit = data.cut(tspan=tspan)
    t_to_fit, sig_to_fit = data_to_fit.grab(signal)
    
    # by plotting, check that everything makes sense
    fig, ax1 = plt.subplots()
    ax1.plot(t, sig, linestyle="", marker="o", markersize="3", markerfacecolor="w",
             markeredgecolor=color, label=signal+ " signal")
    ax1.plot(t_to_fit, sig_to_fit, label="selected")
    ax1.set_xlim(tspan[0]-10,tspan[1]+20)
    ax1.set_xlabel("time / [s]")
    ax1.set_ylabel(ylabel)
    
    sig_fit_params = fit_exp_linear(t_to_fit, sig_to_fit, C=0)
    sig_fit = model_func(t_to_fit, sig_fit_params[0], sig_fit_params[1], C=0)
    
    ax1.plot(t_to_fit, sig_fit, ":", label="fit")
    ax1.legend()
    
    t_half_sig = np.log(2)/-sig_fit_params[1]
    ax1.annotate(f"t_half={t_half_sig:.2f} s", (0.5,0.5), xycoords="subfigure fraction")
    
    print(signal + " t_half at " + str(t_to_fit[0]) + "s = " + str(t_half_sig))
    
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
    # choice of this number is a bit arbitrary and doesn't work reliably
    if gradient_cutoff is None:
        gradient_cutoff = -np.max(m)/25
        
    #plot the gradient vs time to determine the cutoff value. This figure is not saved.
    fig, ax = plt.subplots()
    ax.plot(t, np.gradient(m))
            
    mask1 = np.where(np.gradient(m)<gradient_cutoff)
    # select the range where there is a step in the time because values are removed 
    # by the previous mask
    mask2 = np.where(np.gradient(mask1[0])>1)
    t_list = [t[mask1][0]] + t[mask1][mask2].tolist()
    t_list_clean = t_list[0::2] # necessary because the second mask find 2 times for each step
    
    return t_list_clean


def main():
    if DATA_SOURCE == "raw": # option 1: import both EC and MS data from Zilien data file
        print(DATA_DIRECTORY / ZILIEN_FILENAME)
        full_data = ixdat.Measurement.read(DATA_DIRECTORY / ZILIEN_FILENAME, reader="zilien")
        
        
        axes_a = full_data.plot_measurement(tspan=[0,10000])
        full_plot= axes_a[0].get_figure()
        if SAVE_FIGURES is True:
            full_plot.savefig("./" + EXP_NAME + "full_experiment" + FIGURE_TYPE)
        full_data.export("./" + EXP_NAME + ".csv")
    
    elif DATA_SOURCE == "ixdat":  # option 2: import from ixdat-datafiles
        full_data = ixdat.Measurement.read("./" + EXP_NAME + ".csv", reader="ixdat")
    else:
        raise NameError("DATA_SOURCE not recognized.")
    
    if WHICH_PART == "plot_CVs": # plot the CV part
        cvs = full_data.cut(tspan=[2050,3500])
        cvs.tstamp += 2050
        cvs = cvs.as_cv()
        cvs.redefine_cycle(start_potential=0.9, redox=False)
        axes_b = cvs.plot_measurement()
        cvs_vs_time = axes_b[0].get_figure()
        if SAVE_FIGURES is True:
            cvs_vs_time.savefig("./" + EXP_NAME + "CVs_vs_time" + FIGURE_TYPE)
        # plot one of the CVs (the 4th out of 5) vs potential. 
        axes_c = cvs[3].ec_plotter.plot_vs_potential()
        cvs_ec_vs_pot = axes_c.get_figure()
        if SAVE_FIGURES is True:
            cvs_ec_vs_pot.savefig("./" + EXP_NAME + "CV_vs_potential_EC" + FIGURE_TYPE)
        
        # To instead plot averaged (less noisy) EC data, import biologic file directly 
        cvs_ec_only = ECMeasurement.read(DATA_DIRECTORY / BIOLOGIC_FILENAME, reader="biologic")
        cvs_ec_only = cvs_ec_only.as_cv()
        axes_c = cvs_ec_only[3].plot_vs_potential()
        cvs_ec_vs_pot = axes_c.get_figure()
        if SAVE_FIGURES is True:
            cvs_ec_vs_pot.savefig("./" + EXP_NAME + "CV_vs_potential_EC_biologic" + FIGURE_TYPE)
        
    elif WHICH_PART == "plot+fit_HER": # plot and fit the HER QC
        her = full_data.cut(tspan=[3350,6030])
        her.tstamp += 3355
        axes_her = her.plot_measurement(tspan=[0,3000])
        fig_her = axes_her[0].get_figure()
        fig_her.tight_layout()
        if SAVE_FIGURES is True:
            fig_her.savefig("./" + EXP_NAME + "_HER_CP.png")
        signal = "M2"
        t_list_clean = find_decay_edge(her, signal, gradient_cutoff= -2E-11)  
        # if automatically generated
        # list doesn't catch the right times, cutoff value can be selected manually: 
        # gradient_cutoff=-2E-11 (this value works for HER for the tuning used 
        # in this example dataset.)
        t_half_h2_list = []        
        for time in t_list_clean:
            t_half_h2, fig = exp_fit_signal(her, signal=signal, tspan=[time, time+5])
            if SAVE_FIGURES is True:
                fig[0].savefig("./" + EXP_NAME + "_" + signal + f"decay_at_{time:.0f}s." + FIGURE_TYPE)
            t_half_h2_list.append(t_half_h2)
        np.savetxt("./" + EXP_NAME + "_" + signal + "_decay_times.csv", t_half_h2_list,
               delimiter=", ", fmt='%s')
    
    elif WHICH_PART == "plot+fit_gas_exchange": # plot and fit the gas exchange QC
        gas_exchange = full_data.cut(tspan=[6100, 8250])
        gas_exchange.tstamp += 6100
        axes_gas_ex = gas_exchange.ms_plotter.plot_measurement()
        fig_gas_ex = axes_gas_ex.get_figure()
        if SAVE_FIGURES is True:
            fig_gas_ex.savefig("./" + EXP_NAME + "_gas_exchange" + FIGURE_TYPE)
        times_ar = find_decay_edge(gas_exchange, "M40") #, gradient_cutoff=-1E-10)
        times_he = find_decay_edge(gas_exchange, "M4") #, gradient_cutoff=-1E-10)
        t_half_list = []
        for time in times_ar:
            t_half_ar, fig = exp_fit_signal(gas_exchange, signal="M40", tspan=[time, time+5])
            if SAVE_FIGURES is True:
                fig[0].savefig("./" + EXP_NAME + f"_M40_decay_at_{time:.0f}s" + FIGURE_TYPE)
            t_half_list.append(("M40", t_half_ar))
        for time in times_he:
            t_half_he, fig = exp_fit_signal(gas_exchange, signal="M4", tspan=[time, time+5])
            if SAVE_FIGURES is True:
                fig[0].savefig("./" + EXP_NAME + f"_M4_decay_at_{time:.0f}s" + FIGURE_TYPE)
            t_half_list.append(("M4", t_half_he))
        np.savetxt("./" + EXP_NAME + "_gas_exchange_decay_times.csv", t_half_list,
               delimiter =", ", fmt='%s')
        
    else:
        raise NameError("WHICH_PART not recognized.")
        
    return full_data

        
if __name__ == "__main__":
    main()    
    