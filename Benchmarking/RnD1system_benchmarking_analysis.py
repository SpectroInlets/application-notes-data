# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:30:01 2022

@author: AnnaWiniwarter
"""

import scipy as sp
import numpy as np
import ixdat

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

# define a common file name for all generated data files and plots
exp_name = "RnD1system_benchmark"

if True: # first time/fresh import from zilien
    data_directory = Path.home() / r"C:\Users\AnnaWiniwarter\Dropbox (Spectro Inlets)\Development\Data\New Benchmarking Procedure"
        
    full_data = ixdat.Measurement.read(data_directory / "2022-01-31 11_11_16 test benchmarking 3/2022-01-31 11_11_16 test benchmarking 3.tsv", reader="zilien")
    
    axes_a = full_data.plot_measurement(tspan=[0,10000])
    full_plot= axes_a[0].get_figure()
    full_plot.savefig("./" + exp_name + "full_experiment.png")
    full_data.export("./" + exp_name + ".csv")

if False: # importy from ixdat file
    full_data = ixdat.Measurement.read("./" + exp_name + "csv", reader="ixdat")

if True: # plot the CV part
    cvs = full_data.cut(tspan=[400,1400])
    cvs.tstamp += 480
    cvs = cvs.as_cv()
    axes_b = cvs.plot_measurement()
    cvs_vs_time = axes_b[0].get_figure()
    cvs_vs_time.savefig("./" + exp_name + "CVs_vs_time.png")
    
    # plot one of the CVs vs potential. To only plot EC data, import biologic file directly 
    # TODO: rewrite this to use the EC plotter for zilien data: ixdat.plotters.ECPlotter(measurement=mymeasurement).plot_measurement()
    cvs_ec_only = ECMeasurement.read(data_directory / "2022-01-31 11_11_16 test benchmarking 3/2022-01-31 11_11_16 test benchmarking 3_01_01_CVA_DUSB0_C01.mpt", reader="biologic")
    cvs_ec_only = cvs_ec_only.as_cv()
    axes_c = cvs_ec_only[3].plot_vs_potential()
    cvs_ec_vs_pot = axes_c.get_figure()
    cvs_ec_vs_pot.savefig("./" + exp_name + "CV_vs_potential_EC.png")
    
    
if False: # plot and fit the HER QC
    her = full_data.cut(tspan=[1500,4000])
    her.tstamp += 1589
    her.plot_measurement(tspan=[0,3000])
    
    tm2_her, m2_her = her.grab("M2")
    i_her = her.grab_for_t("I/mA", tm2_her)
    
    # now select the timespan where we want to fit the decay curve
    # the way below is probably not the smartest or most pythonic way, but hopefully
    # it will work
    # select the range where the m2 signal is decreasing drastically (number chosen is a bit arbitrary)
    mask1 = np.where(np.gradient(m2_her)<-1E-11)
    # select the range where there is a step in the time because values are removed 
    # by the previous mask
    mask2 = np.where(np.gradient(mask1[0])>1)
    
    t_list = [tm2_her[mask1][0]] + tm2_her[mask1][mask2].tolist()
    t_list_clean = t_list[0::2]    
    t_half_h2_list = []    
    
    for time in t_list_clean:    
        hydrogen=her.cut(tspan=[time, time+5])
        tm2_dec, m2_dec = hydrogen.grab('M2')
        
        # check that everything makes sense
        fig, ax1 = plt.subplots()
        ax1.plot(tm2_her, m2_her, linestyle="", marker="o", markersize="3", markerfacecolor="w", markeredgecolor="b",label="H2 signal")
        ax1.plot(tm2_dec, m2_dec, label="selected")
        ax1.set_xlim(time-10,time+20)
        ax1.set_xlabel("time / s")
        ax1.set_ylabel("MS signal / A")
        
        h2_fit_params = fit_exp_linear(tm2_dec, m2_dec, C=0)
        h2_fit = model_func(tm2_dec, h2_fit_params[0], h2_fit_params[1], C=0)
        
        ax1.plot(tm2_dec, h2_fit, ":", label="fit")
        ax1.legend()
        
        t_half_h2 = np.log(2)/-h2_fit_params[1]
        t_half_h2_list.append(t_half_h2)
        ax1.annotate(f"t_half={t_half_h2:.2f} s", (0.5,0.5), xycoords="subfigure fraction")
        
        print("Hydrogen T_half at " + str(time) + "s = " + str(t_half_h2))
    
    
if False: # plot and fit the gas exchange QC
    first_ms = MSMeasurement.read(data_directory / "2022-01-31 11_11_16 test benchmarking 3/2022-01-31 11_11_16 test benchmarking 3.tsv", reader="zilien")
    gas_exchange = first_ms.cut(tspan=[5250, 6200])
    gas_exchange.tstamp += 5250
    gas_exchange.plot_measurement()