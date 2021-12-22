# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:27:23 2021

@author: AnnaWiniwarter

This script was used for data treatment and generating figures for Spectro Inlets 
EC-MS application note #11 on CO stripping


"""
import ixdat as ix

from pathlib import Path
from matplotlib import pyplot as plt

from ixdat.techniques.ms import MSMeasurement
from ixdat.techniques.ec_ms import ECMSCalibration


#To import raw data from Zilien or Zilien & BioLogic data files, set the
#following if-clauses to True and enter your datapath and filenames.
#Note, that there is no example data included for this part of the script

#choose the directory of the raw data.
data_directory = Path.cwd() / r"./data"

#data importing
#import from the original data files
    
if False: #option 1: import both EC and MS data from Zilien data file
    full_data = ix.Measurement.read(data_directory / "myzilienfile.tsv")
    #set RHE potential of reference electrode, electrode SA (geometric), 
    #and ohmic resistance
    full_data.calibrate(RE_vs_RHE=0, A_el=0.196, R_Ohm=0)
    #save the important colums as ixdat-datafile
    full_data.export(data_directory / "full_data_COstrip_11-11-21.csv")
    
if False: #option 2: import MS data from Zilien data file and combine with 
    #EC data from BioLogic file 
    zil_data = MSMeasurement.read(data_directory / "myzilienfile.tsv", reader="zilien")
    ec_data = ix.Measurement.read(data_directory / "mybiologicfile.mpt", reader="biologic")
    #set RHE potential of reference electrode, electrode SA (geometric),
    #and ohmic resistance
    ec_data.calibrate(RE_vs_RHE=0, A_el=0.196, R_Ohm=0)
    full_data = zil_data + ec_data
    #save the important colums as ixdat-datafile
    full_data.export(data_directory /"full_data_COstrip_11-11-21.csv")
    
if True:  
    #import from ixdat-datafiles   
    part1 = ix.Measurement.read(data_directory / "data_part1_COstrip_11-11-21.csv", reader="ixdat")    
    part2 = ix.Measurement.read(data_directory / "data_part2_COstrip_11-11-21.csv", reader="ixdat")
    
    full_data = part1 + part2
    
    #calculating the difference between CVs does not (yet) work for ECMSMeasurement
    #objects, therefore import snippets of the pure EC data as directly imported 
    #from 
    ec_part1 = ix.Measurement.read(data_directory / "ecdata_part1_COstrip_11-11-21.csv", reader="ixdat")    
    ec_part2 = ix.Measurement.read(data_directory / "ecdata_part2_COstrip_11-11-21.csv", reader="ixdat")
    
    ec_data = ec_part1 + ec_part2
    
if True: 
    #plot reference measurement and CO strip vs time
    
    #CO strip blank (reference measurement)
    co_blank = full_data.cut(tspan=[8420,10350])  
    co_blank_cv = co_blank.as_cv()
    co_blank_cv.redefine_cycle(start_potential=0.05, redox=False) #this seems not to be optional
    co_blank_cv.tstamp +=8420
    axes_a = co_blank_cv.plot_measurement(mass_lists=[["M2","M32"],["M4","M28","M44"]],logplot=True, legend=False)
    axes_a[0].set_ylim(1.6e-10,1e-9)
    axes_a[3].set_ylim(1e-14,1e-7)
    axes_a[0].set_ylabel("M2, M32 signal / [A]")
    axes_a[3].set_ylabel("M4, M28, M44 signal / [A]")
    coblank_vs_t_fig = axes_a[0].get_figure()
    coblank_vs_t_fig.tight_layout()
    # coblank_vs_t_fig.savefig("co_blank_vs_time_final.svg")

    #CO strip
    co_strip = full_data.cut(tspan=[13525,15438]) 
    co_strip.tstamp +=13525
    co_strip_cv = co_strip.as_cv()
    co_strip_cv.redefine_cycle(start_potential=0.05, redox=False) #this seems not to be optional
    axes_b = co_strip_cv.plot_measurement(mass_lists=[["M4","M28","M44"], ["M2","M32"]],logplot=True, legend=False)
    axes_b[3].set_ylim(1e-10,5e-10)
    axes_b[3].set_ylabel("M2, M32 signal / [A]")
    axes_b[0].set_ylabel("M4, M28, M44 signal / [A]")
    axes_b[1].set_yticks([0,0.5,1,1.5])
    axes_b[2].set_yticks([-0.025, 0, 0.025,0.05])
    axes_b[2].set_yticklabels([-25, 0, 25, 50])
    axes_b[2].set_ylabel("J / [$\mu$A cm$^{-2}$]")
    co_strip_vs_t_fig = axes_b[0].get_figure()
    co_strip_vs_t_fig.tight_layout()
    # co_strip_vs_t_fig.savefig("co_strip_vs_time_final.svg")


if True:
    #plotting CO strip vs potential with either 2nd cycle or blank cycle 
    #(un)comment to select
    stripping_cycle = co_strip.cut(tspan=[500,1286])
    # base_cycle = co_blank_cv[1] #blank cycle
    base_cycle = co_strip.cut(tspan=[1280,1900]) #second cycle of CO strip
    
    axes_c = stripping_cycle.plot_vs_potential(mass_lists=[["M2","M32"],["M44"]], logplot=False, legend=False)
    base_cycle.plot_vs_potential(axes=axes_c, mass_lists=[["M2","M32"],["M44"]], logplot=False, linestyle=":", legend=False)
    axes_c[0].set_ylabel("M2, M32 signal / [A]")
    axes_c[2].set_ylabel("M44 signal / [A]", color="brown")
    axes_c[0].set_xlabel("$U_{RHE}$ / [V]")
    axes_c[1].set_yticks([-0.025, 0, 0.025])
    axes_c[1].set_yticklabels([-25, 0, 25])
    axes_c[1].set_ylabel("J / [$\mu$A cm$^{-2}$]")
    costrip_vs_u_fig = axes_c[0].get_figure()
    costrip_vs_u_fig.tight_layout()
    plotname = "CO_strip+2ndcycle_vs_potential"
    # costrip_vs_u_fig.savefig(plotname + ".png")
    # costrip_vs_u_fig.savefig(plotname + ".svg")
    
if True:
    #integrate the electrochemical CO stripping peak
    ec_co_blank = ec_data.cut(tspan=[7150, 9070])  
    ec_co_blank_cv = ec_co_blank.as_cv()
    ec_co_blank_cv.redefine_cycle(start_potential=0.05, redox=False) #this is NOT optional!
    # ec_co_blank_cv.tstamp +=7150
    
    ec_co_strip = ec_data.cut(tspan=[12250, 14200])  
    ec_co_strip_cv = ec_co_strip.as_cv()
    ec_co_strip_cv.redefine_cycle(start_potential=0.05, redox=False) #this is NOT optional!
    # ec_co_strip_cv.tstamp +=12250
    
    cv_diff_blank = ec_co_strip_cv[1].diff_with(ec_co_blank_cv[1])
    cv_diff = ec_co_strip_cv[1].diff_with(ec_co_strip_cv[2])
    
    #make the difference plot. note:automatically, the color for the first cycle
    #will be set to green. if you want to plot in another color, need to change 
    #that color in your local ixdat ec_plotter.py (and then preferably
    #change it back)
    
    ax_cvdiff = cv_diff.plot()
    ax_cvdiff.set_yticks([-0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04 ])
    ax_cvdiff.set_yticklabels([-40, -30, -20, -10, 0, 10, 20, 30, 40])
    ax_cvdiff.set_ylabel("J / [$\mu$A cm$^{-2}$]")
    cvdiff_fig = ax_cvdiff.get_figure()
    # cvdiff_fig.savefig("CV_diff_EC_integration.svg")
    
    cv_diff.plot_diff()
    
    Q_CO_ox = cv_diff.integrate("raw_current", vspan=[0.4, 0.9]) * 1e-3  # 1e-3 converts mC --> C
    n_CO_ox = Q_CO_ox / (ix.constants.FARADAY_CONSTANT * 2)
    SA = Q_CO_ox*1e6/340 #muC/cm2
    
    Q_CO_ox_blank = cv_diff_blank.integrate("raw_current", vspan=[0.4, 0.9]) * 1e-3  # 1e-3 converts mC --> C
    n_CO_ox_blank = Q_CO_ox_blank / (ix.constants.FARADAY_CONSTANT * 2)
    SA_blank = Q_CO_ox_blank*1e6/340 #muC/cm2

    print(f"vs blank charge passed = {Q_CO_ox_blank*1e6} uC, corresponding to {n_CO_ox_blank*1e9} nmol of CO oxidized and an ECSA of {SA_blank} cm2")
    print(f"vs 2nd cycle charge passed = {Q_CO_ox*1e6} uC, corresponding to {n_CO_ox*1e9} nmol of CO oxidized and an ECSA of {SA} cm2")

if True:
    #integrate the MS part of the CO strip (M44 only)

    axes_co2_blank = co_blank_cv.plot_measurement(mass_list=["M44"], logplot=False, legend=False)
    #integrate and define a background
    co2_int_blank = co_blank_cv[1].integrate_signal('M44', tspan=[750,1100], tspan_bg=[700,750], ax=axes_co2_blank[0])
    #check by plotting
    # axes_co2_blank[0].get_figure().savefig("CO_blank_CV_integrated_vs_time.svg")
    
    axes_co2_strip = co_strip_cv.plot_measurement(mass_list=["M44"], tspan=[500, 2000], logplot=False, legend=False)
    #integrate and define a background
    co2_int_strip_c1 = co_strip_cv[1].integrate_signal('M44', tspan=[700,1050], tspan_bg=[650,700], ax=axes_co2_strip[0])

    # axes_co2_strip[0].get_figure()
    
    co2_int_strip_c2 = co_strip_cv[2].integrate_signal('M44', tspan=[1330,1680], tspan_bg=[1300,1330], ax=axes_co2_strip[0])

    axes_co2_strip[0].set_ylabel("MS signal / [A]")
    axes_co2_strip[2].set_yticks([-0.025, 0, 0.025])
    axes_co2_strip[2].set_yticklabels([-25, 0, 25])
    axes_co2_strip[2].set_ylabel("J / [$\mu$A cm$^{-2}$]")
    # axes_co2_strip[0].get_figure().savefig("CO_strip_CV_integrated_vs_time.svg")
    
    #calculate sensitivity factor F
    F_blank = (co2_int_strip_c1-co2_int_blank)/n_CO_ox_blank
    
    F_2ndcycle = (co2_int_strip_c1-co2_int_strip_c2)/n_CO_ox
    
    print(f"Sensitivity factor F = {F_blank} when subtracting blank cycle or F = {F_2ndcycle} when using 2nd cycle of strip.")
    
    #now use these sensitivity factors for calibration
    from ixdat.techniques.ms import MSCalResult
    
    cal_co2 = MSCalResult(
            name="CO2_M44",
            mol="CO2",
            mass="M44",
            F=F_blank,
            cal_type="internal",
        )
    
if True:
    calibration = cal_co2
    
    co_strip.calibration = ECMSCalibration(ms_cal_results=[calibration], RE_vs_RHE=0, A_el=0.196)

    #plotting calibrated CO strip vs potential 
    strip_cycle = co_strip.cut(tspan=[500,1290])
    # base_cycle = co_blank_cv[1]
    bas_cycle = co_strip.cut(tspan=[1290,1900])
    
    axes_d = strip_cycle.plot_vs_potential(mol_list=["CO2_M44"], logplot=False, legend=False, tspan_bg=[600,700])
    bas_cycle.plot_vs_potential(axes=axes_d, mol_list=["CO2_M44"], logplot=False, linestyle=":", legend=False, tspan_bg=[1700,1800])
    axes_d[0].set_xlabel("$U_{RHE}$ / [V]")
    axes_d[0].set_yticks([0, 1e-12, 2e-12, 3e-12, 4e-12, 5e-12])
    axes_d[0].set_yticklabels([0, 1, 2, 3, 4, 5])
    axes_d[0].set_ylabel("cal. signal / [pmol/s]")
    axes_d[1].set_yticks([-0.025, 0, 0.025])
    axes_d[1].set_yticklabels([-25, 0, 25])
    axes_d[1].set_ylabel("J / [$\mu$A cm$^{-2}$]")
    co_strip_vs_t_fig.tight_layout()
    plotname = "CO_strip+2ndcycle_vs_potential_calibrated"
    # axes_d[0].get_figure().savefig(plotname + ".png")
    # axes_d[0].get_figure().savefig(plotname + ".svg")
        