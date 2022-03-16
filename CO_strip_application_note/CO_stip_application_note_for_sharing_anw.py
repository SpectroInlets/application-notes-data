# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 17:27:23 2021

@author: AnnaWiniwarter

This script was used for data treatment and generating figures for Spectro Inlets 
EC-MS application note #11 on CO stripping


"""
import ixdat

from pathlib import Path
from matplotlib import pyplot as plt

from ixdat.techniques.ms import MSMeasurement, MSCalResult 
from ixdat.techniques.ec_ms import ECMSCalibration


# -----------------------EDIT SETTINGS HERE----------------------------------
#Settings for choosing which part of the script is executed:
DATA_SOURCE = "ixdat"
# DATA_SOURCE can be "raw" (for importing EC and MS data from Zilien .tsv file),
# "raw_biologic" (to import MS data from Zilien .tsv files and EC data from 
# BioLogic .mpt files), or "ixdat" (for importing ixdat .csv files)
WHICH_PART = "integrate+calibrate"
# WHICH_PART can be "vs_time", "vs_potential", "integrate+calibrate"
WHICH_REFERENCE = "2nd_cycle"
# WHICH_REFERENCE determines which cycle is used as a baseline cycle for the 
# CO strip. Can be "reference_cycle" (separate measurement), or "2nd_cycle"
# (cycle directly after CO strip) 
SAVE_FIGURES = True
# SAVE_FIGURES set to True: figures will be saved automatically.
FIGURE_TYPE = ".png"
# FIGURE_TYPE can be any format matplotlib accepts for saving figures, eg. 
# ".png", ".svg", ".pdf" etc.

# choose the directory of the raw data.
# note that the filenames of the files to be imported need to be edited below
data_directory = Path.cwd() / r"./data"

# In addition to the above general settings, when using a diffent dataset, 
# the time of when the different measurements occur in the dataset needs to be 
# adjusted in by changing the "tspan" in each section accordingly.  

#----------------------------------------------------------------------------
def main():
    # -----------------------data importing----------------------------------
    # import from the original data files. Note there is no example data included 
    # for this part of the script
    if DATA_SOURCE == "raw": # option 1: import both EC and MS data from Zilien data file
        full_data = ixdat.Measurement.read(data_directory / "myzilienfile.tsv")
        # set RHE potential of reference electrode, electrode SA (geometric), 
        # and ohmic resistance
        full_data.calibrate(RE_vs_RHE=0, A_el=0.196, R_Ohm=0)
        # save the important colums as ixdat-datafile
        full_data.export(data_directory / "full_data_COstrip_11-11-21.csv")
        
    elif DATA_SOURCE == "raw_biologic": # option 2: import MS data from Zilien data file and combine with 
        # EC data from BioLogic file 
        zil_data = MSMeasurement.read(data_directory / "myzilienfile.tsv", reader="zilien")
        ec_data = ixdat.Measurement.read(data_directory / "mybiologicfile.mpt", reader="biologic")
        # set RHE potential of reference electrode, electrode SA (geometric),
        # and ohmic resistance
        ec_data.calibrate(RE_vs_RHE=0, A_el=0.196, R_Ohm=0)
        full_data = zil_data + ec_data
        # save the important colums as ixdat-datafile
        full_data.export(data_directory /"full_data_COstrip_11-11-21.csv")
             
    elif DATA_SOURCE == "ixdat": # option 3: import from ixdat-datafiles
        part1 = ixdat.Measurement.read(data_directory / "data_part1_COstrip_11-11-21.csv", reader="ixdat")    
        part2 = ixdat.Measurement.read(data_directory / "data_part2_COstrip_11-11-21.csv", reader="ixdat")
        full_data = part1 + part2
        
        # calculating the difference between CVs does not (yet) work for ECMSMeasurement
        # objects, therefore import snippets of the pure EC data as directly imported 
        # from (note this should be fixed in ixdat vs 0.1.7 (?) onwards)
        ec_part1 = ixdat.Measurement.read(data_directory / "ecdata_part1_COstrip_11-11-21.csv", reader="ixdat")    
        ec_part2 = ixdat.Measurement.read(data_directory / "ecdata_part2_COstrip_11-11-21.csv", reader="ixdat")
        ec_data = ec_part1 + ec_part2
        
    else:
        raise NameError("DATA_SOURCE not recognized.")
        
    # ------------------- data treatment & plotting ---------------------------
    # ------------------- cut dataset -----------------------------------------
    # CO strip blank (reference measurement)
    co_blank = full_data.cut(tspan=[8420, 10350])  
    co_blank_cv = co_blank.as_cv()
    co_blank_cv.redefine_cycle(start_potential=0.05, redox=False)
    co_blank_cv.tstamp +=8420
    # CO strip
    co_strip = full_data.cut(tspan=[13525, 15438]) 
    co_strip.tstamp +=13525
    co_strip_cv = co_strip.as_cv()
    co_strip_cv.redefine_cycle(start_potential=0.05, redox=False)
    
    # ------------------- integrate and generate plots --------------------------
    if WHICH_PART == "vs_time":  # plot reference measurement and CO strip vs time
        # CO strip blank (reference measurement)
        axes_a = co_blank_cv.plot_measurement(mass_lists=[["M2", "M32"],["M4", "M28", "M44"]], logplot=True, legend=False)
        axes_a[0].set_ylim(1.6e-10, 1e-9)
        axes_a[3].set_ylim(1e-14, 1e-7)
        axes_a[0].set_ylabel("M2, M32 signal / [A]")
        axes_a[3].set_ylabel("M4, M28, M44 signal / [A]")
        coblank_vs_t_fig = axes_a[0].get_figure()
        coblank_vs_t_fig.tight_layout()
        if SAVE_FIGURES is True:
            coblank_vs_t_fig.savefig("co_blank_vs_time_final" + FIGURE_TYPE)
    
        # CO strip
        axes_b = co_strip_cv.plot_measurement(mass_lists=[["M4", "M28", "M44"], ["M2", "M32"]], logplot=True, legend=False)
        axes_b[3].set_ylim(1e-10, 5e-10)
        axes_b[3].set_ylabel("M2, M32 signal / [A]")
        axes_b[0].set_ylabel("M4, M28, M44 signal / [A]")
        axes_b[1].set_yticks([0,0.5, 1, 1.5])
        axes_b[2].set_yticks([-0.025, 0, 0.025, 0.05])
        axes_b[2].set_yticklabels([-25, 0, 25, 50])
        axes_b[2].set_ylabel("J / [$\mu$A cm$^{-2}$]")
        co_strip_vs_t_fig = axes_b[0].get_figure()
        co_strip_vs_t_fig.tight_layout()
        if SAVE_FIGURES is True:
            co_strip_vs_t_fig.savefig("co_strip_vs_time_final" + FIGURE_TYPE)
    
    
    elif WHICH_PART == "vs_potential":
        # plotting CO strip vs potential with either 2nd_cycle or blank cycle 
        stripping_cycle = co_strip.cut(tspan=[500, 1286])
        
        if WHICH_REFERENCE == "reference_cycle":
            base_cycle = co_blank_cv[1] #blank cycle
        elif WHICH_REFERENCE == "2nd_cycle":
            base_cycle = co_strip.cut(tspan=[1280, 1900]) # second cycle of CO strip
        else:
            raise NameError("WHICH_REFERENCE not recognized.")
            
        axes_c = stripping_cycle.plot_vs_potential(mass_lists=[["M2", "M32"],["M44"]], logplot=False, legend=False)
        base_cycle.plot_vs_potential(axes=axes_c, mass_lists=[["M2", "M32"],["M44"]], logplot=False, linestyle=":", legend=False)
        axes_c[0].set_ylabel("M2, M32 signal / [A]")
        axes_c[2].set_ylabel("M44 signal / [A]", color="brown")
        axes_c[0].set_xlabel("$U_{RHE}$ / [V]")
        axes_c[1].set_yticks([-0.025, 0, 0.025])
        axes_c[1].set_yticklabels([-25, 0, 25])
        axes_c[1].set_ylabel("J / [$\mu$A cm$^{-2}$]")
        costrip_vs_u_fig = axes_c[0].get_figure()
        costrip_vs_u_fig.tight_layout()
        if SAVE_FIGURES is True:
            plotname = "CO_strip+2ndcycle_vs_potential"
            costrip_vs_u_fig.savefig(plotname + FIGURE_TYPE)
        
    elif WHICH_PART == "integrate+calibrate":
        # FIRST integrate the ELECTROCHEMICAL CO stripping peak
        ec_co_strip = ec_data.cut(tspan=[12250, 14200])  
        ec_co_strip_cv = ec_co_strip.as_cv()
        ec_co_strip_cv.redefine_cycle(start_potential=0.05, redox=False)
        # ec_co_strip_cv.tstamp +=12250
        
        if WHICH_REFERENCE == "reference_cycle":
            ec_co_blank = ec_data.cut(tspan=[7150, 9070])  
            ec_co_blank_cv = ec_co_blank.as_cv()
            ec_co_blank_cv.redefine_cycle(start_potential=0.05, redox=False)
            # ec_co_blank_cv.tstamp +=7150
            cv_diff = ec_co_strip_cv[1].diff_with(ec_co_blank_cv[1])
            Q_CO_ox_blank = cv_diff.integrate("raw_current", vspan=[0.4, 0.9]) * 1e-3  # 1e-3 converts mC --> C
            n_CO_ox_blank = Q_CO_ox_blank / (ixdat.constants.FARADAY_CONSTANT * 2)
            SA_blank = Q_CO_ox_blank*1e6/340 #muC/cm2  
            print(f"vs blank charge passed = {Q_CO_ox_blank*1e6} uC, corresponding to {n_CO_ox_blank*1e9} nmol of CO oxidized and an ECSA of {SA_blank} cm2")
        
        elif WHICH_REFERENCE == "2nd_cycle":
            cv_diff = ec_co_strip_cv[1].diff_with(ec_co_strip_cv[2])
            Q_CO_ox = cv_diff.integrate("raw_current", vspan=[0.4, 0.9]) * 1e-3  # 1e-3 converts mC --> C
            n_CO_ox = Q_CO_ox / (ixdat.constants.FARADAY_CONSTANT * 2)
            SA = Q_CO_ox*1e6/340 #muC/cm2
            print(f"vs 2nd cycle charge passed = {Q_CO_ox*1e6} uC, corresponding to {n_CO_ox*1e9} nmol of CO oxidized and an ECSA of {SA} cm2")
        
        else:
            raise NameError("WHICH_REFERENCE not recognized.")
        
        # make the difference plot. note:automatically, the color for the first cycle
        # will be set to green. if you want to plot in another color, need to change 
        # that color in your local ixdat ec_plotter.py (and then preferably
        # change it back)        
        ax_cvdiff = cv_diff.plot()
        ax_cvdiff.set_yticks([-0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04 ])
        ax_cvdiff.set_yticklabels([-40, -30, -20, -10, 0, 10, 20, 30, 40])
        ax_cvdiff.set_ylabel("J / [$\mu$A cm$^{-2}$]")
        cvdiff_fig = ax_cvdiff.get_figure()
        if SAVE_FIGURES is True:
            cvdiff_fig.savefig("CV_diff_EC_integration" + FIGURE_TYPE)
                
      
        # SECOND integrate the MS PART of the CO strip (M44 only)        
        axes_co2_strip = co_strip_cv.plot_measurement(mass_list=["M44"], tspan=[500, 2000], logplot=False, legend=False)
        #integrate and define a background
        co2_int_strip_c1 = co_strip_cv[1].integrate_signal('M44', tspan=[700, 1050], tspan_bg=[650, 700], ax=axes_co2_strip[0])
        
        if WHICH_REFERENCE == "reference_cycle":
            axes_co2_blank = co_blank_cv.plot_measurement(mass_list=["M44"], logplot=False, legend=False)
            # integrate and define a background
            co2_int_blank = co_blank_cv[1].integrate_signal('M44', tspan=[750, 1100], tspan_bg=[700, 750], ax=axes_co2_blank[0])
            # check by plotting
            if SAVE_FIGURES is True:
                axes_co2_blank[0].get_figure().savefig("CO_blank_CV_integrated_vs_time" + FIGURE_TYPE)
            # calculate sensitivity factor F
            f_co2 = (co2_int_strip_c1-co2_int_blank)/n_CO_ox_blank
        elif WHICH_REFERENCE == "2nd_cycle":
            co2_int_strip_c2 = co_strip_cv[2].integrate_signal('M44', tspan=[1330, 1680], tspan_bg=[1300, 1330], ax=axes_co2_strip[0])
            f_co2 = (co2_int_strip_c1-co2_int_strip_c2)/n_CO_ox
        else: #this is technically not required here, but left in anyway
            raise NameError("WHICH_REFERENCE not recognized.")
            
        axes_co2_strip[0].set_ylabel("MS signal / [A]")
        axes_co2_strip[2].set_yticks([-0.025, 0, 0.025])
        axes_co2_strip[2].set_yticklabels([-25, 0, 25])
        axes_co2_strip[2].set_ylabel("J / [$\mu$A cm$^{-2}$]")
        if SAVE_FIGURES is True:
            axes_co2_strip[0].get_figure().savefig("CO_strip_CV_integrated_vs_time" + FIGURE_TYPE)
        
        print(f"Sensitivity factor F = {f_co2} using " + WHICH_REFERENCE + " as baseline.")
        
        # now use these sensitivity factors for calibration
        cal_co2 = MSCalResult(
                name="CO2_M44",
                mol="CO2",
                mass="M44",
                F=f_co2,
                cal_type="internal",
            )
        # and make a plot with calibrated CO2 signal vs potential
        co_strip.calibration = ECMSCalibration(ms_cal_results=[cal_co2], RE_vs_RHE=0, A_el=0.196)
        strip_cycle = co_strip.cut(tspan=[500, 1290])
        if WHICH_REFERENCE == "reference_cycle":
            co_blank_cv.calibration = ECMSCalibration(ms_cal_results=[cal_co2], RE_vs_RHE=0, A_el=0.196)    
            bas_cycle = co_blank_cv[1]
        elif WHICH_REFERENCE == "2nd_cycle":
            bas_cycle = co_strip.cut(tspan=[1290, 1900])
        else: #this is technically not required here, but left in anyway
            raise NameError("WHICH_REFERENCE not recognized.")
        
        axes_d = strip_cycle.plot_vs_potential(mol_list=["CO2_M44"], logplot=False, legend=False, tspan_bg=[600, 700])
        bas_cycle.plot_vs_potential(axes=axes_d, mol_list=["CO2_M44"], logplot=False, linestyle=":", legend=False, tspan_bg=[1700, 1800])
        axes_d[0].set_xlabel("$U_{RHE}$ / [V]")
        axes_d[0].set_yticks([0, 1e-12, 2e-12, 3e-12, 4e-12, 5e-12])
        axes_d[0].set_yticklabels([0, 1, 2, 3, 4, 5])
        axes_d[0].set_ylabel("cal. signal / [pmol/s]")
        axes_d[1].set_yticks([-0.025, 0, 0.025])
        axes_d[1].set_yticklabels([-25, 0, 25])
        axes_d[1].set_ylabel("J / [$\mu$A cm$^{-2}$]")
        if SAVE_FIGURES is True:
            plotname = "CO_strip+2ndcycle_vs_potential_calibrated"
            axes_d[0].get_figure().savefig(plotname + FIGURE_TYPE)
        
    else:
        raise NameError("WHICH_PART not recognized.")
            
if __name__ == "__main__":
    main()