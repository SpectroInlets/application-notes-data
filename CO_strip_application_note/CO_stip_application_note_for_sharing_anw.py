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
# DATA_SOURCE can be "example" (for importing the ixdat .csv files 
# to reproduce the figures from the Spectro Inlets CO-Strip Application Aote), 
# "raw" (for importing EC and MS data from Zilien .tsv file),
# "raw_biologic" (to import MS data from Zilien .tsv files and EC data from 
# BioLogic .mpt files), or "ixdat" (for importing ixdat .csv files exported 
# with the settings "raw" or "raw_biologic")
WHICH_PART = "vs_time"
# WHICH_PART can be "vs_time", "vs_potential", "integrate+calibrate". Note that 
# the combination of DATA_SOURCE="raw" and "integrate+calibrate" is not possible.
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
THIS_DIR = Path(__file__).parent.resolve()
DATA_DIRECTORY = THIS_DIR / "data"
# choose the Zilien file to be imported (example name, data not included)
ZILIEN_FILENAME = "myzilienfile.tsv"
# choose the EC-lab file to be imported (example name, data not included)
ECLAB_FILENAME = "mybiologicfile.mpt"

#choose the directory where the figures will be saved to
FIGURES_DIR = THIS_DIR / "figures"

# In addition to the above general settings, when using a diffent dataset, 
# the time of when the different measurements occur in the dataset needs to be 
# adjusted in by changing the "tspan" in each section accordingly. To this end
# it is recommended to use full_data.plot_measurement(tspan=[t1,t2]) after importing
# and varying t1 and t2 to find the tspans of interest.

#---------------------------END OF EDIT SETTINGS--------------------------


def main():
    # # -----------------------data importing----------------------------------
    # import from the original data files. Note there is no example data included 
    # for this part of the script
    # csv files provided are as close to original datafiles as possible, chosen
    # for the script to run faster
    if DATA_SOURCE == "example": #        
        part1 = ixdat.Measurement.read(DATA_DIRECTORY / "data_part1_COstrip_11-11-21.csv", reader="ixdat")    
        part2 = ixdat.Measurement.read(DATA_DIRECTORY / "data_part2_COstrip_11-11-21.csv", reader="ixdat")
        zil_data = part1 + part2
        # remove the EC columns from the zilien file:
        zil_data.replace_series("Ewe/V", None)
        zil_data.replace_series("I/mA", None)
        # import EC data from csv files of the ECLab data
        ec_part1 = ixdat.Measurement.read(DATA_DIRECTORY / "ecdata_part1_COstrip_11-11-21.csv", reader="ixdat")    
        ec_part2 = ixdat.Measurement.read(DATA_DIRECTORY / "ecdata_part2_COstrip_11-11-21.csv", reader="ixdat")
        # combine MS and EC data
        ec_data = ec_part1 + ec_part2
        full_data = zil_data + ec_data
    
    elif DATA_SOURCE == "raw": # option 1: import both EC and MS data from Zilien data file
        full_data = ixdat.Measurement.read(DATA_DIRECTORY / ZILIEN_FILENAME, reader="zilien")
        # Integrating the EC data of this file will prompt an error!
        print("""WARNING: Data imported from Zilien file only. This will lead 
              to an error when trying to integrate the EC data using 
              cv.diff_with() use DATA_SOURCE=\"raw_biologic\" instead 
              and co-import the EC data directly from the EC-lab file.""")
        # save the important colums as ixdat-datafile
        full_data.export(DATA_DIRECTORY / "full_data_zilien_COstrip_11-11-21.csv")
        
    elif DATA_SOURCE == "raw_biologic": # option 2: import MS data from Zilien data file and combine with 
    # EC data from BioLogic file 
        zil_data = MSMeasurement.read(DATA_DIRECTORY / ZILIEN_FILENAME, reader="zilien")
        # remove the EC columns from the zilien file:
        zil_data.replace_series("Ewe/V", None)
        zil_data.replace_series("I/mA", None)
        ec_data = ixdat.Measurement.read(DATA_DIRECTORY / ECLAB_FILENAME, reader="biologic")
        full_data = zil_data + ec_data
        full_data.export(DATA_DIRECTORY / "full_data_COstrip_11-11-21.csv")
        
    elif DATA_SOURCE == "ixdat": # option 3: import from ixdat-datafiles, both for ec and ms
        # if full data saved as csv using one of the above options
        try:
            full_data = ixdat.Measurement.read(DATA_DIRECTORY / "full_data_COstrip_11-11-21.csv", reader="ixdat")
        except FileNotFoundError:
            full_data = ixdat.Measurement.read(DATA_DIRECTORY / "full_data_zilien_COstrip_11-11-21.csv", reader="ixdat")
            # Integrating the EC data of this file will prompt an error!
            print("""WARNING: Data imported from Zilien file only. This will lead 
                  to an error when trying to integrate the EC data using 
                  cv.diff_with() use DATA_SOURCE=\"raw_biologic\" instead 
                  and co-import the EC data directly from the EC-lab file.""")
    else:
        raise NameError("DATA_SOURCE not recognized.")
    # add an EC calibration
    # and because we'd like a different unit on the EC data but this is not implemented
    # in ixdat yet, we use a trick where we define a different surface area and then
    # manually change the labels in the plots.
    ec_data.calibrate(RE_vs_RHE=0, A_el=0.000196)
    full_data.calibrate(RE_vs_RHE=0, A_el=0.000196) #this is equal to converting from mA/cm2 to muA/cm2
    
    # TODO make sure that everything below can actually run from full_data only.
    # This means in particular all the data treatment of the EC files has to be 
    # adjusted for tspans etc. this is a mess to do but will hopefully make the 
    # code easier to understand. 
    
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
        axes_a = co_blank_cv.plot_measurement(mass_lists=[["M4", "M28", "M44"], ["M2", "M32"]], logplot=True, legend=False)
        axes_a[0].set_ylim(1e-14, 1e-7)
        axes_a[2].set_ylim(1.6e-10, 1e-9)        
        axes_a[0].set_ylabel("M4, M28, M44 signal / [A]")
        axes_a[2].set_ylabel("M2, M32 signal / [A]")
        axes_a[3].set_ylabel("J / [$\mu$A cm$^{-2}$]") #manually change the label for the right current unit
        coblank_vs_t_fig = axes_a[0].get_figure()
        coblank_vs_t_fig.tight_layout()
        if SAVE_FIGURES is True:
            coblank_vs_t_fig.savefig(FIGURES_DIR / ("co_blank_vs_time_final" + FIGURE_TYPE))
    
        # CO strip
        axes_b = co_strip_cv.plot_measurement(mass_lists=[["M4", "M28", "M44"], ["M2", "M32"]], logplot=True, legend=False)
        axes_b[2].set_ylim(1e-10, 5e-10)
        axes_b[2].set_ylabel("M2, M32 signal / [A]")
        axes_b[0].set_ylabel("M4, M28, M44 signal / [A]")
        axes_b[1].set_yticks([0,0.5, 1, 1.5])
        axes_b[3].set_yticks([-25, 0, 25, 50])
        axes_b[3].set_ylabel("J / [$\mu$A cm$^{-2}$]") #manually change the label for the right current unit
        co_strip_vs_t_fig = axes_b[0].get_figure()
        co_strip_vs_t_fig.tight_layout()
        if SAVE_FIGURES is True:
            co_strip_vs_t_fig.savefig(FIGURES_DIR / ("co_strip_vs_time_final" + FIGURE_TYPE))
    
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
        axes_c[1].set_yticks([-25, 0, 25])
        axes_c[1].set_ylabel("J / [$\mu$A cm$^{-2}$]") #manually change the label for the right current unit
        costrip_vs_u_fig = axes_c[0].get_figure()
        costrip_vs_u_fig.tight_layout()
        if SAVE_FIGURES is True:
            plotname = "CO_strip_vs_potential_"
            costrip_vs_u_fig.savefig(FIGURES_DIR / (plotname + WHICH_REFERENCE + FIGURE_TYPE))
        
    elif WHICH_PART == "integrate+calibrate":
        # FIRST integrate the ELECTROCHEMICAL CO stripping peak
        ec_co_strip = ec_data.cut(tspan=[12250, 14200])  
        ec_co_strip_cv = ec_co_strip.as_cv()
        ec_co_strip_cv.redefine_cycle(start_potential=0.05, redox=False)
               
        if WHICH_REFERENCE == "reference_cycle":
            ec_co_blank = ec_data.cut(tspan=[7150, 9070])  
            ec_co_blank_cv = ec_co_blank.as_cv()
            ec_co_blank_cv.redefine_cycle(start_potential=0.05, redox=False)
            cv_diff = ec_co_strip_cv[1].diff_with(ec_co_blank_cv[1])
            #integrate "raw_current" so that it's not affected by SA normalization
            Q_CO_ox_blank = cv_diff.integrate("raw_current", vspan=[0.4, 0.9]) * 1e-3  # 1e-3 converts mC --> C
            n_CO_ox_blank = Q_CO_ox_blank / (ixdat.constants.FARADAY_CONSTANT * 2)
            SA_blank = Q_CO_ox_blank*1e6/340 #muC/cm2  
            print(f"vs blank charge passed = {Q_CO_ox_blank*1e6} uC, corresponding to {n_CO_ox_blank*1e9} nmol of CO oxidized and an ECSA of {SA_blank} cm2")
        
        elif WHICH_REFERENCE == "2nd_cycle":
            cv_diff = ec_co_strip_cv[1].diff_with(ec_co_strip_cv[2])
            #integrate "raw_current" so that it's not affected by SA normalization
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
        ax_cvdiff.set_ylabel("J / [$\mu$A cm$^{-2}$]")
        cvdiff_fig = ax_cvdiff.get_figure()
        if SAVE_FIGURES is True:
            cvdiff_fig.savefig(FIGURES_DIR / ("CV_diff_EC_integration" + FIGURE_TYPE))
                
        # SECOND integrate the MS PART of the CO strip (M44 only)        
        axes_co2_strip = co_strip_cv.plot_measurement(mass_list=["M44"], tspan=[500, 2000], logplot=False, legend=False)
        #integrate and define a background
        co2_int_strip_c1 = co_strip_cv[1].integrate_signal('M44', tspan=[700, 1050], tspan_bg=[650, 700], ax=axes_co2_strip[0])
        
        if WHICH_REFERENCE == "reference_cycle":
            axes_co2_blank = co_blank_cv.plot_measurement(mass_list=["M44"], logplot=False, legend=False)
            axes_co2_blank[3].set_ylabel("J / [$\mu$A cm$^{-2}$]")
            # integrate and define a background + add to plot
            co2_int_blank = co_blank_cv[1].integrate_signal('M44', tspan=[750, 1100], tspan_bg=[700, 750], ax=axes_co2_blank[0])
            if SAVE_FIGURES is True:
                axes_co2_blank[0].get_figure().savefig(FIGURES_DIR / ("CO_blank_CV_integrated_vs_time" + FIGURE_TYPE))
            # calculate sensitivity factor F
            f_co2 = (co2_int_strip_c1-co2_int_blank)/n_CO_ox_blank
        elif WHICH_REFERENCE == "2nd_cycle":
            co2_int_strip_c2 = co_strip_cv[2].integrate_signal('M44', tspan=[1330, 1680], tspan_bg=[1300, 1330], ax=axes_co2_strip[0])
            f_co2 = (co2_int_strip_c1-co2_int_strip_c2)/n_CO_ox
        else: #this is technically not required here, but left in anyway
            raise NameError("WHICH_REFERENCE not recognized.")
            
        axes_co2_strip[0].set_ylabel("MS signal / [A]")
        axes_co2_strip[3].set_ylabel("J / [$\mu$A cm$^{-2}$]")
        if SAVE_FIGURES is True:
            axes_co2_strip[0].get_figure().savefig(FIGURES_DIR / ("CO_strip_CV_integrated_vs_time_" + WHICH_REFERENCE + FIGURE_TYPE))
        
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
        co_strip.calibrate(ms_cal_results=[cal_co2])
        strip_cycle = co_strip.cut(tspan=[500, 1290])
        if WHICH_REFERENCE == "reference_cycle": 
            co_blank_cv.calibrate(ms_cal_results=[cal_co2])
            bas_cycle = co_blank_cv[1]
            tspan_background = [700,750]
        elif WHICH_REFERENCE == "2nd_cycle":
            bas_cycle = co_strip.cut(tspan=[1290, 1900])
            tspan_background = [1700, 1800]
        else: #this is technically not required here, but left in anyway
            raise NameError("WHICH_REFERENCE not recognized.")
        
        axes_d = strip_cycle.plot_vs_potential(mol_list=["CO2"], logplot=False, legend=False, tspan_bg=[600, 700], unit="pmol/s")
        bas_cycle.plot_vs_potential(axes=axes_d, mol_list=["CO2"], logplot=False, linestyle=":", legend=False, tspan_bg=tspan_background, unit="pmol/s")
        axes_d[0].set_xlabel("$U_{RHE}$ / [V]")
        axes_d[0].set_yticks([0, 1, 2, 3, 4, 5])
        axes_d[1].set_ylabel("J / [$\mu$A cm$^{-2}$]")
        if SAVE_FIGURES is True:
            plotname = "CO_strip+2ndcycle_vs_potential_calibrated_using_"
            axes_d[0].get_figure().savefig(FIGURES_DIR / (plotname + WHICH_REFERENCE + FIGURE_TYPE))
        
    else:
        raise NameError("WHICH_PART not recognized.")
    return full_data
            
        
if __name__ == "__main__":
    main()