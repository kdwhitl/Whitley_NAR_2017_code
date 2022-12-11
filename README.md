# Whitley_NAR_2017_code
 
OVERVIEW

This repository contains several scripts and functions associated with Figure 1 of Whitley et al. Nucleic Acids Research (2017), as well as some raw data for testing them.

The purpose of the code overall is to take some raw data containing both fluorescence and optical trap data for fluorescently-labeled oligonucleotides binding under applied force, and to extract various information from this (size of extension change, binding rate, unbinding rate, etc.).

There are several subdirectories here:
- Original_code contains the code as it existed around the time the NAR paper was published (~2017). It is here for the sake of transparency and reproducibility, although it is not especially user-friendly.
- Power_spectra contains the code associated with analyzing the calibration files associated with this data. This is legacy code that dates back to ~2006 or so (written either by Jeff Moffitt or Yann Chemla).
- test_data contains some example data that you can use to test the code and make sure it runs on your machine.
- Updated contains code that has been somewhat modified and cleaned from the original version, to make it somewhat more user-friendly. I recommend using the code in here for analysis.

There are two main scripts/functions that are used for analysis:
- hybrid_step_size_yank (version 3 in the original folder, version 4 in the updated one). This script (v3) or function (v4) loads three data files: a calibration file, an offset file (blank force-extension curve), and a force + trap data file. It identifies events where a fluorescently-labeled oligonucleotide has bound or unbound from a complementary strand under force, and records various information about this (step size, lifetime, etc.). These are outputted in a structure variable AllResults with three fields: binding, unbinding, and lifetimes.
- plot_cut_results (version 3 in original folder, version 4 in the updated one). This script takes the structure variable AllResults, cuts the data based on various things the user chooses (force range, etc.), then does some further analysis and plotting to get various results (extension change, unbinding rate constant, etc.).

INPUTS:

hybrid_step_size_yank4:
- probe: Name of oligonucleotide whose dataset you want to analyze
- DSetsToAnalyze: Which slice of the data trace you want to use. This is associated with the lines in hybridDataFile3 (each slice of data may have its own analysis parameters associated with it).
- beaddiameters: Diameters of the two beads you trapped with this experiment (in nm)
- overview: Whether you want to plot the entire trace (fluorescence, extension, force) with no analysis. Choices are 0 or 1.
- fluor_plots: Whether you wan to plot anything at all (including the overview). Choices are 0 or 1.
- removeyankoffsets: Whether you want to use the offset file to get trap offset values, or use the value from the calibration file (0: use value from calibration file, 1: use offset file).
- step_find: Whether you want to find steps using a sliding t test on the extension trace, or to use the first derivative of the fluorescence (0: use first derivative of fluorescence trace to find steps, 1: use sliding t test on extension trace to find steps).
- step_plots: Whether you want to plot the results of the sliding t test method. Choices are 0 or 1.
- desired_trap_bw: Bandwidth of trap data you want to use for analysis (in Hz).
- avplotfact: Averaging factor you want ot use for plotting fluorescence data (does not affect analysis, only plotting).

plot_cut_results4: This is a script, and so has no direct user inputs. It is run from the Matlab Editor itself, not the command line. Near the top of the file are a series of variables beginning with 'plot'. Depending on what the user wants to plot or analyze, these should be set to 1 (plot this) or 0 (don't plot this). Below this is a series of filters that can optionally be used on the data. Typically these are left off. Below these filters is a variable Tset, where a user chooses which force region of data to analyze/plot.

OUTPUTS:

hybrid_step_size_yank4:
- AllResults: A structure variable with three fields (binding, unbinding, and lifetimes). Each element within each field contains ata on a single event (binding or unbinding).
- YankResults: A cell array of the fluorescence signal across each force range. This is only used to estimate equilibrium constants, so is mostly not needed.

plot_cut_results4:
- plots of various things (lifetime distribution, extension change, etc.)

EXAMPLE USAGE:

hybrid_step_size_yank4:
[AllResults, YankResults] = hybrid_step_size_yank4('10mer0mg', 5:9, [790 880], 1, 1, 1, 0, 0, 5, 20);

plot_cut_results4:
plot_life = 1;
plot_onrate = 1;
(all else set to 0)

Tset = 4

DEMO:

In the subdirectory test_data, there are two things: a folder (120814) and a file (AllResults.mat). The folder contains raw data that can be analyzed by hybrid_step_size_yank4: a calibration file (120814_045), offset files (120814_046), force-extension curve data associated with this DNA tether (120814_49), and a force + fluorescence time trace file (120814_050). This corresponds to the probe '10mer0mg' with lines 5 to 9 (as shown in the 'example usage' section above). The expected output of this data using exactly the line above is a structure variable AllResults, which is also provided in the file AllResults.mat. The file AllResults.mat can be run using plot_cut_results4. The output (if using the parameters described above in 'example usage' should be the plots 221211_10mer0mg_10pN_lifetimes and 22111_10mer0mg_10pN_t_ub, which are provided in the test_data directory.
