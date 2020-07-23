# TectoMagMOR
Toy & analytical model for dike-fault interactions at mid-ocean ridges: predicts the fraction of magmatic plate separation (M) as a function of spreading rate, axial lithospheric thickness, and magmatic input. See Olive & Dublanchet (EPSL, submitted) for details.
MATLAB script M_across_rates.m plots the analytical model, and uses function tectono_magmatic_fcn.m to run toy model simulations. To generate Figs. 1b; 2; 3; 4 and 5, one can run M_across_rates.m by setting variable regime (line 9) to ’thin’ or ’thick’ to select an end-member scenario for lithospheric thickness vs. spreading rate (red and blue curves in Fig. 2a). All scripts must be run from the same folder, which must include the three .txt files containing the M (M data.txt), micro-earthquake depth (depth of microEQs.txt), and AML depth (depth of AMLs.txt) datasets, respectively.
