# AGN-Specfit

**AGN-Specfit** is a modified version of [QSFit](https://github.com/gcalderone/qsfit). 
It is created as a pipeline to perform the analysis of SDSS Type 1 Active Galactic Nuclei optical spectra.

To make use of **AGN-Specfit** pipeline, you have to:
1. Add SDSS spectra files to "data/" directory.
2. Create the list of file names, redshifts, and E(B-V) in "QSO_name.csv".
3. Create directories which are named:
⋅⋅*"output" and "table", to store the calculation result in nested and flattened structure files respectively.
⋅⋅*"plot", to save all the plot data in a ﬁle suitable to be later loaded in GNUPLOT
⋅⋅*"result", to store the merged tables from "table" folder.
4. Start an `IDL` session in your working directory.
5. Compile all the QSFIT procedures and dependencies by calling the `compile.pro` procedure.
5. By using `Python`, run `make_table.py` to combine all tables in "table/" folder into one concatenated table. You can specify which columns to keep by modifying "result/columns_to_keep.txt" files.
6. Edit the scripts if necessary to suit your needs.

Further details about how to use QSFIT can be found in http://qsfit.inaf.it.
The reference paper is https://arxiv.org/abs/1612.01580. 