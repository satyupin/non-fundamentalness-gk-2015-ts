# GK2015 Fundamentalness
Matlab code for a series of exercises on informational sufficiency for SVAR. It applies monetary policy surprise series as in Gertler and Karadi (2015) as an external instrument to identify a monetary policy shock, then proceeds to assess whether the underlying SVAR is fundamental using orthogonality tests as in Forni and Gambetti (2014). Then it compares two versions of an SVAR: an original and an informationally sufficient one by computing the impulse-responses and forecast error variance decompositions.

# Contents

[main.m](main.m): Main file to produce most of the results.

[\functions](/functions): Functions to estimate VAR and IV regressions, compute impulse-responses and other subroutines.

[\data](/data): Data for analysis:

* [fred_md_012022.csv](/data/fred_md_012022.csv): series from FRED-MD database used in factor construction.
* [GK2015_Data.xlsx](/data/GK2015_Data.xlsx): original Gertler and Karadi (2015) data.
* [FRED-MD_updated_appendix.csv](/data/FRED-MD_updated_appendix.csv): an auxiliary file with filenames and transformation codes.



