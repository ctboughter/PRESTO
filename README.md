# PRESTO - PaRsEr of Solved Tcr cOmplexes
An automated pipeline for the analysis of TCR-pMHC structures deposited to the PDB

Newest version (04/25/23) has been updated to be more robust against poorly deposited PDBs. Now able to measure atomic distances
for all class II complex structures and 93% of class I complex structures. These distances are now included to a 6 Angstrom cutoff 
for the PDBs deposited into the TCR3D database thus far. So you don't need to run the parser if you don't want to.

If you'd like to analyze TCR-MHC complexes with your own specifications, you can edit and run the parser_main.py code, which will then 
call some of the functions in tcr_structParse.py . From there you can use the jupyter notebook to analyze the data.

Packages:

````
python2.7
scipy 1.2.1
numpy 1.16.5
pandas 0.24.2
matplotlib 2.2.3  
mdtraj 1.9.3
biopython 1.76
````

mdtraj is the only one of these where I *know* the version matters. Stick with that one

p.s. the name *was* the result of a concerted effort to make the most ridiculous backronym possible

# Reproducing Analysis of Boughter & Meier-Schellersheim 2022
The update on July 26th, 2022 is meant to allow reproduction of structural analysis code in a preprint that is soon to be posted to Biorxiv and submitted for peer review.

The parser_main.py code must be run first, and then the Jupyter Notebook will recreate the figures from that manuscript. Check out that manuscript (link soon) for more details.