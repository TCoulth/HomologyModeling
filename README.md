# HomologyModeling

================

PythoHammering.py

Purpose:        Template and catalytic constraint generation for homology modeling. 

How/What:       Performs Hmmer searches against the PDB and the MCSA
                MCSA hits provide reisudes to be used to generate catalytic constraints for targets with RosettaEnzCM.py
                Filters out PDB and MCSA hits by UniProt ID to provide unique templates only
                
To use:         Need one folder entitled 'Fastas' containing all target fastas
                If using outside of genomecenter cluster, will need new path's for PDB and MCSA databases
                
Python Modules: pandas, os, pypdb


=================
