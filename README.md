General workflow for a CTM and 2 protein binding partners:
Get two proteins with ligands bound (derived from CTM)
Use information on linker length to determine ligand atom restraint cutoffs.
Also include protein restraints (see below).

*binding site residues for protein-protein restraints based on ligand sites and SASA
$SCHRODINGER/run ~/piper-scripts/find_binding_site_residues_sasa.py -r receptor-protein-with-small-molecule.pdb -l ligand-protein-with-small-molecule.pdb

*get JSON restraints file from text file
$SCHRODINGER/run ~/PiperProcessing/GetPiperJson.py -i distances_for_restraints.txt
(defaults are set but you can check the number of restraints you have with --combo)

*run piper
receptor=$1
ligand=$2
het1=6Z3
het2=JQ1
restraints=distance_restraints.json
 
run_piper --rec $receptor --lig $ligand --restraint-set $restraints --vdw-het $het1 $het2 -np 12 --mpi 

(((( If it complains add these flags: )))
--piper-license ~/piper-1.1.0/piper_package/piper_license  --pip er-base ~/piper-1.1.0/piper_package

*process output structure files for analysis in Maestro
Go into the subdirectory with model*pdb files and a results.txt file. You want
to process the model*min*pdb files. Add these to a list.txt file then run:
$SCHRODINGER/run ~/PiperProcessing/processPiper.py -l list.txt -r report.txt

This makes a maefiles/ directory, go in there.
Make another list file of the maefiles (sorry, should combine these scripts in
the future. Pass in reference.
$SCHRODINGER/run ~/PiperProcessing/calculate_piper_rmsd.py -l list -r Brd4_5T35_wligand.pdb 

This makes a rmsd-maefiles/ directory. You can load these files into maestro.






MORE DETAIL ON WORKFLOW BELOW:

For protein binding site restraints first run this to get binding site
residues with high(er) SASA: 
----------------------------------------
$SCHRODINGER/run ~/piper-scripts/find_binding_site_residues_sasa.py -r receptor-protein-with-small-molecule.pdb -l ligand-protein-with-small-molecule.pdb

Run -h to get tuneable options.
Outputs all binding site residues with SASA:
rec_sasa50_binding_site_residues.txt
lig_sasa50_binding_site_residues.txt

Also outputs file of above residues selected by SASA to be used for distance restraints:
distances_for_restraints.txt
The columns here are not restricted to each other within each line, but rather
al combinations of these pairs are added to possible restraints. This is for
more prospective work.


PLEASE CHECK THAT THESE RESIDUES ARE GOOD FOR YOUR PURPOSE. ITERATIONS MAY BE
NECESSARY TO GET THE BEST RESIDUES.
edit the distances_for_restraints.txt to add/remove residues


Make distance restraints file JSON:
------------------------------------------

In vim you can use this to visualize JSON file:
:%!python -m json.tool
