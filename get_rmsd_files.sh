 

for x in `seq 0 2`; do grep "model\.00${x}" rmsd.txt | sort -nk1 >> ../5T35-Brd4-BRD2_JQpr__5LLI-VHL-V02methyl_model${x}_rmsd.txt; done
