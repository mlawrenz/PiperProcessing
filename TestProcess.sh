rec_atom="C1"
lig_atom="N2"
lig_residue='JQ1'
rec_residue='V32'
chain='A'

# tes tthat all these work
$SCHRODINGER/run processPiperNew.py
$SCHRODINGER/run processPiperNew.py --rmsd-ref $ref --chain A
$SCHRODINGER/run processPiperNew.py --rec_residue $rec_residue --lig_residue $lig_residue --lig_atom $lig_atom --rec_atom $rec_atom 
$SCHRODINGER/run processPiperNew.py --rec_residue $rec_residue --lig_residue $lig_residue --lig_atom $lig_atom --rec_atom $rec_atom  --rmsd-ref $ref --chain ${chain}


