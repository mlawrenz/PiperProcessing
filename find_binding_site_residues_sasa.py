__doc__="""
Finding residues near ligand binding site that have high SASA, to use in PIPER
docking
"""

import os
import sys
import argparse
from schrodinger import structure
from schrodinger.utils import fileutils
from schrodinger.infra import mm
from schrodinger.structutils import analyze

def write_prime_inputs(asl, mini=False):
    if mini==True:
        ohandle=open('prime_mini_custom.inp', 'w')
        ohandle.write('''STRUCT_FILE     XXX.mae
JOB_TYPE        REFINE
PRIME_TYPE      REAL_MIN
SELECT  asl = %s
MINIM_NITER     2
MINIM_NSTEP     65
MINIM_RMSG      0.010000
MINIM_METHOD    au
USE_CRYSTAL_SYMMETRY    no
USE_RANDOM_SEED no
SEED    0
OPLS_VERSION OPLS3
EXT_DIEL        80.00
USE_MEMBRANE    no''' % asl)
        ohandle.close()
    else:
        ohandle=open('prime_hybridmc_custom.inp', 'w')
        ohandle.write('''STRUCT_FILE     XXX.mae
JOB_TYPE        REFINE
PRIME_TYPE      MC
SELECT  asl = %s
NSTEPS  25
PROB_SIDEMC     0.000000
PROB_RIGIDMC    0.000000
PROB_HMC        1.000000
TEMP_HMC        900.000000
FIND_BEST_STRUCTURE     yes
NUM_OUTPUT_STRUCT       1
USE_CRYSTAL_SYMMETRY    no
USE_RANDOM_SEED no
SEED    0
OPLS_VERSION OPLS3
EXT_DIEL        80.00
USE_MEMBRANE    no''' % asl)    
        ohandle.close()
    return




def filter_residues_by_sasa(struct, st, all_chains_dict, sasa_cutoff, distance, binding_site_residues, writepdb=False):

    from schrodinger.application.bioluminate import protein
    calculator = protein.PropertyCalculator(st, "my_calculator_jobname")
    all_sasa_dict=dict()
    all_names_dict=dict()
    print "Selecting residues for %s with SASA > %s, Recommend that you double check the selection" % (struct, sasa_cutoff)

    ohandle=open('%s_sasa%s_%s_binding_site_residues.txt' % (struct, sasa_cutoff, distance), 'w')

    if writepdb==True:
        writer=structure.PDBWriter('%s_attract.pdb' % struct)
    binding_sasa_dict=dict()
    delete_indices=[]
    n=0
    for res in structure.get_residues_unsorted(st):
    #for res in structure._ResidueIterator(st, atoms='all'):
    #for res in st._getStructureResidueIterator():
        sasa=calculator.getResidueSASA(res, sidechain=False)
        resnum=res._resnum
        resname=res.pdbres
        chain=res.chain
        if n==0:
            all_chains_dict[chain]=[]
            prev_chain=chain
            n+=1
        if prev_chain!=chain:
            print "******************************************"
            print "******************************************"
            print "chain number changes for %s, make the same!" % struct
            print "DID NOT PRINT DISTANCES"
            print "******************************************"
            print "******************************************"
            sys.exit()
            
        all_sasa_dict[resnum]=sasa
        all_names_dict[resnum]=resname
        if resnum in binding_site_residues:
            print "%s%s SASA %s" % ( resname, resnum, round(sasa, 1))
            if sasa > sasa_cutoff:
                #writer.append(res.extractStructure())
                all_chains_dict[chain].append(resnum)
                ohandle.write("%s%s SASA %s\n" % ( resname, resnum, round(sasa, 1)))
            else:
                indices=res.getAtomList()
                for i in indices: delete_indices.append(i)
        else:
            indices=res.getAtomList()
            for i in indices: delete_indices.append(i)
    st.deleteAtoms(delete_indices)
    if writepdb==True:
        writer.append(st)
        writer.close()
    ohandle.close()
    #at._setAtomLabelUserText('%0.2f' % sasa)
    #print "SASA for %s %s charge %s is %0.2f" % (at.index, at.element, at.formal_charge, sasa)

    maxval=max(all_sasa_dict.values())
    maxval_res=all_sasa_dict.keys()[ all_sasa_dict.values().index(maxval)]
    minval=min(all_sasa_dict.values())
    minval_res=all_sasa_dict.keys()[ all_sasa_dict.values().index(minval)]
    #if verbose==True:
    print "MAX SASA IS %s Res %s%s" % (round(maxval,1), all_names_dict[maxval_res], maxval_res)
    print "MIN SASA IS %s Res %s%s" % (round(minval,1), all_names_dict[minval_res], minval_res)
    print "-------------------------"
    return chain, all_chains_dict


def find_binding_site_residues(input_structure, distance_cutoff):
    st = structure.StructureReader(input_structure).next()
    asl_searcher = analyze.AslLigandSearcher()
    ligands = asl_searcher.search(st)
    if len(ligands)==0:
        print "NO LIGAND IN %s, RESTART" % input_structure
        sys.exit()
        #import pdb
        #pdb.set_trace()
    for lig in ligands:
        print "ligand is called %s" % lig.pdbres
        binding_site = analyze.evaluate_asl(st,"(fillres within %s %s) and (not %s)" % (distance_cutoff, lig.ligand_asl, lig.ligand_asl))
    return st, binding_site


def main(args):
    if args.debug==True:
        import pdb
        pdb.set_trace()
    structures=dict()
    structures['rec'] = args.rec
    if args.lig:
        structures['lig'] = args.lig
    chains=dict()
    # Print parameters
    chain_reference=dict()
    print "Distance cutoff = %s" % args.distance_cutoff
    n=0
    for struct in structures.keys():
        print "on %s %s" % (struct, structures[struct]) 
        st, binding_site_indices=find_binding_site_residues(structures[struct], args.distance_cutoff)
        binding_site_residues=[st.atom[atom].resnum for atom in binding_site_indices]
        if n==0:
            all_chains_dict=dict()
            n+=1
        chain, all_chains_dict=filter_residues_by_sasa(struct, st,
all_chains_dict, args.sasa_cutoff, args.distance_cutoff, binding_site_residues)
        chain_reference[struct]=chain
    rec_chain=chain_reference['rec']
    if args.lig:
        lig_chain=chain_reference['lig']
        # NEED THIS SPECIFIC OUTPUT REC CHAIN REC ID LIG CHAIN LIG ID
        difference=len(all_chains_dict[rec_chain])-len(all_chains_dict[lig_chain])
        if difference > 0:
            all_chains_dict[lig_chain]=all_chains_dict[lig_chain] + ['NA']*abs(difference)
        else:
            all_chains_dict[rec_chain]=all_chains_dict[rec_chain] + ['NA']*abs(difference)
        ohandle=open('distances_for_restraints.txt', 'w')
        for (x,y) in zip(all_chains_dict[rec_chain], all_chains_dict[lig_chain]):
            ohandle.write('%s-%s\t%s-%s\n' % (rec_chain,x, lig_chain, y))
        ohandle.close()
        if args.prime==True:
            rec_num=[]
            for (n, x) in enumerate(all_chains_dict[rec_chain]):
                if x!='NA':
                    rec_num.append(str(x))
            lig_num=[]
            for (n, x) in enumerate(all_chains_dict[lig_chain]):
                if x!='NA':
                    lig_num.append(str(x))
            rec_num=','.join(rec_num)
            lig_num=','.join(lig_num)
            asl='((chain %s) AND (res.num %s)) OR ((chain %s) AND (res.num %s))' % (rec_chain, rec_num, lig_chain, lig_num)
            write_prime_inputs(asl, mini=True)
            write_prime_inputs(asl, mini=False)
    return

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='''Find binding site residues accounting for SASA of the atom.''')
    parser.add_argument('-r', dest="rec", help="PDB receptor file for docking, with ligand")
    parser.add_argument('-l', dest="lig", help="PDB ligand file for docking, with ligand. Can leave this off to just compute SASA list for receptor protein.")

    parser.add_argument('-s', dest='sasa_cutoff', type=float, default=50, help="SASA cutoff for binding site residue selection. Default is 50 Angstrom^2")
    parser.add_argument('-d', dest='distance_cutoff', type=float, default=4, help="Receptor-ligand distance cutoff for residue selection. Default is 4.0 Anstrom")
    parser.add_argument('--debug',action='store_true',help="Debug.")
    parser.add_argument('--prime',action='store_true',help="Write prime outputs with ASL from find SASA")

    args = parser.parse_args()
    main(args)

