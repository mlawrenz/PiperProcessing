_doc__ = """
Script to calculate the RMSD of the CA atoms between the reference structure and
a list of PDB files.

Copyright Schrodinger, LLC. All rights reserved.
"""

###############################################################################
# Globals and constants

import os
import csv
import sys
import argparse

from schrodinger import structure
from schrodinger.utils import cmdline, fileutils
from schrodinger.structutils import rmsd


###############################################################################
def main(reference_file, asl=None, listfile=None, pdb_file=None, chain=None, postprocess=None, writermsd=False):



    """
    Main body of the script.
    """

    if listfile:
        listfile = open(listfile, 'r')
        filenames = listfile.read().splitlines()
    else:
        filenames=[pdb_file]


    if not asl:
        asl = '(atom.ptype " CA ")'
        outfile='rmsd.txt'
        if chain:
            print("specified chain")
            asl = '((chain.name %s)) AND ((atom.ptype " CA "))' % chain
    if postprocess:
        #asl = '((chain.name %s)) AND (backbone)' % args.chain
        #asl = '((((( backbone ) ) AND NOT ((res.ptype "ACE "))) AND NOT((res.ptype "NMA "))) AND ((chain.name %s))) AND NOT ((atom.ele H))' % chain
        asl='(((((( backbone ) ) AND NOT ((res.ptype "ACE "))) AND NOT ((res.ptype "NMA "))) AND NOT ((atom.ele H))) AND NOT ((atom.ptype "OXT")))'
        outfile='new-rmsd.txt'




        
    #basename = fileutils.get_basename(cmd_args.mobile_pdb_file[0])
    #outfile = basename + '-rmsd.mae'

    ref_st = next(structure.StructureReader(reference_file))
    #writer.append(ref_st)
    outfile='rmsd.txt'
    ohandle=open(outfile, 'w')
    for pdb_file in filenames:
        basename=os.path.basename(pdb_file)
        if writermsd:
            writer = structure.StructureWriter('%s-rmsd.mae' % basename.split('.mae')[0].split('.pdb')[0])

        for pdb_st in structure.StructureReader(pdb_file):
            try:
                from schrodinger import structutils
                ref_atoms=  structutils.analyze.evaluate_asl(ref_st, asl)
                model_atoms= structutils.analyze.evaluate_asl(pdb_st, asl)

                conf_rmsd = rmsd.ConformerRmsd(ref_st, pdb_st, asl_expr=asl)
                ca_rmsd = conf_rmsd.calculate()
                pdb_st.property['r_user_RMSD'] = ca_rmsd
                print(pdb_st.title, ca_rmsd)
                if pdb_st.title:
                    ohandle.write('%s\t%0.2f\n' % (pdb_st.title, ca_rmsd))
                else:
                    ohandle.write('%s\t%0.2f\n' % (pdb_file, ca_rmsd))
                if writermsd:
                    writer.append(pdb_st)
                continue
            except RuntimeError:
                print('%s and %s have different number of CA atoms. Skipping.' % (reference_file, pdb_file))
                pass

        if writermsd:
            writer.close()
    ohandle.close()
    if not os.path.exists('rmsd-maefiles'):
        os.mkdir('rmsd-maefiles')
    if writermsd:
        os.system('mv *rmsd*.mae rmsd-maefiles/')
    return


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='compute rmsd for ligand protein to reference for set of maefiles')
    parser.add_argument('-r', dest='reference_file', help='Reference structure file.')

    parser.add_argument('-s', dest='pdb_file', help='PDB file if only single calc')
    parser.add_argument('-c', dest='chain', help='chain name to use for calculation in PDB file if only single calc')

    parser.add_argument('-l','--listfile', dest='listfile', help='file list with names of mae files that will be get cluster property added to')
    parser.add_argument('--writermsd',  action="store_true", dest='writermsd', help='write out maefiles with rmsd.')

    parser.add_argument('--debug',  action="store_true", dest='debug' )
    args = parser.parse_args()
    if args.debug:
        import pdb
        pdb.set_trace()
    main(reference_file=args.reference_file, listfile=args.listfile, pdb_file=args.pdb_file, chain=args.chain, writermsd=args.writermsd)

