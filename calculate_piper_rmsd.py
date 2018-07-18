_doc__ = """
Script to calculate the RMSD of the CA atoms between the reference structure and
a list of PDB files.

Copyright Schrodinger, LLC. All rights reserved.
"""

###############################################################################
# Globals and constants

import os
import csv
import argparse

from schrodinger import structure
from schrodinger.utils import cmdline, fileutils
from schrodinger.structutils import rmsd


###############################################################################
def main(reference_file, listfile=None, pdb_file=None):
    """
    Main body of the script.
    """

    if listfile:
        listfile = open(listfile, 'r')
        filenames = listfile.read().splitlines()
    else:
        filenames=[pdb_file]


    ca_asl = '(atom.ptype " CA ")'
    #basename = fileutils.get_basename(cmd_args.mobile_pdb_file[0])
    #outfile = basename + '-rmsd.mae'

    ref_st = structure.StructureReader(reference_file).next()
    #writer.append(ref_st)

    ohandle=open('rmsd.txt', 'w')
    for pdb_file in filenames:
        basename=os.path.basename(pdb_file)
        writer = structure.StructureWriter('%s-rmsd.mae' % basename.split('.mae')[0])
        for pdb_st in structure.StructureReader(pdb_file):
            try:
                conf_rmsd = rmsd.ConformerRmsd(ref_st, pdb_st, asl_expr=ca_asl)
                ca_rmsd = conf_rmsd.calculate()
                pdb_st.property['r_user_CA_RMSD'] = ca_rmsd
                print pdb_st.title, ca_rmsd
                ohandle.write('%s\t%0.2f\n' % (pdb_st.title, ca_rmsd))
                writer.append(pdb_st)
                continue
            except RuntimeError:
                print '%s and %s have different number of CA atoms. Skipping.' % (reference_file, pdb_file)
                pass

        writer.close()
    ohandle.close()
    os.mkdir('rmsd-maefiles')
    os.system('mv *rmsd*.mae rmsd-maefiles/')


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='compute rmsd for ligand protein to reference for set of maefiles')
    parser.add_argument('-r', dest='reference_file', help='Reference structure file.')

    parser.add_argument('-s', dest='pdb_file', help='PDB file if only single calc')


    parser.add_argument('-l','--listfile', dest='listfile', help='file list with names of mae files that will be get cluster property added to')

    parser.add_argument('--debug',  action="store_true", dest='debug' )
    args = parser.parse_args()
    if args.debug:
        import pdb
        pdb.set_trace()
    main(args.reference_file, args.listfile, args.pdb_file)

