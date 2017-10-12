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


##############################################################################
def parse_args():
    """
    Parse the command line options.

    @return:  All script arguments and options
    @rtype:  class:`argparse.Namespace`
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-r', dest='reference_file', help='Reference structure file.')

    parser.add_argument('-s', dest='pdb_file', help='PDB file')

    parser.add_argument('-l','--listfile', dest='listfile', help='file list with names of mae files that will be get cluster property added to')

    args = parser.parse_args()

    if not os.path.isfile(args.reference_file):
        parser.error('Input reference file %s not found.' % args.reference_file)

    if args.listfile: 
        if not os.path.isfile(args.listfile):
            parser.error('listfile %s not found.' % args.listfile)
    else:
        if not os.path.isfile(args.pdb_file):
            parser.error('PDB file %s not found.' % pdb_file)

    return args


###############################################################################
def main():
    """
    Main body of the script.
    """

    #import pdb
    #pdb.set_trace()
    cmd_args = parse_args()
    if cmd_args.listfile:
        listfile = open(cmd_args.listfile, 'r')
        filenames = listfile.read().splitlines()
    else:
        filenames=[cmd_args.pdb_file]


    ca_asl = '(atom.ptype " CA ")'
    #basename = fileutils.get_basename(cmd_args.mobile_pdb_file[0])
    #outfile = basename + '-rmsd.mae'

    ref_st = structure.StructureReader(cmd_args.reference_file).next()
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
                print '%s and %s have different number of CA atoms. Skipping.' % (cmd_args.reference_file, pdb_file)
                pass

        writer.close()
    ohandle.close()
    os.mkdir('rmsd-maefiles')
    os.system('mv *rmsd*.mae rmsd-maefiles/')


if __name__ == '__main__':
    cmdline.main_wrapper(main)
