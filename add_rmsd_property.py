from schrodinger.structure import StructureReader
import os
import argparse
import glob
import sys
from schrodinger.infra import propedit


def main(args):
    rmsdfile=args.rmsd
    fhandle=open(rmsdfile)
    for line in fhandle.readlines():
        name=line.split()[0].split('.maegz')[0]
        rmsd=float(line.split()[1])
        file='%s.maegz' % name
        st_reader = StructureReader(file)
        for (num, st) in enumerate(st_reader):
            #st.property['r_slow_dE_Local']='%0.2f' % float(docked_num)
            propedit.addProperty(st, 'r_user_CA_RMSD', float(rmsd), add=True)
            st.write('rmsd-%s' % file)
    return


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='add rmsd to complex structures. result of pain in the ass aobut splitting the protein complexes.')
    parser.add_argument('--rmsd', dest='rmsd', help='rmsd.txt with names and rmsd vals')
    args = parser.parse_args()
    main(args)


