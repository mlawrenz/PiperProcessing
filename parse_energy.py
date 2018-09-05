import os
import csv
import sys
import argparse
import glob
from schrodinger import structure
from schrodinger.utils import cmdline, fileutils
from schrodinger.structutils import rmsd

def get_energy(inputfiles, oname):
    monomer=[]
    complex=[]
    energy=dict()
    for file in inputfiles:
        if 'chain' in file:
            monomer.append(file)
        else:
            complex.append(file)
    if len(complex) > 1:
        print("more than one complex")
        print(complex)
        sys.exit()
    if len(monomer) > 2:
        print("more than two monomer files")
        print(monomer)
        sys.exit()
    fhandle=open(complex[0])
    for line in fhandle.readlines():
        #continue until last line with TOTALE, optimized structure
        if 'TOTALE' in line: 
            energy['complex']=float(line.split()[1])
    fhandle.close()
    for (n,file) in enumerate(monomer):
        fhandle=open(file)
        for line in fhandle.readlines():
            if 'TOTALE' in line: 
                energy[n]=float(line.split()[1])
                break
        fhandle.close()
    return energy 

def main(listfile=None, oname='prime_mini', inputfile_basename=None):
    
    if listfile:
        listfile = open(listfile, 'r')
        list_basename_files = listfile.read().splitlines()
        ohandle=open('%s_i_energy_model.000.summary.txt' % oname, 'w')
        for inputfile_basename in list_basename_files:
            if 'log' not in inputfile_basename:
                print("NEED LOG FILES")
                sys.exit()
            inputfiles=glob.glob('*%s*log' % inputfile_basename.rstrip('.log'))
            energy=get_energy(inputfiles, oname)
            interaction_energy=energy['complex']-(energy[0]+energy[1])
            print(inputfile_basename, interaction_energy)
            ohandle.write('%s\t%0.2f\n' % (inputfile_basename, interaction_energy))
        ohandle.close()
    else:
        inputfiles=glob.glob('%s*log' % inputfile_basename.rstrip('.log'))
        energy=get_energy(inputfiles)
        interaction_energy=energy['complex']-(energy[0]+energy[1])
        print(interaction_energy)
        ohandle=open('%s.energy.txt' % inputfile_basename, 'w')
        ohandle.write('%0.2f' % interaction_energy)
        ohandle.close()

    return
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='compute rmsd for ligand protein to reference for set of maefiles')
    parser.add_argument('-i', dest='inputfile_basename', help='basename for the log files with energies to parse.')

    parser.add_argument('-l','--listfile', dest='listfile', help='file list basenames to work on (do not include chain splits')
    parser.add_argument('-o','--oname', dest='oname', help='output file names', default='prime_mini')

    parser.add_argument('--debug',  action="store_true", dest='debug' )
    args = parser.parse_args()
    if args.debug:
        import pdb
        pdb.set_trace()
    main(args.listfile, args.oname, args.inputfile_basename)



