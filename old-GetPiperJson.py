import os
import glob
import json
import numpy
import sys
import argparse
import pickle





def write_json(resdata,  outdir, ngroups=1):
    groups=[resdata,] # later can add more to have more groups
    dataorg=dict()
    dataorg["groups"]=[]
    dataorg["required"]=ngroups
    for n in range(0, ngroups):
        groupdata=dict()
        nrestraints=len(groups[n].keys())
        # default require N-1
        groupdata["required"]=nrestraints-1
        groupdata["restraints"]=[]
        for m in sorted(groups[n].keys()):
            groupdata["restraints"].append(groups[n][m])
        dataorg["groups"].append(groupdata)
    outfile=open("%s/distance_restraints.json" % outdir, 'w')
    json.dump(dataorg, outfile)
    return


def parse_residues(infile):
    resdata=dict()
    fhandle=open(infile)
    for (n, line) in enumerate(fhandle.readlines()):
        receptor=line.split()[0] 
        ligand=line.split()[1]
        residue=dict()
        residue['rec_chain']=receptor.split('-')[0]
        residue['rec_resi']=receptor.split('-')[1]
        residue['lig_chain']=ligand.split('-')[0]
        residue['lig_resi']=ligand.split('-')[1]
        try: 
            residue['dmin']=line.split()[2]
        except IndexError:
             residue['dmin']=2
        try:
            residue['dmax']=line.split()[3]
        except IndexError:
             residue['dmax']=4 
        residue['type']='residue'
        resdata[n]=residue
    return resdata    


def main(args):
    resdata=parse_residues(args.infile)
    outdir=os.path.dirname(os.path.abspath(args.infile))
    print resdata
    write_json(resdata, outdir)
    return

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='make distance restraint json file for piper')
    parser.add_argument('-i', dest="infile", help="file with residue pairs for distance restraints")
    args = parser.parse_args()
    main(args)

