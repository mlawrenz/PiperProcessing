import os
import glob
import json
import numpy
import sys
import argparse
import pickle


def write_combo_json(groups,  outdir, min_contacts):
    dataorg=dict()
    ngroups=len(groups)
    dataorg["groups"]=[]
    if min_contacts=='all':
        dataorg["required"]=ngroups
    else:
        dataorg["required"]=int(min_contacts)
    print "number of groups : %s" % ngroups
    print "required groups %s" % dataorg["required"]
    for n in range(0, ngroups):
        groupdata=dict()
        # require 1 from the group satisfied
        groupdata["required"]=1
        groupdata["restraints"]=[]
        for m in groups[n].keys():
            groupdata["restraints"].append(groups[n][m])
        dataorg["groups"].append(groupdata)
    outfile=open("%s/distance_restraints.json" % outdir, 'w')
    json.dump(dataorg, outfile)
    return

def write_json(resdata,  outdir, nrestraints=None, ngroups=1):
    groups=[resdata,] # later can add more to have more groups
    dataorg=dict()
    dataorg["groups"]=[]
    dataorg["required"]=ngroups
    for n in range(0, ngroups):
        groupdata=dict()
        if nrestraints==None:
            nrestraints=len(groups[n].keys())
        groupdata["required"]=int(nrestraints)
        groupdata["restraints"]=[]
        for m in sorted(groups[n].keys()):
            groupdata["restraints"].append(groups[n][m])
        dataorg["groups"].append(groupdata)
    outfile=open("%s/distance_restraints.json" % outdir, 'w')
    json.dump(dataorg, outfile)
    return

def create_residue_dict(res, lig, dmax):
    residue=dict()
    residue['rec_chain']=res.split('-')[0]
    residue['rec_resi']=res.split('-')[1]
    residue['lig_chain']=lig.split('-')[0]
    residue['lig_resi']=lig.split('-')[1]
    residue['dmin']=float(2)
    residue['dmax']=float(dmax) 
    residue['type']='residue'
    return residue

def parse_combo_residues(infile, dmax):
    groups=[]
    resdata=dict()
    fhandle=open(infile)
    receptor=numpy.loadtxt(infile, usecols=(0,), dtype=str)
    ligand=numpy.loadtxt(infile, usecols=(1,), dtype=str)
    n=0
    for i in receptor:
        if 'NA' in i:
            receptor=numpy.delete(receptor, n)
        else:
            n+=1
    n=0
    for i in ligand:
        if 'NA' in i:
            ligand=numpy.delete(ligand, n)
        else:
            n+=1
    total_restraints=len(receptor)+len(ligand)
    print "total restraints %s" % total_restraints
    for res in receptor:
        # make multiple restraints to ligand, will be in one group, with 1
        # required, do the same for ligand
        for (n, lig) in enumerate(ligand):
            residue=create_residue_dict(res, lig, dmax)
            resdata[n]=residue
        groups.append(resdata)
    for lig in ligand:
        for (n, res) in enumerate(receptor): 
            residue=create_residue_dict(res, lig, dmax)
            resdata[n]=residue
        groups.append(resdata)
    return groups

def parse_ligand_atoms(ligfile, ldmin, ldmax):
    ligdata=dict()
    fhandle=open(ligfile)
    for (n, line) in enumerate(fhandle.readlines()):
        lig_restraint=dict()
        receptor=line.split()[0] 
        ligand=line.split()[1]
        receptor_lig=dict()
        receptor_lig['chain']=receptor.split('-')[0]
        receptor_lig['resi']=receptor.split('-')[1]
        receptor_lig['atom_name']=receptor.split('-')[2]
        ligand_lig=dict()
        ligand_lig['chain']=ligand.split('-')[0]
        ligand_lig['resi']=ligand.split('-')[1]
        ligand_lig['atom_name']=ligand.split('-')[2]
        lig_restraint["rec_atom"]=receptor_lig
        lig_restraint["lig_atom"]=ligand_lig
        lig_restraint["dmin"]=float(ldmin)
        lig_restraint["dmax"]=float(ldmax)
        lig_restraint["rec_type"]="atom" 
        lig_restraint["lig_type"]="atom" 
        ligdata[n]=lig_restraint
    return ligdata    


def parse_residues(infile, dmax):
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
             residue['dmax']=float(dmax) 
        residue['type']='residue'
        resdata[n]=residue
    return resdata    


def main(args):
    if args.debug:
        import pdb
        pdb.set_trace()
    if args.protein_infile:
        outdir=os.path.dirname(os.path.abspath(args.protein_infile))
    elif args.ligand_infile:
        outdir=os.path.dirname(os.path.abspath(args.ligand_infile))
    else:
        print "PLEASE PROVIDE A LIGAND OR PROTEIN FILE"
        sys.exit()
# DOING LIGAND ONLY
    if args.ligand_infile:
        print "reading one-to-one ligand atom contraints"
        if not args.ldmax:
            print "DID NOT SPECIFY MAX DIST FOR LIGAND RESTRAINT"
            print "DEFAULT IS 10"
        ligdata=parse_ligand_atoms(args.ligand_infile, args.ldmin, args.ldmax)            
        if args.protein_infile:
# DOING LIGAND + PROTEIN
            if args.specific==True:
                print "reading one-to-one constraints from protein input file"
                resdata=parse_residues(args.protein_infile, args.dmax)
                groups=[]
                groups.append(resdata)
                groups.append(ligdata)
            else:
                print "creating combindations of constraints from input file"
                groups=parse_combo_residues(args.protein_infile, args.dmax)
                groups.append(ligdata)
            write_combo_json(groups, outdir, args.min_contacts)
        else:
            groups=[]
            groups.append(ligdata)
            if args.specific:
                write_combo_json(groups, outdir, min_contacts='all')
            else:
                print "ADDING SPECIFIC HERE SINCE JUST LIGAND. IF OTHERWISE CHECK SCRIPT"
                write_combo_json(groups, outdir, min_contacts='all')
    else:
# DOING PROTEIN ONLY
        if args.specific==True:
            print "reading one-to-one constraints from protein input file"
            resdata=parse_residues(args.protein_infile, args.dmax)
            print resdata
            write_json(resdata, outdir, args.min_contacts)
        else:
            print "creating combindations of constraints from input file"
            groups=parse_combo_residues(args.protein_infile, args.dmax)
            write_combo_json(groups, outdir, args.min_contacts)
    return

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='''
make distance restraint json file for piper. Default will create unbiased
combinations for your distances.txt file. Add --specific and restraints for the
pairs of residues will be made''')
    parser.add_argument('-proteinfile', dest="protein_infile", help="file with residue pairs for distance restraints")
    parser.add_argument('-ligandfile', dest="ligand_infile", help="file with ligand atom pairs for distance restraints. First receptor protein ligand then ligand protein ligand")
    parser.add_argument('-d', dest="dmax", default=4.5, help="dmax for protein distance restraints, default is 4.5, can probably leave this unless you want less tight criteria for restraints")
    parser.add_argument('--ldmin', dest="ldmin", default=4.0, help="dmin for ligand atom distance restraints, default is 4.0")
    parser.add_argument('--ldmax', dest="ldmax", default=10.0, help="dmax for ligand atom distance restraints, default is 10.0")
    parser.add_argument('--min-contacts', dest="min_contacts", default=3,
help="minimum number of contacts for ternary to have within selection of binding site residues. The default is 3. If you pass in (all) then this will be set to all specified.")
    parser.add_argument('--specific', action="store_true", dest="specific", help="use specific pairs of restraints from file. Default makes unbiased combinations.")
    parser.add_argument('--debug',  action="store_true", dest='debug', help='DEBUG FLAG')
    args = parser.parse_args()
    main(args)

