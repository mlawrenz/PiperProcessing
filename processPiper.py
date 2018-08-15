#
import glob
import argparse, os, sys
from schrodinger import structure
from schrodinger.structure import *
from schrodinger.infra import propedit
import calculate_piper_rmsd

PIPER_BIN=os.environ['PIPER_BIN']

def coeff(val):
    coeffs=dict()
    coeffs[0]='Balanced'
    coeffs[1]='Electrostatics'
    coeffs[2]='Hydrophobics'
    coeffs[3]='No Potential'
    return coeffs[val]


def measure_distance(st, rec_residue, rec_atom, lig_residue, lig_atom):
    st=st.next()
    from schrodinger import structutils
    asl_expr='(res.ptype "%s ") AND (atom. "%s")' % (lig_residue, lig_atom)
    asl_searcher = structutils.analyze.AslLigandSearcher()
    ligands = asl_searcher.search(st)

    atom1=structutils.analyze.evaluate_asl(st, asl_expr)
    asl_expr='(res.ptype "%s ") AND (atom. "%s")' % (rec_residue, rec_atom)
    atom2=structutils.analyze.evaluate_asl(st, asl_expr)
    if len(atom2) > 1:
        atom2=atom2[1]
    else:
        atom2=atom2[0]
    atom1=atom1[0]
    distance=st.measure(atom1, atom2)
    return distance


def read_report(rfile):
    data=dict()
    for line in rfile.readlines():
        if 'ENERGY' in line:
            cluster=True
            if 'SIZE' in line:
                cluster=False
            pass
        else:
            molname='complex-%s' % line.split()[0].split('.pdb')[0]
            if cluster==True:
                charmm_energy=float(line.split()[1])
                cluster_size='N/A'
                coefficient='N/A'
            else:
                cluster_size=line.split()[1]
                charmm_energy=float(line.split()[2])
                val=molname.split('.')[1][2]
                coefficient=coeff(int(val))
            data[molname]=dict()
            data[molname]['cluster_size']=cluster_size
            data[molname]['charmm_energy']=charmm_energy
            data[molname]['potential']=coefficient
    return data

def get_rmsd_data(data):
    if not os.path.exists('rmsd.txt'):
        print "MISSING RMSD FILE. CHECK SCRIPT"
        sys.exit()
    rmsdfile=open('rmsd.txt')
    for line in rmsdfile.readlines():
        molname=line.split()[0].split('.pdb')[0].split('.mae')[0]
        rmsd=float(line.split()[1])
        if molname not in data.keys():
            print "missing %s from data" % molname  
            import pdb
            pdb.set_trace()
        else:
            data[molname]['rmsd']=rmsd
    return data
    

def get_maefiles(rfile, filenames, reference_file=None, chain=None, lig_residue=None, rec_residue=None, lig_atom=None, rec_atom=None):


    data=read_report(rfile)
    for fname in filenames :
        # make sure start is one complex
        # this is important and screwed me for a few days
        os.system('sed "/HEADER/d" < %s | sed "s/END/TER/g" > complex-%s' % (fname, fname))

        #### 
    if reference_file:
        if not chain:
            print "REQUIRES PROTEIN CHAIN NAME"
            sys.exit()
        ohandle=open('triagelist.txt', 'w')
        for fname in filenames :
            new_structure_file='complex-%s' % fname
            st=structure.StructureReader(new_structure_file)
            if rec_atom:
                distance=measure_distance(st, rec_residue, rec_atom, lig_residue, lig_atom)
                if distance > 20:
                    print "DISCARDING POSE %s" % new_structure_file
                    print distance
                    continue
                else:
                    ohandle.write('%s\n' % new_structure_file)
            else:
                    ohandle.write('%s\n' % new_structure_file)
        ohandle.close()
        calculate_piper_rmsd.main(reference_file, listfile='triagelist.txt', chain=chain)
        # add rmsd to maefile data
        data=get_rmsd_data(data)

    else:
        ohandle=open('triagelist.txt', 'w')
        for fname in filenames :
            new_structure_file='complex-%s' % fname
            st=structure.StructureReader(new_structure_file)
            if rec_atom:
                distance=measure_distance(st, rec_residue, rec_atom, lig_residue, lig_atom)
                if distance > 20:
                    print "DISCARDING POSE %s" % new_structure_file
                    print distance
                    continue
                else:
                    ohandle.write('%s\n' % new_structure_file)
            else:
                    ohandle.write('%s\n' % new_structure_file)
        ohandle.close()
    
    newlist=open('triagelist.txt')
    for new_structure_file in newlist.readlines() :
        if reference_file:
            outdir='%s/rmsd-maefiles/' % os.getcwd()
        else:
            outdir='%s/maefiles/' % os.getcwd()
        molname=new_structure_file.split('.pdb')[0]
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        writer = structure.StructureWriter('%s/%s.mae' % (outdir, molname))
        print "writing %s.mae" % molname
        for st in structure.StructureReader(new_structure_file.rstrip()):
            st.title=molname
            propedit.addProperty(st, 'r_piper_charmm_energy', data[molname]['charmm_energy'], add=True)
            propedit.addProperty(st, 's_piper_cluster_size', data[molname]['cluster_size'], add=True)
            propedit.addProperty(st, 's_piper_coefficient', data[molname]['potential'], add=True)
            if reference_file:
                propedit.addProperty(st, 'r_user_CA_RMSD', data[molname]['rmsd'],  add=True)
            writer.append(st)
        writer.close()
    return

def main(args): 

    if args.debug:
        import pdb
        pdb.set_trace()


    if args.listfile:
        listfile = open(args.listfile, 'r')
        filenames = listfile.read().splitlines()
    else:
        filenames=glob.glob('model*min*.pdb')
    
    if not os.path.exists('report.txt'):
        print "NEED TO BE IN DIRECTORY WITH report.txt"
        sys.exit()
    rfile=open('report.txt', 'r')
    get_maefiles(rfile, filenames, args.reference_file, args.chain, args.lig_residue, args.rec_residue, args.lig_atom, args.rec_atom)

    excess=glob.glob('complex*.pdb')
    for f in excess:
        os.remove(f)
    return



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='make maefiles out of piper results, adding data. Pass RMSD ref to compute RMSD and add that data. Add filter based on ligand atom distances by passing in rec residue and atom, ligand residue and atom for distance cutoff. Current distance cutoff is 20 A, could add customization.' )

    parser.add_argument('-l','--listfile', dest='listfile', help='file list with names of files that will be processed for maestro with cluster and energy properties. If you do not pass this, will run on *min*.pdb')
    parser.add_argument('--chain', dest='chain', help='chain in complex for ligand to compute RMSD')
    parser.add_argument('--rmsd-ref', dest='reference_file', help='reference to report RMSD for all structures;')
    parser.add_argument('--rec_residue', dest='rec_residue', help='receptor small molecule residue name')
    parser.add_argument('--lig_residue', dest='lig_residue', help='ligand small molecule residue name')
    parser.add_argument('--lig_atom', dest='lig_atom', help='ligand small molecule atom name')
    parser.add_argument('--rec_atom', dest='rec_atom', help='receptor small molecule atom name')
    parser.add_argument('--debug',  action="store_true", dest='debug' )

    
    args = parser.parse_args()
    main(args)





 

