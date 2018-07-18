
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


def get_maefiles(rfile, filenames):
    global PIPER_BIN
    split_piper='%s/split_piper' % PIPER_BIN
    data=dict()
    for line in rfile.readlines():
        if 'ENERGY' in line:
            cluster=True
            if 'SIZE' in line:
                cluster=False
            pass
        else:
            molname=line.split()[0].split('.pdb')[0]
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

    ohandle=open('maefile_list.txt', 'w')
    for fname in filenames :
        os.system('%s %s' % (split_piper, fname))
        molname=os.path.basename(fname).split('.pdb')[0]
        outdir='%s/maefiles/' % os.getcwd()
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(fname):
            print "missing file %s" % fname
            continue
        if cluster==True:
            writer = structure.StructureWriter('%s/%s.mae' % (outdir, molname))
            print "writing %s.mae" % molname
            ohandle.write('%s/%s.mae\n' % (outdir, molname))
            for st in structure.StructureReader(fname):
                st.title=molname
                propedit.addProperty(st, 'r_piper_charmm_energy', data[molname]['charmm_energy'], add=True)
                writer.append(st)
        else:
            writer = structure.StructureWriter('%s/%s.mae' % (outdir, molname))
            print "writing %s.mae" % molname
            ohandle.write('%s/%s.mae\n' % (outdir, molname))
            rec='%s.rec.pdb' % fname.split('.pdb')[0]
            lig='%s.lig.pdb' % fname.split('.pdb')[0]
            for st in structure.StructureReader(rec):
                st.title='%s_rec' % molname
                propedit.addProperty(st, 's_piper_cluster_size', data[molname]['cluster_size'], add=True)
                propedit.addProperty(st, 'r_piper_charmm_energy', data[molname]['charmm_energy'], add=True)
                propedit.addProperty(st, 's_piper_coefficient', data[molname]['potential'], add=True)
                writer.append(st)
            for st in structure.StructureReader(lig):
                st.title='%s_lig' % molname
                propedit.addProperty(st, 's_piper_cluster_size', data[molname]['cluster_size'], add=True)
                propedit.addProperty(st, 'r_piper_charmm_energy', data[molname]['charmm_energy'], add=True)
                propedit.addProperty(st, 's_piper_coefficient', data[molname]['potential'], add=True)
                writer.append(st)
            writer.close()
    ohandle.close()
    return

def get_rmsd(reference_file, listfile):
    calculate_piper_rmsd.main(reference_file, listfile)
    return


def main(args): 

    if args.debug:
        import pdb
        pdb.set_trace()


    listfile = open(args.listfile, 'r')
    filenames = listfile.read().splitlines()
    
    rfile=open(args.report, 'r')
    # do piper splitting to get maefiles
    get_maefiles(rfile, filenames)
    excess=glob.glob('*rec.pdb')
    for f in excess:
        os.remove(f)
    excess=glob.glob('*lig.pdb')
    for f in excess:
        os.remove(f)
    # run RMSD to reference
    if args.rmsd_reference:
        get_rmsd(args.rmsd_reference, 'maefile_list.txt')
        #import shutil
        #shutil.rmtree('./maefiles/')
        os.remove('./maefile_list.txt')
    return



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='make maefiles out of piper results, adding data')
    parser.add_argument('-l','--listfile', dest='listfile', help='file list with names of model*min*pdb files that will be processed for maestro with cluster and enerfy properties')
    parser.add_argument('-r','--report', dest='report', help='report file with piper data for these files')
    parser.add_argument('--rmsd_ref', dest='rmsd_reference', help='reference to report RMSD for all structures;')
    parser.add_argument('--debug',  action="store_true", dest='debug' )

    
    args = parser.parse_args()
    main(args)





 
