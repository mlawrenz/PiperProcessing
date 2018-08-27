from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy as sp
import argparse
from openeye.oechem import *
import sys
import multiprocessing
import os
from sklearn.metrics.pairwise import pairwise_distances


# In[2]:


def rigid_transform_3D(A, B):
    assert len(A) == len(B)

    N = A.shape[0]; # total points

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    
    # centre the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = np.transpose(AA) * BB

    U, S, Vt = np.linalg.svd(H)

    R = Vt.T * U.T

    # special reflection case
    if np.linalg.det(R) < 0:
       #print "Reflection detected"
       Vt[2,:] *= -1
       R = Vt.T * U.T

    t = -R*centroid_A.T + centroid_B.T

    #print t

    return R, t


# In[3]:


def LibGen(libgen, ofs, unique, isomeric):
    smiflag = OESMILESFlag_DEFAULT  # Canonical|AtomMaps|Rgroup
    if isomeric:
        smiflag |= OESMILESFlag_ISOMERIC
    # access products
    uniqproducts = []
    for mol in libgen.GetProducts():
        smiles = OECreateSmiString(mol, smiflag)
        if not unique or smiles not in uniqproducts:
            uniqproducts.append(smiles)
            OEWriteMolecule(ofs, mol)
    return


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--whA_linkers_file')
    parser.add_argument('--whA_bindingConf_file')
    parser.add_argument('--whB_linkerAtoms_file')
    parser.add_argument('--whB_bindingConf_file')
    parser.add_argument('--whA_prot_file')
    parser.add_argument('--whB_prot_file')
    parser.add_argument('--smirksFile')
    parser.add_argument('--outMolsFile')
    parser.add_argument('--outComplexFile')
    parser.add_argument('--outProteinFile')
    parser.add_argument('--outPDBdir')
    parser.add_argument('--nproc')
    parser.add_argument('--maxProtLigClashes', type=int, default=30)
    parser.add_argument('--maxProtProtClashes', type=int, default=30)
    parser.add_argument('--whAconfs_maxConfs', default=1000, type=int)
    parser.add_argument('--whBconfs_rms', default=0.01, type=float)
    args = parser.parse_args()

    
    
    whA_linkers_file = args.whA_linkers_file
    whA_bindingConf_file = args.whA_bindingConf_file
    whB_bindingConf_file = args.whB_bindingConf_file
    whB_3linkerAtoms_file = args.whB_linkerAtoms_file
    whB_prot_file = args.whB_prot_file
    whA_prot_file = args.whA_prot_file
    smirksFile = args.smirksFile
    outMolsFile = args.outMolsFile
    outComplexesFile = args.outComplexFile
    outproteinsFile = args.outProteinFile
    pdbOutDir = args.outPDBdir



    ##other inputs
    nproc = args.nproc
    nprotligclashes = args.maxProtLigClashes
    nprotprotClashes = args.maxProtProtClashes
    whB_omega_rms = args.whBconfs_rms
    maxConfs_whA_linker = args.whAconfs_maxConfs

    ##notes on inputs
    ## 1) in both linker files, used C14 for terminal carbon
    ## 2) In whA_linkers_file, sd property called linker gives name to linker
    ## 3) In whB_3linkerAtoms_file, sd property called structureForLinkers gives names of whA+linkers which this will be used for
    ## 4) in smirks file, linker smarts first, then warhead smarts
    ## 4) in smirks file, make the left to right atom order of the whA and whB smarts overlap/have them be in proper order


    # In[22]:


    whBprotFile = oemolistream(whB_prot_file)

    for mol in whBprotFile.GetOEGraphMols():
        allProtCoords = []
        coordDict = mol.GetCoords()
        for i, atom in enumerate(mol.GetAtoms()):
            allProtCoords.append(coordDict[i])
        #print 'numatoms = ', i+1
        npAllwhBProtCoords = np.mat(allProtCoords)
        protMol = OEGraphMol(mol)
    #print len(npAllwhBProtCoords)
    #print npAllwhBProtCoords


    whAprotFile = oemolistream(whA_prot_file)

    for mol in whAprotFile.GetOEGraphMols():
        whAprotMol = OEGraphMol(mol)


    whA_linkers_confs_file = whA_linkers_file.split('.sdf')[0]+'_confs.sdf'
    print 'whA+linkers confs file =', whA_linkers_confs_file

    if nproc == 'auto':
        nproc = multiprocessing.cpu_count()
    #print nproc

    print "/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -umatch false -mpi_np {} -flipper true -prefix {} -maxconfs {} -searchff mmff94s -dielectric_constant 78> {}".format(whA_linkers_file, whA_linkers_confs_file, whA_bindingConf_file, str(nproc), 'whA_linkers_omega', maxConfs_whA_linker, 'whA_omega_stdout.txt')
    os.system("/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -flipper true -prefix {} -maxconfs {} > {}".format(whA_linkers_file, whA_linkers_confs_file, whA_bindingConf_file, 'whA_linkers_omega', maxConfs_whA_linker, 'whA_linkers_omega_stdout.txt'))

    whB_3linkerAtoms_confs_file = whB_3linkerAtoms_file.split('.sdf')[0]+'_confs.sdf'
    print 'whB+linkerAtoms confs file = ', whB_3linkerAtoms_confs_file

    print "/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -mpi_np {} -flipper true -rms {} -prefix {} > {}".format(whB_3linkerAtoms_file, whB_3linkerAtoms_confs_file, whB_bindingConf_file, str(nproc), str(whB_omega_rms), 'whB_omega', 'whB_omega_stdout.txt')
    os.system("/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -mpi_np {} -flipper true -rms {} -prefix {} > {}".format(whB_3linkerAtoms_file, whB_3linkerAtoms_confs_file, whB_bindingConf_file, str(nproc), str(whB_omega_rms), 'whB_omega', 'whB_omega_stdout.txt'))


    smartsWHADict = {}
    smartsWHBDict = {}
    smirksDict = {}
    for aline in open(smirksFile, 'r').read().strip().split("\n"):
        smirks, name = aline.split("\t")
        smirksDict[name] = smirks
        reactants = smirks.split('>>')[0]
        linkerSmarts, warheadSmarts = reactants.split('.')
        print name, linkerSmarts, warheadSmarts
        smartsWHADict[name] = linkerSmarts
        smartsWHBDict[name] = warheadSmarts

    whBConfFile = oemolistream(whB_3linkerAtoms_confs_file)

    whBConfCoordDict = {}
       
    allWarheadBCoords = []
    for i,mol in enumerate(oemolistream(whB_3linkerAtoms_confs_file).GetOEGraphMols()):
        
        print 'numAtoms in whB+3atoms = ', mol.GetMaxAtomIdx()
        
        allWarheadBCoords = []
        coordDict = mol.GetCoords()
        for j, atom in enumerate(mol.GetAtoms()):
            allWarheadBCoords.append(coordDict[j])
        print 'numatoms = ', j+1
        npAllWarheadBCoords = np.mat(allWarheadBCoords)
        
        #print 'startLoop'    
        
        #print OECreateIsoSmiString(mol)
        relevantLinkers = OEGetSDData(mol, 'structureForLinkers').split("\n")
        #print relevantLinkers
        
        for alinker in relevantLinkers:
            if alinker not in whBConfCoordDict.keys():
                whBConfCoordDict[alinker] = {}
            #print smartsWHBDict[alinker]
            ss = OESubSearch(smartsWHBDict[alinker])
            OEPrepareSearch(mol, ss)
            
            matchingAtomIndices = []
            for match in ss.Match(mol):
                print 'matchingAtoms'
                for ma in match.GetAtoms():
                    #print ma.target.GetIdx()
                    matchingAtomIndices.append(ma.target.GetIdx())
            #print matchingAtomIndices
            threeAtomCoords_whB = np.mat([coordDict[matchingAtomIndices[0]], coordDict[matchingAtomIndices[1]], coordDict[matchingAtomIndices[2]]])
            #print threeAtomCoords_whB
            
            #print i, alinker
            whBConfCoordDict[alinker][i] = [threeAtomCoords_whB, npAllWarheadBCoords, OEGraphMol(mol)]
            
            
        #print "\n"

    #for akey in whBConfCoordDict.keys():
    #    print akey, whBConfCoordDict[akey]



    if not os.path.isdir(pdbOutDir):
        os.mkdir(pdbOutDir)

    outmols = oemolostream(outMolsFile)
    outComplexes = oemolostream(outComplexesFile)
    outproteins = oemolostream(outproteinsFile)

    complexNum = 0

    for i, mol in enumerate(oemolistream(whA_linkers_confs_file).GetOEGraphMols()):
        #print i
        thisLinker = OEGetSDData(mol, 'linker')

        ss = OESubSearch(smartsWHADict[thisLinker])
        OEPrepareSearch(mol, ss)

        matchingAtomIndices = []
        for match in ss.Match(mol):
            for ma in match.GetAtoms():
                matchingAtomIndices.append(ma.target.GetIdx())
        coordDict = mol.GetCoords()
        threeAtomCoords_whA = np.mat([coordDict[matchingAtomIndices[0]], coordDict[matchingAtomIndices[1]], coordDict[matchingAtomIndices[2]]])

        for aWHBconf in whBConfCoordDict[thisLinker].keys():
            
            R, t = rigid_transform_3D(whBConfCoordDict[thisLinker][aWHBconf][0], threeAtomCoords_whA)

            n = len(whBConfCoordDict[thisLinker][aWHBconf][1])
            whBcoordMat = whBConfCoordDict[thisLinker][aWHBconf][1]
            C = R*whBcoordMat.T + np.tile(t, (1, n))
            C = C.T
            D = np.asarray(C).reshape(-1)


            whbMol = whBConfCoordDict[thisLinker][aWHBconf][2]
            whbMol.SetCoords(D)


            reaction = OEQMol()
            smirks = smirksDict[thisLinker]
            OEParseSmirks(reaction, smirks)
            libgen = OELibraryGen()
            libgen.Init(reaction)
            libgen.AddStartingMaterial(mol, 0, False)
            libgen.AddStartingMaterial(whbMol, 1, False)

            outComplex = OEGraphMol()

            for prodmol in libgen.GetProducts():
                prodMol = OEGraphMol(prodmol)
                OEAddMols(outComplex, prodMol)
                OEWriteMolecule(outmols, prodmol)



            nprot = len(npAllwhBProtCoords)
            E = R*npAllwhBProtCoords.T + np.tile(t, (1,nprot))
            E = E.T
            F = np.asarray(E).reshape(-1)
            protMol.SetCoords(F)


            combinedProtMol = OEGraphMol()
            OEAddMols(combinedProtMol, whAprotMol)
            OEAddMols(combinedProtMol, protMol)

            interactionHintComplex = OEInteractionHintContainer()

            interactionHintComplex.AddMolecule(prodMol, OELigandInteractionHintComponent())
            interactionHintComplex.AddMolecule(combinedProtMol, OEProteinInteractionHintComponent())


            OEPerceiveInteractionHints(interactionHintComplex)
            nrClashintras = len([c for c in interactionHintComplex.GetInteractions(OEIsClashInteractionHint())])
            nrContactIntras = len([c for c in interactionHintComplex.GetInteractions(OEIsContactInteractionHint())])
            nrUnpairedProteinIntras = len([c for c in interactionHintComplex.GetInteractions(OEIsUnpairedProteinInteractionHint())])

            #print nrClashintras

            if nrClashintras > nprotligclashes:
                continue

            protACoords = []
            protBCoords = []
            protACoordDict = whAprotMol.GetCoords()
            protBCoordDict = protMol.GetCoords()
            for a, atom in enumerate(whAprotMol.GetAtoms()):
                protACoords.append(protACoordDict[a])
            npProtACoords = np.mat(protACoords)

            for a, atom in enumerate(protMol.GetAtoms()):
                protBCoords.append(protBCoordDict[a])
            npProtBCoords = np.mat(protBCoords)

            distMat = pairwise_distances(npProtACoords, npProtBCoords)#, n_jobs=str(nproc))
            flatDistMat = distMat.flatten()
            numProtProtClashes = len(np.where(flatDistMat < 1)[0])

            #print numProtProtClashes
            
            if numProtProtClashes > nprotprotClashes:
                continue


            
            print 'whAconf = ', i
            print 'WHBconf = ', aWHBconf
            print 'number ligand-protein clashes = ', nrClashintras
            print 'numProtProtClashes = ', numProtProtClashes
            print 'complexNum = ', complexNum
            
            OEAddMols(outComplex, whAprotMol)
            OEAddMols(outComplex, protMol)
            OEWriteMolecule(outComplexes, outComplex)
            OEWriteMolecule(outproteins, combinedProtMol)
            
            pdbOutFileName = pdbOutDir+'complex{}.pdb'.format(str(complexNum))
            pdbOutFile = oemolostream(pdbOutFileName)
            OEWriteMolecule(pdbOutFile, outComplex)
            pdbOutFile.flush()
            pdbOutFile.close()
            complexNum += 1
            
            reformatFile = open(pdbOutFileName, 'r')
            text = reformatFile.read().strip("\n")
            reformatFile.close()
            firstChunk, secondChunk = text.split("\nENDMDL")
            replacedSecondChunk = secondChunk.replace('ATOM  ', 'HETATM')
            outText = "".join([firstChunk, replacedSecondChunk])
            reformatFile = open(pdbOutFileName, 'w')
            reformatFile.write(outText)
            reformatFile.flush()
            reformatFile.close()
            

            print "\n"
                

    outmols.flush()
    outmols.close()
    outComplexes.flush()
    outComplexes.close()
    outproteins.flush()
    outproteins.close()


if __name__ == '__main__':
    main()


