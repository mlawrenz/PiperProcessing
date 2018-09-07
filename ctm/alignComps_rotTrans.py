from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy as sp
import argparse
from openeye.oechem import *
import sys
import multiprocessing as mp
import os
from sklearn.metrics.pairwise import pairwise_distances
from functools import partial
import time


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
       Vt[2,:] *= -1
       R = Vt.T * U.T

    t = -R*centroid_A.T + centroid_B.T

    return R, t

def LibGen(libgen, ofs, unique, isomeric):
    smiflag = OESMILESFlag_DEFAULT  
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


def detectInternalCTMclashes(ctmMol):
    ctmCoords = []
    ctmCoordDict = ctmMol.GetCoords()
    for a, atom in enumerate(ctmMol.GetAtoms()):
        ctmCoords.append(ctmCoordDict[a])
    npCTMCoords = np.mat(ctmCoords)
    distMat = pairwise_distances(npCTMCoords, npCTMCoords)
    flatDistMat = distMat.flatten()
    numInternalCTMClashes = len(np.where(flatDistMat < 1)[0])
    numAtoms = ctmMol.GetMaxAtomIdx()
    numInternalCTMClashes = numInternalCTMClashes - numAtoms

    return numInternalCTMClashes


def buildComplexes(whAconf, whBConfCoordDict, smartsWHADict, smirksDict, npAllwhBProtCoords, protMol, whAprotMol, nprotligclashes, nprotprotClashes):#whAprot, whBprot, smartsWHADict, whBConfCoordDict, whBConfCoordDict, smirksDict, npAllwhBProtCoords, protMol):
    
    results = []
    inputIndex = whAconf[0]
    mol = whAconf[1]
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
        ##              0       1               2               3               4                   5               6           7       8
        ##results = [ctmMol, complexMol, complexWCTMMol, protProtPassBool, protLigPassBool, numInternalCTMClashes, linkerName, whAconf, whBconf]
        result = [None, None, None, None, None, None, None, None, None]
        result[6] = thisLinker
        result[7] = inputIndex
        result[8] = aWHBconf
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
            prodMol.SetTitle(thisLinker+'_'+str(inputIndex)+'_'+str(aWHBconf))
            OESetSDData(prodMol, 'LinkerName', thisLinker)
            OESetSDData(prodMol, 'whAconf', str(inputIndex))
            OESetSDData(prodMol, 'whBconf', str(aWHBconf))
            OEAddMols(outComplex, prodMol)

        numInternalCTMclashes = detectInternalCTMclashes(prodMol)
        

        result[5] = numInternalCTMclashes

        nprot = len(npAllwhBProtCoords)
        E = R*npAllwhBProtCoords.T + np.tile(t, (1,nprot))
        E = E.T
        F = np.asarray(E).reshape(-1)
        protMol.SetCoords(F)

        combinedProtMol = OEGraphMol()
        OEAddMols(combinedProtMol, whAprotMol)
        OEAddMols(combinedProtMol, protMol)
        OEAddMols(outComplex, whAprotMol)
        OEAddMols(outComplex, protMol)

        combinedProtMol.SetTitle(thisLinker+'_'+str(inputIndex)+'_'+str(aWHBconf))
        OESetSDData(combinedProtMol, 'LinkerName', thisLinker)
        OESetSDData(combinedProtMol, 'whAconf', str(inputIndex))
        OESetSDData(combinedProtMol, 'whBconf', str(aWHBconf))

        outComplex.SetTitle(thisLinker+'_'+str(inputIndex)+'_'+str(aWHBconf))
        OESetSDData(outComplex, 'LinkerName', thisLinker)
        OESetSDData(outComplex, 'whAconf', str(inputIndex))
        OESetSDData(outComplex, 'whBconf', str(aWHBconf))

        interactionHintComplex = OEInteractionHintContainer()
        interactionHintComplex.AddMolecule(prodMol, OELigandInteractionHintComponent())
        interactionHintComplex.AddMolecule(combinedProtMol, OEProteinInteractionHintComponent())

        OEPerceiveInteractionHints(interactionHintComplex)
        nrClashintras = len([c for c in interactionHintComplex.GetInteractions(OEIsClashInteractionHint())])
        nrContactIntras = len([c for c in interactionHintComplex.GetInteractions(OEIsContactInteractionHint())])
        nrUnpairedProteinIntras = len([c for c in interactionHintComplex.GetInteractions(OEIsUnpairedProteinInteractionHint())])
        
        if nrClashintras > nprotligclashes:
            protProtClashBool = 0
            protLigClashBool = 0
            result[3] = protProtClashBool
            result[4] = protLigClashBool
            
        else:
            protProtClashBool = 1
            result[3] = protProtClashBool
        
        if protProtClashBool:
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

            distMat = pairwise_distances(npProtACoords, npProtBCoords)
            flatDistMat = distMat.flatten()
            numProtProtClashes = len(np.where(flatDistMat < 1)[0])

            if numProtProtClashes > nprotprotClashes:
                protLigClashBool = 0
            else:
                protLigClashBool = 1

            result[4] = protLigClashBool

        OESetSDData(prodMol, 'PassProtProtClashFilter', str(bool(protProtClashBool)))
        OESetSDData(prodMol, 'PassProtLigClashFilter', str(bool(protLigClashBool)))
        OESetSDData(prodMol, 'numInternalCTMclashes', str(numInternalCTMclashes))
        result[0] = OEGraphMol(prodMol)
        OESetSDData(combinedProtMol, 'PassProtProtClashFilter', str(bool(protProtClashBool)))
        OESetSDData(combinedProtMol, 'PassProtLigClashFilter', str(bool(protLigClashBool)))
        OESetSDData(combinedProtMol, 'numInternalCTMclashes', str(numInternalCTMclashes))
        result[1] = OEGraphMol(combinedProtMol)
        OESetSDData(outComplex, 'PassProtProtClashFilter', str(bool(protProtClashBool)))
        OESetSDData(outComplex, 'PassProtLigClashFilter', str(bool(protLigClashBool)))
        OESetSDData(outComplex, 'numInternalCTMclashes', str(numInternalCTMclashes))
        result[2] = OEGraphMol(outComplex)

        results.append(result)
    return results


def writeOutPDBsForComplexes(acomplex, pdbOutDir):
        complexWCTMMol, complexNum = acomplex
        pdbOutFileName = pdbOutDir+'complex{}.pdb'.format(str(complexNum))
        pdbOutFile = oemolostream(pdbOutFileName)
        OEWriteMolecule(pdbOutFile, complexWCTMMol)
        pdbOutFile.flush()
        pdbOutFile.close()
        
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
    parser.add_argument('--outMolsFile', default='outmols.sdf')
    parser.add_argument('--outComplexFile', default='outcomplexes.oeb')
    parser.add_argument('--outProteinFile', default='outproteins.oeb')
    parser.add_argument('--outPDBdir', default='outPDBs/')
    parser.add_argument('--outDataFile', default='outdata.txt')
    parser.add_argument('--rocsQuery', default=False)
    parser.add_argument('--nproc', default='auto')
    parser.add_argument('--maxProtLigClashes', type=int, default=30)
    parser.add_argument('--maxProtProtClashes', type=int, default=30)
    parser.add_argument('--maxInternalCTMclashes', type=int, default=5)
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
    outDataFile = args.outDataFile

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

    whBprotFile = oemolistream(whB_prot_file)

    for mol in whBprotFile.GetOEGraphMols():
        allProtCoords = []
        coordDict = mol.GetCoords()
        for i, atom in enumerate(mol.GetAtoms()):
            allProtCoords.append(coordDict[i])
        npAllwhBProtCoords = np.mat(allProtCoords)
        protMol = OEGraphMol(mol)

    whAprotFile = oemolistream(whA_prot_file)

    for mol in whAprotFile.GetOEGraphMols():
        whAprotMol = OEGraphMol(mol)


    whA_linkers_confs_file = whA_linkers_file.split('.sdf')[0]+'_confs.sdf'
    print 'whA+linkers confs file =', whA_linkers_confs_file

    if nproc == 'auto':
        nproc = mp.cpu_count()
    else:
        nproc = int(nproc)

    print "/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -umatch false -maxmatch 5 -mpi_np {} -flipper true -prefix {} -maxconfs {} -searchff mmff94s -dielectric 78 &> {}".format(whA_linkers_file, whA_linkers_confs_file, whA_bindingConf_file, str(nproc), 'whA_linkers_omega', maxConfs_whA_linker, 'whA_linkers_omega_stdout.txt')
    os.system("/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -fixrms 0.3 -umatch false -maxmatch 5 -mpi_np {} -flipper true -prefix {} -maxconfs {} -searchff mmff94s -dielectric 78 &> {}".format(whA_linkers_file, whA_linkers_confs_file, whA_bindingConf_file, str(nproc), 'whA_linkers_omega', maxConfs_whA_linker, 'whA_linkers_omega_stdout.txt'))

    whB_3linkerAtoms_confs_file = whB_3linkerAtoms_file.split('.sdf')[0]+'_confs.sdf'
    print 'whB+linkerAtoms confs file = ', whB_3linkerAtoms_confs_file

    print "/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -mpi_np {} -flipper true -rms {} -prefix {} &> {}".format(whB_3linkerAtoms_file, whB_3linkerAtoms_confs_file, whB_bindingConf_file, str(nproc), str(whB_omega_rms), 'whB_omega', 'whB_omega_stdout.txt')
    os.system("/home/pnovick/installers/openeye/bin/omega2 -in {} -out {} -fixfile {} -fixrms 0.3 -mpi_np {} -flipper true -rms {} -prefix {} &> {}".format(whB_3linkerAtoms_file, whB_3linkerAtoms_confs_file, whB_bindingConf_file, str(nproc), str(whB_omega_rms), 'whB_omega', 'whB_omega_stdout.txt'))

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
        
        #print 'numAtoms in whB+3atoms = ', mol.GetMaxAtomIdx()

        print mol.GetTitle()
        allWarheadBCoords = []
        coordDict = mol.GetCoords()
        for j, atom in enumerate(mol.GetAtoms()):
            allWarheadBCoords.append(coordDict[j])
        #print 'numatoms = ', j+1
        npAllWarheadBCoords = np.mat(allWarheadBCoords)
        
        relevantLinkers = OEGetSDData(mol, 'structureForLinkers').split("\n")
        print relevantLinkers
        
        for alinker in relevantLinkers:
            print alinker
            if alinker not in whBConfCoordDict.keys():
                whBConfCoordDict[alinker] = {}
            ss = OESubSearch(smartsWHBDict[alinker])
            OEPrepareSearch(mol, ss)
            
            matchingAtomIndices = []
            for match in ss.Match(mol):
                #print 'matchingAtoms'
                for ma in match.GetAtoms():
                    matchingAtomIndices.append(ma.target.GetIdx())
            threeAtomCoords_whB = np.mat([coordDict[matchingAtomIndices[0]], coordDict[matchingAtomIndices[1]], coordDict[matchingAtomIndices[2]]])

            whBConfCoordDict[alinker][i] = [threeAtomCoords_whB, npAllWarheadBCoords, OEGraphMol(mol)]
            
    if not os.path.isdir(pdbOutDir):
        os.mkdir(pdbOutDir)

    chunksize = 10
    print 'original chunksize nproc', chunksize, nproc
    molfile = oemolistream(whA_linkers_confs_file)
    molList = []
    for i,mol in enumerate(molfile.GetOEGraphMols()):
        molList.append((i, OEGraphMol(mol)))
    numConfs = len(molList)
    
    if numConfs < nproc:
        nproc=numConfs

    if numConfs/float(nproc) < chunksize:
        chunksize = int(math.ceil(numConfs/float(nproc)))

    print 'final chunksize nproc', chunksize, nproc
    
    Pool = mp.Pool(processes=nproc)
    poolPart = partial(buildComplexes, whBConfCoordDict=whBConfCoordDict, smartsWHADict=smartsWHADict, smirksDict=smirksDict, npAllwhBProtCoords=npAllwhBProtCoords, protMol=protMol, whAprotMol=whAprotMol, nprotligclashes=nprotligclashes, nprotprotClashes=nprotprotClashes)

    results = Pool.map_async(poolPart, iterable=molList, chunksize=chunksize)
    allResults = results.get()

    outmols = oemolostream(outMolsFile)
    outComplexes = oemolostream(outComplexesFile)
    outproteins = oemolostream(outproteinsFile)
    outResults = open(outDataFile, 'w')
    outResults.write("complexNum\tlinkerName\twhAconf\twhBconf\tnumInternalCTMclashes\tprotProtPassBool\tprotLigPassBool\n")
    print outMolsFile, outComplexesFile, outproteinsFile, outDataFile
    print len(allResults)

    ##              0       1               2               3               4                   5               6           7       8
    ##results = [ctmMol, complexMol, complexWCTMMol, protProtPassBool, protLigPassBool, numInternalCTMclashes, linkerName, whAconf, whBconf]    

    complexNum = 0
    ##writeListForCTMpdbs = [ [ mol, complexNum], ...]
    writeListForCTMpdbs = []
    for oneWorker in allResults:
        for oneResult in oneWorker:
            ctmMol, complexMol, complexWCTMMol, protProtPassBool, protLigPassBool, numInternalCTMclashes, linkerName, whAconf, whBconf = oneResult
            outResults.write("\t".join(['complex'+str(complexNum), linkerName, str(whAconf), str(whBconf), str(numInternalCTMclashes), str(protProtPassBool), str(protLigPassBool)]))
            outResults.write("\n")
            print whAconf, whBconf#, numInternalCTMclashes, args.maxInternalCTMclashes
            if numInternalCTMclashes <= args.maxInternalCTMclashes:
                #print 'begin ctm mol write', time.time()
                OEWriteMolecule(outmols, ctmMol)
                #print 'end ctm mol write', time.time()
                
            if protProtPassBool:
                #print 'being prot prot mol write', time.time()
                OEWriteMolecule(outproteins, complexMol)
                #print 'end prot prot mol write', time.time()

            if (numInternalCTMclashes <= args.maxInternalCTMclashes) and protProtPassBool and protLigPassBool:
                OEWriteMolecule(outComplexes, complexWCTMMol)
                writeListForCTMpdbs.append([OEGraphMol(complexWCTMMol), complexNum])

                
##                print 'begin whole ctm initial PDB write', time.time()
##                pdbOutFileName = pdbOutDir+'complex{}.pdb'.format(str(complexNum))
##                pdbOutFile = oemolostream(pdbOutFileName)
##                OEWriteMolecule(pdbOutFile, complexWCTMMol)
##                pdbOutFile.flush()
##                pdbOutFile.close()
##                print 'end whole ctm initial PDB write', time.time()
##                
##
##                print 'begin reformat', time.time()
##                reformatFile = open(pdbOutFileName, 'r')
##                text = reformatFile.read().strip("\n")
##                reformatFile.close()
##                firstChunk, secondChunk = text.split("\nENDMDL")
##                replacedSecondChunk = secondChunk.replace('ATOM  ', 'HETATM')
##                outText = "".join([firstChunk, replacedSecondChunk])
##                reformatFile = open(pdbOutFileName, 'w')
##                reformatFile.write(outText)
##                reformatFile.flush()
##                reformatFile.close()
##                print 'end reformat', time.time()
            complexNum += 1

    writePoolPart = partial(writeOutPDBsForComplexes, pdbOutDir=pdbOutDir)
    writepool = mp.Pool(processes=nproc)
    writeresults = writepool.map_async(writePoolPart, iterable=writeListForCTMpdbs)#, chunksize=chunksize)
    writeresults.get()

    outmols.flush()
    outmols.close()
    outComplexes.flush()
    outComplexes.close()
    outproteins.flush()
    outproteins.close()
    outResults.flush()
    outResults.close()


    if args.rocsQuery:
        print 'Performing ROCS run using query', args.rocsQuery
        print "/home/pnovick/installers/openeye/bin/rocs -query {} -dbase {} -besthits 1 -oformat sdf  -mpi_np {} -prefix rocs_vs_{} &> {}".format(args.rocsQuery, outMolsFile, str(nproc), args.rocsQuery.split('.')[0], 'rocs_stdout.txt')
        os.system("/home/pnovick/installers/openeye/bin/rocs -query {} -dbase {} -besthits 1 -oformat sdf  -mpi_np {} -prefix rocs_vs_{} &> {}".format(args.rocsQuery, outMolsFile, str(nproc), args.rocsQuery.split('.')[0], 'rocs_stdout.txt'))
        

    
if __name__ == '__main__':
    main()
