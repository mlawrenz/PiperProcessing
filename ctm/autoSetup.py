from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import Descriptors
import argparse
import sys
import multiprocessing as mp
import os

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputCTMs', default=False)
    parser.add_argument('--whA_bindingConf_file')
    parser.add_argument('--whB_bindingConf_file')
    parser.add_argument('--outprefix', default='autoSetupOut_')
    args = parser.parse_args()


    whAbindingConf = args.whA_bindingConf_file
    whBbindingConf = args.whB_bindingConf_file
    ctmFile = args.inputCTMs

    whA = Chem.MolFromMolFile(whAbindingConf)
    whB = Chem.MolFromMolFile(whBbindingConf)

    for mol in Chem.SDMolSupplier(whBbindingConf):#, removeHs=False):
        #mol.assignStereochemistryFrom3D()
        print Chem.MolToSmiles(mol, isomericSmiles=True)
        whB = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
    
    print 'whA smiles', Chem.MolToSmiles(whA, isomericSmiles=True)
    print 'whB_smiles', Chem.MolToSmiles(whB, isomericSmiles=True)

    atomNumToIsotopeDict = { 6:'[14CH3]', 7:'[15NH2]', 8:'[17OH1]', 16:'[33SH1]' }
    atomNumToElementDict = { 6:['C', 'c'], 7:['N', 'n'], 8:['O', 'o'], 16:['S', 's']}

    ligPrefix = None

    outWhALinkers = Chem.SDWriter(args.outprefix+'whA_linkers.sdf')
    outWhBAttach = Chem.SDWriter(args.outprefix+'whB_attach.sdf')
    outSmirks = open(args.outprefix+'smirks.smk', 'w')
    

    if ctmFile:
        ctmFileRoot = ctmFile.split('/')[-1]
        for i, ctm in enumerate(Chem.SDMolSupplier(ctmFile)):
            print Chem.MolToSmiles(ctm, True)
            if i == 0:
                ligTitle = ctm.GetProp('_Name') 
                if ligTitle == ctmFileRoot:
                    ligPrefix = 'linkerNum'

            if ligPrefix:
                ligName = ligPrefix+str(i)
            else:
                ligName = ctm.GetProp('_Name')
            
            ## first, separate the ctm into whA+linker and whB.  
            ## Dummy atom attachment points are left

            whA_plus_linker = Chem.ReplaceCore(ctm, whB)
            whB_plus_attach = Chem.ReplaceSidechains(ctm, whB)
            #print Chem.MolToSmiles(whA_plus_linker, isomericSmiles=True)
            #print Chem.MolToSmiles(whB_plus_attach, isomericSmiles=True)

            ## we are going to replace the dummy atoms with the appropriate isotope
            ## first, find out the atom on whB bonded to attachment point
            [x.GetAtomicNum() for x in whB_plus_attach.GetAtomWithIdx(whB_plus_attach.GetSubstructMatch(Chem.MolFromSmarts('[1*]'))[0]).GetNeighbors()]
            #print [x.GetAtomicNum() for x in whB_plus_attach.GetAtomWithIdx(whB_plus_attach.GetSubstructMatch(Chem.MolFromSmarts('[1*]'))[0]).GetNeighbors()]
            #neighborAtomNum = whB_plus_attach.GetAtomWithIdx(0).GetNeighbors()[0].GetAtomicNum()
            neighborAtomNum = whB_plus_attach.GetAtomWithIdx(whB_plus_attach.GetSubstructMatch(Chem.MolFromSmarts('[1*]'))[0]).GetNeighbors()[0].GetAtomicNum()
            neighborAtom_isAro = whB_plus_attach.GetAtomWithIdx(whB_plus_attach.GetSubstructMatch(Chem.MolFromSmarts('[1*]'))[0]).GetNeighbors()[0].GetIsAromatic()
            neighborAtom_implicitVal = whB_plus_attach.GetAtomWithIdx(whB_plus_attach.GetSubstructMatch(Chem.MolFromSmarts('[1*]'))[0]).GetNeighbors()[0].GetImplicitValence()
            #print neighborAtomNum, neighborAtom_isAro

            ## then, replace the dummy atom on whA+linker with the right isotope
            whA_plus_linker_C14 = Chem.MolFromSmiles(Chem.MolToSmiles(whA_plus_linker, isomericSmiles=True).replace('[1*]', atomNumToIsotopeDict[neighborAtomNum]))

            ## now, figure out what the two atoms are that we will put on whB
            firstAtom_idx = [x.GetIdx() for x in whA_plus_linker.GetAtomWithIdx(whA_plus_linker.GetSubstructMatch(Chem.MolFromSmarts('[1*]'))[0]).GetNeighbors()][0]
            firstAtom_atomNum = [x.GetAtomicNum() for x in whA_plus_linker.GetAtomWithIdx(whA_plus_linker.GetSubstructMatch(Chem.MolFromSmarts('[1*]'))[0]).GetNeighbors()][0]
            firstAtom_isAro = whA_plus_linker.GetAtomWithIdx(firstAtom_idx).GetIsAromatic()
            firstAtom_implicitVal = whA_plus_linker.GetAtomWithIdx(firstAtom_idx).GetImplicitValence()

            nextAtoms_idx = [x.GetIdx() for x in whA_plus_linker.GetAtomWithIdx(firstAtom_idx).GetNeighbors()]
            nextAtoms_atomNum = [x.GetAtomicNum() for x in whA_plus_linker.GetAtomWithIdx(firstAtom_idx).GetNeighbors()]

            for i, atom in enumerate(nextAtoms_atomNum):
                if atom == 6 or atom == 7:
                    secondAtom_idx = nextAtoms_idx[i]
                    secondAtom_atomNum = atom
                    break

            secondAtom_isAro = whA_plus_linker.GetAtomWithIdx(secondAtom_idx).GetIsAromatic()
            secondAtom_implicitVal = whA_plus_linker.GetAtomWithIdx(secondAtom_idx).GetImplicitValence()

            #whB_plus_attach_C14 = AllChem.ReplaceSubstructs(whB_plus_attach, Chem.MolFromSmiles('[1*]'), Chem.MolFromSmiles('CC'))[0]
            replaceString = atomNumToIsotopeDict[secondAtom_atomNum]+atomNumToElementDict[firstAtom_atomNum][int(firstAtom_isAro)]
            whB_plus_attach_C14 = Chem.MolFromSmiles( Chem.MolToSmiles(whB_plus_attach, isomericSmiles=True).replace('[1*]', replaceString) )

            ## make smirks
            whA_linker_smirks = atomNumToIsotopeDict[neighborAtomNum] + '[' + atomNumToElementDict[firstAtom_atomNum][int(firstAtom_isAro)] + 'H{}:1]['.format(firstAtom_implicitVal) + atomNumToElementDict[secondAtom_atomNum][int(secondAtom_isAro)] + 'H{}:2]'.format(secondAtom_implicitVal)
            whB_linker_smirks = '[' + atomNumToElementDict[neighborAtomNum][neighborAtom_isAro] + 'H{}:3]'.format(neighborAtom_implicitVal) + atomNumToElementDict[firstAtom_atomNum][int(firstAtom_isAro)] + atomNumToIsotopeDict[secondAtom_atomNum]
            productSmirks = '[' + atomNumToElementDict[neighborAtomNum][neighborAtom_isAro] + ':3]' + '[' + atomNumToElementDict[firstAtom_atomNum][int(firstAtom_isAro)] + ':1][' + atomNumToElementDict[secondAtom_atomNum][int(secondAtom_isAro)] + ':2]'
            smirksText = whA_linker_smirks + '.' + whB_linker_smirks + '>>' + productSmirks

            whA_plus_linker_C14.SetProp('linker', ligName)
            AllChem.Compute2DCoords(whA_plus_linker_C14)
            #print Chem.MolToSmiles(whA_plus_linker_C14)
            outWhALinkers.write(whA_plus_linker_C14)

            whB_plus_attach_C14.SetProp('structureForLinkers', ligName)
            AllChem.Compute2DCoords(whB_plus_attach_C14)
            #print Chem.MolToSmiles(whB_plus_attach_C14)
            outWhBAttach.write(whB_plus_attach_C14)

            outSmirks.write("{}\t{}\n".format(smirksText, ligName))


    outSmirks.flush()
    outSmirks.close()
    outWhALinkers.flush()
    outWhALinkers.close()
    outWhBAttach.flush()
    outWhBAttach.close()
            
    

if __name__=='__main__':
    main()
