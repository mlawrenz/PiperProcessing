import os
import pickle
import argparse
import glob
import pandas as pd
import pylab, numpy



def coeffs(val):
    if val==0:
        label='Balanced'
    if val==1:
        label='Electrostatics'
    if val==2:
        label='Hydrophobics'
    return label




def get_native_contacts(df):
    contacts=[]
    boolean=pd.notnull(df['Closest'])
    for (res1, res2) in zip(df['Residue'].loc[boolean],
df['Closest'].loc[boolean]):
        res2=res2.replace('\n', ' ')
        if 'V32' in res1:
            continue
        if len(res2.split()) > 1:
            for x in res2.split():
                if (x, res1) not in contacts:
                    if 'V32' not in x:
                        contacts.append((res1, x))
                else:
                    pass
        else:
            res2=res2.split()[0]
            if (res2, res1) not in contacts:
                if 'V32' not in res2:
                    contacts.append((res1, res2))
            else:
                pass
    return contacts

#============================================================
def main(reference, listfile=None):
    df=pd.read_csv(reference, skipinitialspace=True)
    native=get_native_contacts(df)

    for contact in native:
        print(contact)
    total_native=len(native)
    print("%s total native %s contacts" % (reference, len(native)))
    
    fraction_native_contacts=[]

    if listfile:
        listfile = open(listfile, 'r')
        interaction_files = listfile.read().splitlines()
    else:
        interaction_files = sorted(glob.glob('*interaction.csv'))

    for file in interaction_files:
        #print("on %s" % file)
        df=pd.read_csv(file, skipinitialspace=True)
        contacts=get_native_contacts(df)
        count=0
        for contact in contacts:
            if contact in native:
                count+=1
            elif (contact[1], contact[0]) in native:
                #print("reverse counted")
                count+=1
            else:
                pass
        #fraction_native_contacts.append(('prep-%s' % file.split('prep-')[1].split('-out')[0], float(count)/total_native)) 
        fraction_native_contacts.append((file, float(count)/total_native)) 

        dir=os.path.dirname(file)
        if dir=='':
            dir='./'
        numpy.savetxt('%s/contacts_model.000.summary.txt' % dir, fraction_native_contacts, fmt='%s')

    surface_comp=[]
    total_hb=[]
    total_sb=[]
    total_buried_SASA=[]
    for file in interaction_files:
        df=pd.read_csv(file, skipinitialspace=True)
        #surface_comp.append(('prep-%s' % file.split('prep-')[1].split('-out')[0], sum(df['Surface\nComplementarity'])))
        surface_comp.append(('prep-%s' % file, sum(df['Surface\nComplementarity'])))
        #total_hb.append(('prep-%s' % file.split('prep-')[1].split('-out')[0], sum(df['# HB'])))
        total_hb.append(('prep-%s' % file, sum(df['# HB'])))
        #total_sb.append(('prep-%s' % file.split('prep-')[1].split('-out')[0],sum(df['# Salt\nBridges'])))
        total_sb.append(('prep-%s' % file,sum(df['# Salt\nBridges'])))
        #total_buried_SASA.append(('prep-%s' % file.split('prep-')[1].split('-out')[0], sum([float(i.split('%')[0]) for i in df['Buried\nSASA']])))
        total_buried_SASA.append(('prep-%s' % file, sum([float(i.split('%')[0]) for i in df['Buried\nSASA']])))
    numpy.savetxt('%s/surface_comp_model.000.summary.txt' % dir, surface_comp, fmt='%s')
    numpy.savetxt('%s/total_hb_model.000.summary.txt' % dir, total_hb, fmt='%s')
    numpy.savetxt('%s/total_sb_model.000.summary.txt' % dir, total_sb, fmt='%s')
    numpy.savetxt('%s/total_buried_SASA_model.000.summary.txt' % dir, total_buried_SASA, fmt='%s')


    return

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='parse interaction csv files and use native structure to compute % native contacts')
    parser.add_argument('-r', dest='reference', help='INTERACTION CSV from the reference (native) structure')

    parser.add_argument('-l','--listfile', dest='listfile', help='file list of interaction.csv files. If you do not pass this, will run on all *interaction.csv files in directory.')
    parser.add_argument('--debug',  action="store_true", dest='debug' )
    args = parser.parse_args()
    if args.debug:
        import pdb
        pdb.set_trace()
    main(args.reference, args.listfile)

