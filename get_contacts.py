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
def main(listfile, reference):
    df=pd.read_csv(reference, skipinitialspace=True)
    native=get_native_contacts(df)

    for contact in native:
        print contact
    print "%s total %s contacts" % (reference, len(native))
    total_native=len(native)    
    print "total native %s" % total_native
    
    fraction_native_contacts=[]

    listfile = open(listfile, 'r')
    interaction_files = listfile.read().splitlines()

    for file in interaction_files:
        print file
        df=pd.read_csv(file, skipinitialspace=True)
        contacts=get_native_contacts(df)
        count=0
        for contact in contacts:
            if contact in native:
                count+=1
            elif (contact[1], contact[0]) in native:
                print "reverse counted"
                count+=1
            else:
                pass
        #print count
        fraction_native_contacts.append(('prep-%s' % file.split('prep-')[1].split('-out')[0], float(count)/total_native)) 

        dir=os.path.dirname(file)
        if dir=='':
            dir='./'
        numpy.savetxt('%s/contacts_model.000.summary.txt' % dir, fraction_native_contacts, fmt='%s')

    return

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='parse interaction csv files and use native structure to compute % native contacts')
    parser.add_argument('-r', dest='reference', help='reference interaction csv')
    parser.add_argument('-l','--listfile', dest='listfile', help='file list of interaction.csv files')
    parser.add_argument('--debug',  action="store_true", dest='debug' )
    args = parser.parse_args()
    if args.debug:
        import pdb
        pdb.set_trace()
    main(args.listfile,args.reference)

