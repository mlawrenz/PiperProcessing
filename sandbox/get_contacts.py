import os
import glob
import pandas as pd
import statsmodels.formula.api as sm
import pylab, numpy

topdir='/home/mlawrenz/PIPER_WORK_2018/validation_2018/vhl_bromodomain/ligand_restraints_20_only/mmgbsa/'
global_name='vhl_bromodomain_ligand_restraints_20_only'
#topdir='/home/mlawrenz/PIPER_WORK_2018/validation_2018/vhl_bromodomain/protein_restraints_only/mmgbsa/'
#global_name='vhl_bromodomain_protein_restraints_only'
#native='5T35-Brd4-BRD2_JQpr'
#topdir='/home/mlawrenz/PIPER_WORK_2018/validation_2018/cbln_bromodomain/ligand_restraints_20_only/mmgbsa'
#global_name='cbln_bromodomain_ligand_restraints_20_only'

os.chdir(topdir)
print os.getcwd()
# RUN PROPLISTER TO GET DATAFILES AND THEN LOAD HERE. ITS JUST EASIER.
# RUN SPLIT ON MINIMIZED COMPLEXES AND THEN RUN PRIME CALC ENERGY TO GET LOG
# FILES, PARSE THESE
# RUN PROTEIN INTERACTION SCRIPT AND PARSE THESE
#============================================================


# COMPUTE INTERACTION ENERGIES
dir='5T35-Brd4-BRD2_JQpr__clean-prep-aligned_vhl_5T35_V02methyl_rmsd-maefiles'
os.chdir(dir)
minimized_interaction_energy=numpy.loadtxt('all.energy.txt', dtype=str)
print minimized_interaction_energy
#============================================================

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

reference='/home/mlawrenz/PIPER_WORK_2018/validation_2018/vhl_bromodomain/Complex_aligned/ref-clean-prep-aligned-5T35-Brd4-BRD2_5T35.interaction.csv'
df=pd.read_csv(reference, skipinitialspace=True)
native=get_native_contacts(df)
for contact in native:
    print contact
print "%s total %s contacts" % (reference, len(native))
total_native=len(native)    
    
#============================================================

interaction_files=glob.glob('./*model.000.*interaction.csv')
fraction_native_contacts=[]
for file in interaction_files:
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
    fraction_native_contacts.append((file, float(count)/total_native)) 

#============================================================

rmsddata=numpy.loadtxt('rmsd_model.000.summary.txt', dtype=str)

#============================================================

data=dict()
for entry in rmsddata:
    new='model.%s' % entry[0].split('model.')[1].split('.min')[0]
    rmsd=float(entry[1])
    data[new]=dict()
    data[new]['rmsd']=rmsd
for entry in minimized_interaction_energy:
    new='model.%s' % entry[0].split('model.')[1].split('.min')[0]
    data[new]['i_energy']=float(entry[1])
for entry in fraction_native_contacts:
    new='model.%s' % entry[0].split('model.')[1].split('.min')[0]
    data[new]['f_native']=float(entry[1])
df=pd.DataFrame.from_dict(data, orient='index')
print df.sort_values('i_energy', ascending=True)

#============================================================

result=sm.ols(formula="f_native ~ rmsd", data=df).fit()
print result.rsquared
fig=pylab.figure()
df.plot.scatter(x='rmsd', y='f_native', marker='o', color='k')
pylab.plot(df['rmsd'].values, result.predict(), 'r-', label='%0.2f' %
result.rsquared)
pylab.title('RMSD vs Native Contact')
pylab.legend(loc=2)

result=sm.ols(formula="f_native ~ i_energy", data=df).fit()
print result.rsquared
fig=pylab.figure()
df.plot.scatter(x='i_energy', y='f_native', marker='o', color='k')
pylab.plot(df['i_energy'].values, result.predict(), 'r-', label='%0.2f' %
result.rsquared)
pylab.title('I_energy vs Native Contact')
pylab.legend(loc=2)


pylab.show()
