import glob
import pickle
import os
import numpy
def coeffs(val):
    if val==0:
        label='Balanced'
    if val==1:
        label='Electrostatics'
    if val==2:
        label='Hydrophobics'
    return label

keys=['rmsd', 'charmm', 'cluster', 'prime_mini', 'prime_hybridmc']
#keys=['rmsd', 'prime_mini', 'prime_hybridmc']

if os.path.exists('data.pickle'):
    ohandle=open('data.pickle', 'rb')
    data=pickle.load(ohandle)
    ohandle.close
else:
    data=dict()
    for dir in glob.glob('./*/'):
        name=dir.split('/')[1].split('rmsd_')[0]
        data[name]=dict()
        for key in keys:
            summary_files=glob.glob('%s/%s_model.*.summary.txt' % (dir, key))
            for file in summary_files:
                m=os.path.basename(file)
                if 'model.000' in m:
                    model=coeffs(0)
                if 'model.001' in m:
                    model=coeffs(1)
                if 'model.002' in m:
                    model=coeffs(2)
                if 'model.003' in m:
                    print m, "skip these results"
                    continue
                if model not in data[name].keys():
                    data[name][model]=dict()
                if key not in data[name][model].keys():
                    data[name][model][key]=dict()
                names=numpy.loadtxt(file, usecols=(0,), dtype=str)
                vals=numpy.loadtxt(file, usecols=(1,), dtype=float)
                for (n,v) in zip(names,vals):
                    if key=='cluster':
                        if v < 0:
                            import pdb
                            pdb.set_trace()
                    data[name][model][key][n]=v
    ohandle=open('data.pickle', 'wb')
    pickle.dump(data, ohandle)
    ohandle.close()


##########################################
import pandas as pd
import statsmodels.formula.api as sm
import pylab

metadata=dict()

models=['Balanced',] # 'Electrostatics', 'Hydrophobics']

for model in models:
    metadata[model]=dict()

corr_keys=['rmsd', 'charmm', 'cluster', 'prime_mini', 'prime_hybridmc']
#for key in metadata.keys():
#    for ckey in corr_keys:
#        metadata[key][ckey]=dict()



for model in metadata.keys():
    for name in data.keys():  #names of proteins docked
        if model in data[name].keys():
            if name not in metadata[model].keys():
                metadata[model][name]=dict()
            for ckey in data[name][model].keys():
                if ckey not in metadata[model][name].keys():
                    metadata[model][name][ckey]=dict()
                for entry in sorted(data[name][model][ckey].keys()):
                    val=data[name][model][ckey][entry]
                    metadata[model][name][ckey][entry]=val
        else:
            pass

results=dict()

for model in metadata.keys():
    results[model]=dict()
    for name in metadata['Balanced'].keys():
        results[model][name]=dict()
        df=pd.DataFrame.from_dict(metadata[model][name]) #, orient='index')
        print df.head()
        for ckey in corr_keys:
            if ckey=='rmsd':
                results[model][name]['rmsdmin']=min(df['rmsd'])
                results[model][name]['rmsdavg']=round(numpy.mean(df['rmsd']),2)
                results[model][name]['rmsdmax']=numpy.max(df['rmsd'])
            else:
                result=sm.ols(formula="%s ~ rmsd" % ckey, data=df).fit()
                print result.rsquared, ckey
                results[model][name][ckey]=round(result.rsquared,2)

df=pd.DataFrame.from_dict(results['Balanced'], orient='index')
print df

