import os
import glob
import pickle
import numpy
import pandas as pd
import statsmodels.formula.api as sm
import pylab

def coeffs(val):
    if val==0:
        label='Balanced'
    if val==1:
        label='Electrostatics'
    if val==2:
        label='Hydrophobics'
    return label


def get_data(dir):
    keys=['contacts', 'rmsd', 'cluster' ] #'charmm', 'prime_mini_i_energy',] # 'prime_hybridmc']

    data=dict()
    for subdir in glob.glob('./*mae*/'):
        name='%s-%s' % (dir.rstrip('/'), subdir.split('/')[1].split('rmsd_')[0])
        data[name]=dict()
        for key in keys:
            summary_files=glob.glob('%s/%s_model.*.summary.txt' % (subdir, key))
            if len(summary_files)==0:
                print "MISSING SUMMARY TXT FILES for %s %s. RUN PROPLISTER" % (subdir, key)


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
                print file
                names=numpy.loadtxt(file, usecols=(0,), dtype=str)
                vals=numpy.loadtxt(file, usecols=(1,), dtype=float)
                for (n,v) in zip(names,vals):
                    if key=='cluster':
                        if v < 0:
                            import pdb
                            pdb.set_trace()
                    #print n, file
                    #print n
                    n='model%s' % n.split('.mae')[0].split('model')[1].split('.interaction')[0].split('.pdb')[0]
                    data[name][model][key][n]=v
    
    return data

def get_df(data, global_name):
    metadata=dict()
    models=['Balanced',] # 'Electrostatics', 'Hydrophobics']
    
    for model in models:
        metadata[model]=dict()
    corr_keys=['contacts', 'rmsd', 'cluster'] #'charmm', 'prime_mini_i_energy',] # 'prime_hybridmc']

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
    pylab.show(False)
    for model in metadata.keys():
        results[model]=dict()
        for name in metadata['Balanced'].keys():
            results[model][name]=dict()
            df=pd.DataFrame.from_dict(metadata[model][name]) #, orient='index')
            for ckey in corr_keys:
                if ckey not in df.columns:
                    continue
                if ckey=='rmsd':
                    results[model][name]['rmsdmin']=min(df['rmsd'])
                    results[model][name]['rmsdavg']=round(numpy.mean(df['rmsd']),2)
                elif ckey=='contacts':
                    results[model][name]['contactsmax']=max(df['contacts'])
                    results[model][name]['contactsavg']=round(numpy.mean(df['contacts']),2)
                    native= df['contacts']>0
                    results[model][name]['contacts_above_0.0']=len(df.loc[native])
                    native= df['contacts']>0.45
                    results[model][name]['contacts_above_0.5']=len(df.loc[native])
                    ranks=[5,10,15, 20]
                    #energy_df=df.sort_values('i_energy', ascending=True)
                    for rank in ranks:
                        test=df.iloc[:rank].loc[native]
                        if test['contacts'].any():
                            results[model][name]['contacts_0.5_rank%s' % rank]=int(len(test))
                            break
                else:
                    result=sm.ols(formula="%s ~ contacts" % ckey, data=df).fit()
                    print "correlation", result.rsquared, ckey
                    #results[model][name][ckey]=round(result.rsquared,2)
                    #fig=pylab.figure()
                    #df.plot.scatter(x='contacts', y=ckey, marker='o', color='k')
                    #pylab.plot(df['contacts'].values, result.predict(), 'r-', label='%0.2f' % result.rsquared)
                    #pylab.title('%s-%s' % (ckey, name))
                    #pylab.legend(loc=2)
                    #pylab.savefig('%s-%s.png' % (ckey, name), dpi=300)

    df=pd.DataFrame.from_dict(results['Balanced'], orient='index')
    return df

def main():
    df_list=[]
    for dir in numpy.loadtxt('dirs', dtype=str):
        print "on %s" % dir
        topdir='/home/mlawrenz/PIPER_WORK_2018/validation_2018/vhl_bromodomain/%s/results/' % dir.rstrip('/')
        global_name='vhl_bromodomain_%s' % dir.rstrip('/')
        os.chdir(topdir)
        print os.getcwd()
        data=get_data(dir)
        df=get_df(data, global_name)
        df_list.append(df)
    all_df=pd.concat(df_list) 
    print all_df
    all_df.to_csv('vhl_bromodomain_tests.noenergy.csv')
    return

if __name__=="__main__":
    main()
