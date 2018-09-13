import numpy
data=dict()
import glob
files=glob.glob('rmsd*txt')
for file in files:
    name=file.split('members-')[1].split('.txt')[0]
    subnames=numpy.loadtxt(file, usecols=(0,), dtype=str)
    vals=numpy.loadtxt(file, usecols=(1,), dtype=float)
    for (subname, val) in zip(subnames, vals):
        fullname='%s-%s' % (name, subname.split('.pdb')[0])
        data[fullname]=dict()
        data[fullname]['rmsd']=val

files=glob.glob('roc*txt')
for file in files:
    name=file.split('members-')[1].split('.summary')[0]
    subnames=numpy.loadtxt(file, usecols=(0,), dtype=str)
    vals=numpy.loadtxt(file, usecols=(1,), dtype=float)
    for (subname, val) in zip(subnames, vals):
        subname=subname.split('outmos_mol-')[1].split('_VSP')[0]
        fullname='%s-%s' % (name, subname)
        if fullname not in data.keys():
            print "missing %s" % fullname
        else:
            data[fullname]['rocs']=val
import pandas as pd
df=pd.DataFrame.from_dict(data, orient='index')
print df.describe()
print df.sort_values('rocs', ascending=False)[:20]
corr=numpy.corrcoef(df['rocs'], df['rmsd'])[0][1]
print corr
sorted_df=df.sort_values('rocs', ascending=False)
ohandle=open('top50_matches.txt', 'w')
for name in sorted_df.index[0:50]:
    model=name.split('-complex')[0]
    cluster=name.split('-complex.')[1]
    file='cluster-members-%s/top-outmos_mol-complex.%s_VSPiper_hits_1.sdf'  % (model, cluster)
    ohandle.write('%s\n' % file)
ohandle.close()


