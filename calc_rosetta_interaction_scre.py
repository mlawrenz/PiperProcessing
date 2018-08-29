import numpy
import sys

def check_title_of_array(array, title):
    # making sure you are loading in what you think from score file
    if title in array[0]:
        array=array[1:]
    else:
        print "column name problem for %s" % title
    return array

def get_raw_data(property_keys, infile):
    raw_data=dict()
    for key in sorted(property_keys.keys()):
        with open(infile) as f: values = [line.strip('\n').split()[key] for line in f.readlines()] 
        raw_data[key]=check_title_of_array(values, property_keys[key])
    return raw_data

def main():
    property_keys=dict()
    property_keys[1]='total_score'
    property_keys[4]='fa_atr'
    property_keys[10]='fa_sol'
    property_keys[11]='hbond_bb_sc'
    property_keys[12]='hbond_lr_bb'
    property_keys[13]='hbond_sc'
    property_keys[25]='description'
    infile='score_only.sc'
    interaction_data=dict()    
    raw_data=get_raw_data(property_keys, infile)
    names=raw_data[25]
    raw_data.pop(25)
    property_keys.pop(25)
    for key in sorted(property_keys.keys()):
        prevname=0
        for (i, name) in enumerate(names):
            if 'chain' not in name:
                complex=True
                if name not in interaction_data.keys():
                    interaction_data[name]=dict()
                else:
                    pass
                if prevname:
                    if name.split('_split')[0]!=prevname:
                        # new complex evaluated
                        interaction_data[prevname]['%s_i_energy' % property_keys[key]]=complex_score-(monomer_score[0]+monomer_score[1])
                        prevname=name
                else:
                    prevname=name
            else:
                complex=False
            if complex==True:
                complex_score=float(raw_data[key][i])
                monomer_score=dict()
                n=0
            else:
                if name.split('.split')[0]==prevname.split('_split')[0]:
                    monomer_score[n]=float(raw_data[key][i])
                    n+=1
                else:
                    print("name issue %s %s" % (name.split('.split')[0], prevname))
                    sys.exit()
        # deal with last entry
        interaction_data[prevname]['%s_i_energy' % property_keys[key]]=complex_score-(monomer_score[0]+monomer_score[1])
        #
    for key in property_keys.keys():
        ohandle=open('%s_i_energy_model.000.summary.txt' % property_keys[key], 'w')
        for name in sorted(interaction_data.keys()):
            title=name.split('_')[0].split('complex-')[1]
            ohandle.write('%s\t%s\n' % (title, interaction_data[name]['%s_i_energy' % property_keys[key]]))
        ohandle.close()
    return



if __name__=="__main__":
    main()
