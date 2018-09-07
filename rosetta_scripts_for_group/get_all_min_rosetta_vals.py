import os
import sys
import pickle
import argparse
import glob
import pandas as pd
import pylab, numpy


def get_interaction_energies(infile, columns):
    i_e_data=dict()
    df=pd.read_csv(infile, delimiter=" ", skipinitialspace=True, header=1)
    if 'chain' in df['description'][0]:
        print("check order of files in score file", infile)
        print(df['description'][0],  df['description'][1], df['description'][2])

        sys.exit()
    if 'chain' not in df['description'][1]:
        print("check order of files in score file", infile)
        print(df['description'][0],  df['description'][1], df['description'][2])

        sys.exit()
    if 'chain' not in df['description'][2]:
        print("check order of files in score file", infile)
        print(df['description'][0],  df['description'][1], df['description'][2])

        sys.exit()
    i_e_data=dict()
    for col in df.columns:
        if 'SCORE' in col:
            continue
        if 'desc' in col:
            continue
        if 'Unname' in col:
            continue
        i_e_data[col]=float(df[col][0])-(float(df[col][1])+float(df[col][2]))
    return i_e_data


def main():
    columns=['total_score', 'fa_atr', 'fa_dun',  'fa_elec', 'fa_rep',   'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb'   , 'hbond_sc', 'lk_ball_wtd']
    dirs=glob.glob('./com*/')
    #ohandle=open('%s-rosetta_interaction_energy_data.txt' % outname, 'w')
    line=['name',] + sorted(columns) + ['\n']
    min_i_e_data=dict()
    for dir in sorted(dirs):
        min_i_e_data[dir]=dict()
        min_i_e_data[dir]['total_score']=1000
        for infile in glob.glob('%s/rosetta*.sc' % dir):
            i_e_data=get_interaction_energies(infile, columns)
            if i_e_data['total_score'] <= min_i_e_data[dir]['total_score']:
                for col in columns:
                    min_i_e_data[dir][col]= i_e_data[col]
            else:
                pass 
    ohandles=dict()
    for key in sorted(min_i_e_data.keys()):
        for col in columns:
            if col not in ohandles.keys():
                ohandles[col]=open('rosetta_%s_model.000.summary.txt' % col, 'w')
            ohandles[col].write('%s\t%s\n' % (key.rstrip('/'), (min_i_e_data[key][col])))


    for key in ohandles.keys():
        ohandles[key].close()
    #line=[outname,] + [str(round(i_e_data[col],2)) for col in sorted(columns)] + ['\n']
    return



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='parse rosetta score file from interaction energy calculation')
    #parser.add_argument('-i', dest='infile', help='input file score_only.sc with only complex and two monomer energies')
    #parser.add_argument('-o', dest='outname', help='output name prefix for rosetta interaction energy data')
    parser.add_argument('--debug',  action="store_true", dest='debug' )
    args = parser.parse_args()

    if args.debug:
        import pdb
        pdb.set_trace()
    main()
