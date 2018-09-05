import os
import sys
import pickle
import argparse
import glob
import pandas as pd
import pylab, numpy




def main():
    with open('dirs') as f: dirs = [line.strip('/').split()[0] for line in f.readlines()]
    columns=['total_score', 'fa_atr', 'fa_dun',  'fa_elec', 'fa_rep',   'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb'   , 'hbond_sc', 'lk_ball_wtd']
    i_e_data=dict()
    ohandle=open('rosetta_interaction_energy_data.txt', 'w')
    line=['name',] + sorted(columns) + ['\n']
    ohandle.write(','.join(line))
    for dir in dirs:
        infile='%s/score_only.sc' % dir 
        df=pd.read_csv(infile, delimiter=" ", skipinitialspace=True, header=1)
        if 'chain' in df['description'][0]:
            print "check order of files in score file"
            sys.exit()
        if 'chain' not in df['description'][1]:
            print "check order of files in score file"
            sys.exit()
        if 'chain' not in df['description'][2]:
            print "check order of files in score file"
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
        print i_e_data
        line=[dir,] + [str(round(i_e_data[col],2)) for col in sorted(columns)] + ['\n']
        ohandle.write(','.join(line))
    return

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='parse rosetta score file from interaction energy calculation')
    parser.add_argument('--debug',  action="store_true", dest='debug' )
    #parser.add_argument('-r', dest='reference', help='INTERACTION CSV from the reference (native) structure')



    args = parser.parse_args()
   

    if args.debug:
        import pdb
        pdb.set_trace()
    main()



