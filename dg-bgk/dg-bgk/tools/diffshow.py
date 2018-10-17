#!/usr/bin/python3
''' 
    Simple tool to visualize the relative and absolute difference in 
    data from RESULTS1.RES and RESULTS2.RES.
    Takes two directory names as arguments which must contain these two files, 
    plus the mesh data file and the metis partition file. See visualize.py.
'''
import visualize as vs
from sys import argv
pd_df, nxy_df = vs.geom_info(argv[1])
newr1s = vs.parse_results_all(argv[1])
refr1s = vs.parse_results_all(argv[2])
d_df , rd_df = vs.diff(newr1s[10],refr1s[10])

part_dfs = [ pd_df.loc[ pd_df['P'] == i ] for i in range(int(pd_df.values.min()),int(pd_df.values.max()+1))]
parts = [ part_df.join(nxy_df)[['x','y']].values for part_df in part_dfs ]


def show(col, pm = False) :
    vs.plt.figure(1)
    vs.plot_stuff(nxy_df,d_df,col)
    vs.plt.title("Delta, " + col)
    if pm:
        for part in parts:
            vs.plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
    vs.plt.figure(2)
    
    vs.plot_stuff(nxy_df,rd_df,col)
    vs.plt.title("Relative Delta, " + col)
    if pm: 
        for part in parts:
            vs.plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
    vs.plt.show()


show('ND',True)
show('U',True)
show('V',True)
show('RHO',True)
show('PS',True)
show('TEMP',True)
