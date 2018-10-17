#!/usr/bin/python3
''' 
A simple tool to show the spatial configuration of partitions.
Looks into directories and looks for the first file that matches the pattern
'reentry.con.npart.*', and for reentry.plt.
See visualize.py. 
'''

import visualize as vs
from sys import argv


def get_parts(dirname):
    pd_df, nxy_df = vs.geom_info(dirname)
    
    part_dfs = [ pd_df.loc[ pd_df['P'] == i ] for i in range(int(pd_df.values.min()),int(pd_df.values.max()+1))]
    parts = [ part_df.join(nxy_df)[['x','y']].values for part_df in part_dfs ]
    return parts
    
parts1 = get_parts(argv[1])
# showing partitions for first directory
for part in parts1:
    vs.plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
vs.plt.title(argv[1])

do_second_file = True
try:
    parts2 = get_parts(argv[2])
except:
    print("Second filename not given, or problems in opening it")
    do_second_file = False

if do_second_file:
    # comparison plot
    vs.plt.figure()
    for part in parts1:
        vs.plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
    
    yshift = 0.01
    xshift = 0.01
    
    for part in parts2:
        vs.plt.plot(part[:,0]+xshift,part[:,1]+yshift,linestyle='None',marker='.')
    
    vs.plt.title('"'+argv[1]+ '"(+) vs "' + argv[2] +'"(.)')

    # showing partitions for second directory
    vs.plt.figure()
    for part in parts2:
        vs.plt.plot(part[:,0],part[:,1],linestyle='None',marker='.')
    vs.plt.title(argv[2])
    

vs.plt.show()
