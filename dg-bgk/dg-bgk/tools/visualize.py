#!/usr/bin/python3
''' 
   Small library to visualize partitions and physical results.
'''
from sys import argv
from os import path
import numpy as np
import pandas as pd
import glob
from scipy.interpolate import griddata
from matplotlib import pyplot as plt


# MESH INFO
def geom_info(dirname):
    '''
    directory name as input, must contain mesh file and metis partition file
    '''
    f = open(path.join(dirname,'reentry.plt'))
    
    nxy_data= f.readlines()
    f.close()
    start = [ i for i,t in enumerate(nxy_data) if 'coordinates' in t ][0]+1
    end = [ i for i,t in enumerate(nxy_data) if 'boundary sides' in t ][0]
    nxy_data = '\n'.join(nxy_data[start:end])
    nxy_data = np.fromstring(nxy_data,sep=' ').reshape((-1,5))[:,0:3]
    nxy_df = pd.DataFrame(data=nxy_data[:,1:],index = nxy_data[:,0].astype(np.int),columns=['x','y'])
    #PARTITIONING INFO
    part_data = np.loadtxt(glob.glob(path.join(dirname,'reentry.con.npart.*'))[0])
    pd_df = pd.DataFrame(data=part_data,index = np.arange(1,part_data.size+1), columns = ['P'])
    
    return pd_df,nxy_df


# RESULTS 
def parse_results(file_name,columns):
    f = open(file_name)
    r_data = f.readlines()
    f.close()
    idxs = [i for i,t in enumerate(r_data) if 'RESULTS' in t] + [len(r_data)]
    starts = idxs[:-1]
    ends = idxs[1:]
    results = []
    for s,e in zip(starts,ends):
        data = '\n'.join(r_data[s+1:e])
        data = np.fromstring(data,sep=' ').reshape((-1,4))
        results.append(data)
    res_dfs = [ pd.DataFrame(data=r_data[:,1:],
                        index = r_data[:,0].astype(np.int),
                        columns=columns) for r_data in results ] 

    return res_dfs

def parse_results_all(dirname):
    r1file_name = path.join(dirname,'RESULTS1.RES')
    r2file_name = path.join(dirname,'RESULTS2.RES')
    r1dfs = parse_results(r1file_name,['ND','U','V'])
    r2dfs = parse_results(r2file_name,['RHO','PS','TEMP'])

    return [r1.join(r2) for r1,r2 in zip(r1dfs,r2dfs)]
        
    


def diff(df1,df2):
    if not np.all(df1.columns == df2.columns):
        print("Don't compare apples with oranges")
        print(df1.columns)
        print(df2.columns)
        return None

    d  = df1.values - df2.values
    s  = np.hypot(df1.values,df2.values)
    rd = d/s
    d_df = pd.DataFrame(data=d,index = df1.index, columns = df1.columns)
    rd_df = pd.DataFrame(data=rd,index = df1.index, columns = df1.columns)
    return d_df, rd_df


def plot_stuff(nxy_df,ndata_df,col):
   
    xyz_df = nxy_df.join(ndata_df)[['x','y',col]]
    x = xyz_df['x'].values.flatten()
    y = xyz_df['y'].values.flatten()
    z = xyz_df[col].values.flatten()
    
    xi = np.linspace(x.min(),x.max(),800)
    yi = np.linspace(y.min(),y.max(),800)
 
    zi = griddata((x,y),z, (xi[None,:], yi[:,None]), method='nearest')
    plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
    plt.colorbar()
 
   
if __name__ == "__main__":
    try: 
        dirname = argv[1]
    except:
        dirname = '.'
    try: 
        res_frame = int(argv[2])
    except:
        res_frame = -1
    pd_df,nxy_df = geom_info(dirname)
    part_dfs = [ pd_df.loc[ pd_df['P'] == i ] for i in range(int(pd_df.values.min()),int(pd_df.values.max()+1))]
    parts = [ part_df.join(nxy_df)[['x','y']].values for part_df in part_dfs ]
     
    r1_dfs = parse_results(path.join(dirname,'RESULTS1.RES'),['ND','U','V'])
    r2_dfs = parse_results(path.join(dirname,'RESULTS2.RES'),['RHO','PS','TEMP'])
    
    all_data = [nxy_df.join(r1_df) for r1_df in r1_dfs ]
    all_data = [d.join(r2_df) for d,r2_df in zip(all_data,r2_dfs) ]
    
    x = all_data[0][['x']].values.flatten()
    y = all_data[0][['y']].values.flatten()
    
    xi = np.linspace(x.min(),x.max(),800)
    yi = np.linspace(y.min(),y.max(),800)
    
    for qname in ['ND','U','V','RHO','PS','TEMP']:
    #for qname in ['TEMP']:
        qs = [ data[[qname]].values for data in all_data ] 
        plt.figure()
        plt.title(qname)
        for part in parts:
            plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
        zi = griddata((x,y),qs[res_frame].flatten(), (xi[None,:], yi[:,None]), method='nearest')
        plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
        plt.colorbar()
    
    plt.show()
