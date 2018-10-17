#!/usr/bin/python3
'''
This script takes a file name as argument. The file must have 4 columns,
in this order:
0. no-of-position-space-partitions 
1. no-of-velocity-space-partitions 
2. no-of-nodes 
3. time-per-iteration
and does the math and many plots to show how good scaling is.

If the script is given a second file with the same format it will plot 
the data from the two files together (only in some plots).

'''
import numpy as np
from sys import argv
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
times1 = np.loadtxt(argv[1])
comparison = False
try:
    times2 = np.loadtxt(argv[2])
    comparison = True
except:
    print("Second file not provided or error in reading it.")
    print("Not printing comparisons")



#speedup vs total number of working ranks

# grouping per p-space partitions
plt.figure()
t_per_p= [ (times1[times1[:,0]==n_p_part ][:,[1,3]],n_p_part) for n_p_part in set(times1[:,0]) ]
if comparison:
    t_per_p2= [ (times2[times2[:,0]==n_p_part ][:,[1,3]],n_p_part) for n_p_part in set(times2[:,0]) ]
for arr,n_p_part in t_per_p:
    plt.plot(arr[:,0],arr[0,1]/arr[:,1]/arr[:,0],label = "Pp="+str(n_p_part))
      
plt.xlabel("V-partitions")
plt.ylim([0,2])
plt.ylabel("Deviation from exp.scaling (wrt 1 rank)")
plt.legend()

plt.figure()
for arr,n_p_part in t_per_p:
    nranks = arr[:,0] * n_p_part
    speedup = times1[0,3]/arr[:,1]
    plt.plot(nranks,speedup,label = "Pp="+str(n_p_part), linestyle = "None", marker = '+')
if comparison:
    for arr,n_p_part in t_per_p2:
       nranks = arr[:,0] * n_p_part
       speedup = times1[0,3]/arr[:,1]
       plt.plot(nranks,speedup, linestyle = "None", marker = '.')
     
plt.xlabel("total no of working ranks")
plt.ylabel("speedup (wrt 1 rank)")
plt.legend()


plt.figure()
for arr,n_p_part in t_per_p:
    nranks = arr[:,0] * n_p_part
    time_needed = arr[:,1]
    plt.plot(1.0/nranks,time_needed,label = "Pp="+str(n_p_part), linestyle = "None", marker = '+')
if comparison:
   for arr,n_p_part in t_per_p2:
       nranks = arr[:,0] * n_p_part
       time_needed = arr[:,1]
       plt.plot(1.0/nranks,time_needed, linestyle = "None", marker = '.')
      
plt.xlabel("inverse of no of working ranks")
plt.ylabel("time per step")
plt.legend()



# grouping per v-space partitions
plt.figure()
t_per_v= [ (times1[times1[:,1]==n_v_part ][:,[0,3]],n_v_part) for n_v_part in set(times1[:,1])]
if comparison:
    t_per_v2= [ (times2[times2[:,1]==n_v_part ][:,[0,3]],n_v_part) for n_v_part in set(times2[:,1])]

for arr,n_v_part in t_per_v:
    plt.plot(arr[:,0],arr[0,1]/arr[:,1]/arr[:,0],label = "Pv="+str(n_v_part))
    
plt.xlabel("P-partitions")
plt.ylabel("Deviation from exp.scaling (wrt 1 rank)")
plt.ylim([0,2])
plt.legend()

plt.figure()
for arr,n_v_part in t_per_v:
    nranks = arr[:,0] * n_v_part
    speedup = times1[0,3]/arr[:,1]
    plt.plot(nranks,speedup,label = "Pv="+str(n_v_part), linestyle = "None", marker = '+')
if comparison:
    for arr,n_v_part in t_per_v2:
        nranks = arr[:,0] * n_v_part
        speedup = times1[0,3]/arr[:,1]
        plt.plot(nranks,speedup, linestyle = "None", marker = '.')


plt.xlabel("no of working ranks")
plt.ylabel("speedup (wrt 1 rank)")
plt.legend()

plt.figure()
for arr,n_v_part in t_per_v:
    nranks = arr[:,0] * n_v_part
    time_needed = arr[:,1]
    plt.plot(1.0/nranks,time_needed,label = "Pv="+str(n_v_part), linestyle = "None", marker = '+')
if comparison:
    for arr,n_v_part in t_per_v2:
        nranks = arr[:,0] * n_v_part
        time_needed = arr[:,1]
        plt.plot(1.0/nranks,time_needed,linestyle = "None", marker = '.')

plt.xlabel("inverse of no of working ranks")
plt.ylabel("time")
plt.legend()



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
z = times1[:,0]*times1[:,1]*times1[:,3]
ax.set_xlabel('Pp')
ax.set_ylabel('Pv')
ax.set_zlabel('Deviation from exp.scaling (wrt 1 rank)')
ax.plot_trisurf(times1[:,0], times1[:,1],z[0]/z)

plt.show()
