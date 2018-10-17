#!/usr/bin/env python
'''
   Simple tool to compare numerical differences between ascii files.
'''


from sys import argv, exit
from math import fabs

def to_number(s):
    try:
        return float(s)
    except ValueError:
        pass


f1 = open(argv[1],'r')
f2 = open(argv[2],'r')

line1s,line2s = f1.readlines(), f2.readlines()

f1.close()
f2.close()

if len(line1s) != len(line2s):
    print "The files have different number of lines, %d vs %d" \
        % (len(line1s),len(line2s))
    exit(1)
else:
    print "Lines: %d" % len(line1s)

l2_difference = 0
l1_difference = 0
sup_difference = 0

l2r_difference = 0
l1r_difference = 0
supr_difference = 0



for line1, line2 in zip(line1s,line2s):
    f1s = line1.split()
    f2s = line2.split()

    f1sn = [ to_number(s) for s in f1s if to_number(s) is not None  ]
    f2sn = [ to_number(s) for s in f2s if to_number(s) is not None  ]
    
    if len(f1sn) != len(f2sn):
        print "The files lines have different number of numeric fields"
    
    if len(f1sn) != 0:
        for f1n, f2n in zip(f1sn,f2sn):
            l2_difference += (f1n - f2n)**2
            l1_difference += fabs(f1n - f2n)
            difference = fabs(f1n - f2n)
            if difference > sup_difference :
                sup_difference = difference
            denominator = max(fabs(f1n),fabs(f2n)) 
            if denominator > 0 :
                l2r_difference += ((f1n - f2n)/denominator)**2
                l1r_difference += fabs((f1n - f2n)/denominator)
                rdifference = fabs((f1n - f2n)/denominator)
                if rdifference > supr_difference :
                    supr_difference = rdifference
     
    
print "l2: ", l2_difference
print "l1: ", l1_difference
print "sup: ", sup_difference

print "relative l2: ", l2r_difference
print "relative l1: ", l1r_difference
print "relative sup: ", supr_difference





    

