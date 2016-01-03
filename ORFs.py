# Vicky Torrance
#
# Library: VickyModules
# Module: ORFs
#
# Description: Module for ORF-related scripts

import re

def find_ATG(sequence):
    ATG = [(m.start()) for m in re.finditer('ATG', sequence)]
    return ATG
    
def find_STOP(sequence):
    STOP=[]
    for m in re.finditer('TGA', sequence):
        STOP.append(m.start())
    for m in re.finditer('TAG', sequence):
        STOP.append(m.start())   
    for m in re.finditer('TAA', sequence):
        STOP.append(m.start())
    STOP.sort(key=int)    
    return STOP             

def find_orf(seq):
    starts=find_ATG(seq)
    stops=find_STOP(seq)
    orfs=[]
    for x in starts:
        for y in stops:
            if y>x and (y-x)%3==0 and x>=4:
                orfs.append(x)
                break
    return orfs 
