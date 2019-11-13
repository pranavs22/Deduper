# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 17:10:22 2019

@author: Pranav
"""
umi_file="C:/Bi624/deduper/Deduper/STL96.txt"
umis={}
with open(umi_file) as u:
    while True:
        
        randomers=u.readline().strip()
        if len(randomers)==0:
            break
        else:
            umis[randomers]=[]
    print(umis)
        
sam="C:/Bi624/sorted_Dataset.sam"
def parse_sam(sam):
    x=0
    with open (sam) as fh:
        read=''
        while True:
            read=fh.readline().strip()
            if len(read)==0:
                break
            if not read.startswith("@"):
                info=read.split()
                if int(info[2])==2:
                    x+=1
                    
    return x
parse_sam(sam)
