# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 10:28:20 2019

@author: Pranav
"""
#SETUP
__author__="Pranav Sahasrabudhe"
__email__="pranavs@uoregon.edu"

print("\n---------------------------------------------------------")
print("Sahasrabudhe_deduper.py - A tool that identifies and removes PCR dedupliates")
print('Author:',__author__)
print("Email:",__email__,)
print("Github: https://github.com/pranavs22/Deduper")
print("---------------------------------------------------------\n\n")

#Imports and declarations
import re
import sys
import argparse
import time
start=time.clock()


#check Python version, exit if python <=2.7
if float((sys.version[:3]))<= 2.7:
    print('This script only works for Python version > 3.5.0')
    sys.exit(1)


#Functions
def get_args():
    ''' This function creates arguments to be passed to the program'''
    parser=argparse.ArgumentParser(description="Please provide path to your folder as well the desired text")
    parser.add_argument("-f","--FILE",help="File path of input file",required=True,type =str)
    parser.add_argument("-p","--PAIRED",help="If using paired-end library (note- this option is not supported currently)",required=False,type =str)
    parser.add_argument("-u","--UMI",help="List of files containing UMIs",required=True,type =str)
    parser.add_argument("-H","--HELP",help="Show this help message and quit",required=False,type =help)
    parser.add_argument("-o","--OUTFILE",help="Specify Output File name",required=True,type =str)
    
    return parser.parse_args()

def parse_sam(FILE,OUTFILE):
    '''
    This function takes in a SAM file identifies and removes dupplicate reads 
    and writes unique reads to a new user specified file.
    Duplicates are identifies based on -
    1) Chromosome position
    2) Strand
    3) UMI
    4) CIGAR string
    '''
    print("\nDeduplicating reads...")
    pattern=re.compile('([DMSNI])')
    unknown_umi=with_umi=total=fwd=rev=dup=uniq=0                #initialize all counter variables to zero
    read=chrom=''
    rec_set=set()

    with open (FILE) as fh,open(OUTFILE,"w")as o:
        while True:
            read=fh.readline().strip()
            if read=='':
                break
            else:
                if read.startswith("@"):                                #check headers
                    o.write(read + "\n")                                #write header to new file
                else:
                    info=read.split()
                    total+=1
                    umi=info[0].split(':')[-1]
                    if umi in UMIs:                  #check if UMI is present in known list of UMIs
                        with_umi+=1
                        chrom_h=info[2]
                        if chrom_h==chrom:                         #check chromosome
######FWD######################                        
                            split_cigar=''                         #initialize CIGAR
                            
                            pos=int(info[3])                       #initialize POS
                            cigar=info[5]                   
                            chrom_h=info[2]                        #initialize chrom
                            flag=int(info[1])                      #initialize flag
###STRANDEDNESS
                            if flag & 16 !=16:                          #check strand                   #Forward  Strand
                                fwd+=1                                  #increment counter for forward strand
                                split_cigar=pattern.split(cigar)
                                if split_cigar[1]=="S":                 # check if 'S' in the initial position in CIGAR
                                    pos=pos-int(split_cigar[0])         # adjust position
                                else:
                                      pass     
                                rec='_'.join([umi,chrom_h,str(pos),str(flag)]) #create a joint string of umi,chrom,pos and flag

                                if rec in rec_set:                       #if record is not seeen
                                    dup+=1                               #increment the duplicate counter
    
                                else:
                                    o.write('\t'.join(info)+"\n")        #if it is already seen
                                    uniq+=1                              #increment counter for unique reads   
                                    rec_set.add(rec)
                            else:                                  #Read is on reverse sraand
                                rev+=1                             #increment counter for reverse reads      
                                split_cigar=pattern.split(cigar)

                                for i in range(len(split_cigar)//2): #Iteraate through cigar string
                                    c=split_cigar[2*i+1]
                                    s=int(split_cigar[2*i])
#CIGAR operations
                                    if c=="S" and i!=1:
                                        pos+=s
                                    if c=="M":
                                        pos+=s
                                    if c=="D":
                                        pos+=s
                                    if c=="N":
                                        pos+=s
                                rec='_'.join([umi,chrom_h,str(pos),str(flag)])

                                if rec in rec_set:
                                    dup+=1                               #increment the duplicate counter
                                else:
                                     rec_set.add(rec)
                                     o.write('\t'.join(info)+"\n")
                                    
    
###########################
                        else:
                            chrom=chrom_h                                      #change value of chromosome
                            print("Chromosome {} finished".format(chrom))
                            
                    else:                                                 #UMMI not in the list of know UMIs
                        unknown_umi+=1
                        
    print('\nFinished Deduplicating reads!')
    print("\n---------------------------------------------------------")
    print("Summary")
    print("---------------------------------------------------------")
    print("Total \t\t",total)
    print("With umi \t",with_umi)
    print("without umi \t",unknown_umi)
    print("Forward \t",fwd)
    print("Reverse \t",rev)
    print("duplicate \t",dup)
    print("unique \t\t",uniq)

    return "\nRecords written successfully in '{}' file".format(OUTFILE)

#Get user arguments
args=get_args()
FILE=args.FILE
OUTFILE=args.OUTFILE + '_deduped.sam'
UMI=args.UMI
PAIRED=args.PAIRED

print("\nRunning command python Sahasrabudhe_deduper.py -f {} -o {} -u {}\n".format(FILE,OUTFILE,UMI))

#############Conditions before execution


#check if paired-end option is given, if yes print an error message and quit
if PAIRED:
    print( 'Sorry, the current version of script does not support paired end files!. Exiting Program ')
    sys.exit(1)

if UMI:
    UMIs=set()
    with open(UMI) as u:
        for i in u:
            UMIs.add(i.strip())
    print("UMIs added...")

else:
    print('This script requires a list of UMIs')
    sys.exit(1)

#Run the script
if __name__=="__main__":
    print(parse_sam(FILE,OUTFILE))
    print("\nFinished in {0:5f} seconds".format(time.clock()-start))