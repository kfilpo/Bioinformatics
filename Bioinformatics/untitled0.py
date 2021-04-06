#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 21:20:41 2020

@author: katherinefilpolopez
"""
import re

def transponMethod(seq, transSeq):
    n = len(transSeq)
    n = int(n/2)
    trans1 = transSeq[:n]
    trans2 = transSeq[n+1:]
    

    tracker = 0
    for i in range(n-2):
        for j in range(n-2):
            if(trans1[i] == trans2[j]):
                tracker = tracker +1
            j= j+1
        i = i +1
        
        if tracker > 3:
            print("It's a transpon")
            regexObj = re.compile(transSeq)
            matchObj = regexObj.search(seq)
            
            transTable = str.maketrans("ATCG", "TAGC")
            comp = seq.translate(transTable)
            reverse = comp[::-1]
            
            regexObj2 = re.compile(transSeq)
            matchObj2 = regexObj2.search(reverse)
            
            if matchObj or matchObj2:
                print("The transpon has been inserted")
            else: 
                print("The transpon has not been inserted")   
    if tracker < 3:
        raise Exception("It's not a transpon")         
    return seq
    


seq1 = 'ATGCAGTTTTTTT'
seq2 = 'ATATGGGGGATAT'

seq3 = 'ATATGGGGG'
seq4 = 'ATGCAGTTTTTTTATATGGGGGATAT'


print(transponMethod(seq1,seq2))
print(transponMethod(seq4,seq2))
