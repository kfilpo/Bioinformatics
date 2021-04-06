#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 20:19:07 2020

@author: katherinefilpolopez
"""

STANDARD_GENETIC_CODE = { 
          'UUU':'Phe', 'UUC':'Phe', 'UCU':'Ser', 'UCC':'Ser',
          'UAU':'Tyr', 'UAC':'Tyr', 'UGU':'Cys', 'UGC':'Cys',
          'UUA':'Leu', 'UCA':'Ser', 'UAA':'Stop','UGA':'Stop',
          'UUG':'Leu', 'UCG':'Ser', 'UAG':'Stop','UGG':'Trp',
          'CUU':'Leu', 'CUC':'Leu', 'CCU':'Pro', 'CCC':'Pro',
          'CAU':'His', 'CAC':'His', 'CGU':'Arg', 'CGC':'Arg',
          'CUA':'Leu', 'CUG':'Leu', 'CCA':'Pro', 'CCG':'Pro',
          'CAA':'Gln', 'CAG':'Gln', 'CGA':'Arg', 'CGG':'Arg',
          'AUU':'Ile', 'AUC':'Ile', 'ACU':'Thr', 'ACC':'Thr',
          'AAU':'Asn', 'AAC':'Asn', 'AGU':'Ser', 'AGC':'Ser',
          'AUA':'Ile', 'ACA':'Thr', 'AAA':'Lys', 'AGA':'Arg',
          'AUG':'Met', 'ACG':'Thr', 'AAG':'Lys', 'AGG':'Arg',
          'GUU':'Val', 'GUC':'Val', 'GCU':'Ala', 'GCC':'Ala',
          'GAU':'Asp', 'GAC':'Asp', 'GGU':'Gly', 'GGC':'Gly',
          'GUA':'Val', 'GUG':'Val', 'GCA':'Ala', 'GCG':'Ala', 
          'GAA':'Glu', 'GAG':'Glu', 'GGA':'Gly', 'GGG':'Gly'}
        
AMINO_ACID_CODE = {
          'Ala':'A', 'Cys':'C', 'Asp':'D', 'Glu':'E',
          'Phe':'F', 'Gly':'G', 'His':'H', 'Ile':'I',
          'Lys':'K', 'Leu':'L', 'Met':'M','Asn':'N',
          'Pro':'P', 'Gln':'Q', 'Arg':'R','Ser':'S',
          'Thr':'T', 'Val':'V', 'Trp':'W', 'Tyr':'Y',
          'Stop':'_', 
        }
 
class Sequence:
    def __init__(self, DNAfile):
        file = open(DNAfile)
        sequence = ''
        
        for line in file:
            sequence += line.strip()
            
        validation = ['A','C','T','G', 'N']
        
        for el in sequence:
           if el not in validation:
               raise Exception('Not a valid DNA sequence')
                   
        self.dna = sequence
        
 
    def translateDNA(self, rfP):
        rfs = []
        #starts from 0
        for i in range(1,8):
            if i>0 and i<4:
                seqRep = self.dna
                seqRep = self.transcribe()
                rfs.append(self.translateRNAtoAA(seqRep[i:]))
            else:
                seqRep = self.reverseTrans()
                rfs.append(self.translateRNAtoAA(seqRep[i:]))
       
        if rfP == 0:
            for j in range(1,7):
                print("Reading Frame", j, rfs[j])
                return rfs[j]
            
        elif rfP in range(1,7):
            print("Reading Frame", rfP, rfs[rfP])
            return rfs[rfP]
           
        else:
            raise Exception("RF not in 1-6")
        
    def predictRF(self):

        for i in range (1,7):
            readFrame = []
            seqIn = self.translateDNA(i)
            stop='1'

            for el in seqIn:
                tracker = []
                tracker.append(self.helper(tracker))
                for el in tracker:
                    if el == '_':
                        stop = tracker.index(el)
                        readFrame.append(tracker[:stop])
            if stop == None:
                print("No start and/or stop codon in reading frame ", i)
            else:
                readFrame.append(seqIn[:stop])       
    
        return readFrame
    
    def helper(self, seq):
    	helperReturn = []
    	for el in seq:
    		if el =='M':
    			helperReturn = helperReturn.append(seq[el:])
    			return seq.index('M') + helperReturn
    
    def longestRE(self):
        possibleFrames = []
        for i in range(len(self.dna)):
            possibleFrames.append(self.predictRF)
        for el in possibleFrames:
            largest = ''
            if len(largest) > 0 and max(el,key = len) > len(largest):
                largest = el
        return largest                   
    
    def reverseTrans(self):    
        convert = self.dna.replace("T", "U")
 
        transTable = str.maketrans("AUCG", "UAGC")
        comp = convert.translate(transTable)
        reverse = comp[::-1]
        
        return reverse
    
    def transcribe(self):
        convert = self.dna.replace("T", "U")
        return convert
    
    def translateRNAtoAA(self, sequence):
        proteinSeq = []
        i = 0
        while i +2 < len(sequence):
            codon = sequence[i:i+3]
            aminoAcid = STANDARD_GENETIC_CODE[codon]
            AAProtein = AMINO_ACID_CODE[aminoAcid]
            
            proteinSeq.append(AAProtein)        
            i += 3
        return proteinSeq
 
        
        
seq = Sequence("RFdata.txt")
#print(seq.translateDNA(0))
#seq.predictRF()

seq2 = Sequence("RFdata3.txt")
#print(seq2.transcribe())
#print(seq2.predictRF())
print(seq2.longestRE())


