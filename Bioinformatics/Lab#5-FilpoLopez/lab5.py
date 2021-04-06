#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 20:01:23 2020

@author: katherinefilpolopez
"""
import random
import re
from MultipleSeqAlignment import DNA_2, sequenceAlign, calcSeqSimilarity, consensusMultipleAlign

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



class Mutation:
    def __init__(self, sequence):
        self.seq = sequence
#        nucleotides = ['A','T','G','C']
#        for i in range(0, length):
#          self.seq += random.choice(nucleotides)
#          
    def getSeq(self):
        return self.seq
    
    def getNucleotides(self):
        return self.nucleotides
      
    def pointMutation(self):
        sequenceList = list(self.getSeq())
        mutationSite = random.randint(0,len(sequenceList)-1)
        choices = set(['A','T','G','C']) - set(sequenceList[mutationSite])
        sequenceList[mutationSite] = random.choice(list(choices))
        return ''.join(sequenceList)
        
    def frameshiftInsertion(self):
        sequenceList = list(self.getSeq())
        mutationSite = random.randint(0,len(sequenceList)-1)
        randomChoice = random.choice(['A','T','G','C'])
        sequenceList.insert(mutationSite,randomChoice)  
        return ''.join(sequenceList)
          
    def frameshiftDeletion(self):
        sequenceList = list(self.getSeq())
        mutationSite = random.randint(0,len(sequenceList)-1)
        sequenceList.remove(sequenceList[mutationSite])
        return ''.join(sequenceList)

    def NPointMutation(self, n):
        sequenceList = list(self.getSeq())
        if n > len(sequenceList):   
            return "N is longer than sequence. Please pick another n!"
        while n > 0:
            mutationSite = random.randint(0,len(sequenceList)-1)
            choices = set(['A','T','G','C']) - set(sequenceList[mutationSite])
            sequenceList[mutationSite] = random.choice(list(choices))
            n = n-1
        return ''.join(sequenceList)
    
    def crossoverMutation(self, seq2):
        sequenceList1 = list(self.getSeq())
        sequenceList2 = list(seq2)
        placeholderlist = sequenceList1
        mutationSite = random.randint(0,len(sequenceList1)-1)
        sequenceList1 = sequenceList1[:mutationSite] + sequenceList2[mutationSite+1:]
        sequenceList2 = sequenceList2[:mutationSite] + placeholderlist[mutationSite+1:]
        return ''.join(sequenceList1), ''.join(sequenceList2)
    
    def motifAddition(self, motif):
        sequenceList = list(self.getSeq())
        ret = re.compile(motif)
        matchObj = ret.search(self.seq)
        ref = matchObj.start()
        motifList = list(motif)
        if ref:
            sequenceList = sequenceList[:ref] + motifList + sequenceList[ref:]
        return ''.join(sequenceList)
    
    def motifDeleter(self, motif):
        sequenceList = list(self.getSeq())
        ret = re.compile(motif)
        matchObj = ret.search(self.seq)
        ref = matchObj.start()
        if ref:
            sequenceList = sequenceList[:ref] + sequenceList[ref+4:]
        return ''.join(sequenceList)
    
    def transposon(self, motif, tseq):
        sequenceList = list(self.getSeq())
        tseqList = list(tseq)
        ret = re.compile(motif)
        matchObj = ret.search(self.seq)
        ref = matchObj.start()
        if ref:
            sequenceList = sequenceList[:ref] + tseqList + sequenceList[ref:]
        return ''.join(sequenceList)
        
    
class SequenceVariationDetecion:
    def __init__(self, originalSeq, mutatedSeq):
        self.originalSeq = originalSeq
        self.mutatedSeq = mutatedSeq
    
    def mutationDetection(self, originalSeq, seq2, tseq= []):
        if self.originalSeq== seq2:
            print("There is no mutation")
        elif len(self.originalSeq)==len(seq2):
            print ("There is a point or crossover mutation.")
            self.pointMutationDetection(self.originalSeq,seq2)
        elif len(self.originalSeq) < len(seq2):
            print ("There is a frame shift insertion.")
        elif len(self.originalSeq) > len(seq2):
            self.transposonDetection(seq2, tseq)
        else:
            print ("There is a frame shift deletion.")
        return


    def pointMutationDetection(self, originalSeq, seq2):
         if self.crossoverDetection(originalSeq, seq2):
             return
         i = 0
         baseRNA = self.originalSeq.replace('T', 'U')
         RNA2 = seq2.replace('T', 'U')
         while i +2 < len(originalSeq):
                codon1 = baseRNA[i:i+3]
                codon2 = RNA2[i:i+3]
                baseAA = STANDARD_GENETIC_CODE[codon1]
                seq2AA = STANDARD_GENETIC_CODE[codon2]
                if seq2AA == 'Stop' and baseAA != 'Stop':
                    print("There is a nonsense mutation.")
                    return
                if baseAA != seq2AA:
                    print("There is a missense mutation.")
                    return
                i += 3
         print("There is a silent mutation.") 
         return
     
       
    def transposonDetection(self, seq, tseq):
        ret = re.compile(tseq)
        matchObj = ret.search(self.seq)
        ref = matchObj.start()
        if ref:
            print("There is a transposon mutation")
        return
   
    def crossoverDetection(self, seq, seq2):
        seqIn = list(seq)
        seq2In = list(seq2)
        i = len(seqIn)
        threshold = 0
        while i > 0: 
            if seqIn[i] != seq2In[i]:
                i= i-1
                threshold+1
                if threshold > 1:
                      print("There is a crossover mutation")
                      return
         

seq1 = 'GCACGTATTGATTGGCCTGTACCTA'
test = Mutation(seq1)
motifseq = 'ACGT'

seq2 = test.NPointMutation(15)

seq3 = test.crossoverMutation(seq2)[0]
seq4 = test.crossoverMutation(seq2)[1]

seq5 = test.motifAddition(motifseq)
seq6 = test.motifDeleter(motifseq)

tseq = 'ACGTGGTTGCACGT'
seq7 = test.transposon(motifseq, tseq)


arraySeq =[seq1, seq2, seq3, seq4, seq5, seq6, seq7]
scoreArray = []
for el in arraySeq:
    scoreArray.append(calcSeqSimilarity(seq1, el, DNA_2))

print(scoreArray)
print(consensusMultipleAlign(arraySeq, 0.25, DNA_2))









