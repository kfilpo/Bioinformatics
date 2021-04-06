# -*- coding: utf-8 -*-
"""
Lab Report #1 DNA and RNA 
by Katherine Filpo 
September 13th 2020

"""

test = "ACAGGCATTTTgG"
testRNA = "GCAUCUUgG"

#Tests if the sequence given is DNA(returns 1) or RNA (returns 0)
def isDNA(sequence):
    for el in sequence:
        if el =='U':
            return 0
    else: return 1

#Combines two sequences
def ligase(sequence1,sequence2):
    return (sequence1+sequence2)

#Cuts part of a sequence
def nuclease(sequence, index):
    return sequence[:index]

#Makes seuences from DNA to RNA
def transcription(sequence):
    return sequence.replace("T", "U")

#Replicates the DNA by giving the resverse complement
def replication(sequence):
    table = str.maketrans("ATGC", "TACG")
    sequence = sequence.translate(table)
    return sequence
    #this is commented out because while it is technically right, 
    #I prefer to see the sequence printed the other way, so it simulates 
    #the antiparallel strand
    #reverseSeq = sequence[::-1]
    #return reverseSeq

#Replicates both DNA and RNA by giving the reverse complement
def reverseTranscription(sequence):
    if (isDNA(sequence) == 0):
        if (sequence.isupper() != 'True'):
            return sequence.upper() + "\n" + sequence.upper().replace("U", "T")
        else: return sequence.replace("U", "T")
    elif (isDNA(sequence) == 1):
        return sequence.upper() + "\n" + replication(sequence.upper())
   

print(reverseTranscription(testRNA))
#print(testRNA.upper())
#print(ligase(test,test))
#print(nuclease(test, 4))
#print(transcription(test))
#print(replication(test))
