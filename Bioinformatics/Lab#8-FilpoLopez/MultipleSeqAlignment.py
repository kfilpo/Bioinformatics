# -*- coding: utf-8 -*-
"""
Consensus Multiple Alignment
Created on Fri Feb  2 10:08:02 2018

@author: Kerri Norton
"""

DNA_2 = {'G': { 'G': 1, 'C':-3, 'A':-3, 'T':-3, 'N':0 },
         'C': { 'G':-3, 'C': 1, 'A':-3, 'T':-3, 'N':0 },
         'A': { 'G':-3, 'C':-3, 'A': 1, 'T':-3, 'N':0 },
         'T': { 'G':-3, 'C':-3, 'A':-3, 'T': 1, 'N':0 },
         'N': { 'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N':0 }}  

BLOSUM62 = {'A':{'A': 4,'R':-1,'N':-2,'D':-2,'C': 0,'Q':-1,'E':-1,'G': 0,'H':-2,'I':-1,
                 'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 0,'W':-3,'Y':-2,'V': 0,'X':0},
            'R':{'A':-1,'R': 5,'N': 0,'D':-2,'C':-3,'Q': 1,'E': 0,'G':-2,'H': 0,'I':-3,
                 'L':-2,'K': 2,'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3,'X':0},
            'N':{'A':-2,'R': 0,'N': 6,'D': 1,'C':-3,'Q': 0,'E': 0,'G': 0,'H': 1,'I':-3,
                 'L':-3,'K': 0,'M':-2,'F':-3,'P':-2,'S': 1,'T': 0,'W':-4,'Y':-2,'V':-3,'X':0},
            'D':{'A':-2,'R':-2,'N': 1,'D': 6,'C':-3,'Q': 0,'E': 2,'G':-1,'H':-1,'I':-3,
                 'L':-4,'K':-1,'M':-3,'F':-3,'P':-1,'S': 0,'T':-1,'W':-4,'Y':-3,'V':-3,'X':0},
            'C':{'A': 0,'R':-3,'N':-3,'D':-3,'C': 9,'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,
                 'L':-1,'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1,'X':0},
            'Q':{'A':-1,'R': 1,'N': 0,'D': 0,'C':-3,'Q': 5,'E': 2,'G':-2,'H': 0,'I':-3,
                 'L':-2,'K': 1,'M': 0,'F':-3,'P':-1,'S': 0,'T':-1,'W':-2,'Y':-1,'V':-2,'X':0},
            'E':{'A':-1,'R': 0,'N': 0,'D': 2,'C':-4,'Q': 2,'E': 5,'G':-2,'H': 0,'I':-3,
                 'L':-3,'K': 1,'M':-2,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
            'G':{'A': 0,'R':-2,'N': 0,'D':-1,'C':-3,'Q':-2,'E':-2,'G': 6,'H':-2,'I':-4,
                 'L':-4,'K':-2,'M':-3,'F':-3,'P':-2,'S': 0,'T':-2,'W':-2,'Y':-3,'V':-3,'X':0},
            'H':{'A':-2,'R': 0,'N': 1,'D':-1,'C':-3,'Q': 0,'E': 0,'G':-2,'H': 8,'I':-3,
                 'L':-3,'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y': 2,'V':-3,'X':0},
            'I':{'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I': 4,
                 'L': 2,'K':-3,'M': 1,'F': 0,'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V': 3,'X':0},
            'L':{'A':-1,'R':-2,'N':-3,'D':-4,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-3,'I': 2,
                 'L': 4,'K':-2,'M': 2,'F': 0,'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V': 1,'X':0},
            'K':{'A':-1,'R': 2,'N': 0,'D':-1,'C':-3,'Q': 1,'E': 1,'G':-2,'H':-1,'I':-3,
                 'L':-2,'K': 5,'M':-1,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
            'M':{'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q': 0,'E':-2,'G':-3,'H':-2,'I': 1,
                 'L': 2,'K':-1,'M': 5,'F': 0,'P':-2,'S':-1,'T':-1,'W':-1,'Y':-1,'V': 1,'X':0},
            'F':{'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I': 0,
                 'L': 0,'K':-3,'M': 0,'F': 6,'P':-4,'S':-2,'T':-2,'W': 1,'Y': 3,'V':-1,'X':0},
            'P':{'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,
                 'L':-3,'K':-1,'M':-2,'F':-4,'P': 7,'S':-1,'T':-1,'W':-4,'Y':-3,'V':-2,'X':0},
            'S':{'A': 1,'R':-1,'N': 1,'D': 0,'C':-1,'Q': 0,'E': 0,'G': 0,'H':-1,'I':-2,
                 'L':-2,'K': 0,'M':-1,'F':-2,'P':-1,'S': 4,'T': 1,'W':-3,'Y':-2,'V':-2,'X':0},
            'T':{'A': 0,'R':-1,'N': 0,'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,
                 'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 5,'W':-2,'Y':-2,'V': 0,'X':0},
            'W':{'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,
                 'L':-2,'K':-3,'M':-1,'F': 1,'P':-4,'S':-3,'T':-2,'W':11,'Y': 2,'V':-3,'X':0},
            'Y':{'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H': 2,'I':-1,
                 'L':-1,'K':-2,'M':-1,'F': 3,'P':-3,'S':-2,'T':-2,'W': 2,'Y': 7,'V':-1,'X':0},
            'V':{'A': 0,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I': 3,
                 'L': 1,'K':-2,'M': 1,'F':-1,'P':-2,'S':-2,'T': 0,'W':-3,'Y':-1,'V': 4,'X':0},
            'X':{'A': 0,'R': 0,'N': 0,'D': 0,'C': 0,'Q': 0,'E': 0,'G': 0,'H': 0,'I': 0,
                 'L': 0,'K': 0,'M': 0,'F': 0,'P': 0,'S': 0,'T': 0,'W': 0,'Y': 0,'V': 0,'X':0}}

#Assumes sequences are aligned

def calcSeqIdentity(seqA, seqB):
    #find the shortest length
    slen = min(len(seqA), len(seqB))
    #define score variable
    score = 0.0
    
    #loop through sequence
    for i in range(slen):
        #determine if the sequences are the same
        if(seqA[i] == seqB[i]):
            #increment score
            score += 1
    
    #calculate %        
    return 100.0*score/slen

if __name__ == '__main__':

  seq1 = 'ALIGNMENTS'
  seq2 = 'ALIGDVENTS'
  seq3 = 'ALIGDPVENTS'
  seq4 = 'ALIGN-MENTS'

 # print(calcSeqIdentity(seq1, seq2)) # 80.0%
 # print(calcSeqIdentity(seq1, seq3)) # 40.0%
 # print(calcSeqIdentity(seq4, seq3)) # 72.7%
  
def calcSeqSimilarity(seqA, seqB, seqMatrix):
   #find shortest sequence length
   slen = min(len(seqA), len(seqB))
   score = 0.0
   
   #loop through sequence and find each pos 
   for i in range(slen):
       val1 = seqA[i]
       val2 = seqB[i]
       #find score
       score += seqMatrix[val1][val2]
       
   return score
       
 # DNA example
#print(calcSeqSimilarity('AGCATCGCTCT', 'AGCATCGTTTT', DNA_2))
#print(calcSeqSimilarity('ALIGNMENT', 'AYIPNVENT', BLOSUM62))

seq1 = 'ALIGNMENTS'
seq2 = 'ALIGDVENTS'
seq3 = 'ALIGDPVENTS'
seq4 = 'ALIGNXMENTS'

#print(calcSeqSimilarity(seq1, seq2, BLOSUM62)) # 80.0%
#print(calcSeqSimilarity(seq1, seq3, BLOSUM62)) # 40.0%
#print(calcSeqSimilarity(seq4, seq3, BLOSUM62)) # 72.7%


  
def pairAlignScore(alignA, alignB, simMatrix, insert =8, extend = 4):
    #variable for score
    score = 0.0;
    
    #find the shortest sequence length
    n = min(len(alignA), len(alignB))
    
    #loop through sequence
    for i in range(n):
        #find each residue in the sequence
        res1 = alignA[i]
        res2 = alignB[i]
        
        #compare each residue
        #first check whether either is a -
        if '-' not in (res1, res2):
            #add to score
            score += simMatrix[res1][res2]
        #else if '-' not before either of them    
        elif(i > 0) and ('-' in (alignA[i-1], alignB[i-1]) ): 
            #subtract extend score 
            score -= extend  #pentaly for second insertion of dash
            
        else:
            #subtract insert score
            score -= insert  #penalty for first dash
            
            
    return score  
  
# Test
#print(pairAlignScore('GCAATC', 'GGAA-C', DNA_2))
#print(pairAlignScore('GGAATC', 'GCA--C', DNA_2)) # -12
#print(pairAlignScore('ALIGDPPVENTS', '--ALIGNMENTS', BLOSUM62)) # -3
          
def sequenceAlign(seqA, seqB, simMatrix, insert=8, extend=4):
    #
  #numI = len(seqA) + 1 
  numI = len(seqA) + 1
  numJ = len(seqB) + 1

  #
  scoreMatrix = [[0] * numJ for x in range(numI)]
  routeMatrix = [[0] * numJ for x in range(numI)]
  
  #
  for i in range(1, numI):
    routeMatrix[i][0] = 1
  #
  for j in range(1, numJ):
    routeMatrix[0][j] = 2
  
    #
  for i in range(1, numI):
    for j in range(1, numJ):
    
      penalty1 = insert
      penalty2 = insert
      
      #
      if routeMatrix[i-1][j] == 1:
        penalty1 = extend
        
      elif routeMatrix[i][j-1] == 2:
        penalty2 = extend
      
        #
      similarity = simMatrix[ seqA[i-1] ][ seqB[j-1] ]
      
      #
      paths = [scoreMatrix[i-1][j-1] + similarity, # 
               scoreMatrix[i-1][j] - penalty1, # 
               scoreMatrix[i][j-1] - penalty2] #                      
      
      #
      best = max(paths)
      route = paths.index(best)           
      
      #
      scoreMatrix[i][j] = best
      routeMatrix[i][j] = route
      
  #
  alignA = []
  alignB = []
  
  #
  i = numI-1
  j = numJ-1
  #
  score = scoreMatrix[i][j]

  #  
  while i > 0 or j > 0:
      #
    route = routeMatrix[i][j] 
    #

    #
    if route == 0: # 
        #
      alignA.append( seqA[i-1] )
      alignB.append( seqB[j-1] )
      i -= 1
      j -= 1

    elif route == 1: # 
      alignA.append( seqA[i-1] )
      alignB.append( '-' )
      i -= 1      

    elif route == 2: # 
      alignA.append( '-' )
      alignB.append( seqB[j-1] ) 
      j -= 1
  
  # 
  alignA.reverse()
  alignB.reverse()
  #
  alignA = ''.join(alignA)
  alignB = ''.join(alignB)

  return score, alignA, alignB 
  

#Test function
seqA = 'WFSEPAIST'
seqB = 'FSRPAGVIST'
seq1 = 'ATCG'
seq2 = 'ATG'

#score, alignA, alignB = sequenceAlign(seqA, seqB, BLOSUM62)
score, alignA, alignB = sequenceAlign(seq1, seq2, BLOSUM62)

#print(score)  # 22
#print(alignA) # WFSEPE--IST
#print(alignB) # -FSRPAVVIST

def consensus(alignment, threshold=0.25):

  #
  n = len(alignment[0])
  #
  nSeq = float(len(alignment))
  consensus = ''

  #
  for i in range(n):
    #
    counts = {}

    #
    for seq in alignment:
      #
      letter = seq[i]
      #
      if letter == '-':
        continue
    #
    #
      counts[letter] = counts.get(letter, 0) + 1
    
    #new list
    fractions = []
    # 
    for letter in counts:
        #
      frac = counts[letter]/nSeq
      #
      #
      fractions.append([frac, letter])
      
    #
    fractions.sort()
    #
    bestFraction, bestLetter = fractions[-1]
    
    #
    #
    if bestFraction <=  threshold:
      consensus += 'N'
    #
    else:
      consensus += bestLetter

  return consensus

#i
def consensusMultipleAlign(seqs, threshold, simMatrix):
  #
  n = len(seqs)
  #
  insert = 2;
  extend = 4;
  multipleAlign = []
  
  i = 0
  #
  #
  for j in range(i+1,n):
    #
    seqB = seqs[j]

    #
    if not multipleAlign:
      seqA = seqs[i]
      #
      #
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix,insert, extend)
      #
      multipleAlign.append(alignA)
      multipleAlign.append(alignB)
      
    
    else:
      # 
      seqA = consensus(multipleAlign, threshold)
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix,insert, extend)
      #
      gaps = []
      #
      for k, letter in enumerate(alignA):
        # 
        if letter == '-':
          #
          #
          gaps.append(k)
      
      #
      for k, seq in enumerate(multipleAlign):
        #
        for gap in gaps:
          # 
          seq = seq[:gap] + '-' + seq[gap:]
        
        # 
        multipleAlign[k] = seq
      
      #
      multipleAlign.append(alignB)
      
  
  #    
  for k, seq in enumerate(multipleAlign):
    print(k, seq)


# if __name__ == '__main__':

#   print('\nConsensus sequence')
  
#   alignment = ['SRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK',
#                'TVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK',
#                'HHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIK',
#                'HEFSELVIALIIFGVMAGVIGTILFISYGSRRLIK']

#   print(consensus(alignment))
# #  
# #  seqs = ['SRPAPVVLIILCVMAGVIGTILLISYGIRLLIK',
# #          'TVPAPVVIILIILCVMAGIIGTILLLIISYTIRRLIK',
# #          'HHFSEPEITLIIFGVMAGVIGTILLLIISYGIRLIK',
# #          'HFSELVIALIIFGVMAGVIGTILFISYGSRLIK']


# #  print('\nConsensus paired alignment')
# #
# #  consensusMultipleAlign(seqs, 0.25, BLOSUM62)
  
#   seqsn = ['ATTGGC', 'ATTCGC', 'ATTGAC', 'ATTGC']
#   consensusMultipleAlign(seqsn, 0.25, DNA_2)

#   seqsn = ['ATTGGCTTC', 'ATTCGCC', 'ATTGACTT', 'ATTGTTC']
#   consensusMultipleAlign(seqsn, 0.25, DNA_2)