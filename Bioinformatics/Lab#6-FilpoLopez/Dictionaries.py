#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 20:02:34 2020

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


DNA_1 = {'G': { 'G':1, 'C':0, 'A':0, 'T':0 },
         'C': { 'G':0, 'C':1, 'A':0, 'T':0 },
         'A': { 'G':0, 'C':0, 'A':1, 'T':0 },
         'T': { 'G':0, 'C':0, 'A':0, 'T':1 }}


