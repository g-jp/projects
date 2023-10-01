import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import geneddit as ge
from collections import Counter

def readFASTA(filename):
    """This function reads a FASTA format file and
    returns a pair of strings
    with the header and the sequence
    """
    with open(filename) as a:
        lines = [line.strip() for line in a]
    lines = [line for line in lines if len(line) > 0]
    
    if lines[0].startswith('>'):
        return lines[0], ''.join(lines[1:])
    else:
        return '', ''.join(lines)

# A função 'open' lê cada linha como um elemento individual do ficheiro.
# É quase como se cada linha fosse um elemento da lista gerada com a função.
# A função retorna uma string de acordo com as operações mencionadas.

basesDNA = 'ATGC'
basesRNA = 'AUGC'

aa_residues   = "ACDEFGHIKLMNPQRSTVWY"

complementDNA = { 'A':'T', 'T':'A', 'G':'C', 'C':'G'}
complementRNA = { 'A':'U', 'T':'A', 'G':'C', 'C':'G'}

gencode = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
     'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y',
     'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L',
     'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P',
     'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
     'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I',
     'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'ACU': 'T', 'ACC': 'T',
     'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K',
     'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A',
     'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D',
     'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
     'GGG': 'G', 'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'}

trans123 = {'A': 'Ala', 'C': 'Cys', 'E': 'Glu', 'D': 'Asp', 'G': 'Gly',
            'F': 'Phe', 'I': 'Ile', 'H': 'His', 'K': 'Lys', 'M': 'Met',
            'L': 'Leu', 'N': 'Asn', 'Q': 'Gln', 'P': 'Pro', 'S': 'Ser',
            'R': 'Arg', 'T': 'Thr', 'W': 'Trp', 'V': 'Val', 'Y': 'Tyr'}

aa_masses = {'A': 71.0788, 
             'C': 103.1448, 
             'E': 129.1155, 
             'D': 115.0886, 
             'G': 57.0519, 
             'F': 147.1766, 
             'I': 113.1594, 
             'H': 137.1411, 
             'K': 128.1741, 
             'M': 131.1926, 
             'L': 113.1594, 
             'N': 114.1038, 
             'Q': 128.1307, 
             'P': 97.11667, 
             'S': 87.0782, 
             'R': 156.1875, 
             'T': 101.1051, 
             'W': 186.2132, 
             'V': 99.1326, 
             'Y': 163.176
             }

aa_classes = { 'small'       : 'PCAGVTDSN',
    'tiny'        : 'AGCS',
    'aliphatic'   : 'ILV',
    'aromatic'    : 'FYWH',
    'positive'    : 'KHR',
    'negative'    : 'DE',
    'charged'     : 'KHRDE',
    'hydrophobic' : 'CAGIVLTMHYWF',
    'polar'       : "CSTNDYWHKQRDE",
    'proline'     : 'P'}


def complement(seq):
    complementair = []
    for n in seq:
        complementair.append(complementDNA[n])
    comp = ''.join(complementair)
    return comp

def prot_noncod(seq): #From the noncoding chain it produces the associated protein
    comp = complement(seq)
    mRNA = []
    for n in comp:
        mRNA.append(complementRNA[n])
    mRNA = ''.join(mRNA)

    protein = []
    cod = [mRNA[i] + mRNA[i+1] + mRNA[i+2] for i in range(0, len(mRNA), 3)]
    for n in cod:
        protein.append(gencode[n])
    protein = '-'.join(protein)
    return protein

def prot_cod(seq): #From the noncoding chain it produces the associated protein
    mRNA = []
    for n in seq:
        mRNA.append(complementRNA[n])
    mRNA = ''.join(mRNA)

    protein = []
    cod = [mRNA[i] + mRNA[i+1] + mRNA[i+2] for i in range(0, len(mRNA), 3)]
    for n in cod:
        protein.append(gencode[n])
    protein = '-'.join(protein)
    return protein
    
def BLAST(seq, seq1): #For a sequence (seq), this function analyses if another, smaller sequence (seq1) is contained inside the first one
    complementair = [] #To find the complementair sequence
    for n in seq:
        complementair.append(complementDNA[n])
    comp = ''.join(complementair)
#________________________________________________
    anneal = str(seq1) 
#________________________________________________
    if anneal in seq:
        return print(f'{anneal} appears in the requested sequence')
    elif anneal in comp:
        return print(f'{anneal} appears in the sequence complementary to the requested one')
    else:
        return print(f'{anneal} does not appear in any sequence') 

def gene_size(gene, primerfw, primerrv):
    if primerfw in gene:
        a = gene.split(primerfw)[1]
        a = primerfw + a
        try:
            a = a.split(complement(primerrv))[0]
            a = a + complement(primerrv)
            size = len(a)
            return print(f'The requested amplicon has {size} bp')
        except:
            print('Primer reverse not found in gene')
    else: print('Primer foward not found in gene')

def amplicon(gene, primerfw, primerrv):
    if primerfw in gene:
        try:
            a = gene.split(primerfw)[1]
            a = primerfw + a
            try:
                a = a.split(complement(primerrv))[0]
                a = a + complement(primerrv)
                return a
            except: 
                print('Primer reverse not found in gene')
        except:
            print('Primer fowrd not found in gene')

def reverse(x):
    return x[::-1]

import math as m
from matplotlib import pyplot as plt

