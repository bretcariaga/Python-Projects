#!/usr/bin/env python3
# Name: Bret Cariaga (bcariaga)
# Group Members: None

class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    nucleotide = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0}
    aa2mw = {
             'A':80.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C':121.158, 
             'H': 155.155, 'N': 132.118, 'T':119.119, 'D': 133.103, 'I': 131.173,
             'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q':146.145,
             'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y':181.189
             }
    
    def __init__ (self, inString=''):
        # takes the rnacodon dictionary and switches keys and values. Then it makes the values for the new dictionary 0
        self.aaComp = {values:keys for keys, values in self.rnaCodonTable.items()}
        self.aaComp = {x:0 for x in self.aaComp}
        
        # makes values of rnacodontable 0 for each key for a new dictionary
        self.codonComp = {key:0 for key in self.rnaCodonTable.keys}
        
        # just taking the dictionary from above
        self.nucComp = self.nucleotide
        
    def addSequence (self, sequence):
        
        sequence = ''.join(self.split()).upper()
        
        for nuc in sequence:
            if nuc in self.nucleotide.keys():  
                self.nucComp[nuc] += 1  
        
        rnaSeq = sequence.replace('T', 'U')  

        for num in range(0, len(rnaSeq), 3):
            codon = rnaSeq[num: num + 3]
            if codon in self.rnaCodonTable.keys():
                self.codonComp[codon] += 1  
                aa = self.rnaCodonTable[codon]
                if aa != '-':
                    self.aaComp[aa] += 1
                    
    def aaComposition(self):
        return self.aaComp
    
    def nucComposition(self):
        return self.nucComp
    
    def codonComposition(self):
        return self.codonComp
    
    def nucCount(self):
        totalNuc = 0
        for nucs in self.nucComp.values():
            totalNuc += nucs
        return sum(self.nucComp)

'''
 The FastAreader class was created by Professor David Bernick to read fasta files. It is very useful to
 iterate through each header and sequence of a fasta file.
'''
import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
