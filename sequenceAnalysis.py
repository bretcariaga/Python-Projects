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
ProteinParam is the same as last lab so no changes were made.
'''
#I used this numpy package to create a float range between 0 and 14
# I learned about this package from this website:
# https://pynative.com/python-range-for-float-numbers/
import numpy
class ProteinParam(str):

    #method that saves an attribute that is just the input string
    def __init__(self, protein):
        self.myProtein = protein
    
    # counts all single letter aa abbreviations in the input string
    def aaCount(self):
        upperSeq = self.upper()
        aaSum = 0
        aa2mw = {
             'A':80.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C':121.158, 
             'H': 155.155, 'N': 132.118, 'T':119.119, 'D': 133.103, 'I': 131.173,
             'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q':146.145,
             'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y':181.189
             }
        # whenever a key in the aa2mw dictionary is found in the input string, we add 1 amino acid
        for aa in upperSeq:
            if aa in aa2mw.keys():
                aaSum += 1
        return(aaSum)
    
    # creates a list in the pH range 0-14 with a step of .01
    # The numpy.arange function was taken from the internet in:
    # https://pynative.com/python-range-for-float-numbers/
    # it is also cited above where I imported numpy
    def pI(self):
        pHList = list(numpy.arange(0,14,.01))
        # iterates the numbers between 0 and 14 and finds the pH where the charge of the sequence is close to 0
        # calls the _charge_ method
        for pHvalues in pHList:
            if self._charge_(pHvalues) <= 0.01:
                return(pHvalues)
        
    # creates a dictionary that gives the counts of each amino acid in a string
    #the values are called in the main method that was provided so the dictionary is the only thing returned
    def aaComposition(self):
        upperSeq = self.upper()
        aaComp = {
             'A': upperSeq.count('A'), 'G': upperSeq.count('G'), 'M': upperSeq.count('M'), 'S': upperSeq.count('S'), 'C':upperSeq.count('C'), 
             'H': upperSeq.count('H'), 'N': upperSeq.count('N'), 'T':upperSeq.count('T'), 'D': upperSeq.count('D'), 'I': upperSeq.count('I'),
             'P': upperSeq.count('P'), 'V': upperSeq.count('V'), 'E': upperSeq.count('E'), 'K': upperSeq.count('K'), 'Q':upperSeq.count('Q'),
             'W': upperSeq.count('W'), 'F': upperSeq.count('F'), 'L': upperSeq.count('L'), 'R': upperSeq.count('R'), 'Y':upperSeq.count('Y')
             }
        return(aaComp)
        
    # the pH parameter is provided by the pI method
    # essentially follows the equation given to find charge at any given pH
    def _charge_(self, pH):
        upperSeq = self.upper()
        aa2chargePos={'K':10.5,'R':12.4,'H':6}
        aa2chargeNeg={'D':3.86,'E':4.25,'C':8.33,'Y':10}
        aaNterm=9.69
        aaCterm=2.34
        # initializes the sum of the positive and negative charge as
        # the P-terminus and N-terminus effects on the charge since there is always a P and N-terminus
        sumPosCharge = (10**aaNterm)/((10**aaNterm)+10**pH)
        sumNegCharge = (10**pH)/((10**aaCterm)+10**pH)
        # differentiates between positive and negative charged aa's
        # adds the charge by following the equation with different pH's
        for aa in upperSeq:
            if aa in aa2chargePos.keys():
                aaPosCharge = 10**(aa2chargePos.get(aa))/((10**aa2chargePos.get(aa))+(10**pH))
                sumPosCharge += aaPosCharge
            elif aa in aa2chargeNeg.keys():
                aaNegCharge = 10**pH/((10**aa2chargeNeg.get(aa))+10**pH)
                sumNegCharge += aaNegCharge
        return(sumPosCharge-sumNegCharge)
        
    #finds the amount of light absorbed at a wavelength of 280 nm
    # Has the optional parameter which allows the molar extinction to be
    # calculated under oxidized and reduced conditions
    def molarExtinction(self, Cystine = True):
        aa2abs280 = {'Y':1490, 'W':5500, 'C':125}
        # added a reduced dictionary which is the same as the oxidized without cysteine.
        redAA2abs280 = {'Y':1490, 'W':5500}
        upperSeq = self.upper()
        totalAbs = 0.0
        if Cystine: # when oxidized
            for aa in upperSeq:
                if aa in aa2abs280.keys():
                    aaAbs = aa2abs280.get(aa)
                    totalAbs += aaAbs
        elif not Cystine: #when reduced
            for aa in upperSeq:
                if aa in redAA2abs280.keys():
                    redAAabs = redAA2abs280.get(aa)
                    totalAbs += redAAabs
        return(totalAbs)
        
    # provided in the template
    # added an optional parameter that allows calculation under both oxidized and reduced conditions
    # It was not necessary to change the method with the way I set up the molarExtinction method
    #but I did so because the extra credit said to.
    def massExtinction(self,Cystine = True):
        myMW = self.molecularWeight()
        if Cystine:
            return self.molarExtinction() / myMW if myMW else 0.0
        elif not Cystine:
            return self.molarExtinction() / myMW if myMW else 0.0
        
    #finds the molecular weight of any aa string (ignores other characters that aren't amino acids)
    # follows the given equation but I could not get the same answer as the one provided in the lab PDF
    # When I calculated it by hand I got the same number my program is getting.
    def molecularWeight(self):
        upperSeq = self.upper()
        mwH2O = 18.015
        aa2mw = {
             'A':80.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C':121.158, 
             'H': 155.155, 'N': 132.118, 'T':119.119, 'D': 133.103, 'I': 131.173,
             'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q':146.145,
             'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y':181.189
             }
        totalMW = 0.0
        #iterates the whole input string adding the masses of each amino acid subtracted by the molecular weight of water from each and adding a single water in the end.
        waterWeight =  mwH2O * (self.aaCount()-1)
        for aa in upperSeq:
            if aa in aa2mw.keys():
                aaMW = aa2mw.get(aa)
                totalMW += (aaMW)
        
        return(totalMW-waterWeight)


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