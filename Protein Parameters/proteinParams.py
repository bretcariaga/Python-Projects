#!/usr/bin/env python3
# Name: Bret Cariaga (bcariaga)
#Group Members: None

'''
This program takes in an input amino acid single letter abbreviation string
with any letter, number, or special character, and will only look at the 
single letter amino acid abbreviations. It will then print out the parameters
of the amino acid sequence. It also takes into consideration the aa sequence
in oxidized or reduced conditions.

Resubmission: I found out that my dictionary value was wrong for one of the amino acids.

Input: VLSPADKT%3.NVKAAW5 (amino acid sequence with other characters in it)
Output:
Number of Amino Acids: 14
Molecular Weight: 1472.7
molar Extinction coefficient: 5500.00
mass Extinction coefficient: 3.73
Theoretical pI: 9.88
Amino acid composition:
A = 21.43%
C = 0.00%
D = 7.14%
E = 0.00%
F = 0.00%
G = 0.00%
H = 0.00%
I = 0.00%
K = 14.29%
L = 7.14%
M = 0.00%
N = 7.14%
P = 7.14%
Q = 0.00%
R = 0.00%
S = 7.14%
T = 7.14%
V = 14.29%
W = 7.14%
Y = 0.00%
'''
#I used this numpy package to create a float range between 0 and 14
# I learned about this package from this website:
# https://pynative.com/python-range-for-float-numbers/
import numpy
class ProteinParam(str):
    #method that saves an attribute that is just the input string
    def __init__(self, protein):
        self.myProtein = protein
    
    ''' counts all single letter aa abbreviations in the input string '''
    def aaCount(self):
        upperSeq = self.upper()
        aaSum = 0
        aa2mw = {
             'A':89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C':121.158, 
             'H': 155.155, 'N': 132.118, 'T':119.119, 'D': 133.103, 'I': 131.173,
             'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q':146.145,
             'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y':181.189
             }
        # whenever a key in the aa2mw dictionary is found in the input string, we add 1 amino acid
        for aa in upperSeq:
            if aa in aa2mw.keys():
                aaSum += 1
        return(aaSum)
    '''
     creates a list in the pH range 0-14 with a step of .01
     The numpy.arange function was taken from the internet in:
     https://pynative.com/python-range-for-float-numbers/
     it is also cited above where I imported numpy
    '''
    def pI(self):
        pHList = list(numpy.arange(0,14,.01))
        # iterates the numbers between 0 and 14 and finds the pH where the charge of the sequence is close to 0
        # calls the _charge_ method
        for pHvalues in pHList:
            if self._charge_(pHvalues) <= 0.01:
                return(pHvalues)
        
    ''' creates a dictionary that gives the counts of each amino acid in a string
     the values are called in the main method that was provided so the dictionary is the only thing returned
    '''
    def aaComposition(self):
        upperSeq = self.upper()
        aaComp = {
             'A': upperSeq.count('A'), 'G': upperSeq.count('G'), 'M': upperSeq.count('M'), 'S': upperSeq.count('S'), 'C':upperSeq.count('C'), 
             'H': upperSeq.count('H'), 'N': upperSeq.count('N'), 'T':upperSeq.count('T'), 'D': upperSeq.count('D'), 'I': upperSeq.count('I'),
             'P': upperSeq.count('P'), 'V': upperSeq.count('V'), 'E': upperSeq.count('E'), 'K': upperSeq.count('K'), 'Q':upperSeq.count('Q'),
             'W': upperSeq.count('W'), 'F': upperSeq.count('F'), 'L': upperSeq.count('L'), 'R': upperSeq.count('R'), 'Y':upperSeq.count('Y')
             }
        return(aaComp)
        
    ''' the pH parameter is provided by the pI method
     essentially follows the equation given to find charge at any given pH '''
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
        charge = sumPosCharge-sumNegCharge
        return(charge)
        
    '''
    finds the amount of light absorbed at a wavelength of 280 nm
    Has the optional parameter which allows the molar extinction to be
    calculated under oxidized and reduced conditions
    '''
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
        
    '''
    provided in the template
    added an optional parameter that allows calculation under both oxidized and reduced conditions
    It was not necessary to change the method with the way I set up the molarExtinction method
    but I did so because the extra credit said to.
    '''
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
             'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C':121.158, 
             'H': 155.155, 'N': 132.118, 'T':119.119, 'D': 133.103, 'I': 131.173,
             'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q':146.145,
             'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y':181.189
             }
        totalMW = 0.0
        for aa in upperSeq:
            if aa in aa2mw:
                totalMW += (aa2mw[aa]-mwH2O)
        return(totalMW+mwH2O)
        

# the following was not changed from the template
import sys
def main():
    inString = input('protein sequence?')
    while inString:
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print("Number of Amino Acids: {aaNum}".format(aaNum=myAAnumber))
        print("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction(True)))
        print("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction(False)))
        print("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print("Amino acid composition:")
        myAAcomposition=myParamMaker.aaComposition()
        keys=list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber==0:myAAnumber=1
        for key in keys :
            print("\t{} = {:.2%}".format(key,myAAcomposition[key]/myAAnumber))
        
        inString=input('protein sequence? ')

if __name__ == '__main__':
    main()
