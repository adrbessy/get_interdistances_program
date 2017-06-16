### ----------------------------
### TFBS Interdistances program
### ----------------------------

'''
This program allows to calculate interdistances between transcription factor binding sites.
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks), and unbound sequences).
You can display a plot for several thresholds.
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import re
import time
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
from interdistances_functions import *
from plot_functions import *
from operator import truediv
import operator
from collections import Counter
import matplotlib.patches as mpatches
from matplotlib import pylab
import types
import argparse
import logging
from optparse import OptionParser
from scipy import stats

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type = str, default= "ARF2")
parser.add_argument("--negative_sets", "-neg", nargs='*')
parser.add_argument("--pseudoCount", "-pc",type = float, default = 0.001)
parser.add_argument("--threshold", "-th",nargs='+',type = int, default= [-8, -9, -10])
parser.add_argument("--Interdistance_maxValue", "-maxInter", type = int, default= 20)
parser.add_argument("--sequence_number", "-sequence_number", type = int, default= 11654)
parser.add_argument("--histo", "-histo", type = bool, default= False)
parser.add_argument("--points", "-points", type = bool, default= True)
parser.add_argument("--curve", "-curve", type = bool, default= False)
parser.add_argument("--sum_threshold", "-sum_threshold", type = bool, default= False)

args = parser.parse_args()

#python get_interdistances.py -fac "ARF2" -pc 0.001 -maxInter 30 -th -6 -7 -8 -9 -10 -11 -points True -neg ../../sequences/ARF2_neg1.fas ../../sequences/ARF2_neg2.fas ../../sequences/ARF2_neg3.fas ../../sequences/ARF2_neg4.fas
#python get_interdistances.py -fac "ARF5" -pc 0.001 -maxInter 20 -th -10 -sequence_number 26659 -neg ../sequences/ARF5_neg1.fas ../sequences/ARF5_neg2.fas ../sequences/ARF5_neg3.fas ../sequences/ARF5_neg4.fas -points True
#python get_interdistances.py -fac "ARF5" -pc 0.001 -maxInter 20 -sequence_number 26659 -th -9 -10 -11 -12 -13 -neg ../sequences/ARF5_neg1.fas ../sequences/ARF5_neg2.fas -points True
#python get_interdistances.py -fac "ARF2" -pc 0.001 -maxInter 30 -sequence_number 26659 -th -6 -7 -8 -9 -10 -11 -points True -neg ../../Adrien_Arnaud_Raquel_Francois/Arnaud/fasta/ARF2_neg_1.fas ../../Adrien_Arnaud_Raquel_Francois/Arnaud/fasta/ARF2_neg_2.fas ../../Adrien_Arnaud_Raquel_Francois/Arnaud/fasta/ARF2_neg_3.fas ../../Adrien_Arnaud_Raquel_Francois/Arnaud/fasta/ARF2_neg_4.fas
factorTranscription = args.factor
negative_sets = args.negative_sets
pseudoCount = args.pseudoCount
threshold = args.threshold
Interdistance_maxValue = args.Interdistance_maxValue 
sequence_number = args.sequence_number
histo = args.histo
points = args.points
sum_threshold = args.sum_threshold
    
if histo == True and len(threshold) > 1 : 
	print("Impossible to display an histogram with several thresholds, you have to change a parameter.")
	sys.exit(0)
	
###################Parameters we can change#################

if factorTranscription == "ARF2" :
	FastaFile = "../../Adrien_Arnaud_Raquel_Francois/Arnaud/fasta/ARF2.fas" 
	MatrixFile = "../../matrices/ARF2_OMalley_matrixC.txt" 
	matrixType = "freq" 
	dependencyFile = ""
	
if factorTranscription == "ARF5" :
	FastaFile = "../../sequences/ARF5_bound_sequences.fas" 
	MatrixFile = "../../matrices/ARF5_OMalley_matrixC.txt" 
	matrixType = "freq" 
	dependencyFile = ""
	
if factorTranscription == "LFY_matrix_19nucl" :
	FastaFile = "../sequences/LFY_bound_sequences.fas"
	MatrixFile = "../matrix/LFY_scores_matrix_19nucl.txt" 
	dependencyFile = "../matrix/interdependent_bases_matrix_for_LFY.txt"
	matrixType = "score" 

########################################### About the main matrix #######################

''' The sens of the matrix is important: The positions are on the vertical sens and the bases are on the horizontal sens as described in the example.
separation between numbers can be spaces, tabulation, comas...

                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
                                                                        '''
####################################################################################

################### To capture file names where there are unbound sequences ###################

#FastaFileNnumber = input("\nHow many fasta files with unbound sequences ? ")
#d = {}
#for i in range (1,FastaFileNnumber+1) :
	#d["FastaFileN{0}".format(i)] = raw_input("\nName of fasta file with unbound sequences: ")

###############################################################################################

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read().replace("\r","\n") + "\n"
F.close()

# These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list

num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)
#print("list(reversed(Mdata)) : ",list(reversed(Mdata)))
matScore, lenMotif = get_score_matrix(Mdata,matrixType,pseudoCount)

# The following line allows to produce the reversed matrix
'''if we take the example given before : A T G C
			Position 1:      0.4444  0.155  0.654   0.645
			Position 2:      0.1645  0.1565 0.21614 0.16456
Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
So we can calculate with this reverse matrix, the score of the complementary strand.
'''
matRev = list(reversed(matScore))

	
########## get INTERDISTANCE VALUES for POSITIVE sets:

InterDR, InterER, InterIR = get_interdist(matScore,matRev,FastaFile,threshold,factorTranscription,Interdistance_maxValue,sum_threshold,lenMotif,dependencyFile,sequence_number)

##### Create empty lists to store interdistances occurences for the negative set:

InterDR_N = []
InterER_N = []
InterIR_N = []
lenThr = 0
listThr = []
for a in threshold :
	InterDR_N.append( [0] * (Interdistance_maxValue + 1) )
	InterER_N.append( [0] * (Interdistance_maxValue + 1) )
	InterIR_N.append( [0] * (Interdistance_maxValue + 1) )	
	listThr.append(lenThr)
	lenThr = lenThr + 1

########## get INTERDISTANCE occurences for NEGATIVE sets				
if negative_sets :
	for fastafileN in negative_sets :
		InterDR_N_temp, InterER_N_temp, InterIR_N_temp = get_interdist(matScore,matRev,fastafileN,threshold,factorTranscription,Interdistance_maxValue,sum_threshold,lenMotif,dependencyFile,sequence_number)
		# addition of the occurences of every negative sets
		for a,b,c,d in zip(InterDR_N_temp,InterER_N_temp,InterIR_N_temp,listThr) :
			InterDR_N[d] = [x + y for x, y in zip(InterDR_N[d], a)]
			InterER_N[d] = [x + y for x, y in zip(InterER_N[d], b)]
			InterIR_N[d] = [x + y for x, y in zip(InterIR_N[d], c)]
	# divide by the number of negative sets in order to have a mean
	if len(negative_sets) > 0 :
		for a,b,c,d in zip(InterDR_N,InterER_N,InterIR_N,listThr) :
			InterDR_N[d] = [x / float(len(negative_sets)) for x in a]
			InterER_N[d] = [x / float(len(negative_sets)) for x in b]
			InterIR_N[d] = [x / float(len(negative_sets)) for x in c]
		
interdist_sum = []
interdist_sum_N = []
rate = []
if negative_sets :
	for a,b,c,d,e,f in zip(InterDR,InterER,InterIR,InterDR_N,InterER_N,InterIR_N) :
		interdist_sum.append(sum(a) + sum(b) + sum(c))
		interdist_sum_N.append(sum(d) + sum(e) + sum(f))
		rate.append([divide(sum(a), sum(d)),divide(sum(b) , sum(e)),divide(sum(c) , sum(f))])
	
relative_DR = []
relative_ER = []
relative_IR = []
relative_DR_neg = []
relative_ER_neg = []
relative_IR_neg = []

if negative_sets :
	for a,b,c,d,e in zip(threshold,InterDR,interdist_sum,InterDR_N,interdist_sum_N) :
		relative_DR.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_DR_neg.append( [divide(x , float(e)) for x in d] )
	for a,b,c,d,e in zip(threshold,InterER,interdist_sum,InterER_N,interdist_sum_N) :
		relative_ER.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_ER_neg.append( [divide(x , float(e)) for x in d] )
	for a,b,c,d,e in zip(threshold,InterIR,interdist_sum,InterIR_N,interdist_sum_N) :
		relative_IR.append( [divide(x , float(c)) for x in b] )
		if len(negative_sets) > 0 :
			relative_IR_neg.append( [divide(x , float(e)) for x in d] )
#print("relative_DR : ",relative_DR)
#print("relative_ER : ",relative_ER)
#print("relative_IR : ",relative_IR)
command = sys.argv			
			
if histo == True and negative_sets :
	negative_sets_histo(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,factorTranscription,threshold,rate,command)
	
if histo == False and points == False and negative_sets :
	negative_sets_curve(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,factorTranscription,threshold,rate,command)
	
if points == True and negative_sets :
	negative_sets_points(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,factorTranscription,threshold,rate,command)