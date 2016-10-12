from __future__ import print_function, division
import gRNAtoOligo as grna
from csv import reader
from os.path import splitext


class oligoSequence:
	"""
	An oligo sequence. Contains methods to do the aligment and 
	show the aligment graphically on the CLI. 

	Note, the alignment algoithm assumes that the two sequences
	align perfectly (ie without any sequence breaks or mismatches).
	If you are looking for a different algorithm, don't use this one.
	"""
	def __init__(self,sequence):
		self.sequence = sequence.upper()
	def __str__(self):
		return '<Oligo: %s >' %self.sequence

	def alignment(self, seqSelf, seqComp):
		# Assumes perfect alignment.
		both = [seqSelf, seqComp]
		both.sort(key = len) # Sort based on length, shortest first
		smaller, larger = both[0], both[1] 

		revSmall = grna.reverseComp(smaller) # take reverse comp of smaller oligo, for string matching
		for i in range(len(larger)): # iterate over larger oligo
			if larger[i:i+len(revSmall)] == revSmall: # 
				return i, larger
		return None, larger

	def showAlignment(self,comp):

		seqForward, seqReverse = self.sequence, comp.sequence

		index, oligoLarger = self.alignment(seqForward, seqReverse)
		size = 0

		if index is None:
			raise IOError(
				"Alignment unsuccessful. Please check input.")

		if seqForward is not oligoLarger:
			size = len(seqForward)
			seqForward = (' ' * index) + seqForward
		else:
			size = len(seqReverse)
			seqReverse = (' ' * index) + seqReverse

		output = "Successful alignment: \n" + \
			seqForward + '\n' + \
			' ' * index + '|' * size + '\n' + \
			seqReverse[::-1] + '\n'
		return output


def readCSV(file):
	"""	Read in a csv and take the first 2 gRNA sequences.
	"""
	if not checkIfCSV(file):
		raise IOError("Not a .csv")
	with open(file) as f:
		seqs = []
		read = reader(f, delimiter=',')
		for row in read:
			test = row[1]
			if checkIfDNA(test):
				seqs.append(test)
			if len(seqs) == 2:
				break
		if len(seqs) == 0:
			raise IOError('No valid DNA sequences in .csv')
	return seqs

def checkIfCSV(file):
	"""  Check if the file is a valid csv file.
	"""
	return splitext(file)[1] == '.csv'

def checkIfDNA(string):
	""" Returns true if the string contains elements 
	"""
	if set(string) <= set('ATGCatgc'): 
		return True
	return False


		

if __name__ == "__main__":
	import sys
	try:
		grnas = readCSV(sys.argv[1])
		seq1, seq2 = oligoSequence(grnas[0]), oligoSequence(grnas[1])
	except IOError:
		seq1, seq2 = oligoSequence(sys.argv[1]), oligoSequence(sys.argv[2])
		
	print(seq1.showAlignment(seq2))
