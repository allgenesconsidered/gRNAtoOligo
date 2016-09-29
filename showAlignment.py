from __future__ import print_function, division
import gRNAtoOligo as grna
from csv import reader
from os.path import splitext


class oligoSequence:
	def __init__(self,sequence):
		self.sequence = sequence
	def __str__(self):
		return '<Oligo: %s >' %self.sequence

	def alignment(self, seq1, seq2):
		# Assumes perfect alignment.
		both = [seq1, seq2]
		both.sort(key = len)
		smaller, larger = both[0], both[1]

		revSmall = grna.reverseComp(smaller)
		for i in range(len(larger)):
			if larger[i:i+len(revSmall)] == revSmall:
				return i, larger
		return None, larger

	def showAlignment(self,comp):
		index, oligo1 = self.alignment(self.sequence, comp.sequence)

		if index is None:
			raise IOError(
				"Alignment unsuccessful. Please check input.")

		for seq in [self.sequence, comp.sequence]:
			if seq is not oligo1:
				oligo2 = seq
		output = "Successful alignment: \n" + \
			oligo1 + '\n' + \
			' ' * index + '|' * len(oligo2) + '\n' + \
			' ' * index + oligo2[::-1] + '\n'
		return output


def readCSV(file):
	if not checkIfCSV(file):
		raise IOError("Not a .csv")
	with open(file) as f:
		seqs = []
		read = reader(f, delimiter=',')
		for row in read:
			test = row[1]
			if checkIfDNA(test):
				seqs.append(test)
			if len(seqs) >= 2:
				break
		if len(seqs) == 0:
			raise IOError('No valid DNA sequences in .csv')
	return seqs

def checkIfCSV(file):
	return splitext(file)[1] == '.csv'

def checkIfDNA(s):
	if set(s) <= set('ATGCatgc'): 
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
