from __future__ import print_function, division
import pandas as pd
from os.path import dirname, basename
import sys, csv


def parseFile(doc):
	"""
	Input:  Document of gRNAs (for now, a .csv).
	Output: panda df with gRNA and names
	"""
	if not '.csv' in basename(doc):
		raise IOError(
		"Please convert your output to a .csv with:\
		\n grna name | grna sequence ")

	with open(doc, 'r') as f:
		reader = csv.reader(f)
		grnaList = [row for row in reader]

	if not set(grnaList[0][1]) <= set('ATGCatgc'): # File has a header
		grnaList.remove(grnaList[0])
	return grnaList # File has no header, starts with first gRNA

def reverseComp(sequence):
	"""
	Input:  string of DNA
	Output: string of reverse complement DNA
	"""
	sequence = sequence.upper()
	comp = {'A':"T","T":"A","C":"G","G":"C"}
	out = ''.join([comp[i] for i in sequence]) # Generator function to 
											   # replace all characters in s.
	return out[::-1] # Return inverse string

def addAdapter(grna, backbone):
	"""
	Function to add flanks to gRNAs depending on the backbone.
	Input: gRNA and backbone.
	Output: Two lines of a dataframe coorisponding to the forward and 
	reverse strand of the gRNA with flanks. 
	"""
	flanks = {
		'f_137':('TTGG','GTTTAAGAGC'),
		'r_137':('TTAGCTCTTAAAC','CCAACAAG'),
		'f_330':('CACCG', ''),
		'r_330':('AAAC','C')}

	grna.upper()
	if grna[0] == 'G' and len(grna) == 20: #Too many G's
		grna = grna[1:]

	if backbone in ('p1371','p1372'):
		line1 = flanks['f_137'][0] + grna + flanks['f_137'][1]
		line2 = flanks['r_137'][0] + reverseComp(grna) + flanks['r_137'][1]
	elif backbone == 'px330':
		line1 = flanks['f_330'][0] + grna + flanks['f_330'][1]
		line2 = flanks['r_330'][0] + reverseComp(grna) + flanks['r_330'][1]
	else:
		raise ValueError('Backbone argument must be either:\
		 \n p1371 \
		 \n p1372 \
		 \n px330.')
	return line1, line2

def getFilename(file):
	"""
	Use os.path to grab the filename for use in generating the output csv
	"""
	return dirname(file), basename(file).split('.')[0]


def saveCSV(output, file):
	"""
	Saves the dataframe
	"""
	path, filename = getFilename(file)

	if path == '':
		path = '.'

	save = path + '/' + filename + '_oligo_output.csv'

	with open(save, 'wb') as outfile:
	    csv_writer = csv.writer(outfile)
	    for row in output:
	        csv_writer.writerow(row)
	return

def generateOutput(dat, backbone):
	"""
	Main funciton for gRNAtoOligo.py
	Input: the file and the backbone argument
	Output: A .csv file with adapters for each gRNa in the list.
	"""
	output_csv = [['Name', 'Sequence']]

	for row in dat:
		oligos = addAdapter(row[1], backbone)
		name = row[0]
		output_csv.append([name + '_F', oligos[0]])
		output_csv.append([name + '_R', oligos[1]])

	return output_csv

def main():
	
	dat = parseFile(sys.argv[1])
	backbone = sys.argv[2]

	saveCSV(generateOutput(dat, backbone), sys.argv[1])


if __name__ == "__main__":

	main()

