from __future__ import print_function, division
import pandas as pd
from os.path import dirname, basename
import sys, csv, re, argparse


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
	comp = {'A':"T","T":"A","C":"G","G":"C"}
	out = ''.join([comp[i] for i in sequence]) # Generator function to 
											   # replace all characters in s.
	return out[::-1] # Return inverse string

def formatString(sequence):
	"""
	Format string to remove spaces and have all cases be upper.
	"""
	return re.sub(" ", "", sequence).upper()

def addAdapter(grna, backbone, delete_g):
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
		'r_330':('AAAC','C'),
		'f_1010.1':('CTTG',''),
		'r_1010.1':('','CAAA'),
		'f_1010.2':('TTGG',''),
		'r_1010.2':('','CAAA')}

	if delete_g and grna[0] == 'G' and len(grna) == 20: #Too many G's
		grna = grna[1:]

	if backbone in ('p1371','p1372'):
		line1 = flanks['f_137'][0] + grna + flanks['f_137'][1]
		line2 = flanks['r_137'][0] + reverseComp(grna) + flanks['r_137'][1]
	elif backbone == 'px330':
		line1 = flanks['f_330'][0] + grna + flanks['f_330'][1]
		line2 = flanks['r_330'][0] + reverseComp(grna) + flanks['r_330'][1]
	elif backbone == 'p1010.1':
		line1 = flanks['f_1010.1'][0] + grna + flanks['f_1010.1'][1]
		line2 = flanks['r_1010.1'][0] + reverseComp(grna) + flanks['r_1010.1'][1]
	elif backbone == 'p1010.2':
		line1 = flanks['f_1010.2'][0] + grna + flanks['f_1010.2'][1]
		line2 = flanks['r_1010.2'][0] + reverseComp(grna) + flanks['r_1010.2'][1]
	else:
		raise ValueError('Backbone argument must be either:\
		 \n p1371 \
		 \n p1372 \
		 \n px330 \
		 \n p1010.1\
		 \n px1010.2')
	return line1, line2

def getFilename(file):
	"""
	Use os.path to grab the filename for use in generating the output csv
	Output: The path to the original file and the file name.
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

def generateOutput(dat, backbone, delete_g):
	"""
	Main funciton for gRNAtoOligo.py
	Input: the file and the backbone argument
	Output: A .csv file with adapters for each gRNa in the list.
	"""
	output_csv = [['Name', 'Sequence']]

	for row in dat:
		oligos = addAdapter(formatString(row[1]), backbone, delete_g)
		name = row[0]
		output_csv.append([name + '_F', oligos[0]])
		output_csv.append([name + '_R', oligos[1]])

	return output_csv

def main():

	parser = argparse.ArgumentParser(description='Autogenerate gRNA oligos')
	parser.add_argument('backbone', help='Backbone: either p1371, p1372, px330, or p1010(.1,.2).')
	parser.add_argument('guides', help='csv file of name : guide combos')
	parser.add_argument('-d','--delete_g', help='Delete front G, default False.', default=False)
	args = parser.parse_args()
	
	dat = parseFile(args.guides)

	saveCSV(generateOutput(dat, args.backbone, args.delete_g), args.guides)


if __name__ == "__main__":

	main()

