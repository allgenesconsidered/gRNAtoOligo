from __future__ import print_function, division
import pandas as pd
from os.path import splitext




def parseFile(doc):
	"""
	Input:  Document of gRNAs
	Output: panda df with gRNA and names
	"""
	if not '.csv' in doc:
		raise IOError(
		"Please convert your output to a .csv with:\
		\n grna name | grna sequence ")

	with open(doc, 'r') as f:
		first_line = f.readline()
	if not set(first_line.split(',')[1]) <= set('ATGCatgc'): # File has a header
		return pd.read_csv(doc, header = 1)
	return pd.read_csv(doc, header = None) # File has no header, starts with first gRNA

def reverseComp(s):
	"""
	Input:  string of DNA
	Output: string of reverse complement DNA
	"""
	s = s.upper()
	comp = {'A':"T","T":"A","C":"G","G":"C"}
	out = ''.join([comp[i] for i in s])
	return out[::-1]

def addAdapter(grna, backbone):
	"""
	Function to add flanks to gRNAs depending on the backbone.
	Input: gRNA and backbone.
	Output: Two lines of a dataframe coorisponding to the forward and 
	reverse strand of the gRNA with flanks. 
	"""
	flanks = {
		'f_137':('TTG','GTTTAAGAGC'),
		'r_137':('TTAGCTCTTAAAC','CAACAAG'),
		'f_330':('CACCG', ''),
		'r_330':('AAAC','C')}

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
	name = splitext(file)[0]
	return name.split('/')[-1]


def saveCSV(output, file):
	"""
	Saves the dataframe
	"""

	save = './' + getFilename(file) + '_oligo_output.csv'
	output.to_csv(save, index = False)

	return

def generateOutput(dat, backbone):
	"""
	Main funciton for gRNAtoOligo.py
	Input: the file and the backbone argument
	Output: A .csv file with adapters for each gRNa in the list.
	"""
	output_df = pd.DataFrame(columns=('Name','Sequence'))
	for row in dat.iterrows():
		oligos = addAdapter(row[1][1], backbone)
		(row[1][0])
		new_row = pd.DataFrame([
			[row[1][0] + '_F', oligos[0]],
			[row[1][0] + '_R', oligos[1]]],
			columns=('Name','Sequence'))
		output_df = output_df.append(new_row)
	return output_df


if __name__ == "__main__":
	import sys
	dat = parseFile(sys.argv[1])
	backbone = sys.argv[2]

	saveCSV(generateOutput(dat, backbone), sys.argv[1])

