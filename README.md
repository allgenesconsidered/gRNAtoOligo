# gRNA to Oligos

### Dependencies: 
* python 2.7._
* pandas (a python package)
See the Appendix for installingt these dependencies.

### Usage

This script can be used for generating oligo sequences for gRNAs. You can then use the oligo sequences for ordering through IDT. The functions are pretty basic and take only two arguments:
* A .csv with two columns, the name of the gRNA,  and the gRNA sequence. If the first row contains the column names, they can be ignored.
* The plasmid you are designing your insert for. Either 'p1371', 'p1372', or 'px330'.

### For those who are not use to working with the ternimal

The best usage of this script is using it in the terminal.
* Open your Terminal window (In Applications/Utilities on OSX).
* Navigate to the gRNAtoOligo folder:
```bash
$ cd ./PATH_TO_gRNAtoOligo/
```
if you just downloaded the zip file, you can probably run
```bash
$ cd ~/Downloads/gRNAtoOligo/ 
```

Lets see an example:
```bash
$ python gRNAtoOligo.py test.csv p1371
```
which will output the test_oligo_output.csv file. The first argument must be the path to the .csv file with your gRNAs. The
second argument is the plasmid. 

For further validation, you can run showAlignment.py, which will show you alignment in the command line. The function will return
the first pair of oligos in the lost and show you the alignment pattern for quick visualization.

```bash
$ python showAlignment.py test.csv
Successful alignment:
TTAGCTCTTAAACTTGTTGTCCTCCTGCTCGCCCAACAAG
   |||||||||||||||||||||||||||||||||
   CGAGAATTTGAACAACAGGAGGACGAGCGGGTT
```
