Program function:
The program takes an LRG file, and using the exon coordinates outputs the sequence for each exon for each transcript. 
The script can be modified in future versions to be able to specify one transcript and to add intronic sequences.
The output is in fasta format with the transcript and Exon name in the header and sequence in second line .

To run this file you will need: 
-A PC with python (python version 2.7) installed.
-Python packages installed - xml.etree.ElementTree
		 	   - sys
-LRG file(s) downloaded from http://www.lrg-sequence.org

Instructions:
Save the LRG file in a folder. 
Open a Linux terminal or Windows command prompt within this folder.
Type python parse_LRG.py "filename.xml"
nb. the filename much match exactly including capital letters and including the extension (.xml) and in double quotation marks.
The output is saved in the existing directory and named: filename.fa

Further development:
The code is currently designed to allow retrieval of specific exons from a particular source sequence (eg. genomic sequence, cDNA, protein) by accessing elements of the exon dictionaries. This functionality is currently inaccessible from the terminal. Further development would allow the user to specify these options as further arguments at the terminal.

The set_exon_seq() function which retrieves the exon sequences can accept upstream and downstream flanking ranges, but these are currently defaulted to 0. Further development would enable the user to input desired values via the command line. This functionality is disabled as until further checks on the adjusted positions (i.e they fall within the source sequence) have been implemented (see comments on assert statements).