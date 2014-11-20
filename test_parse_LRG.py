from parse_LRG import *
#from nose.tools import *

# User notes.
#This file has been created to ensure that the file parse_LRG.py processed the file LRG_292.xml as expected.
#This test will only work with this file as static information is required to best test the code's function.
#This will catch any failures in the code designed to retrieve the desired information
#This script can be modified to work with a alternate LRG file but many elements below will need replacing manually.

filename = "LRG_292.xml"

def test_all():
    '''Execute all tests listed...'''
    test_set_id()
    test_exons()
    print "----Testing complete----"


def test_set_id():
    '''Test that LRG id extracted from xml file matches expected value'''
    test_LRG = LRG(filename)            # opens the LRG file and runs it through the LRG class imported from parse_LRG.py
    test_id = test_LRG.set_id()

    print "test_id :", test_id
    assert test_id == "LRG_292", "Error: id != 'LRG_292'"       # ensure the ID atribute has been correctly selected


def test_exons():
    ''' This function checks that the number of exons is correct, that there are the correct number of transcripts in each exon, these are named correctly and the sequences are correct '''
    test_LRG = LRG(filename)
    
    #Correct number of exons for each transcript
    assert len(test_LRG.exons) == 23, "Error: Incorrect number of exons in dictionary" # ensure all exons have been selected

    #Correct number of reference sequences in exon. 
    assert len(test_LRG.exons["1"]) == 2, "Error: Incorrect number of reference sequences for exon 1" # exon 1 is UTR so does not have protein sequence
    assert len(test_LRG.exons["2"]) == 3, "Error: Incorrect number of reference sequences for exon 2"
    assert len(test_LRG.exons["11"]) == 3, "Error: Incorrect number of reference sequences for exon 11"
    assert len(test_LRG.exons["22"]) == 3, "Error: Incorrect number of reference sequences for exon 22"
    assert len(test_LRG.exons["23"]) == 3, "Error: Incorrect number of reference sequences for exon 23"

    #Correct identity of reference sequences in exon
    assert test_LRG.exons["1"].keys() == ["LRG_292","LRG_292t1"], "Error: Incorrect identity of reference sequences in exon 1"
    assert test_LRG.exons["2"].keys() == ["LRG_292","LRG_292t1", "LRG_292p1"], "Error: Incorrect identity of reference sequences in exon 2"
    assert test_LRG.exons["11"].keys() == ["LRG_292","LRG_292t1", "LRG_292p1"], "Error: Incorrect identity of reference sequences in exon 11"
    assert test_LRG.exons["22"].keys() == ["LRG_292","LRG_292t1", "LRG_292p1"], "Error: Incorrect identity of reference sequences in exon 22"
    assert test_LRG.exons["23"].keys() == ["LRG_292","LRG_292t1", "LRG_292p1"], "Error: Incorrect identity of reference sequences in exon 23"
    
    #Correct sequence of exon in genomic DNA. (This was pulled out manually). This can be repeated for cDNA (LRG_292t1) and protein (LRG_292p1)
    assert test_LRG.exons["1"]["LRG_292"]["Sequence"] == "GTACCTTGATTTCGTATTCTGAGAGGCTGCTGCTTAGCGGTAGCCCCTTGGTTTCCGTGGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTGCGCGGCGTGAGCTCGCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAG", "Error: incorrect sequence for exon 1, reference LRG_292"
    assert test_LRG.exons["2"]["LRG_292"]["Sequence"] == "TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTG", "Error: incorrect sequence for exon 2, reference LRG_292"
    assert test_LRG.exons["11"]["LRG_292"]["Sequence"] == "GTGAAGCAGCATCTGGGTGTGAGAGTGAAACAAGCGTCTCTGAAGACTGCTCAGGGCTATCCTCTCAGAGTGACATTTTAACCACTCAG", "Error: incorrect sequence for exon 11, reference LRG_292"
    assert test_LRG.exons["22"]["LRG_292"]["Sequence"] == "GGTGTCCACCCAATTGTGGTTGTGCAGCCAGATGCCTGGACAGAGGACAATGGCTTCCATG", "Error: incorrect sequence for exon 2, reference LRG_292"
    assert test_LRG.exons["23"]["LRG_292"]["Sequence"] == "CAATTGGGCAGATGTGTGAGGCACCTGTGGTGACCCGAGAGTGGGTGTTGGACAGTGTAGCACTCTACCAGTGCCAGGAGCTGGACACCTACCTGATACCCCAGATCCCCCACAGCCACTACTGACTGCAGCCAGCCACAGGTACAGAGCCACAGGACCCCAAGAATGAGCTTACAAAGTGGCCTTTCCAGGCCCTGGGAGCTCCTCTCACTCTTCAGTCCTTCTACTGTCCTGGCTACTAAATATTTTATGTACATCAGCCTGAAAAGGACTTCTGGCTATGCAAGGGTCCCTTAAAGATTTTCTGCTTGAAGTCTCCCTTGGAAATCTGCCATGAGCACAAAATTATGGTAATTTTTCACCTGAGAAGATTTTAAAACCATTTAAACGCCACCAATTGAGCAAGATGCTGATTCATTATTTATCAGCCCTATTCTTTCTATTCAGGCTGTTGTTGGCTTAGGGCTGGAAGCACAGAGTGGCTTGGCCTCAAGAGAATAGCTGGTTTCCCTAAGTTTACTTCTCTAAAACCCTGTGTTCACAAAGGCAGAGAGTCAGACCCTTCAATGGAAGGAGAGTGCTTGGGATCGATTATGTGACTTAAAGTCAGAATAGTCCTTGGGCAGTTCTCAAATGTTGGAGTGGAACATTGGGGAGGAAATTCTGAGGCAGGTATTAGAAATGAAAAGGAAACTTGAAACCTGGGCATGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGCAGATCACTGGAGGTCAGGAGTTCGAAACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAATACAGAAATTAGCCGGTCATGGTGGTGGACACCTGTAATCCCAGCTACTCAGGTGGCTAAGGCAGGAGAATCACTTCAGCCCGGGAGGTGGAGGTTGCAGTGAGCCAAGATCATACCACGGCACTCCAGCCTGGGTGACAGTGAGACTGTGGCTCAAAAAAAAAAAAAAAAAAAGGAAAATGAAACTAGAAGAGATTTCTAAAAGTCTGAGATATATTTGCTAGATTTCTAAAGAATGTGTTCTAAAACAGCAGAAGATTTTCAAGAACCGGTTTCCAAAGACAGTCTTCTAATTCCTCATTAGTAATAAGTAAAATGTTTATTGTTGTAGCTCTGGTATATAATCCATTCCTCTTAAAATATAAGACCTCTGGCATGAATATTTCATATCTATAAAATGACAGATCCCACCAGGAAGGAAGCTGTTGCTTTCTTTGAGGTGATTTTTTTCCTTTGCTCCCTGTTGCTGAAACCATACAGCTTCATAAATAATTTTGCTTGCTGAAGGAAGAAAAAGTGTTTTTCATAAACCCATTATCCAGGACTGTTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTGTACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTCCA", "Error: incorrect sequence for exon 2, reference LRG_292"

    
test_all()
