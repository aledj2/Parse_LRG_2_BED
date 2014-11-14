from parse_LRG import *

def test_all():
    '''Execute all tests listed...'''
    test_set_id()
    
    
    print "----Testing complete----"



def test_set_id():
    '''Test that LRG id extracted from xml file matches expected value'''
    test_LRG = LRG("LRG_292.xml")
    test_id = test_LRG.set_id()

    print "test_id :", test_id
    assert test_id == "LRG_292", "Error: id != 'LRG_292'"

test_all()
