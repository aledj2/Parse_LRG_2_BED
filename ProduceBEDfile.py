'''
Created on 29 May 2015

@author: Aled
'''
import sys
import xml.etree.ElementTree as etree


class parse_LRG:
    def __init__(self):
        pass
    
    def open_XML(self,filename):
        self.xmlfile = filename
        
        self.tree = etree.parse(self.xmlfile)
        self.root = self.tree.getroot()
        
        self.id = self.set_id()                 #eg "LRG_292"
        
        self.genename=self.get_gene_name() # get Gene name
        #print self.genename
        
        self.chromnumber=0
        self.start_coord=0
        self.genomebuild="."
        self.get_chrom_num()
        print self.chromnumber
        print self.start_coord
        print self.genomebuild
        
        self.exons = self.set_exons()           #Nested dictionaries of exon number,reference seq,attribute(start, end, sequence)
        parse_LRG().create_bed_file(self.exons,self.genename,self.chromnumber,self.start_coord)
        
        if self.root.attrib['schema_version'] <> '1.8':
            print 'This LRG file is not the correct version for this script'
            print 'This is designed for v.1.8'
            print 'This file is v.' + self.root.attrib['schema_version']
    
    def set_id(self):
        '''
        Returns LRG id string from LRG xml
        e.g. "LRG_292"
        '''
        LRG_ids = self.root.findall("./fixed_annotation/id")
        assert len(LRG_ids) == 1                #Check that one LRG was found within the file
        for LRG in LRG_ids:
            LRG_id = LRG.text
            return LRG_id
        
    def get_gene_name(self):
        '''
        Returns LRG gene string from LRG xml
        e.g. "LRG_292"
        '''
        LRG_locus = self.root.findall("./updatable_annotation/annotation_set/lrg_locus")
        assert len(LRG_locus) == 1                #Check that one LRG was found within the file
        for locus in LRG_locus:
            genename = locus.text
            return genename
        
    def get_chrom_num(self):
        '''
        Returns LRG gene string from LRG xml
        e.g. "LRG_292"
        '''

        mapping = self.root.findall("./updatable_annotation/annotation_set/mapping")
        #assert len(LRG_locus) == 1                #Check that one LRG was found within the file
        for line in mapping:
            if line.attrib['coord_system'] == "GRCh37.p13":
                #print "hit"
                self.genomebuild= line.attrib['coord_system']
                self.chromnumber= line.attrib['other_name']
                self.start_coord= line.attrib['other_start']
            
        
    def set_exons(self):
        '''Returns dictionary containing exon numbers as keys, dictionary as values.
           Nested dictionary contains reference sequence as keys, dictionary as values
           Second nested dictionary contains "Start","End" as keys,
               coordinates within parent reference seq as values'''
        
        exons = {}
        for exon in self.tree.iter(tag = 'exon'):
            if 'label' in exon.attrib:
                #this_exon = {}
                #print "----------------\n\nExon:" , exon.attrib["label"]
                
                for coords in exon.findall('coordinates'):
                    if coords.attrib['coord_system'][-2:] == "p1":
                        pass
                    elif coords.attrib['coord_system'][-2:] == "t1":
                        pass
                    else:
                        exon_number = exon.attrib["label"]
                        #exon_ref = coords.attrib['coord_system']
                        exon_start = coords.attrib['start']
                        exon_end = coords.attrib['end']
                        #assert exon_start< exon_end          # assert that the exon end is after exon start
                        exons[int(exon_number)] = tuple([int(exon_start)]+[int(exon_end)])
                        #print this_exon
                
        return exons
    
    def create_bed_file(self,exon_dict,genename,chromnumber,startcoord):
        #print exon_dict
        #print len(exon_dict)
        
        for exon in exon_dict:
            exon_number=int(exon)
            exon_start=int(exon_dict[exon][0])+int(startcoord)
            exon_stop=int(exon_dict[exon][1])+int(startcoord)
            print "Chr"+chromnumber+"\t"+str(exon_start)+"\t"+str(exon_stop)+"\t"+genename+"_exon"+str(exon_number)

file2open="C:\Users\Aled\workspace\Parse_LRG_2_BED\LRG_292.xml"
parse_LRG().open_XML(file2open)