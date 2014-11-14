class LRG(object):
            
    def __init__(self, filepath='LRG_292.xml'):
        import xml.etree.ElementTree as etree
        self.xmlfile = filepath
        self.tree = etree.parse(self.xmlfile)
        self.root = self.tree.getroot()
        self.id = self.set_id()
        self.sequences = self.set_sequences()
        self.exons = self.set_exons()
        #self.references = self.set_references()
        self.set_exon_seq()

    def set_id(self):
        '''Returns LRG id string from LRG xml
            e.g. "LRG_292"          '''
        LRG_ids = self.root.findall("./fixed_annotation/id")
        for LRG in LRG_ids:
            LRG_id = LRG.text
        return LRG_id
                

    def set_sequences(self):
        '''Returns dictionary containing sequence ids as keys, sequences as values from LRG xml
        Includes genomic, transcript, translation
        e.g. {"LRG_292":"ATCG....", LRG_292t1:"ATCG....", LRG292p1:"AACE"}'''
        sequences = {}
        sequences.update(self.set_genomic_seq())
        sequences.update(self.set_cDNA())
        sequences.update(self.set_protein())
        return sequences

        
    def set_genomic_seq(self):
        '''Returns dictionary containing sequence ids as keys, sequences as values from LRG xml
        Genomic only
        e.g. {"LRG_292":"ATCG...."}'''
        genomic_id_seq = {}
        for item in (items for items in self.root[0] if items.tag == 'sequence'):
            genomic_seq = item.text
        for item in self.tree.iter(tag = 'id'):
            genomic_id = item.text
        genomic_id_seq[genomic_id]=genomic_seq
        return genomic_id_seq


    def set_cDNA(self):
        '''Returns dictionary containing sequence ids as keys, sequences as values from LRG xml
        Transcripts of genomic sequence only
        e.g. {"LRG_292t1":"ATCG...."}'''
        transcripts = {}
        transcript_elems = self.root.findall("./fixed_annotation/transcript")
        for transcript in transcript_elems:
            transcript_id = self.id + transcript.attrib["name"]
            for sequence in transcript.findall("./cdna/sequence"):
                transcript_seq = sequence.text
            transcripts[transcript_id] = transcript_seq
        return transcripts


    def set_protein(self):
        '''Returns dictionary containing sequence ids as keys, sequences as values from LRG xml
        Proteins from genomic transcripts only
        e.g. {"LRG_292p1":"AACE...."}'''
        proteins = {}
        protein_elems = self.root.findall("./fixed_annotation/transcript/coding_region/translation")
        for protein in protein_elems:
            print protein.tag
            protein_id = self.id + protein.attrib["name"]
            print type(protein_id)
            for sequence in protein.findall("./sequence"):
                protein_seq = sequence.text
            proteins[protein_id] = protein_seq
        return proteins


    def set_exons(self):
        '''Returns dictionary containing exon numbers as keys, dictionary as values.
           Nested dictionary contains reference sequence as keys, dictionary as values
           Second nested dictionary contains "Start","End" as keys,
               coordinates within parent reference seq as values'''
        exons = {}
        for exon in self.tree.iter(tag = 'exon'):
            if 'label' in exon.attrib:
                this_exon = {}
                print "----------------\n\nExon:" , exon.attrib["label"]
                exon_number = exon.attrib["label"]
                for coords in exon.findall('coordinates'):
                    print "\nReference: ", coords.attrib['coord_system']
                    print "Start: ", coords.attrib['start']
                    print "End: ", coords.attrib['end']
                    exon_ref = coords.attrib['coord_system']
                    exon_start = coords.attrib['start']
                    exon_end = coords.attrib['end']
                    
                    this_exon[exon_ref] = {"Start":exon_start,"End":exon_end}
                    print this_exon
                exons[exon_number] = this_exon
        return exons


    def set_references(self):
        '''Returns list of reference sequence names
           for genomic, transcript, translation sequences
           Not currently used as redundant (use self.sequences.keys() instead)'''
    
        for item in self.tree.iter(tag = 'id'):
            genomic_id = item.text

        transcript_ids = []
        for item in self.tree.iter(tag = 'transcript'):
            if 'name' in item.attrib:
                transcript_id = genomic_id + item.attrib["name"]
                transcript_ids.append(transcript_id)

        protein_ids = []
        for item in self.tree.iter(tag = 'translation'):
            if 'name' in item.attrib:
                protein_id = genomic_id + item.attrib["name"]
                protein_ids.append(protein_id)

        references = [genomic_id]
        for ids in transcript_ids, protein_ids:
            references.extend(ids)
        return references

    def get_exon_list(self):
        '''Prints a list of  sorted exon numbers and sequenes they map to within the LRG xml
           E.g. "2 ['LRG_292', 'LRG_292t1', 'LRG_292p1']"    '''

        sorted_exon_list = sorted(map(int, self.exons.keys()))
        for exon in sorted_exon_list:
            exon_refs = self.exons[str(exon)].keys()
            print "Exon:", exon, exon_refs
            
    def set_exon_seq(self, exon_list = None, upstream = 0, downstream = 0):
        '''Retrieves the exon sequence for each exon/reference sequence combination
           Saves the sequence in the second nested dictionary in self.exons
               i.e. in addition to keys "Start", "End"
                 now also contains {Sequence":exon sequence} '''
        output_fasta_file="fasta_fname" # name of the output fasta file can be changed to accept the filename as an argument
        f=open(output_fasta_file,"w") #creates an output file. If any data pre-exists it'll wipe the file
        f.close() # closes the file to prevent it being open while looping through
        if exon_list == None:
            exon_list = self.exons.keys()
            sorted_exon_list = sorted(map(int, exon_list))
            print sorted_exon_list

        for reference in self.sequences.keys():
            print "@@@@@@@@@@@@@@@@@@@@@@@@\n@@@@@@@@@@@@@@@@@@@@@@@@"
            print "\n\nRef:     ", reference, "\n"
            for exon in sorted_exon_list:
                exon = str(exon)
                print "Ref:     ", reference
                print "Exon:    ", exon
                if reference in self.exons[exon].keys():
                    exon_start = self.exons[exon][reference]["Start"]
                    exon_end = self.exons[exon][reference]["End"]
                    print "Start:   ", exon_start
                    print "End:     ", exon_end
                    print "Length:  ", int(exon_end) - int(exon_start) +1
                    adjusted_start = int(exon_start) - upstream
                    adjusted_end = int(exon_end) + downstream
                    
                    #get sequence from reference within this range
                    exon_seq = self.sequences[reference][int(exon_start)-1:int(exon_end)]
                    #save to nested dicts
                    self.exons[exon][reference]["Sequence"] = exon_seq
                        #cache results if time
                    print "Sequence:", exon_seq, "\n"
                    f=open(output_fasta_file,"a") #opens the file created before
                    Fasta_str_seq=">"+reference+" Exon: "+exon+"\n"+exon_seq+"\n\n" ## assembles all the info required for the output file
                    f.write(Fasta_str_seq)
                    f.close()
                else:
                    print "Start:    -"
                    print "End:      -"
                    
    def select_output(reference, exons, upstream, downstream):
        a=1
    
    def write_fasta(dictionary_of_id_seq):
        a=1

        
    
myLRG = LRG()
