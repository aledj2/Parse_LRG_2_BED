class LRG(object):
            
    def __init__(self, filepath='C:/lrg_xml_parser/LRG_292.xml'):
        import xml.etree.ElementTree as etree
        self.xmlfile = filepath
        self.tree = etree.parse(self.xmlfile)
        self.root = self.tree.getroot()
        self.id = self.set_id()
        self.sequences = self.set_sequences()
        self.exons = self.set_exons()
        #self.references = self.set_references()


    def set_id(self):
        LRG_ids = self.root.findall("./fixed_annotation/id")
        for LRG in LRG_ids:
            LRG_id = LRG.text
        return LRG_id
                

    def set_sequences(self):
        sequences = {}
        sequences.update(self.set_genomic_seq())
        sequences.update(self.set_cDNA())
        sequences.update(self.set_protein())
        return sequences

        
    def set_genomic_seq(self):
        genomic_id_seq = {}
        for item in (items for items in self.root[0] if items.tag == 'sequence'):
            genomic_seq = item.text
        for item in self.tree.iter(tag = 'id'):
            genomic_id = item.text
        genomic_id_seq[genomic_id]=genomic_seq
        return genomic_id_seq


    def set_cDNA(self):
        transcripts = {}
        transcript_elems = self.root.findall("./fixed_annotation/transcript")
        for transcript in transcript_elems:
            transcript_id = self.id + transcript.attrib["name"]
            for sequence in transcript.findall("./cdna/sequence"):
                transcript_seq = sequence.text
            transcripts[transcript_id] = transcript_seq
        return transcripts


    def set_protein(self):
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
                    
                    #Yo dawg...I heard you like dictionaries...
                    this_exon[exon_ref] = {"Start":exon_start,"End":exon_end}
                    print this_exon
                exons[exon_number] = this_exon
        return exons


    def set_references(self):

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
        sorted_exon_list = sorted(map(int, self.exons.keys()))
        print sorted_exon_list
        for exon in sorted_exon_list:
            exon_refs = self.exons[str(exon)].keys()
        
            print exon, exon_refs
            
    def get_exon_seq(self, exon_list = None, upstream = 0, downstream = 0):
        if exon_list == None:
            exon_list = self.exons.keys()
            sorted_exon_list = sorted(map(int, exon_list))
            print sorted_exon_list

        for reference in self.sequences.keys():
            print "@@@@@@@@@@@@@@@@@@@@@@@@\n@@@@@@@@@@@@@@@@@@@@@@@@"
            print "\n\nRef:     ", reference
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
                else:
                    print "Start:    -"
                    print "End:      -"
                    
    def select_output(reference, exons, upstream, downstream):
        a=1
    
    def write_fasta(dictionary_of_id_seq):
        a=1

        
    
myLRG = LRG()
