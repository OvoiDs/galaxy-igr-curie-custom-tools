#!/usr/local/bin/python2.7
# encoding: utf-8
'''

@author:     Yannick Boursin

@copyright:  2014 Institut Gustave Roussy. All rights reserved.

@contact:    elipsoid@gmail.com

@version:    stable - 1.2 - IonTorrent Suite 4.2.1 Variant Caller

@deffield    updated: Updated

@commentary: Changed some behavior: now prefered transcripts will be enforced !
'''

from __future__ import with_statement
from warnings import warn
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
import os
import time
from collections import defaultdict
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import xlwt
from subprocess import *
from ParseTab import TabInput as TI
import gzip
import random
import uuid

__all__ = []
__version__ = '1.3.1'
__date__ = '2014-03-05'
__updated__ = '2017-22-08'

SNPSIFT = "/data/softwares/snpEff/4.3/SnpSift.jar"
SNPEFF = "/data/softwares/snpEff/4.3/snpEff.jar"
DBSNP = "/data/db/dbsnp/150_GRCh37p13/All_20170710.vcf.gz"
SNPEFF_CONFIG = "/data/softwares/snpEff/4.3/snpEff.config"
DBNSFP = "/data/db/dbnsfp/2.9_20160107/dbNSFP2.9.txt.gz"
TABIX = "/data/softwares/htslib/1.5/bin/tabix"
REFSEQPOSITIONS = "/data/db/refSeqPositions/15062014/refGene_15062014.sorted.txt.gz"
COSMIC = "/data/db/Cosmic/Cosmic77CodingMuts.vcf.gz"
CANCERPANEL = "/data/db/CancerPanel/CancerPanel.tsv"
_1KG = "/data/db/1000Genome/all.vcf.gz"
ESP = "/data/db/ESP/ESP6500SI-V2-SSA137/all.ucsc.vcf.gz"
CLINICS = "/data/db/BDD_variants/BDD_variants.tsv.gz"
JAVA = "/data/softwares/java/1.8.0_144/bin/java"
tabFileName = ''

RESIDUE_DICO = {"Gly": "G",
                "Ala": "A",
                "Val": "V",
                "Leu": "L",
                "Ile": "I",
                "Met": "M",
                "Phe": "F",
                "Trp": "W",
                "Pro": "P",
                "Ser": "S",
                "Thr": "T",
                "Cys": "C",
                "Tyr": "Y",
                "Asn": "N",
                "Gln": "Q",
                "Asp": "D",
                "Glu": "E",
                "Lys": "K",
                "Arg": "R",
                "His": "H"
}

#####################################################################################################################
#############################################   Class definition   ##################################################
#####################################################################################################################  

def filter_dic(x):
    if (type(x) == dict):
        return x
    elif (type(x) == list):
        return x
    else:
        return str(x)

def removeDups(x):
    x = x.split(",")
    removed = list(set(x))
    return ",".join(removed)
    

class CancerPanel(object):
    def __init__(self, cancerpaneltsv):
        self.dico_panel = defaultdict(bool)
        with open(cancerpaneltsv) as saphirh:
            line = saphirh.readline()
            i = 0
            while (line != ""):
                if (i == 0):
                    line = saphirh.readline()
                    i += 1
                    continue
                else:
                    splitted = line.rstrip("\n").split("\t")
                    gene, ref, type, AZD2014, AZD4547, AZD5363, AZD8931, Selumetinib, Vandetanib, Olaparib = splitted
                    self.dico_panel[gene] = ObjectPanel(gene, ref, type, AZD2014, AZD4547, AZD5363, AZD8931, Selumetinib, Vandetanib, Olaparib)
                line = saphirh.readline()
        temp_NM = []
        for (key,value) in self.dico_panel.iteritems():
            ref = value.ref
            splref = ref.split(".")[0]
            temp_NM.append(splref)
        self.all_NM = temp_NM

class DB(object):
    def __init__(self, sName, barcode, runName):
        self.date = time.strftime("%c")
        self.sName = sName
        self.barcode = barcode
        self.runName = runName
        self.runDate = "None"
        self.histo = "None"
        self.cassette = "None"
        self.conclusion = "None"
        self.project = "None"
        self.batch = "None"
        self.gao = "None"

    def parse(self):
        a = self.runName.split('_')
        b = [x.split('-') for x in a]
        runNameArray = []
        # print b => [['R'], ['2014'], ['11'], ['06'], ['16'], ['16'], ['14'], [''], ['GAO', '469', 'RT00911'], ['SAFIR', 'Batch303'], ['poolC'], ['20141106']]
        b = [ x for item in b for x in item ]
        # print b => ['R', 1: '2014', 2: '11', 3: '06', 4: '16', 5: '16', 6: '14', '', 'GAO', '469', 'RT00911', 'SAFIR', 'Batch303', 'poolC', '20141106']   

        self.runDate = '{0}/{1}/{2} - {3}:{4}:{5}'.format(b[1], b[2], b[3], b[4], b[5], b[6])

        if ('GAO' in b):
            self.gao = b[b.index('GAO') + 1]
        else:
            self.gao = False

        self.project = []
        for el in b:
            if el.startswith("RT") or el.startswith('rt'):
                self.project.append(el)
        if (len(self.project) == 0):
            self.project = 'False'
        elif (len(self.project) == 1) and (type(self.project) == list):
            self.project = ''.join(self.project)
        elif (type(self.project) == list):
            self.project = ';'.join(self.project)
        else: self.project = 'False'
        
        self.batch = []
        for el in b:
            el = el.upper()
            if el.startswith('BATCH'):
                self.batch.append(el.lstrip('BATCH'))
        if (len(self.batch) == 0):
            self.batch = "False"
        elif (len(self.batch) == 1) and (type(self.batch) == list):
            self.batch = ''.join(self.batch)
        elif (type(self.batch) == list):
            self.batch = ';'.join(self.batch)
        else: self.batch = "False"

        #Cassette
        self.cassette = []
        for el in self.sName.split('_'):
            if el.startswith('K') or el.startswith('k'):
                self.cassette.append(el)
        if (len(self.cassette) == 0):
            self.cassette = "False"
        elif (len(self.cassette) == 1) and (type(self.cassette) == list):
            self.cassette = ''.join(self.cassette)
        elif (type(self.cassette) == list):
            self.cassette = ';'.join(self.cassette)
        else: self.cassette = "False"

        
    def addValues(self, histo, conclusion, annotator, cassette):
        self.histo = histo if histo != None else "None"
        self.conclusion = conclusion if conclusion != None else "None"
        self.annotator = annotator if annotator != None else "CLI"
        if (cassette != ""):
            self.cassette = cassette
 
    def __str__(self):
        try:
            tmp = []
            for el in self.project.split(';'):
                tmp.append('\t'.join([self.date, self.sName, self.barcode, self.runName, self.runDate, self.histo, self.cassette, self.conclusion, el, self.batch, self.gao, self.annotator]))
            return '\n'.join(tmp)
        except:
            # print type(self.cassette), type(self.conclusion), type(self.project)
            return 'An error occured'
    def printHeader(self):
        return '\t'.join(["Date", "Sample Name", "Barcode", "runName", "runDate", "N° Histo", "Cassette", "Conclusion", "Project", "Batch", "GAO", "Referent"])


class ObjectPanel(object):
    def __init__(self, gene, ref, type, AZD2014, AZD4547, AZD5363, AZD8931, Selumetinib, Vandetanib, Olaparib):
        self.gene = gene
        self.ref = ref
        self.type = type
        self.AZD201 = AZD2014
        self.AZD454 = AZD4547
        self.AZD536 = AZD5363
        self.AZD893 = AZD8931
        self.Selumetinib = Selumetinib
        self.Vandetanib = Vandetanib
        self.Olaparib = Olaparib

class OutputXls(object):
    def __init__(self,outputfile, cancerpanel, db):
        self.cancer = cancerpanel
        self.output = outputfile
        self.db = db
        self.workbook = xlwt.Workbook()
        self.vc = self.workbook.add_sheet("VC", cell_overwrite_ok=False)
        self.dbs = self.workbook.add_sheet("DB", cell_overwrite_ok=False)
        self.NM_ref = self.workbook.add_sheet("NM_ref", cell_overwrite_ok=False)
        self.columns = ["Global_Conclusion", "Manual_Var_Comment", "Manual_Var_Classif", "Gene_Symbol", "Protein_Change", "Exon",
                        "Variant_Freq", "Position_Cov", "RefSeq_Id", "cDNA_Change", "Codon", "Chr", "Start_Position",
                        "End_Position", "Strand", "Type", "Reference_Seq", "Variant_Seq", "Variant_Cov", "MAF_classification",
                        "ESP_Freq", "by1000G_Freq", "DbSNP_Id", "COSMIC_Id", "SIFT_Prediction", "Polyphen2_Prediction",
                        "Quality", "Strand_Bias", "Amplicon_Ref", "Panel"]
        self.columns_db = ["Date", "Sample Name", "Barcode", "runName", "runDate", "N_Histo", "Cassette", "Conclusion", "Project", "Batch", "GAO", "Referent"]
        
    def column_headers(self):
        i = 0
        for element in self.columns:
            self.vc.write(0,i,element)
            i += 1
    def column_DB(self):
        i = 0
        for element in self.columns_db:
            self.dbs.write(0,i,element)
            i += 1
    
    def write_DB(self):
        i = 0
        k = 1
        tmp = []
        for el in self.db.project.split(';'):
            tmp.append([self.db.date, self.db.sName, self.db.barcode, self.db.runName, self.db.runDate, self.db.histo, self.db.cassette, self.db.conclusion, el, self.db.batch, self.db.gao, self.db.annotator]) 
        for element in tmp:
            for el2 in element:
                self.dbs.write(k,i,el2)
                i += 1
            k += 1

    def write_cancer_panel(self):
        header = ["Gene", "Ref", "Type", "AZD2014", "AZD4547", "AZD5363", "AZD8931", "Selumetinib", "Vandetanib", "Olaparib"]
        i = 0
        for element in header:
            self.NM_ref.write(0,i,element)
            i += 1
        k = 1
        for key,value in self.cancer.dico_panel.iteritems():
            self.NM_ref.write(k,0,value.gene)
            self.NM_ref.write(k,1,value.ref)
            self.NM_ref.write(k,2,value.type)
            self.NM_ref.write(k,3,value.AZD201)
            self.NM_ref.write(k,4,value.AZD454)
            self.NM_ref.write(k,5,value.AZD536)
            self.NM_ref.write(k,6,value.AZD893)
            self.NM_ref.write(k,7,value.Selumetinib)
            self.NM_ref.write(k,8,value.Vandetanib)
            self.NM_ref.write(k,9,value.Olaparib)
            k += 1
    
    def write_TSV(self, outname):
        output = outname
        with open(output, "w") as tab:
            tab.write("\t".join(self.columns) + "\n")
            temp = []
            for variant in self.variants:
                temp_temp = []
                for value in variant:
                    if (type(value) == list): 
                        if (len(value) == 1): 
                            value = dict(value[0])
                        else:
                            temp_value = []
                            for el in value:
                                temp_value.append(el["id"])
                            value = ", ".join(temp_value)
                            temp_temp.append(value)
                    if (type(value) == dict):
                                temp_temp.append(value["id"])
                    else:
                        temp_temp.append(value)
                temp.append(temp_temp)
            for el in temp:
                el = [ str(x) for x in el ]
                tab.write("\t".join(el) + "\n")
            
    def write_info(self):
        i = 1
        style = xlwt.XFStyle()
        
        fnt = xlwt.Font()
        fnt.name = 'Arial'
        
        borders = xlwt.Borders()
        
        pattern = xlwt.Pattern()
        pattern.pattern = xlwt.Pattern.SOLID_PATTERN
        pattern.pattern_fore_colour = 0x0A
        
        style = xlwt.XFStyle()
        style.font = fnt
        style.borders = borders
        style.pattern = pattern
        
        self.vc.col(19).width = 256*32
        self.vc.col(23).width = 256*32
        self.vc.col(8).width = 256*16
        self.vc.col(28).width = 256*32
        self.vc.col(25).width = 256*32
        for variant in self.variants:
            k = 0
            for value in variant:
                if (type(value) == list): 
                    if (len(value) == 1): 
                        value = dict(value[0])
                    else:
                        temp_value = []
                        for el in value:
                            temp_value.append(el["id"])
                        value = ", ".join(temp_value)
                    
                if (type(value) == dict):
                    if ("colorize" in value.keys()):
                        if ("hyperlink" in value.keys()):
                            self.vc.write(i,k, xlwt.Formula('HYPERLINK("{0}";"{1}")'.format(value["hyperlink"], value["id"])), style)
                        else:   
                            self.vc.write(i,k, value["id"], style)
                    else:
                        if ("hyperlink" in value.keys()):
                            self.vc.write(i,k, xlwt.Formula('HYPERLINK("{0}";"{1}")'.format(value["hyperlink"], value["id"])))
                        else: 
                            self.vc.write(i,k, value["id"])
                else:
                    self.vc.write(i,k, value)
                k += 1
            i += 1
    
    def get_variant_info(self, dico_variants):
        
        def map_position_on_transcript(exons, position, cds):
                    temp_exon = []
                    for x,y in exons:
                        temp_exon.append((int(x),int(y)))
                    exons = temp_exon
                    position = int(position)
                    intra_exon_len = [ (y - x) for x,y in exons ]
                    relative_position = 0
                    exon_nb = 1
                    diff_flag = False
                    for exon in exons:
                        if (position > exon[0]) and (position < exon[1]):
                            relative_position += position - exon[0]
                            exon_nb += 1
                            break
                        elif (position > exon[0]):
                            diff_cds = exon[1] - int(cds[0])
                            if (diff_cds <= 0):
                                exon_nb += 1
                                continue
                            else:
                                if not (diff_flag):
                                    relative_position += diff_cds
                                    diff_flag = True
                                else:
                                    relative_position += intra_exon_len[exon_nb - 1]
                            exon_nb += 1                            
                        else:
                            continue
                    return relative_position + 1
                
        def map_effect_to_maf(effect, type):
            effect_dico={
                         "coding_sequence_variant": "",
                         "chromosome": "",
                         "coding_sequence_variant": "Missense_Mutation",
                         "inframe_insertion": "In_Frame_Ins",
                         "disruptive_inframe_insertion": "In_Frame_Ins",
                         "inframe_deletion": "In_Frame_Del",
                         "disruptive_inframe_deletion": "In_Frame_Del",
                         "downstream_gene_variant": "",
                         "exon_variant": "",
                         "exon_loss_variant": "In_Frame_Del",
                         "frameshift_variant": "Frame_Shift",
                         "gene_variant": "",
                         "intergenic_region": "IGR",
                         "conserved_intergenic_variant": "IGR",
                         "intragenic_variant": "",
                         "intron_variant": "Intron",
                         "conserved_intron_variant": "Intron",
                         "miRNA": "RNA",
                         "non_coding_exon_variant": "Silent",
                         "missense_variant": "Missense_Mutation",
                         "initiator_codon_variant": "Silent",
                         "stop_retained_variant": "Silent",
                         "rare_amino_acid_variant": "Missense_Mutation",
                         "splice_acceptor_variant": "Splice_Site",
                         "splice_donor_variant": "Splice_Site",
                         "splice_region_variant": "Splice_Site",
                         "splice_region_variant": "Splice_Site",
                         "splice_region_variant": "Splice_Site",
                         "stop_lost": "Nonstop_Mutation",
                         "5_prime_UTR_premature_start_codon_gain_variant": "Translation_Start_Site",
                         "start_lost": "Translation_Start_Site",
                         "stop_gained": "Nonsense_Mutation",
                         "synonymous_variant": "Silent",
                         "start_retained": "Silent",
                         "stop_retained_variant": "Silent",
                         "transcript_variant": "RNA",
                         "regulatory_region_variant": "RNA",
                         "upstream_gene_variant": "",
                         "3_prime_UTR_variant": "3'UTR",
                         "3_prime_UTR_truncation": "3'UTR",
                         "5_prime_UTR_variant": "5'UTR",
                         "5_prime_UTR_truncation": "5'UTR" }
            try:
                translated_effect = effect_dico[effect]
            except:
                try:
                    print 'Unrecognized effect: {0}'.format(effect)
                    uEffect = effect.split("+")
                    print uEffect
                    try:
                        uEffect_translated = [effect_dico[x] for x in uEffect]
                        translated_effect = "+".join(uEffect_translated)
                    except:
                        print uEffect
                        try:
                            print [effect_dico[x] for x in uEffect]
                        except: 
                            translated_effect = "{0}".format(effect)
                except:
                    translated_effect = "{0}".format(effect)
            if (translated_effect == "Frame_Shift"):
                if (type == "DEL"):
                    translated_effect = "Frame_Shift_Del"
                elif (type == "INS"):
                    translated_effect = "Frame_Shift_Ins"
                else: translated_effect = "{0}".format(effect)
            return translated_effect
                
        def filter(element, i, freq):
            globbool = True
            globbool = globbool and (float(freq) > float(0))
            if (re.search(r"regulatory|non_coding|intron|UTR_v|UTR_t|upstream|downstream|intergenic|synonymous_variant|start_retained|stop_retained_variant|intragenic_variant", element.highest_id)
                and not re.search(r"stop_gained|start_lost|splice|rare_amino_acid_variant|missense_variant|coding_sequence_variant|insertion|deletion|exon_loss_variant|frameshift", element.highest_id)):
                 globbool = False
            return globbool
        
        self.variants = []
        counter = 0
        glob_counter = 0
        print('DISCARDED VARIANTS')
        print('["Concl.Bio", "Comm.Bio", "Var.Class", "Gene", "p.", "exon", "Var.freq", "Pos.Cov.", "Ref.NM", "c.", "Codon", "Chr.", "Start_Position", "End_Position", "Strand", "type", "ref.seq", "var.seq", "var.cov", "Variant_Class", "ESP Freq", "1000G freq", "DbSNP ID", "COSMIC ID", "SIFT prediction", "Polyphen-2 prediction", "Quality", "Strand Bias", "Amplicon"]')
        for element in dico_variants:
            howMany = len(element.genotype[0]["AO"].split(","))
            i = 0
            while (i < howMany):
                print element.highest_id
                do_not_write = 0
                glob_counter += 1
                #transcript_start = element.info[element.transcriptID.split(".")[0]]
                #try:
                #    exon_start = transcript_start.split('|')[0].split('-')[0].split(',')[0:-1]
                #    exon_stop = transcript_start.split('|')[0].split('-')[1].split(',')[0:-1]
                #    cds = transcript_start.split('|')[1].split('~')
                #except:
                #    print transcript_start
                #exons = zip(exon_start, exon_stop)
                
                #position_on_transcript = map_position_on_transcript(exons, element.position, cds)
                temp_list = []
                temp_list.append(tabFileName) #0
                temp_list.append("") #1
                annot_clinics = element.info['clinics'] if element.info['clinics'] != "" else "/"
                temp_list.append(annot_clinics) #2
                temp_list.append(element.geneName) #3
                try:
                    p = element.HGVSp if element.HGVSp != "" else "/"
                    c = element.HGVSc if element.HGVSc != "" else "/"
                except:
                    print "ERROR while parsing p & c"
                try:
			if (p != "/"):
                        	regexp = re.compile(r"p\.([0-9]+)?(Gly|Ala|Val|Leu|Ile|Met|Phe|Trp|Pro|Ser|Thr|Cys|Tyr|Asn|Gln|Asp|Glu|Lys|Arg|His)(_?[0-9]+)(Gly|Ala|Val|Leu|Ile|Met|Phe|Trp|Pro|Ser|Thr|Cys|Tyr|Asn|Gln|Asp|Glu|Lys|Arg|His)(fs|del)?")
                        	match = regexp.search(p)
				aa1 = match.group(2)
                        	aa2 = match.group(4)
                        	new_aa1 = RESIDUE_DICO[aa1]
                        	new_aa2 = RESIDUE_DICO[aa2]
				group1 = match.group(1) if match.group(1) != None else ""
				group3 = match.group(3) if match.group(3) != None else ""
				group5 = match.group(5) if match.group(5) != None else ""
                        	new_codon_prot = "p.{0}{1}{2}{3}{4}".format(group1, new_aa1, group3, new_aa2, group5)
                        	p = new_codon_prot
		except:
			pass
                temp_list.append(p) #4
                temp_list.append(element.exon) #5
                try:
                    fao = element.info["FAO"]
                    cov = element.coverage
                    gt = Decimal(100) * Decimal(Decimal(fao)/Decimal(cov)).quantize(Decimal('0.0001'), ROUND_HALF_UP)
                    temp_list.append(gt) #6
                except:
                    do_not_write = 1
                    fao = ""
                    temp_list.append("ERROR") #6
                try:
                    fdp = Decimal(element.coverage)
                except:
                    fdp = ""
                temp_list.append(fdp) #7

                #temp_list.append(element.transcriptID)
                if (element.transcriptID.split(".")[0] not in self.cancer.all_NM): #8
                    temp_list.append({"id": element.transcriptID,
                                      "colorize": 1,
                                      "hyperlink": "http://www.ncbi.nlm.nih.gov/nuccore/" + element.transcriptID})
                else:
                    temp_list.append({"id": element.transcriptID,
                                      "hyperlink": "http://www.ncbi.nlm.nih.gov/nuccore/" + element.transcriptID})
                regex_codon2 = re.search(r"([0-9]+)", p)
                codon2 = regex_codon2.group(0) if regex_codon2 != None else "/"
		codon = c
		if (codon != "/"):
                	temp_list.append({"id": c,
				"hyperlink": "http://localhost:10000/show?request={0}:{1}".format(element.transcriptID, c)})
                else: temp_list.append("/")
		temp_list.append(codon2)
                temp_list.append(element.chrName) #11
                pos = {"id": element.good_pos,
                       "hyperlink": "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=" + element.chrName + ":" + str(int(element.position) - 50) + "-" + str(int(element.position) + 50)}
                temp_list.append(pos) #12 http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr5%3A145508588-145508788
                offset = len(element.refNucleo) - 1 #13
                temp_list.append(str(int(element.good_pos) + offset)) #14
		try:
                	temp_list.append(str(element.info["STRAND[{0}]".format(element.transcriptID.split('.')[0])])) #15
		except:
			temp_list.append(str(element.info["STRAND"]))
                temp_list.append(element.type.split(',')[i].upper()) #16
                temp_list.append(element.refNucleo) #17
                temp_list.append(element.altNucleo.split(",")[i]) #18
                temp_list.append(fao) #19
                temp_list.append(map_effect_to_maf(element.highest_id, element.info["TYPE"].split(',')[i].upper())) #20
                try:
                    esp = Decimal( element.esp.split(',')[2] ).quantize(Decimal('0.0001'), ROUND_HALF_UP)
                except:
                    esp = "/"
                temp_list.append(esp) #21
                try:
                    _1000g = Decimal(100) * Decimal(element.al_freq_alt).quantize(Decimal('0.0001'), ROUND_HALF_UP)
                except:
                    _1000g = "/"
                temp_list.append(_1000g) #22
                id = element.dbId.split(";")
                if (id != ['']):
                    rs = []
                    junk = []
                    for el in id:
                        if (re.match(r"^rs", el)):
                            rs.append(el)
                        else:
                            junk.append(el)
                    rs = removeDups(",".join(rs))
                    hyperlink_rs = [] 
                    if (rs.split(",")[0] != ''):
                        for el in rs.split(","):
                            hyperlink_rs.append({"id": el,
                                                 "hyperlink": "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=" + el.lstrip("rs")})
                    else: hyperlink_rs = "/"
                    cosm = element.info["COSMIC[" + element.altNucleo.split(",")[i] + "]"]
                    hyperlink_cosm = []
                    if (cosm != ""):
                        hyperlink_cosm.append({"id": cosm,
                                                   "hyperlink": "http://cancer.sanger.ac.uk/cosmic/search?q=" + cosm})
                    else:
                        hyperlink_cosm = ("/")
                    temp_list.append(hyperlink_rs) #23
                    temp_list.append(hyperlink_cosm) #24
                else:
                    temp_list.append("/")
                    temp_list.append("/")
                prot = element.uniprot_acc.split(",")
                sift = element.sift.split(',')
                polyphen2 = element.polyphen2.split(',')
                sift = [ x for x in sift if x != '.' ]
                polyphen2 = [ x for x in polyphen2 if x != '.']
                if (len(prot) == len(sift)) and (sift != ['']) and (prot != ['']):
                    a_sift = [ "{0}({1})".format(x,y) for x,y in zip(sift, prot) ]
                else:
                    a_sift = sift
                if (len(prot) == len(polyphen2)) and (polyphen2 != ['']) and (prot != ['']):
                    a_polyphen = [ "{0}({1})".format(x,y) for x,y in zip(polyphen2, prot) ]
                else:
                    a_polyphen = polyphen2
                a_sift = ", ".join(a_sift)
                a_polyphen = ", ".join(a_polyphen)
                temp_list.append(a_sift) #25
                temp_list.append(a_polyphen) #26
                temp_list.append(element.quality) #27
                temp_list.append(element.std_bias) #28
                temp_list.append(element.amp) #29
                temp_list.append(element.info['PANEL']) # Panel ?
                boolean = filter(element, i, gt)
                #converted_list = [ filter_dic(x) for x in temp_list ]
                converted_list = temp_list
                if (do_not_write == 0 and boolean == True):
                    counter += 1
                    self.variants.append(converted_list)
                else:
                    print(converted_list)
                i += 1
        print "Nombre de variants dans le rapport final: {0}".format(counter)
        print "Nombre de variants en sortie du VariantCaller: {0}".format(glob_counter)
        print "Ratio: {0}%".format(float(float(counter)/float(glob_counter)) * 100)
            
    def save(self):
        self.workbook.save(self.output)
       
class Annotation(object):
    def __init__(self, annotations, type_of_analysis, cancer):
        self.annotations = annotations
        self.type = type_of_analysis
        self.cancer = cancer
        
    #Here are defined the ranks for the different analysis paradigm
    def rank(self, annotation, prioritize):
        add = 1
        if (prioritize is True):
            add = 100000000000000
        if (self.type == "exonseq"):
            Very_High = re.compile(r"missense_variant|gained|lost|frameshift|exon_loss_variant|disruptive_inframe|coding_sequence_variant|inframe_insertion|inframe_deletion")
            High = re.compile(r"splice|initiator_codon_variant|stop_retained_variant|protein_protein_contact|structural_interaction_variant")
            Low = re.compile(r"regulatory|non_coding|intron|UTR|upstream|downstream|intergenic")
            Medium = re.compile(r"synonymous_variant|start_retained|stop_retained_variant|intragenic_variant")
            #How informative exactly is that annotation ?? Well, in exome-seq, let's give it a Medium-Low weight, so we will give an advantage to Synonymous coding
            Unknown = re.compile(r"exon_variant")
        elif (self.type == "cancer"):
            Very_High = re.compile(r"missense_variant|gained|lost|frameshift|exon_loss_variant|disruptive_inframe|coding_sequence_variant|inframe_insertion|inframe_deletion")
            High = re.compile(r"splice|initiator_codon_variant|stop_retained_variant|protein_protein_contact|structural_interaction_variant")
            Low = re.compile(r"regulatory|non_coding|intron|UTR|upstream|downstream|intergenic")
            Medium = re.compile(r"synonymous_variant|start_retained|stop_retained_variant|intragenic_variant")
            Unknown = re.compile(r"exon_variant")
        
        if (Very_High.search(annotation)): return add * 10000000000
        elif (High.search(annotation)): return add * 1000000000
        elif (Low.search(annotation)): return add * 1
        elif (Medium.search(annotation)): return add * 50000
        elif (Unknown.search(annotation)): return add * 25000
    
    #This method parses the annotations and ranks them using the rank method
    def parse(self):
        howManyAnnotations = len(self.annotations.split(','))
        self.annot = defaultdict(list)
        if (howManyAnnotations > 1): #If more than one
            i = 0
            Highest_Length = 0
            if (self.type == "cancer"):
                to_look_in = self.cancer.all_NM
            else:
                to_look_in = []
            while i < howManyAnnotations:
                annotation_details = self.annotations.split(",")[i].split("(")[1].split("|")
                if (annotation_details[4] == ''): annotation_length = 1
                else: annotation_length = int(annotation_details[4])
                if (annotation_length > Highest_Length): Highest_Length = annotation_length
                i += 1
            i = 0
            while i < howManyAnnotations:
                annotation_type = self.annotations.split(',')[i].split('(')[0]
                annotation_details = self.annotations.split(",")[i].split("(")[1].split("|")
                if (annotation_details[4] == ''): annotation_length = 1
                else: annotation_length = int(annotation_details[4])
                #We mitigate ranks by length.
                try:
                    current_annotation = [annotation_type, annotation_details, annotation_length,i]
                    if (current_annotation[1][8].split('.')[0] in to_look_in): 
                        prioritize = True
                    else: prioritize = False
                    current_rank = Decimal(self.rank(annotation_type, prioritize)) * Decimal(Decimal(annotation_length) / Decimal(Highest_Length))
                    print current_rank
                    current_rank = int(current_rank.quantize(Decimal('1'), rounding=ROUND_HALF_UP))
                except:
                    print "Debug"
                    print "Annotation details"
                    print annotation_details
                    print "Annotation type"
                    print annotation_type
                    print "Annotation Rank"
                    print self.rank(annotation_type, True)
                    print "Protein Length"
                    print annotation_length
                    print "Highest Length"
                    print Highest_Length
                    raise
                if (len(self.annot[current_rank]) == 0):
                    self.annot[current_rank].append(current_annotation)
                else:
                    self.annot[current_rank].append(current_annotation)
                i += 1
        else:   #If there is only one annotation, we take it no matter what and give it a standard weight of 1000
            annotation_type = self.annotations.split(',')[0].split('(')[0]
            annotation_details = self.annotations.split(",")[0].split("(")[1].split("|")
            if (annotation_details[4] == ''): annotation_length = 1
            else: annotation_length = int(annotation_details[4])
            Highest_Length = annotation_length
            
            current_annotation = [annotation_type, annotation_details, annotation_length, 0]
            
            self.annot[1000].append(current_annotation)

        score_array = self.annot.keys()
        best_score = sorted(score_array)[-1]
        work_on = self.annot[best_score]
        #This function must return a 3-tuple containing these informations
        return (work_on[0][1], work_on[0][0], work_on[0][2], work_on[0][3])
                
          
class Variant(object):
    def __init__(self, cancer_panel):
        #standard VCF informations
        self.header = defaultdict(list)
        self.chrName = ""
        self.position = ""
        self.refNucleo = ""
        self.altNucleo = ""
        self.dbId = "."
        self.quality = ""
        self.info = defaultdict(str)
        self.genotype = defaultdict(lambda: defaultdict(str))
        self.line = ""
        self.info_dico = defaultdict(list)
        self.sample_names = []
        self.tag = ""
        self.genoraw = []
        self.cancer = cancer_panel

    def annotate(self):
        #additionnal informations present in annotated VCF files   
        #Format EFF : Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | 
        #Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | 
        #Transcript_ID | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )'
        self.coverage = self.info["COVERAGE"]
        self.effect = self.info["EFF"]
        self.uniprot_acc = self.info["dbNSFP_Uniprot_acc"]
        self.esp = self.info["ESPMAF"]
	self.amp = self.info["AMPLICON"].rstrip("'").lstrip("'")
        self.std_bias = self.info["STD_BIAS"]
        self.type = self.info["TYPE"]
        self.polyphen2 = self.info["dbNSFP_Polyphen2_HVAR_pred"]
        self.al_freq_alt = self.info["1000G"]
        self.al_freq_ref = ""
        self.cosm = False
        self.good_pos = self.info["GOOD_POS"]
        if (self.info["COSMIC"]): self.cosm = True
        self.sift = self.info["dbNSFP_SIFT_pred"]
        
        if (self.effect != ""):
            Ann = Annotation(self.effect,"cancer", self.cancer)
            (splitted_effect,self.highest_id,self.length,annot_id) = Ann.parse()
            self.geneName = splitted_effect[5]
            self.functionalClass = splitted_effect[1].replace("_", " ")
            self.codonChange = splitted_effect[2]
            self.aaChange = splitted_effect[3]
            splitted_aaCh = self.aaChange.split("/")
            self.HGVSp = ""
            self.HGVSc = ""
            for elem in splitted_aaCh:
                if (elem.startswith("p.")):
                    self.HGVSp = elem
                if (elem.startswith("c.")):
                    self.HGVSc = elem
            self.transcriptID = splitted_effect[8]
            self.exon = splitted_effect[9]
        else:
            self.length, self.highest_id, self.geneName, self.functionalClass, self.codonChange = "", "", "", "", ""
            self.aaChange, self.transcriptID, self.exon = "", "", ""
            self.HGVSc = ""
            self.HGVSp = ""
        
        #CAF : An ordered, comma delimited list of allele frequencies based on 1000Genomes, 
        #starting with the reference allele followed by alternate alleles as ordered in the ALT column. 
        #Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the 
        #ALT column.  The minor allele is the second largest value in the list, and was previously reported in VCF 
        #as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter
        
        self.info_dico = [self.chrName, self.position, self.refNucleo, self.altNucleo,
                          self.dbId.lstrip("rs"), self.quality, self.geneName,
                          self.functionalClass, self.codonChange, self.aaChange,
                          self.transcriptID, self.exon, self.uniprot_acc, self.esp,
                          self.polyphen2, self.al_freq_ref, self.cosm, self.length,
                          self.highest_id.replace("_", " "), self.coverage, self.al_freq_alt, self.sift]
    
            
class VCF():
    def __init__(self, VCF_file, data_vcf, global_samples, cancer_panel):
        self.VCF_file = VCF_file
        self.VCF_handle = open(VCF_file, "r")
        self.VCF = data_vcf
        self.global_samples = global_samples
        self.cancer = cancer_panel
    def parse(self):
        tags = []
        a = 0
        line = self.VCF_handle.readline()
        while (line != ''):
            c_VCF = Variant(self.cancer)
            if line.startswith("#"):
                if (line.startswith("#CHROM")):
                    newsamplelist = []
                    Samples = list(line.rstrip('\n').split('\t')[9:])
                    HowMany = len(Samples)
                    for sample in Samples:
                        self.global_samples.append(sample)
                if line.startswith("##"):
                    splitted = line.split("=")
                    if (len(splitted) == 2):
                        c_VCF.header[splitted[0]].append(splitted[1])
                line = self.VCF_handle.readline()
                continue
            

            # Informations basiques communes pour une ligne à un ou plusieurs échantillons
            field = line.split("\t")
            #On compte le nombre d'échantillons du VCF
            c_VCF.sample_names = Samples # noms des échantillons (array)
            c_VCF.nb_samples = HowMany # nombre d'échantillons
            offset = HowMany - 1
            
            #On crée une empreinte de l'objet et on la stocke.
            #Elle est ensuite comparée avec les autres empreintes.
            #Si l'objet existe déjà, on passe en mode "addition"
            #sinon on le crée
            
            c_VCF.tag = field[0] + field[1] + field[3] + field[4]
            tags.append(c_VCF.tag)
            c_VCF.chrName = field[0]
            c_VCF.position = field[1]
            c_VCF.dbId = field[2]
            c_VCF.refNucleo = field[3]
            c_VCF.altNucleo = field[4]
            c_VCF.quality = field[5]
            filter = False
            g = line.split("\t")[:9]
            c_VCF.line = '\\t'.join(g) #Ligne brute. On échappe les whitespaces car javascript reparsera la chaine
        
            #Champ INFO : on indexe les informations contenues dans ce champ dans un dictionnaire
            for value in field[7].split(";"):
                splitted_value = value.split("=")
                if (len(splitted_value) == 1):
                    c_VCF.info[splitted_value[0]] = True
                elif(len(splitted_value) == 2):
                    c_VCF.info[splitted_value[0]] = splitted_value[1]
                else:
                    print("Error {0}".format(splitted_value))
                    continue
            infs = ""
            for k, v in c_VCF.info.iteritems():
                if k.startswith('CLINICS'):
                    if k.startswith('CLINICS[{0}]'.format(c_VCF.altNucleo)):
                        if infs != "":
                            infs = infs + ", " + v
                        else: infs = v
            c_VCF.info['clinics'] = infs
            print c_VCF.tag, c_VCF.info['clinics']
            #Samples. On crée un tableau contenant pour chaque échantillon les informations dans un dictionnaire
            id_sample = c_VCF.nb_samples - HowMany
            offset = 0
            
            for sample in Samples: #field[9:] => données sur les échantillons (VCF multi-sample)
                selected_field = field[9 + offset]
                c_VCF.genotype[id_sample]["raw"] = ''.join(field[9 + offset]).rstrip("\n")
                value_counter = 0
                temp_genotype_info = selected_field.split(":")
                for value in field[8].split(":"): #Champ FORMAT qui indique le format des données contenues dans SAMPLE
                    #c_VCF.genotype est un dictionnaire de dictionnaire du format suivant :
                    # defaultdict { sample_number : { field_name : value ... } ... }
                    try:
                        c_VCF.genotype[id_sample][value] = temp_genotype_info[value_counter]
                    except:
                        continue
                    value_counter += 1
                c_VCF.genotype[id_sample]["name"] = sample.replace('-', '_')
                offset += 1
                id_sample += 1
                #print "New Length = {0}".format(sample_counter)
                
            #c_VCF, un objet "Variant", contient donc les informations de la ligne parsée
            line = self.VCF_handle.readline() #On passe à la ligne suivante
            
            #On annote le fichier VCF, permettant de produire le tableau de résultat indispensable à la méthode json puis on ajoute
            #les informations retournées par la méthode json dans un tableau (tableau de données passées à javascript).
            #Ensuite, on ajoute l'objet Variant crée (un objet variant par ligne de fichier VCF) dans un tableau pour
            #une hypothétique utilisation future.
            c_VCF.annotate()
            if not filter: self.VCF.append(c_VCF)
        return (self.VCF, self.global_samples)
    
    
        
        
#####################################################################################################################
################################# Argument parsing ##################################################################
#####################################################################################################################

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

######################################################################################################################
#################################################   main   ###########################################################
######################################################################################################################

def main(argv=None):
    global tabFileName
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
    
%s

Created by Yannick Boursin on %s.
Copyright 2014 Institut Gustave Roussy. All rights reserved.


USAGE
''' % (program_version_message, program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-i", "--input", dest="input", type=str, required=True, help="Please input the TSV file outputed by IonTorrent VariantCaller")
        parser.add_argument("-o", "--output", dest="output", help="output report name", required=True, type=str)
        parser.add_argument("-t", "--tabular", dest="tab", help="output tabular file report name", required=True, type=str)
        parser.add_argument("-n", "--name", dest="name", required=True, type=str)
        parser.add_argument("--histo", dest="histo", type=str, default="/")
        parser.add_argument("--conclusion", dest="conclusion", type=str, default="/")
        parser.add_argument("--cassette", dest="cassette", type=str, required=False, default="/")
        parser.add_argument("--annotator", dest="annotator", type=str, default="/")
        parser.add_argument("--panel", dest='panel', type=str, default="/")
	# Process arguments
        args = parser.parse_args()
        tabFile = args.input
        output = args.output
        cancerpanel = CANCERPANEL
        histo, conclusion = args.histo, args.conclusion
        annotator = args.annotator
	tabFile2 = args.tab
        tabFileName = args.name
        cassette = args.cassette
        panel = args.panel

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2
    

    #Initialisation des variables utilisées ultérieurement
    data_vcf, global_samples = [], []
    cancer_panel = CancerPanel(cancerpanel)

    #Nouvelle étape: parsing du fichier tabulé seul et génération d'un fichier VCF    
    #On parse le fichier tabulé
    tab_data = TI(tabFile)
    tab_data.parse()
    
    vcf_base = "temp" + str(uuid.uuid4())

    #On crée un VCF avec une annotation par ligne
    vcf_data = []
    sName, barcode, runName = False, False, False
    for line in tab_data.getLineIter():
        temp = []
        if not sName:
            sName = line[tab_data.getIDByColname("Sample Name")]
        if not barcode:
            barcode = line[tab_data.getIDByColname("Barcode")]
        if not runName:
            runName = line[tab_data.getIDByColname("Run Name")]
        chrom = line[tab_data.getIDByColname("Chrom")]
        pos = line[tab_data.getIDByColname("VCF Position")]
        id = line[tab_data.getIDByColname("Allele Name")] if line[11].startswith("COSM") else "."
        ref = line[tab_data.getIDByColname("VCF Ref")]
        alt = line[tab_data.getIDByColname("VCF Variant")]
        qual = line[tab_data.getIDByColname("Quality")]
        filter = line[tab_data.getIDByColname("Allele Call")]
        amplicon = line[tab_data.getIDByColname("Region Name")]
        bias = line[tab_data.getIDByColname("Strand Bias")]
        allele_freq = line[tab_data.getIDByColname("Frequency")]
        cov = line[tab_data.getIDByColname("Coverage")]
        type = line[tab_data.getIDByColname("Type")]
        good_pos = line[tab_data.getIDByColname("Position")]
        fao = line[tab_data.getIDByColname("Allele Cov")]
        info = "AMPLICON='{0}';STD_BIAS={1};TYPE={2};COVERAGE={3};GOOD_POS={4};FAO={5}".format(amplicon, bias, type, cov, good_pos, fao)
        gt_info = "GT"
        gt = allele_freq
        temp.append(str(chrom))
        temp.append(str(pos))
        temp.append(str(id))
        temp.append(str(ref))
        temp.append(str(alt))
        temp.append(str(qual))
        temp.append(str(filter))
        temp.append(str(info))
        temp.append(str(gt_info))
        temp.append(str(gt))
        if not ((filter == "Absent") or (filter == "No Call")):
            line = "\t".join(temp)
            vcf_data.append(line)
    newDBEntry = DB(sName, barcode, runName)
    newDBEntry.parse()
    print newDBEntry.printHeader()
    newDBEntry.addValues(histo, conclusion, annotator, cassette)
    print newDBEntry

    vcf_header = []
    vcf_header.append("##fileformat=VCFv4.1")
    vcf_header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tXXX\n")
    str_header = "\n".join(vcf_header)
    with open(vcf_base + ".vcf","w") as tada:
        tada.write(str_header)
        tada.write("\n".join(vcf_data))

    vcf_input=vcf_base + ".vcf"

    #Pour chaque fichier:
    print "File : {0}".format(vcf_input)
    print "Annotation with SnpEFF & SnpSift"
    
    #Begin Annotation Step
    filename1 = vcf_base + ".dbsnp.vcf"
    filename2 = vcf_base + ".dbsnp.eff.vcf"
    filename3 = vcf_base + ".dbsnp.eff.dbnsfp.vcf"
    filename4 = vcf_base + ".annotated.vcf"
    filename5 = filename4.rstrip(".dbsnp.eff.dbnsfp.vcf")
    report = vcf_input + ".snpEff.html"
    
    #Creating File Objects in order to gather annotations
    filename1_h = open(filename1, "w")
    filename2_h = open(filename2, "w")
    filename3_h = open(filename3, "w")
    
    #Annotation
    SnpSift1 = [JAVA, "-jar", "-Xmx8000m", SNPSIFT, "annotate", "-v", DBSNP, vcf_input]
    SnpEff1 = [JAVA, "-jar", "-Xmx8000m", SNPEFF, "eff", "-c", SNPEFF_CONFIG, "-v","-formatEff", "-t", "-noStats", "-noMotif", "-noNextProt", "hg19", filename1]
    SnpSift2 = [JAVA, "-jar", "-Xmx8000m", SNPSIFT, "dbnsfp", "-v", "-db", DBNSFP, filename2]
    print " ".join(SnpSift1)
    print " ".join(SnpEff1)
    print " ".join(SnpSift2)

    #Execution des programmes d'annotation
    
    out1 = Popen(SnpSift1, bufsize=-1, stdout=filename1_h, stderr=PIPE)
    output1, err1 = out1.communicate()
    out1.wait()
    print err1
    filename1_h.close()
        
    out2 = Popen(SnpEff1, bufsize=-1, stdout=filename2_h, stderr=PIPE)
    output1, err2 = out2.communicate()
    out2.wait()
    print err2
    filename2_h.close()
        
    out3 = Popen(SnpSift2, bufsize=-1, stdout=filename3_h, stderr=PIPE)
    output1, err3 = out3.communicate()
    out3.wait()
    print err3
    filename3_h.close()
    
    #Annotating strand orientation using tabix and dbNSFP
    i = 0
    with open(filename3, "r") as filename3_h:
        with open(filename4, "w") as filename4_h:
            line = filename3_h.readline()
            while (line != ""):
                k = False
                if (line.startswith("#")):
                    filename4_h.write(line)
                    pass
                else:
                    splitted = line.split('\t')
                    tabix1 = "{0}:{1}-{1}".format(splitted[0].lstrip("chr"), splitted[1])
                    tabix2 = "{0}:{1}-{1}".format(splitted[0], splitted[1])
                    ref = splitted[3]
                    alt = splitted[4].split(",")
                    cmd1 = [TABIX, COSMIC, tabix1]
                    cmd2 = [TABIX, REFSEQPOSITIONS, tabix2]
                    cmd3 = [TABIX, _1KG, tabix1]
                    cmd4 = [TABIX, ESP, tabix1]
                    cmd5 = [TABIX, CLINICS, tabix2]
                    out4 = Popen(cmd1, bufsize=-1, stdout=PIPE, stderr=PIPE)
                    output_1, error_1 = out4.communicate()
                    out4.wait()
                    out5 = Popen(cmd2, bufsize=-1, stdout=PIPE, stderr=PIPE)
                    output_2, error_2 = out5.communicate()
                    out5.wait()
                    out6 = Popen(cmd3, bufsize=-1, stdout=PIPE, stderr=PIPE)
                    output_3, error_3 = out6.communicate()
                    out6.wait()
                    out7 = Popen(cmd4, bufsize=-1, stdout=PIPE, stderr=PIPE)
                    output_4, error_4 = out7.communicate()
                    out7.wait()
                    out8 = Popen(cmd5, bufsize=-1, stdout=PIPE, stderr=PIPE)
                    output_5, error_5 = out8.communicate()
                    out8.wait()
                    cosmic_output = output_1.split("\n")
                    cosmic_output = [ x.split('\t') for x in cosmic_output ]
                    refseq_output = output_2.split("\n")
                    refseq_output = [ x.split('\t') for x in refseq_output ]
                    _1KG_output = output_3.split('\n')
                    _1KG_output = [x.split('\t') for x in _1KG_output]
                    ESP_output = output_4.split('\n')
                    ESP_output = [x.split('\t') for x in ESP_output]
                    clinics_output = output_5.split('\n')
                    clinics_output = [x.split('\t') for x in clinics_output]
                    print clinics_output
                    for nuc in alt:
                        for cel in cosmic_output:
                            if (len(cel) == 8) and (nuc == cel[4]):
                                splitted[7] = splitted[7] + ";COSMIC[{0}]={1}".format(nuc, cel[2])
                        for el in refseq_output:
                            if (len(el) == 16):
                                splitted[7] = splitted[7] + ";" + "STRAND=" + el[3] + ";" + "STRAND[{0}]=".format(el[1]) + el[3] + ";" + el[1] + "=" + el[9] + "-" + el[10] + "|" + el[6] + "~" + el[7]
                        for el in _1KG_output:
                            if (len(el) == 8):
                                splitted[7] = splitted[7] + ";1000G=" + el[7].split(';')[1].split('=')[1]
                        for el in ESP_output:
                            if (len(el) == 8):
                                splitted[7] = splitted[7] + ";ESPMAF=" + el[7].split(';')[4].split('=')[1]
                        cnt=1
                        for el in clinics_output:
                            if len(el) == 8 and el[3] == nuc:
                                toadd = el[6]
                                splitted[7] = splitted[7] + ";CLINICS[{0}]{2}={1}".format(nuc, toadd,cnt)
                                cnt += 1 
                    splitted[7] = splitted[7] + ";PANEL={}".format(panel) + ";ANNOTATEUR={}".format(annotator)
                    to_write = "\t".join(splitted)
                    filename4_h.write(to_write)

                line = filename3_h.readline()
            
    print "Now parsing annotated VCF"
    
    #On initialise et on parse le VCF, en prennant garde à bien garder les variables
    init_data = VCF(filename4, data_vcf, global_samples, cancer_panel)
    data_vcf, global_samples = init_data.parse()
    #Liste des échantillons
    list_of_samples = '|'.join(global_samples)
    list_of_samples = list_of_samples.replace('-', '_')
    book = OutputXls(output, cancer_panel, newDBEntry)
    book.column_headers()
    book.get_variant_info(data_vcf)
    book.write_info()
    book.write_cancer_panel()
    book.write_TSV(tabFile2)
    book.column_DB()
    book.write_DB()
    book.save()

    with open('/data/db/BDD_variants/Annotated.db', 'a') as dbh:
        print "Adding entry to DB"
        dbh.write(str(newDBEntry) + "\n")
    
if __name__ == "__main__":
    sys.exit(main())
