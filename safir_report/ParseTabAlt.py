#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Remplir plus tard

@author:     Yannick Boursin

@copyright:  2014 Institut Gustave Roussy. All rights reserved.

@contact:    elipsoid@gmail.com

@version:    stable - beta 1

@deffield    updated: Updated
'''
import xlrd
from warnings import warn
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
import os
from math import log10
from decimal import Decimal, ROUND_HALF_UP
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import magic

__all__ = []
__version__ = 1.0
__date__ = '2014-03-05'
__updated__ = '2014-23-05'

class TabInput(object):
    def __init__(self, path_to_tab):
        self.map_things = {}
        self.dico = []
        self.path = path_to_tab
        self.type = self.sniff()
        if (self.type == "xls"):
            self.tabh = xlrd.open_workbook(self.path)
        elif (self.type == "tsv"):
            self.tabh = open(path_to_tab, "r")
        else:
            self.tabh = None
        
        
    def sniff(self):
        m = magic.Magic(mime=True)
        mime = m.from_file(self.path)
        xls_mimes = ["application/CDFV2-corrupt", "application/vnd.ms-excel", "application/msexcel",
                      "application/x-msexcel", "application/x-ms-excel", "application/x-excel",
                      "application/x-dos_ms_excel", "application/xls", "application/x-xls", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                      "application/vnd.ms-office", "CDF V2 Document, No summary info"]
        tsv_mimes = ["text/plain"]
        if (mime in xls_mimes):
            return "xls"
        elif (mime in tsv_mimes):
            return "tsv"
        else:
            print mime 
            raise TypeError(mime)
    
    def parse(self, sheet_nb=0):
        if (self.type == "tsv"):
            line_nb = 0
            for line in self.tabh:
                if (line_nb == 0):
                    header = line.split("\t")
                    index = 0
                    line_nb += 1
                    for value in header:
                        self.map_things[value.rstrip('\n')] = index
                        index += 1
                    self.check_length = index
                else:
                    line = line.rstrip("\n").split("\t")
                    index = 0
                    temp_dico = {}
                    for value in line:
                        temp_dico[index] = value 
                        index += 1
                    if (index == self.check_length):
                        self.dico.append(temp_dico)
                    line_nb += 1
            return self
        elif (self.type == "xls"):
            sheets = self.tabh.sheet_names()
            current_sheet = self.tabh.sheet_by_name(sheets[sheet_nb])
            c = current_sheet
            
            nrows = c.nrows
            i = 0
            while (i < nrows):
                if (i == 0):
                    header = c.row_values(0,0)
                    index = 0
                    i += 1
                    for value in header:
                        self.map_things[value.rstrip('\t').rstrip('\n')] = index
                        index += 1
                    self.check_length = index
                else:
                    line = c.row_values(i,0)
                    index = 0
                    temp_dico = {}
                    for value in line:
                        if (index in [1, 14]):
                            temp_dico[index] = value
                        else:
                            temp_dico[index] = value
                        index += 1
                    if (index == self.check_length):
                        self.dico.append(temp_dico)
                    i += 1
            return self
        else: raise TypeError('input file is neither text nor xls {0}'.format(self.type))
    
    def getIDByColname(self, Colname):
        #print self.map_things
        return self.map_things[Colname]
    def getInfoByLineAndColname(self, line_number, Colname):
        return self.dico[line_number][self.getIDByColname(Colname)]
    def getLine(self, line_number):
        return self.dico[line_number]
    def getLineIter(self):
        for line in self.dico:
            yield line
    def getLineNumberByPosition(self, position):
        i = 0
        for line in self.dico:
            pos = line[self.getIDByColname("VCF Position")]
            if (pos == position):
                return i
            i += 1
        return "Not Found"
        
    def getPath(self):
        return self.path
