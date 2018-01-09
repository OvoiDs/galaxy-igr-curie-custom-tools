#!/usr/bin/env python2
# encoding: utf-8
# This program treats LRT Database and extracts variants
# It will also treat annotated variant files in order to add them to DB

import ParseTabAlt as TI
from argparse import ArgumentParser
from subprocess import Popen, PIPE
import re
from operator import itemgetter
from collections import defaultdict
import time
import sys

nuc = re.compile('[ATCG]+')

import fcntl
pid_file = '/data/galaxy/genome_db/BDD.lock'
fp = open(pid_file, 'w')
try:
    fcntl.lockf(fp, fcntl.LOCK_EX | fcntl.LOCK_NB)
except IOError:
    # another instance is running
    raise "Another instance is running. Please run the program again when all other job have finished"

print >>sys.stdout, "Command line:", sys.argv
parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input", type=str, required=False, help="Please input the either nothing, or an annotated file outputed by safir_create.py")
parser.add_argument('-d', '--database', dest='database', type=str, required=True, help='Please input the database tsv.gz file')
parser.add_argument('-m', '--mode', dest='mode', type=str, required=True, help='Please choose between "append" and "prepare" modes')
parser.add_argument('-o', '--output', dest='output', type=str, required=False, default=None, help="If you want a dump, please specify a path")
args = parser.parse_args()

# If "prepare" mode, only -d arg is to be provided
# If "append" mode, please provide input xls and database

# Gathered from: https://github.com/slowkow/pytabix
# Modified to suit our business

def bgzip(filename):
    """Call bgzip to compress a file."""
    Popen(['bgzip', '-f', filename]).wait()
    

def tabix_index(filename, chrom=1, start=2, end=2, skip=0, comment="#"):
    """Call tabix to create an index for a bgzip-compressed file."""
    chrom, start, end, skip, comment = str(chrom), str(start), str(end), str(skip), str(comment)
    Popen(['tabix', '-s', chrom, '-b', start, '-e', end, filename]).wait()

def tabix_query(filename, chrom, start):
    """Call tabix and generate an array of strings for each line it returns."""
    query = '{}:{}-{}'.format(chrom, start, start)
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
    for line in process.stdout:
        yield line.strip().split()

def read_database(filename):
    process = Popen(['zcat', filename], stdout=PIPE)
    output, error = process.communicate()
    process.wait()
    return output

def database_as_list(database):
    with open(database, "r") as dbh:
        dump = dbh.readlines()
    dump2 = [x.rstrip('\n') for x in dump]
    dump3 = [x.split('\t') for x in dump2]
    dump4 = [x[0:8] for x in dump3]
    return dump4

def pruneForward(pos, ref, alt):
    pos = int(pos)
    stop = False
    tref, talt = [], []
    for n,m in zip(ref, alt):
        if n == m and not stop:
            pos = pos + 1
            continue
        else:
            stop = True
            tref.append(n)
            talt.append(m)
    ref = ''.join(tref)
    alt = ''.join(talt)
    return pos, ref, alt

def pruneBackward(pos, ref, alt):
    stop = False
    counter = 0
    tref, talt = [], []
    for n, m in zip(reversed(ref), reversed(alt)):
        if n == m and not stop:
            counter += 1
            continue
        else:
            #print "count: ",counter, "ref: " ,n,"alt: ", m
            stop = True
            tref.append(n)
            talt.append(m)
    ref = ''.join(reversed(tref))
    alt = ''.join(reversed(talt))
    return pos, ref, alt

def validate(chr, pos, ref, alt, panel, nd, vp, date):
    if not chr.startswith('chr'):
        print 'Will not validate: chr does not start with chr'
        return None
    try:
        int(chr.split('chr')[1])
    except:
        print 'Will not validate: no-int chr'
        return None
    try:
        int(pos)
    except:
        print 'Will not validate: position cannot be turned to integer'
        return None
    if not nuc.search(ref.strip()):
        print 'Will not validate: reference nucleotide is not matching regex'
        return None
    if not nuc.search(alt.strip()):
        print 'Will not validate: alternate nucleotide is not matching regex'
        return None
    if alt.strip().strip().strip() == '':
        print 'Will not validate: non valid stripping'
        return None
    if ref.strip() == alt.strip():
        print 'Will not validate: alt and ref nucs are identical'
        return None
    if len(date.split('/')) != 3:
        print 'Will not validate: date is invalid'
        return None
    return chr, pos, ref.strip(), alt.strip(), panel, nd, vp, date

def first_cure(variants):
    cured = []
    tmp = []
    #banana = None
    for el in variants:
        #print len(el)
        if len(el) != 8: print el
        try:
            chr, pos, ref, alt, panel, nd, vp, date = el
        except:
            print "Invalid: {}".format(el)
            continue
        if validate(chr, pos, ref, alt, panel, nd, vp, date) is not None:
            tmp.append(el)
        else:
            print 'Could not validate: {}'.format(el)

    for el in tmp:
        chr, pos, ref, alt, panel, nd, vp, date = el
        if pos == "Position":
            print "Skipping header"
            continue
        if len(ref) == 1 or len(alt) == 1:
            cured.append([int(chr.lstrip('chr')), int(pos), ref, alt, panel, nd, vp, date])
        else:
            pos, ref, alt = pruneBackward(pos, ref.strip(), alt.strip())
            pos, ref, alt = pruneForward(pos, ref.strip(), alt.strip())
            cured.append([int(chr.lstrip('chr')), int(pos), ref, alt, panel, nd, vp, date])
    return sorted(cured, key=itemgetter(0, 1))

def create_var_dico(scured):
    varDic = defaultdict(lambda: defaultdict(list))
    for el in scured:
        chr = el[0]
        pos = el[1]
        ref, alt = el[2], el[3]
        panel, ann = el[4], el[5]
        pathos = el[6].upper()
        date = el[7]
        hash = '{0}-{1}-{2}-{3}'.format(chr, pos, ref, alt)
        if len(varDic[hash][pathos]) == 0:
            varDic[hash][pathos] = el
        else:
            elseDate = [int(x) for x in varDic[hash][pathos][7].split('/')]
            thisDate = [int(x) for x in date.split('/')]
            #Compare dates
            if thisDate[2] == elseDate[2] and thisDate[1] == elseDate[1] and thisDate[0] == elseDate[0]:
                pass 
            elif thisDate[2] > elseDate[2] or thisDate[1] > elseDate[1] or thisDate[0] > elseDate[0]:
                varDic[hash][pathos] = el
            else: print "Not recording this variant since {0} < {1}".format('/'.join([str(x) for x in elseDate]), '/'.join([str(x) for x in thisDate]))
    return varDic

def output_var_dico(varDic, output):
    toOutput = []
    for k, v in varDic.iteritems():
        for k2, v2 in v.iteritems():
            toOutput.append(v2)
    toOutput = sorted(toOutput, key=itemgetter(0, 1))
    with open(output, "w") as bddh:
        for el in toOutput:
            #print '\t'.join(['chr' + str(el[0])] + [str(x) for x in el[1:]])
            print >>bddh, '\t'.join(['chr' + str(el[0])] + [str(x) for x in el[1:]])
    regenerate_index(output)
    return toOutput

def regenerate_index(output):
    bgzip(output)
    tabix_index(output + ".gz")

if args.mode == "prepare":
    variants = database_as_list(args.input)

    scured = first_cure(variants)

    varDic = create_var_dico(scured)
    toOutput = output_var_dico(varDic, args.database)

    print "Finished preparing database. {0} variants imported.".format(len(toOutput))

elif args.mode == "append":
    # Input is some annotated file. We will use ParseTab so we do not care about tsv or xls filetype
    ti = TI.TabInput(args.input).parse()
    # Get index of columns we are interested in
    chr = ti.getIDByColname('Chr')
    pos = ti.getIDByColname('End_Position')
    ref = ti.getIDByColname('Reference_Seq')
    alt = ti.getIDByColname('Variant_Seq')
    try:
        panel = ti.getIDByColname('Panel')
    except:
        panel = 31
    try:
        ann = ti.getIDByColname('Annotator')
    except:
        ann = 30
    pathos = ti.getIDByColname('Manual_Var_Classif')
    date = time.strftime('%d/%m/%Y')

    l_extract = lambda x: [str(x[chr]), str(int(x[pos]) - len(x[ref]) + 1), str(x[ref]), str(x[alt]), str(x[panel]), str(x[ann]), str(x[pathos]), str(date) + '\n']

    database = read_database(args.database).split('\n')
    #print len(database)
    for line in ti.getLineIter():
        print "appending: {}".format(line)
        toAppend = l_extract(line)
        if toAppend[6] == '/':
            print 'No annotation provided. Skipping.'
        else:
            database.append('\t'.join(toAppend))
    #print database[-30:]
    #print len(database)
    
    # replace old database with new one
    with open(args.database + ".tmp", "w") as bddh:
        for line in database:
            if line.strip() != "": 
                print >>bddh, line.strip()
    variants = database_as_list(args.database + ".tmp")
    print len(variants)
    #for el in variants:
    #    print el
    #print variants
    scured = first_cure(variants)
    print len(scured)
    varDic = create_var_dico(scured) 
    toOutput = output_var_dico(varDic, args.database.rstrip('.gz'))

    print "Finished preparing database. {0} variants imported.".format(len(toOutput))

    print "Success"
    if args.output:
        a = read_database(args.database)
        with open(args.output, "w") as outputh:
            for line in a:
                outputh.write(line)
else:
    pass
