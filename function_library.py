#!/usr/bin/env python3

# by Robert R Butler III
# Library of python functions

# error logging function plus config
def error_handling():
    import logging
    import sys
    logging.basicConfig(level=logging.WARNING)
    
    return ' {}. {}, line: {}'.format(sys.exc_info()[0],
                                           sys.exc_info()[1],
                                           sys.exc_info()[2].tb_lineno)


# opens gzip text files
def myopen(fn):
    import gzip
    try:
        h = gzip.open(fn)
        ln = h.read(2) # read arbitrary bytes so check if @param fn is a gzipped file
    except:
        # cannot read in gzip format
        return open(fn)
    h.close()
    return gzip.open(fn, 'rt')

# looks up rsids from RS_MERGE and RS_HISTORY see rsupdate scripts for details
def lookup_rsid(rsid): # needs 
    import re
    if not re.match(r'^\d+$', rsid):
        return rsid
    elif rsid not in RS_MERGE: # rs number not appear in RS_MERGE -> there is no merge on this rs
        return rsid
    elif rsid in RS_MERGE: # lift rs number
        rsLow, rsCurrent = RS_MERGE[rsid]
        if rsCurrent not in RS_HISTORY:
            print("replacing %s with %s" % (rsid, rsCurrent))
            return rsCurrent
        else:
            print("%s has been updated to %s which is retired, last known rsid was %s" % (rsid, rsCurrent, rsLow))
            return rsLow
    else:
        return rsid

# updates columns in table files using pandas,re 
def update_csv(file, up_column):
    import pandas as pd
    import re
    print("reading %s" % file)
    csv_tbl = pd.read_csv(file, sep='\t', dtype=object)
    print("updating %s in %s" % (up_column, file))
    csv_tbl[up_column] = csv_tbl[up_column].map(function_to_update_stuff) # modify this column to change function
    outfile = re.sub(r'current_filename', r'updated_filename', file)
    print("saving %s" % outfile)
    csv_tbl.to_csv(outfile, sep='\t', index=False, na_rep='NA')
    return

# looks up info from dbSNP based on chromosome and position, returns dict of rsid key and allele values (not ordered)
def read_dbsnp(CHR, CHRPOS):
    from Bio import Entrez
    import re

    Entrez.email = "no@no.org"
    search = Entrez.read(Entrez.esearch(db="snp", term="%s[CHR] AND %s[CHRPOS]" % (CHR, CHRPOS)))
    snp_hits = search["IdList"]
    snp_str = ",".join(snp_hits)
    local_snps = Entrez.read(Entrez.esummary(db="snp", id=snp_str))
    SNP_MATCHES = dict()
    for snp in local_snps:
        match = re.search(r"[A|G|C|T]{5}\[(\S+)\][A|G|C|T]{5}", snp['DOCSUM'])
        match_list = match.group(1).split("/")
        ref = snp['SNP_ID']
        SNP_MATCHES[ref] = match_list
        print(ref, SNP_MATCHES[ref])
    return SNP_MATCHES

# sifts through files matching a pattern in a dir and runs function
def file_list(mydir, pattern):
    import glob, os
    dir_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(mydir)
    for filename in glob.glob(pattern):
        function_of_choice(args) # replace with desired function
    os.chdir(dir_path)
    return

# loads *.coverage file with no headers into dict. For cov files made with "samtools depth" function.
# Assumes filenames are research_ids for samples 
def load_coverages_dict(mydir, filename):
    import re
    covs_key = re.match(r'(\d+).coverage', filename)
    print("reading %s" % filename)
    coverages[covs_key.group(1)] = {}
    file = open('%s%s' % (mydir, filename))
    for ln in file:
        f_row = ln.strip().split('\t')
        chrpos = '%s:%s' % (f_row[0].lstrip('chr'), f_row[1])
        site_cov = f_row[2]
        coverages[covs_key.group(1)][chrpos] = site_cov
    file.close()
    return

# loads csv as dataframe and performs function to manipulate table
def func_on_csv(file):
    import pandas as pd
    import re
    print("reading %s" % file)
    csv_tbl = pd.read_csv(file, sep='\t', dtype=str)
    print("running function on %s" % file)
    csv_tbl['plusminus3'] = csv_tbl.apply(lambda x: plusminus3(x['Chr'], x['Position'], x['ResearchID']), axis=1) # modify this column to change function
    outfile = re.sub(r'current_filename', r'updated_filename', file)
    print("saving %s" % outfile)
    csv_tbl.to_csv(outfile, sep='\t', index=False, na_rep='NA')
    return

# Identifies snps in vcf +- three positions from reference vcf, requires ref_tbl to be loaded globally
def plusminus3(Chr, Position, ResearchID):
    import pandas as pd
    pos = int(Position)
    pos_list = []
    for n in range(pos - 3, pos + 4):
        pos_list.append(str(n))
    hits_tbl = ref_tbl[(ref_tbl['Chr'] == Chr) & (ref_tbl['ResearchID'] == ResearchID) & ref_tbl['Position'].isin(pos_list)]
    if len(hits_tbl.index) != 0:
        print('%s:%s_%s has match(es)' % (Chr, Position, ResearchID))
        hits_tbl.to_csv('plusminus3_matches/%s:%s_%s.vcf' % (Chr, Position, ResearchID), sep='\t', index=False)
        return True
    else:
        return False
