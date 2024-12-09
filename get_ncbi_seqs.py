#!/usr/env/bin python3
"""
get_ncbi_seqs.py

Uses Entrez to retrieve sequences from NCBI
Takes in text file with list of accession IDs
"""

from pathlib import Path
from Bio import Entrez
import textwrap
import argparse

Entrez.email = "yhsieh@alaska.edu"

## COMMON FUNCTIONS

def readaccids(infile):
    """Reads in file with accession IDs, returns a list"""
    accids = []
    with open(infile, 'r') as f:
        for line in f:
            accids.append(line.strip())
    return accids

def check_moltype(somestr):
    if somestr.lower() in 'proteinaaminoacid':
        return 'protein'
    elif somestr.lower() in 'nucleotidednarna':
        return 'nucleotide'

def batch_get(listofaccids, moltype='protein'):
    """Submits Entrez batch request for prot (default)
    or nucl seqs, returns handle.
    Use when you have relatively high confidence 
    your accids have existing entries
    
    :param listofaccids: list of str
    :param moltype: str
    
    :returns: http.client.HTTPResponse"""
    #TODO: add error handling
    nonemptylist = [item for item in listofaccids if item]
    moltype = check_moltype(moltype)
    request = Entrez.epost(moltype, id=",".join(nonemptylist))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db=moltype,\
                           retmode="xml",\
                               webenv=webEnv,\
                                   query_key=queryKey)
    return handle

def single_get(accid, moltype='protein'):
    """Submits single Entrez request for prot(default)
    or nucl seqs, returns handle
    
    :param accid: str
    :param moltype: str
    
    :returns: http.client.HTTPResponse"""
    moltype = check_moltype(moltype)
    try:
        handle = Entrez.efetch(db=moltype, id=accid, retmode="xml")
    except Exception as e:
        print(e)
        return 1
    else:
        return handle

def parse_entrez_handle_to_list(handleobj):
    """Parses Entrez handle containing one or more records
    into dict object
    
    :param handleobj: http.client.HTTPResponse
    
    :returns: list"""
    try:
        recordslist = list(Entrez.parse(handleobj))
    except Exception as e:
        print(e)
        return 1
    else:
        return recordslist

def find_organism(entrez_record):
    """Takes in record object, finds organism name"""
    feature_table_list = entrez_record['GBSeq_feature-table']
    for item in feature_table_list:
        idtag = item['GBFeature_key']
        featquals = item['GBFeature_quals']
        if idtag == 'source':
            for item in featquals:
                if item['GBQualifier_name'] == 'organism':
                    return item['GBQualifier_value']
                
## PROT FUNCTIONS

def find_protein(entrez_record):
    """Takes in record object, finds protein name"""
    feature_table_list = entrez_record['GBSeq_feature-table']
    for item in feature_table_list:
        idtag = item['GBFeature_key']
        featquals = item['GBFeature_quals']
        if idtag == 'Protein':
            for item in featquals:
                if item['GBQualifier_name'] == 'product':
                    return item['GBQualifier_value']
                
def find_gene_accid(entrez_record):
    """Takes in record object, returns GenBank accession"""
    tag1 = 'GBSeq_feature-table'
    tag2 = 'GBFeature_quals'
    tag3 = 'GBQualifier_name'
    tag4 = 'GBQualifier_value'
    try:
        infolist = entrez_record[tag1][-1][tag2]
    except Exception:
        return 1
    else:
        for gbdict in infolist:
            if gbdict[tag3] == 'coded_by':
                accid = gbdict[tag4].split(':')[0]
                if '(' in accid:
                    return accid.split('(')[1]
                else:
                    return accid
    
def get_gene_accid(handleobj):
    """Takes in handle, iterates find gene accid, returns list"""
    records = Entrez.parse(handleobj)
    geneaccid_list = []
    for record in records:
        geneaccid_list.append(find_gene_accid(record)) # could append 0/None
    return geneaccid_list

def write_gene_accid_to_file(accidlist, outfile="accids.txt"):
    """Writes out accids to a file"""
    with open(outfile, 'w+') as outf:
        for item in accidlist:
            if item: # skip None
                outf.write(f'{item}\n')

## FASTA FUNCTIONS

def extract_fasta_info(entrez_record, moltype='protein'):
    """Takes in a record object from Entrez handle,
    returns a dict with header and fasta"""
    moltype=check_moltype(moltype)
    organism = find_organism(entrez_record)
    seq = entrez_record['GBSeq_sequence']
    if moltype=='protein':
        accid = entrez_record['GBSeq_primary-accession']
        protein = find_protein(entrez_record)
        header = f'{accid}|{protein}|{organism}'
    elif moltype=='nucleotide':
        accid = entrez_record['GBSeq_locus']
        primeaccid = entrez_record['GBSeq_primary-accession']
        seqdef = entrez_record['GBSeq_definition']
        header = f'{accid}|{primeaccid}|{seqdef}'
        
        #### TODO: extract protid and append to header if it exists
        # for featitem in entrez_record['GBSeq_feature-table']:
        #     if not 'GBFeature_quals' in featitem.keys():
        #         protid = ''
        #     else:
        #         ...
        
        # if entrez_record['GBSeq_feature-table'][-1]['GBFeature_quals']:
        #     for entry in entrez_record['GBSeq_feature-table'][-1]['GBFeature_quals']:
        #         if entry['GBQualifier_name'] == 'protein_id':
        #             protid = entry['GBQualifier_value']
        #         else:
        #             protid = ''
        # else:
        #     protid = ''
        # header = f'{accid}|{primeaccid}|{protid}|{seqdef}'
        ####
        
    return {header:seq}

def write_fasta(fasta_dict, out_path,flag='w'):
    with open(out_path, flag) as f:
        for k, v in fasta_dict.items():
            f.write(''.join(['>', k, '\n']))
            seq = ''.join([v, '\n'])
            seq = seq.upper()
            wrappedseq = textwrap.wrap(seq, width=80)
            f.write('\n'.join(wrappedseq))
        f.write('\n')
        
def get_fastas_from_record(entrezrecordlist):
    pass

def write_fastas_to_file(somefastadict):
    pass

def get_fastas_from_entrez(handleobj, outputpath, seqtype='protein'):
    """Takes an Entrez handle object 
    and saves the seqs in fasta files"""
    records = Entrez.parse(handleobj)
    for record in records:
        outfile = f"{record['GBSeq_primary-accession']}.fasta"
        outpath = outputpath.joinpath(outfile)
        print(f"Writing to {outfile}")
        fastadict = extract_fasta_info(record, seqtype)
        write_fasta(fastadict, outpath)
        print('---> written')

def argumentsparser():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] accidfile outpath",
                                     description="Retrieves seq from NCBI",)
    parser.add_argument("accidfile", help="path to accession ID file", type=str)
    parser.add_argument("outpath", help="path to output files", type=str)
    args = parser.parse_args()
    return args

# if __name__=="__main__":
#     # args = argumentsparser()
    # accidfilepath = Path(args.accidfile)
    # outpath = Path(args.outpath)
    # idlist = readaccids(accidfilepath)
    # protrecs = batch_get_protein(idlist)
    # genelist = get_gene_accid(protrecs)
    # generecs = batch_get_nucleotide(genelist)
    # get_fastas_from_entrez(generecs, outpath, 'nucleotide')
    
    # accidfilepath=Path('accidsmall')
    # idlist=readaccids(accidfilepath)
    # protrecs = batch_get_protein(idlist)
    # genelist = get_gene_accid(protrecs)
    # generecs = batch_get_nucleotide(genelist)

        
# protid = 'KAI4832452.1'

# protrec = batch_get_protein([protid])

# geneid = get_gene_accid(protrec)

# generec = batch_get_nucleotide(geneid)

# recs = Entrez.parse(generec)
