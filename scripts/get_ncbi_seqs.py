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

def readaccids(infile):
    """Reads in file with accession IDs, returns a list"""
    accids = []
    with open(infile, 'r') as f:
        for line in f:
            accids.append(line.strip())
    return accids

def batch_get(listofaccids):
    """Submits Entrez batch request, returns handle"""
    request = Entrez.epost("protein", id=",".join(listofaccids))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="protein", retmode="xml", webenv=webEnv, query_key=queryKey)
    return handle


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
    except Exception as e:
        return 0
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
        geneaccid_list.append(find_gene_accid(record))
    return geneaccid_list

def write_gene_accid_to_file(accidlist, outfile="accids.txt"):
    """Writes out accids to a file"""
    with open(outfile, 'w+') as outf:
        for item in accidlist:
            if item:
                outf.write(f'{item}\n')

def extract_fasta_info(entrez_record):
    """Takes in a record object from Entrez handle,
    returns a dict with header and fasta"""
    accid = entrez_record['GBSeq_primary-accession']
    organism = find_organism(entrez_record)
    protein = find_protein(entrez_record)
    seq = entrez_record['GBSeq_sequence']
    header = f'{accid}|{protein}|{organism}'
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


def get_fastas_from_entrez(handleobj, outputpath):
    """Takes an Entrez handle object 
    and saves the seqs in fasta files"""
    records = Entrez.parse(handleobj)
    for record in records:
        outfile = f"{record['GBSeq_primary-accession']}.fasta"
        outpath = outputpath.joinpath(outfile)
        print(f"Writing to {outfile}")
        fastadict = extract_fasta_info(record)
        write_fasta(fastadict, outpath)
        print('---> written')

    
def argumentsparser():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] accidfile outpath",
                                     description="Retrieves seq from NCBI",)
    parser.add_argument("accidfile", help="path to accession ID file", type=str)
    parser.add_argument("outpath", help="path to output files", type=str)
    args = parser.parse_args()
    return args

if __name__=="__main__":
    # args = argumentsparser()
    # accidfilepath = Path(args.accidfile)
    # outpath = Path(args.outpath)
    # idlist = readaccids(accidfilepath)
    # testrecs = batch_get(idlist)
    # get_fastas_from_entrez(testrecs, outpath)
    
    accidfilepath=Path("acc_ids_boxseqs.txt")
    idlist=readaccids(accidfilepath)
    testrecs=batch_get(idlist)

        
