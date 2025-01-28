#!/usr/env/bin python3
"""
get_ncbi_taxonomy.py

Uses Entrez to retrieve taxonomic info
Takes in list of taxids or common names or scientific names
"""

from pathlib import Path
from Bio import Entrez
import textwrap
import argparse

Entrez.email = "yhsieh@alaska.edu"

def readtaxids(infile):
    """Reads in file with taxonomic IDs, returns a list"""
    taxids = []
    with open(infile, 'r') as f:
        for line in f:
            taxids.append(line.strip())
    return taxids

def batch_get(listoftaxids):
    """Submits Entrez batch request, returns handle"""
    request = Entrez.epost("taxonomy", id=",".join(listoftaxids))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="taxonomy", retmode="xml", webenv=webEnv, query_key=queryKey)
    return handle

def get_lineage(entrez_record):
    key1 = 'LineageEx'
    key2 = 'ScientificName'
    position = 1 # to get kingdom
    return entrez_record[key1][position][key2]

def get_info_pair(entrez_record, mode=1):
    key1 = 'TaxId'
    key2 = 'ScientificName'
    if mode:
        return (entrez_record[key1], entrez_record[key2])
    else:
        kingdom_name = get_lineage(entrez_record)
        return (entrez_record[key1], f'{entrez_record[key2]}({kingdom_name})')

def get_info_tuples(handleobj, mode=1):
    """Takes in handle, iterates find scientific names, returns of tuples"""
    records = Entrez.parse(handleobj)
    taxinfolist = []
    for record in records:
        taxinfolist.append(get_info_pair(record, mode))
    return taxinfolist

def write_info_tuples_to_file(taxinfolist, outfile="taxid_pair.txt"):
    """Writes out taxinfo tuples to a file"""
    with open(outfile, 'w+') as outf:
        for taxid, name in taxinfolist:
            if name:
                name = '_'.join(name.split())
                outf.write(f'{name}^{taxid}\n')
            else:
                outf.write(taxid)

def get_taxinfofile_from_entrez(handleobj, outputpath):
    """Takes in a list of taxids and saves the taxid-name tuple
    into a file"""
    taxinfolist = get_info_tuples(handleobj)
    print(f"Writing to {outputpath}")
    write_info_tuples_to_file(taxinfolist, outputpath)
    print('---> written')

    
def argumentsparser():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] taxidfile outpath",
                                     description="Retrieves taxid and scientific name from NCBI",)
    parser.add_argument("taxidfile", help="path to taxonomic ID file", type=str)
    parser.add_argument("outpath", help="path to output file", type=str)
    args = parser.parse_args()
    return args

if __name__=="__main__":
    
    args = argumentsparser()
    taxidfile = Path(args.taxidfile)
    outpath = Path(args.outpath)

    taxidlist = readtaxids(taxidfile)
    handle = batch_get(taxidlist)
    get_taxinfofile_from_entrez(handle, outpath)
