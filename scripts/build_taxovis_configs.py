#!/usr/env/bin python3
"""
build_taxovis_configs.py

Takes in a CSV of format: accid,protein,species
Makes the following config files for ProTaxoVis, taxovis.py:
    proteinlist.txt
    limits.txt
    heatmap_config.txt
"""

from pathlib import Path
import argparse
import pandas as pd
import re

def read_in_csv(pathtocsv):
    columns = ['accid', 'protein', 'species']
    return pd.read_csv(pathtocsv, names=columns)

def remove_whites(somestr):
    return re.sub(r'\s+', '_', somestr)

def make_protein_list(dataframe):
    #formatrow_names = lambda row: f'{row["protein"].replace(" ","_")}_\
     #                         {row["species"].replace(" ","_")}'
    formatrow_names = lambda row: f'{remove_whites(row["protein"])}_{remove_whites(row["species"])}'
    formatrow_fastas = lambda row: f'{row["accid"]}.fasta'
    formattednames = dataframe.apply(formatrow_names, axis=1)
    formattedfasta = dataframe.apply(formatrow_fastas, axis=1)
    return pd.concat([formattednames, formattedfasta], axis=1)

def write_protein_list(dataframe):
    dataframe.to_csv('proteinlist.txt', sep='\t', index=False)

def make_limits():
    pass

def make_heatmap_config():
    pass


if __name__ == "__main__":
    #parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] csvfile",
    #                                 description="makes taxovis config files",)
    #parser.add_argument("csvfile", help="path to csv file format: accids,protein,species", type=str)
    #args = parser.parse_args()
    df = read_in_csv("box_table.txt")
    df1 = make_protein_list(df)
    write_protein_list(df1)
