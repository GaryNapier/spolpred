#!/usr/bin/env python

import argparse
import json
import os
from tqdm import tqdm
import csv
from python_scripts.utils import *


def main(args):

    # FILES

    # in
    tbprofiler_results_dir = args.tbprofiler_results_dir
    # out
    lineage_file = args.lineage_file

    # VARIABLES
    suffix = args.suffix
    files = os.listdir(tbprofiler_results_dir)
    samples = {x.replace(suffix, '') for x in files}

    # Append the directory
    files = [tbprofiler_results_dir + file for file in files]

    # Loop over tbprofiler
    lineage_dict = {}
    # for file in tqdm([tbprofiler_results_dir + 'SRR5709951.results.json']):
    for file in tqdm(files):
        data = json.load(open(file))

        # Skip mixed samps
        if ';' in data['sublin']: continue

        # Pull all the lineage data
        lineage = data['lineage']
        lowest_level = lineage[len(lineage)-1] # STUPID!

        lineage_dict[data['id']] = {
            'id': data['id'],
            'main_lin': data['main_lin'],
            'sublin': data['sublin'],
            'family': lowest_level['family'],
            'spoligotype': lowest_level['spoligotype'],
            'rd': lowest_level['rd']
        }


    # Write out data

    with open(lineage_file, 'w') as f:
            writer = csv.DictWriter(f, fieldnames = list(get_embedded_keys(lineage_dict)))
            writer.writeheader()
            for row in lineage_dict:
                writer.writerow(lineage_dict[row])


parser = argparse.ArgumentParser(description='get lineage data from all TB-profiler samples',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# in
parser.add_argument('--tbprofiler-results-dir', type = str, help = 'directory of TB-profiler json files')
parser.add_argument('--suffix', default = '.results.json', type = str, help = 'suffix of json files in tbprofiler_results_dir')

# out
parser.add_argument('--lineage-file', default = '', type = str, help = 'name of output file - outputs metadata of lineages for all samples in TB-profiler')


parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)