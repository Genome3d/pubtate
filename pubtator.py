import requests
import io
import json
import sys
import os
#from bioc import biocjson
import pandas as pd
import argparse
import time

import logger

def fetch_annot(input_fp, form, bioconcept):
    inp_json = {}
    with io.open(input_fp,'r',encoding="utf-8") as file_input:
        inp_json = {"pmids": [pmid.strip() for pmid in file_input.readlines()]}

    if bioconcept != 'all': 
        inp_json["concepts"] = bioconcept
    url = 'https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/'
    r = requests.post(url + form, json = inp_json)
    if r.status_code != 200 :
            print("[Error]: HTTP code "+ str(r.status_code))
    else:
        df = []
        res = r.text
        for i in res.split('\n'):
            try:
                pub = json.loads(i.encode('utf-8'))
                for passages in pub['passages']:
                    for annotation in passages['annotations']:
                        annot_type = annotation['infons']['type']
                        if annot_type in ['Disease', 'Gene']:
                            df.append([pub['pmid'],
                                       annot_type,
                                       annotation['infons']['identifier'],
                                       annotation['text']])
            except:
                pass
        df = pd.DataFrame(df, columns = ['pmid', 'type', 'identifier', 'text'])
        return df

    
def parse_args():
    parser = argparse.ArgumentParser(
        description='Retrieve Pubtator annotations of articles.')
    parser.add_argument(
        '-i', '--in-dir', required=True,
        help='Directory containing files each having PUBMED IDs of interest.')
    parser.add_argument(
        '-f', '--format', default='biocjson',
        choices=['pubtator', 'biocxml', 'biocjson'],
        help='Return type format: pubtator(PubTator), biocxml (BioC-XML), and biocjson (JSON-XML).')
    parser.add_argument(
        '-b', '--bioconcept', nargs='+', default='all',
        choices=['all', 'gene', 'disease', 'chemical',
                 'species', 'proteinmutation', 'dnamutation', 'snp', 'cellline'],
        help='''Type of annotation of interest: gene, disease, chemical, species, 
        proteinmutation, dnamutation, snp, and cellline. Default: all''')
    parser.add_argument(
        '-o', '--output', required=True, help='Filepath to write results.') 

    return parser.parse_args()
    
if __name__ == "__main__":
    pd.options.mode.chained_assignment = None
    args = parse_args()
    start_time = time.time()
    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)
    global logger
    logger = logger.Logger(logfile=os.path.join(out_dir, 'pubtator.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    logger.write('\n')
    df = []
    for fp in sorted([fp for fp in os.listdir(args.in_dir) if fp.endswith('.txt')]):
        print(fp)
        input_fp = os.path.join(args.in_dir, fp)
        df.append(fetch_annot(input_fp, args.format, args.bioconcept))
    df = pd.concat(df)
    df.to_csv(args.output, sep='\t', index=False)

    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
