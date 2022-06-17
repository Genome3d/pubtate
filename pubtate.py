#! /usr/bin/env python

import pandas as pd
import istarmap
import argparse
import requests
import sys
import os
import time
import multiprocessing as mp
from tqdm import tqdm
from itertools import repeat
import json

import logger
import pubtator

def write_results(df, output_fp, logger):
    out_dir = os.path.dirname(output_fp)
    os.makedirs(out_dir, exist_ok=True)
    #logger.write('Writing output...')
    df.to_csv(output_fp, sep='\t', index=False)


def parse_genes(genes_args, logger):
    logger.write('Parsing gene input...')
    df = pd.DataFrame()
    if os.path.isfile(genes_args[0]):
        df = pd.read_csv(genes_args[0], sep='\t')
        if not 'gene' in df.columns:
            sys.exit('Input file has no gene column.')
        df = df[df.gene.str.lower() != 'gene']
        df = df[['gene']].drop_duplicates()
    else:
        df = pd.DataFrame({'gene': [i.upper() for i in genes_args]})
        df = df[['gene']].drop_duplicates()
    return df

def join_path(*args):
    fp = ''
    for arg in args:
        fp = os.path.join(fp, arg)
    return fp
            
    
def read_gene2pub(fp, logger):
    logger.write('Loading gene PubMed entries...')
    # TODO: Make SQL db.
    return pd.read_csv(fp, compression='gzip', sep = '\t', header=None,
                       names=['pubmed_id', 'type', 'concept_id', 'mentions', 'resource'])

def get_gene_pubmed(genes, gene2pub, logger):
    logger.write('Identifying articles related to gene input...')
    cond = (gene2pub.mentions.str.contains(r'\b' + r'|'.join(genes)
                                           .replace('|', r'\b|\b') + r'\b', regex=True, na=False))
    pubmed_ids = gene2pub[cond]
    pubmed_ids = (pubmed_ids[['pubmed_id', 'mentions']]
                  .drop_duplicates()
                  .rename(columns={'pubmed_id': 'pmid'})
                  )
    gene_pubmed = []
    for gene in genes:
        df = pubmed_ids[pubmed_ids.mentions.str.contains(r'\b' + gene + r'\b', regex=True, na=False)]
        if not df.empty:
            df.loc[:, 'gene'] = gene
            gene_pubmed.append(df)
    gene_pubmed = pd.concat(gene_pubmed)
    return gene_pubmed[['pmid', 'gene']].drop_duplicates()

def chunk_pmids(pmids, chunksize, out_dir):
    pmids = pmids.unique()
    chunks = [pmids[i:i+chunksize] for i in range(0, len(pmids), chunksize)]
    for i, chunk in enumerate(chunks):
        pd.Series(chunk).to_csv(os.path.join(out_dir, f'{i}.txt'), sep='\t',
                                index=False, header=None)

def get_mesh(df, mesh_type, logger):
    logger.write(f'Interpreting {mesh_type} annotations...')
    res = []
    df = df.query(f'type == "{mesh_type}"')
    df.loc[:, 'id'] = df['identifier'].str[5:]
    for meshid in df.id.unique():
        try:
            r = requests.post(f'https://id.nlm.nih.gov/mesh/{meshid}.json')
            term = json.loads(r.text)
            res.append([term['identifier'], term['label']['@value']])
        except:
            pass
    res = pd.DataFrame(res, columns=['id', 'label'])
    df = df.merge(res, how='left', on=['id'])
    return df


def parse_args():
    parser = argparse.ArgumentParser(
        description='Identify diseases associated with genes based on literature.')
    parser.add_argument(
        '-g', '--genes', required=True, nargs='+',
        help='''Space-separated gene names or a file with gene names in a column.''')
    parser.add_argument(
        '-o', '--output-dir', required=True,
        help='Directory to write results.')
    parser.add_argument(
        '--gene2pub', default=os.path.join(os.path.dirname(__file__),
                                           'data/gene2pubtatorcentral.gz'),
        help='''Path to Pubtator's gene2pubtatorcentral.gz''')
    parser.add_argument(
        '-f', '--format', default='biocjson',
        choices=['putator', 'biocxml', 'biocjson'],
        help='Return type format: putator(PubTator), biocxml (BioC-XML), and biocjson (JSON-XML).')
    parser.add_argument(
        '-b', '--bioconcept', nargs='+', default='all',
        choices=['all', 'gene', 'disease', 'chemical',
                 'species', 'proteinmutation', 'dnamutation', 'snp', 'cellline'],
        help='''Type of annotation of interest: gene, disease, chemical, species,
        proteinmutation, dnamutation, snp, and cellline. Default: all''')
    
    return parser.parse_args()

if __name__=='__main__':
    pd.options.mode.chained_assignment = None
    args = parse_args()

    start_time = time.time()
    os.makedirs(args.output_dir, exist_ok=True)
    pmid_dir = os.path.join(args.output_dir, 'pmids')
    os.makedirs(pmid_dir, exist_ok=True)
    global logger
    logger = logger.Logger(logfile=os.path.join(args.output_dir, 'pubtate.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    logger.write('\n')
    genes = parse_genes(args.genes, logger)
    gene2pub = read_gene2pub(args.gene2pub, logger)
    gene_pubmed = get_gene_pubmed(genes['gene'], gene2pub, logger)
    write_results(gene_pubmed, join_path(args.output_dir, 'gene_pubmed.txt'), logger)
    chunk_pmids(gene_pubmed['pmid'], 1000, pmid_dir)
    logger.write('Retrieving annotations in relevant articles...')
    annotations = []
    for fp in [fp for fp in os.listdir(pmid_dir) if fp.endswith('.txt')]:
        batch_fp = os.path.join(pmid_dir, fp)
        annotations.append(
            pubtator.fetch_annot(batch_fp, args.format, args.bioconcept))
    annotations = pd.concat(annotations)
    write_results(annotations, join_path(args.output_dir, 'annotations.txt'), logger)
    mesh = get_mesh(annotations, 'Disease', logger)
    write_results(mesh, join_path(args.output_dir, 'mesh.txt'), logger)
    '''
    pubscore = (
        mesh
        .query('type == "Disease"')
        .merge(gene_pubmed, on=['pmid'], how='inner')
        .drop_duplicates(subset=['pmid', 'gene', 'label'])
        .groupby(['gene', 'label'])['pmid'].agg('nunique').reset_index()
        .assign(total = lambda x: x.groupby(['gene'])['pmid'].transform(sum),
                perc = lambda x: round(x['pmid']/x['total'], 2))
        .sort_values(by=['pmid'], ascending=False)
    )
    write_results(pubscore, join_path(args.output_dir, 'pubscore.txt'), logger)
    print(pubscore)
    '''
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
