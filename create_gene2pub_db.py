import requests
import io
import json
import sys
import os
#from bioc import biocjson
import pandas as pd
from sqlalchemy import create_engine
import argparse
import time

import logger
    

def create_db(db_url, fp):
    with create_engine('sqlite:///{}'.format(db_url), echo=False).connect() as db:
        names = ['pmid', 'type', 'concept_id', 'mentions', 'resource']
        i = 0
        mode = 'replace'
        for batch in pd.read_csv(fp, sep='\t', compression='gzip', chunksize=20000,
                                 header=None, names=names):
            batch = batch[~batch['mentions'].isna()]
            batch['mentions'] = batch['mentions'].astype(str)
            batch['idx'] = batch['pmid'].astype(str) + '__' + batch['mentions'].astype(str)
            df = (pd.DataFrame(batch.mentions.str.split('|').values.tolist(), index=batch['idx'])
                  .stack()
                  .reset_index()[[0, 'idx']]
                  .rename(columns={0: 'gene'})
            )
            df = (batch[['pmid', 'idx']]
                  .merge(df, on=['idx'], how='inner')[['pmid', 'gene']]
                  .drop_duplicates()
            )
            df['gene'] = df['gene'].str.strip()
            print(df)
            if i > 0:
                mode = 'append' 
            df.to_sql('lookup', db, if_exists=mode, index=False)
            i += 1
        db.execute('''CREATE INDEX idx_pmid on lookup (pmid)''')
        db.execute('''CREATE INDEX idx_gene on lookup (gene)''')


    
def parse_args():
    parser = argparse.ArgumentParser(
        description='Retrieve Pubtator annotations of articles.')
    parser.add_argument(
        '--gene2pub',
        default=os.path.join(os.path.dirname(__file__),'data/gene2pubtatorcentral.gz'),
        help='gene2pubtatorcentral.gz from PubTator sit.')
    parser.add_argument(
        '--db-url', default=os.path.join(os.path.dirname(__file__), 'data/gene2pubtatorcentral.db'),
        help='''Filepath to SQL DB''')

    return parser.parse_args()
    
if __name__ == "__main__":
    pd.options.mode.chained_assignment = None
    args = parse_args()
    start_time = time.time()
    out_dir = os.path.dirname(args.db_url)
    os.makedirs(out_dir, exist_ok=True)
    create_db(args.db_url, args.gene2pub)
    '''global logger
    logger = logger.Logger(logfile=os.path.join(out_dir, 'mesh.log'))
    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    logger.write('\n')

    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
    '''
