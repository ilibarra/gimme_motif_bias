'''
Created on 7/22/2018

DESCRIPTION

@author: ignacio
'''

from .DataFrameAnalyzer import DataFrameAnalyzer
from .path_functions import *
import numpy as np

class TerminalRepressors():
    @staticmethod
    def get_tss(species='human', extend=0, **kwargs):
        df = None
        assert species in {'human', 'mouse'}
        datadir = kwargs.get('datadir',
                             '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data_2018-07-09_00-00-00/data')
        if species == 'mouse':
            df = DataFrameAnalyzer.read_tsv(join(datadir, 'mm10_refseq_n_ensembl_n_pos.tsv'))
            df['start'] = np.where(df['pos'] == 1, df['tss'].astype(int), (df['tss'] - extend).astype(int))
            df['end'] = np.where(df['pos'] == 1, (df['tss'] + extend).astype(int), df['tss'].astype(int))
        elif species == 'human':
            bed_tss_hg19 = None
            bed_tss_hg19 = join(datadir, 'hg19_refseq_n_ensembl_n_pos.tsv')
            df = DataFrameAnalyzer.read_tsv(bed_tss_hg19)
            df['start'] = df['tss'].astype(int) - extend
            df['end'] = (df['start'] + extend * 2).astype(int)
        df = df.sort_values(['chr', 'start'], ascending=[True, True])
        return df
