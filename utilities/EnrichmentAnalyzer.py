'''
Created on Nov 21, 2016

@author: ignacio
'''

from scipy.stats import fisher_exact
import sys
from .path_functions import *
from itertools import product
from scipy.stats import ttest_ind
from .DataFrameAnalyzer import DataFrameAnalyzer
import numpy as np
import pandas as pd
from numpy import median

class EnrichmentAnalyzer(object):
    @staticmethod
    def get_motif_enrichments_by_pairwise_grouping(genes_by_grouping, motif_hits, label=None, column_gene='ensembl',
                                                   stopat=None, key_c1=None):
        enrichments = []

        grouping_names = list(genes_by_grouping.keys())
        # print genes_df.keys()
        queries = [[c1, c2] for c1, c2 in product(grouping_names, repeat=2)]
        print('# queries', len(queries))

        motif_hits_by_ont = {ont: motif_hits[motif_hits[column_gene].isin(genes_by_grouping[ont])] for ont in grouping_names}
        for qi, q in enumerate(queries):

            c1, c2 = q
            if qi % 100 == 0:
                print(qi, 'pairwise queries out of', len(queries))

            if c1 == c2:
                t = [c1, c2, None, None, None, None, label, None, None, 1, 0.0, 1.0, 0, 1.0]
                enrichments.append(t)
                continue
            if key_c1 is not None and not key_c1 in c1:
                continue


            genes2 = set(genes_by_grouping[c2])
            genes1 = set(genes_by_grouping[c1])  # - genes2

            # genes1_mus = {engmus_by_enghuman[g] for g in genes1 if g in engmus_by_enghuman}
            # genes2_mus = {engmus_by_enghuman[g] for g in genes2 if g in engmus_by_enghuman}
            # print grp.head()

            # print grp.shape[0]
            # print grp.head()
            # check in mouse and human genes, respectively
            # sel1 = motif_hits[motif_hits['ensembl'].isin(genes1_mus)]
            # if sel1.shape[0] == 0:
            sel1 = motif_hits_by_ont[c1]

            # sel2 = motif_hits[motif_hits['ensembl'].isin(genes2_mus)]
            # if sel2.shape[0] == 0:
            sel2 = motif_hits_by_ont[c2]

            # print c1, len(genes1), sel1.shape[0], c2, len(genes2), sel2.shape[0]
            hits1_by_gene = {}
            for ensg_id, grp2 in sel1.groupby(column_gene):
                hits1_by_gene[ensg_id] = grp2.shape[0]
            for gene_id in genes1:
                if not gene_id in hits1_by_gene:
                    hits1_by_gene[gene_id] = 0
            hits2_by_gene = {}
            for ensg_id, grp2 in sel2.groupby(column_gene):
                hits2_by_gene[ensg_id] = grp2.shape[0]
            for gene_id in genes1:
                if not gene_id in hits2_by_gene:
                    hits2_by_gene[gene_id] = 0

            t_stat, pval_ttest = ttest_ind(list(hits1_by_gene.values()),
                                           list(hits2_by_gene.values()))

            scores1, scores2 = sel1['score'], sel2['score']
            a, b = len(set(sel1[column_gene])), len(genes1)
            c, d = len(set(sel2[column_gene])), len(genes2)

            odd_ratio, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
            log2FC = np.log2(odd_ratio)

            t = [c1, c2, a, len(genes1), c, len(genes2), label, \
                 median(scores1), median(scores2), odd_ratio,
                 -1 if log2FC == float("-inf") else (1 if log2FC == float("inf") else log2FC), pval, t_stat,
                 pval_ttest]

            # print odd_ratio, pval
            enrichments.append(t)

            if stopat is not None and len(enrichments) >= stopat:
                break
            # print enrichments[-1]

        enrichments = pd.DataFrame(enrichments, columns=['a', 'b', 'n.a.w.motif', 'n.a', 'n.b.w.motif', 'n.b',
                                                         'motif', 'median.a', 'median.b', 'odds.ratio', 'log2FC',
                                                         'p.val',
                                                         'ttest.ind', 'p.val.ttest'])
        return enrichments

if __name__ == '__main__':
    from os.path import exists
    # test our script by loading some examples
    path = '../armando_snps/data/BRD2-CTCF.bed.gz'
    print(exists(path))
    enrichment_analyzer = EnrichmentAnalyzer()
    enrichment_analyzer.get_enrichment_from_bed(path, path, path)
    sys.exit()
    matrix = [[1000, 100], [1000, 120]]
    print(enrichment_analyzer.get_enrichment_from_matrix(matrix))