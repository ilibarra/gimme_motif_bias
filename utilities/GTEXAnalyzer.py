'''
Created on 3/11/2018

DESCRIPTION

@author: ignacio
'''
from .path_functions import *
from .DataFrameAnalyzer import DataFrameAnalyzer
from .FastaAnalyzer import FastaAnalyzer
from .SequenceMethods import SequenceMethods
import numpy as np
from numpy import average
import pandas as pd
from .MyGeneAnalyzer import MyGeneAnalyzer

class GTEXAnalyzer():
    @staticmethod
    def get_gene_tss():
        all_ids_path = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/1000_genes_gtex_tissues/all_ids.tsv.gz'
        gtex_dir = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/1000_genes_gtex_tissues'
        if not exists(all_ids_path):
            ids = set()
            for f in listdir(gtex_dir):
                p = join(gtex_dir, f)
                df = DataFrameAnalyzer.read_tsv_gz(p)
                for s in df['Gene stable ID']:
                    ids.add(s)
            print(len(ids))
            df = pd.DataFrame(list(ids), columns=['ensembl'])
            print(df.head())
            DataFrameAnalyzer.to_tsv(df, '/tmp/input_rscript.tsv')
        return DataFrameAnalyzer.read_tsv_gz(all_ids_path)

    @staticmethod
    def get_sg_scores():
        p = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/zscores_by_motif_bkp_cisbp_GTEx.tsv.gz'
        df = DataFrameAnalyzer.read_tsv_gz(p)
        df['sg'] = np.sqrt(((df['median.z.score.motif']) ** 2) + (df['z.score.expr'] ** 2))
        return df

    @staticmethod
    def get_tss_fa_path():
        tss_df = GTEXAnalyzer.get_gene_tss()
        tss_df['promoter.start'] = np.where(tss_df['strand'] == 1, tss_df['transcription_start_site'] - 1000,
                                            tss_df['transcription_start_site'])
        tss_df['promoter.end'] = tss_df['promoter.start'] + 1000

        print(tss_df.shape)
        tss_df = tss_df.drop_duplicates('ensembl_gene_id')
        tss_df['chromosome_name'] = "chr" + tss_df['chromosome_name']
        fa_path = '../../data/sequences_gtex_human_tss/all.fa'
        if not exists(fa_path):
            bed_path = fa_path.replace(".fa", '.bed')
            DataFrameAnalyzer.to_tsv(tss_df[['chromosome_name', 'promoter.start', 'promoter.end', 'ensembl_gene_id']],
                                     bed_path, header=None)
            FastaAnalyzer.convert_bed_to_fasta(bed_path, fa_path, genome='hg19')
            tss_df = SequenceMethods.parse_range2coordinate(tss_df,
                                                            ['chromosome_name', 'promoter.start', 'promoter.end'])
            ensembl_by_range = DataFrameAnalyzer.get_dict(tss_df, 'range', 'ensembl_gene_id')
            headers, seqs = list(zip(*FastaAnalyzer.get_fastas(fa_path, uppercase=True)))
            headers = [ensembl_by_range[h] for h in headers]
            print(len(headers), headers[:5])
            FastaAnalyzer.write_fasta_from_sequences([[h, s] for h, s in zip(headers, seqs)], fa_path)
        return fa_path

    @staticmethod
    def get_medians_by_gene():
        medians_bkp = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/gtex_tpm_by_tissue.tsv.gz'
        medians_bkp = "../../data/gtex_tpm_by_tissue.tsv.gz"
        if not exists(medians_bkp):
            gtex_dir = '/g/scb2/zaugg/rio/data/GTEx/v7'
            gtex_path = join(gtex_dir, 'GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz')

            samples = DataFrameAnalyzer.read_tsv(join(gtex_dir,
                                                      'GTEx_v7_Annotations_SampleAttributesDS.txt'))

            smts_by_id = pd.Series(samples['SMTSD'].values,
                                   index=samples['SAMPID'].values).to_dict()

            df = DataFrameAnalyzer.read_tsv_gz(gtex_path,
                                               sep='\t', skiprows=2)
            df.columns = [c + ":" + smts_by_id[c] if c in smts_by_id else c for c in
                          df.columns]

            # for c in df.columns[2:]:
            # print c, c.split(':')[1]
            tissues = {c.split(':')[1] for c in df.columns[2:]}
            print(tissues)
            print(df.head())

            all = []

            for t in tissues:
                sel = pd.concat([df[['Name', 'Description']],
                                 df[[c for c in df.columns[2:] if c.split(":")[1] == t]]],
                                axis=1)
                print(t, sel.shape)
                sel['median.' + t] = sel[sel.columns[2:]].median(axis=1)
                sel.index = sel['Name']
                all.append(sel)

            print(all[0].head())

            medians_df = pd.concat([sel[[c for c in sel.columns if 'median.' in c]]
                                    for sel in all], axis=1)
            medians_df.index = all[0].index

            print(medians_df.head())
            DataFrameAnalyzer.to_tsv_gz(medians_df, medians_bkp, index=True)

        medians_df = DataFrameAnalyzer.read_tsv_gz(medians_bkp, index_col=0)

        medians_df = DataFrameAnalyzer.read_tsv_gz(medians_bkp, index_col=0)
        return medians_df


    @staticmethod
    def get_zscores_by_gene(overwrite=False):
        zscores_bkp = "/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/zscores_gtex.tsv.gz"
        print('getting GTEx Z-scores')
        if not exists(zscores_bkp) or overwrite:
            print('calculating medians per gene...')
            medians_df = GTEXAnalyzer.get_medians_by_gene()
            zscores = []
            for ri, r in medians_df.iterrows():
                mu = average(r.values)
                std = np.std(r.values)
                next_zscores = (r.values - mu) / std
                zscores.append(next_zscores)

            zscores = pd.DataFrame(zscores, columns=medians_df.columns, index=medians_df.index)
            DataFrameAnalyzer.to_tsv_gz(zscores, zscores_bkp, index=True)

        print('reading')
        print(zscores_bkp)
        return DataFrameAnalyzer.read_tsv_gz(zscores_bkp, index_col=0)

    @staticmethod
    def get_genes_by_tissue(return_all=False):
        genes_by_tissue = '../../data/gene_subsets_five_tissues_1000_genes.tsv'
        genes_df = DataFrameAnalyzer.read_tsv(genes_by_tissue)
        genes_df = {k: set(genes_df[k]) for k in genes_df.columns}

        # add custom names from Butte's data
        butte_dir = '../../data/gene_sets_custom'
        for f in listdir(butte_dir):
            p = join(butte_dir, f)
            ids = {s.strip() for s in open(p)}
            # print f, len(ids), ids
            genes_df[f.replace(".txt", '')] = ids

        # load IDs for Pancreas and Heart from GTEx directly
        if not return_all:
            for tissue_id in ['Pancreas', 'Heart - Left Ventricle', 'Muscle - Skeletal']:
                p = '../../data/1000_genes_gtex_tissues/$1.tsv.gz'
                next_df = DataFrameAnalyzer.read_tsv_gz(p.replace("$1", tissue_id))
                genes_df[tissue_id.replace('Pancreas', 'pancreas').replace("Heart - Left Ventricle", 'Heart')] = set(
                    next_df['mouse.ensembl'])
        else:
            for f in listdir('../../data/1000_genes_gtex_tissues/'):
                if 'all_ids.tsv.gz' in f:
                    continue
                p = join('../../data/1000_genes_gtex_tissues/', f)
                next_df = DataFrameAnalyzer.read_tsv_gz(p)
                genes_df[f.replace(".tsv.gz", '')] = set(next_df['mouse.ensembl'])

        return genes_df

    @staticmethod
    def get_updated_ensembl_identifiers(**kwargs):
        ensembl_mygene_by_symbol_bkp = join(kwargs.get('datadir',
                                                       '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data'),
                                                       'ensembl_mygene_by_symbol.pkl')
        bkp_names_gtex = join(kwargs.get('datadir',
                                         '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data'),
                                         'gtex_ensembl_by_name.tsv.gz')
        df = DataFrameAnalyzer.read_tsv_gz(bkp_names_gtex)
        # convert old to new IDs and remove duplicates
        df['ensembl.gtex'] = df['Name'].str.split(".").str[0]
        df = df.drop_duplicates('ensembl.gtex')


        if not exists(ensembl_mygene_by_symbol_bkp):
            # convert ensembl IDs to the newest versions according to what MyGeneAnalyzer maps
            ensembl_mygene_by_symbol = MyGeneAnalyzer.get_ensembl_by_symbol(df['Description'])
            DataFrameAnalyzer.to_pickle(ensembl_mygene_by_symbol, ensembl_mygene_by_symbol_bkp)

        ensembl_mygene_by_symbol = DataFrameAnalyzer.read_pickle(ensembl_mygene_by_symbol_bkp)
        df['ensembl.mygene'] = df['Description'].map(ensembl_mygene_by_symbol)

        df = df.append(pd.DataFrame([['ENSG00000273439.1', 'ZNF8', 'ENSG00000273439', 'ENSG00000278129']],
                                    columns=['Name', 'Description', 'ensembl.gtex', 'ensembl.mygene'])).reset_index(drop=True)
        # df = df.drop_duplicates('ensembl.mygene')
        return df

