'''
Created on 1/31/2018, 2018

@author: Ignacio Ibarra Del Rio

Description:
'''

import mygene
mg = mygene.MyGeneInfo()

class MyGeneAnalyzer:

    @staticmethod
    def get_ensembl_by_symbol(names, species='human', returnall=False, fields='ensembl.gene'):
        # get the IDs and then convert them to ENSG
        print('querying genes...')
        if len(names) == 0:
            print('list emtpy. Return none.')
            return {}
        out = mg.querymany(list(names), scopes='symbol', fields=fields, species=species, as_dataframe=True,
                           returnall=returnall)
        if not 'ensembl.gene' in out:
            print('nothing found..')
            return None
        return {ri: (r['ensembl']['gene'] if isinstance(r['ensembl'], dict) else
                     (r['ensembl'][0]['gene'] if isinstance(r['ensembl'], list) else (r['ensembl.gene'] if 'ensembl.gene' in out else np.nan)))
                for ri, r in out.iterrows()}

    @staticmethod
    def get_symbol_by_ensembl_series(series, species='human', **kwargs):
        symbol_by_ensembl = MyGeneAnalyzer.get_symbol_by_ensembl(set(series), species=species, **kwargs)
        return series.map(symbol_by_ensembl)

    @staticmethod
    def get_gene_tss(species='hg19', extend=0, **kwargs):
        from lib.TerminalRepressors import TerminalRepressors
        return TerminalRepressors.get_tss(species, extend, **kwargs)

    @staticmethod
    def get_ensemblgene_by_ensemblprotein(ensembl_protein, species='human'):
        # get the IDs and then convert them to ENSG
        print('querying protein IDs...')
        out = mg.querymany(list(ensembl_protein), scopes='ensembl.protein', fields='ensembl.gene', species=species,
                           as_dataframe=True)

        if not 'ensembl' in out:
            print('nothing found..')
            return None
        return {ri: (r['ensembl']['gene'] if isinstance(r['ensembl'], dict) else
                     r['ensembl'][0]['gene'] if isinstance(r['ensembl'], list) else r['ensembl'])
                for ri, r in out.iterrows()}

    @staticmethod
    def get_ensemblprotein_by_ensemblgene(ensembl_gene, species='human'):
        # get the IDs and then convert them to ENSG
        print('querying protein IDs...')
        out = mg.querymany(list(ensembl_gene), fields='ensembl.protein', scopes='ensembl.gene', species=species,
                           as_dataframe=True)
        if not 'ensembl' in out:
            print('nothing found..')
            return None
        return {ri: (r['ensembl']['protein'] if isinstance(r['ensembl'], dict) else
                     r['ensembl'][0]['protein'] if isinstance(r['ensembl'], list) else r['ensembl'])
                for ri, r in out.iterrows()}

    @staticmethod
    def get_ensembl_by_refseq(names, species='human'):
        # get the IDs and then convert them to ENSG
        print('querying genes...')
        out = mg.querymany(list(names), scopes='refseq', fields='ensembl.gene', species=species, as_dataframe=True)

        if not 'ensembl' in out:
            print('nothing found..')
            return None
        return {ri: (r['ensembl']['gene'] if isinstance(r['ensembl'], dict) else
                     r['ensembl'][0]['gene'] if isinstance(r['ensembl'], list) else r['ensembl'])
                for ri, r in out.iterrows()}


    @staticmethod
    def get_refseq_by_ensembl(names, species='human'):
        # get the IDs and then convert them to ENSG
        print(('querying genes with species (%s)...' % species))
        out = mg.querymany(list(names), fields='refseq', scopes='ensembl.gene', species=species, as_dataframe=True)

        print((out.head()))
        if not 'refseq.genomic' in out:
            print('nothing found..')
            return None

        return {ri: (r['refseq.genomic'][0] if (not isinstance(r['refseq.genomic'], float) and isinstance(r['refseq.genomic'], list)) else
                     r['refseq.genomic'] if not isinstance(r['refseq.genomic'], float) else None)
                for ri, r in out.iterrows()}

    @staticmethod
    def get_name_by_ensembl(ensembl, species='human'):
        # get the IDs and then convert them to ENSG
        print('querying genes...')
        out = mg.querymany(list(ensembl), scopes='ensembl.gene',
                           fields='refseq,name', species=species, as_dataframe=True)

        print((out.head()))
        return {ri: r['name'] for ri, r in out.iterrows()}

    @staticmethod
    def get_symbol_by_ensembl(ensembl, species='human'):
        # get the IDs and then convert them to ENSG
        print('querying genes...')
        out = mg.querymany(list(ensembl), scopes='ensembl.gene', fields='symbol',
                           species=species, as_dataframe=True)
        print('removing repeated genes...')
        return DataFrameAnalyzer.get_dict(out, None, 'symbol')

    @staticmethod
    def get_ensemblgene_by_ensemblprotein_obsolete(ensembl, species='mouse'):
        # get the IDs and then convert them to ENSG
        print('querying genes...')
        out = mg.querymany(list(ensembl), scopes='ensembl.protein', fields='ensembl.gene',
                           species=species, as_dataframe=True)

        return {ri: (r['ensembl']['gene'] if isinstance(r['ensembl'], dict) else None) for ri, r in out.iterrows()}


    @staticmethod
    def get_gene_ensg_by_entrez(entrez, species='human'):
        # get the IDs and then convert them to ENSG
        print('querying genes...')
        out = mg.getgenes(list(entrez), scopes='entrezgene',
                          fields='ensembl.gene', as_dataframe=True,returnall=True,
                          species=species)

        return {vi: v['ensembl']['gene'] if isinstance(v['ensembl'], dict) else
                (v['ensembl'][0]['gene'] if not isinstance(v['ensembl'], float) else None)
                for vi, v in out.iterrows()}


    @staticmethod
    def get_ensembl_by_uniprot(uniprot, species='human', returnall=False):
        out = mg.querymany(uniprot, scopes = 'uniprot',
                           fields = 'ensembl.gene', species = species,
                           as_dataframe = True, returnall=returnall)
        d = {}
        for k, r in out.iterrows():
            if isinstance(r['ensembl'], float):
                d[k] = None
            elif isinstance(r['ensembl'], list):
                d[k] = r['ensembl'][0]['gene']
            else:
                d[k] = r['ensembl']['gene']
        return d

    @staticmethod
    def get_uniprot_by_ensembl(ensembl, species='human', entry_type='Swiss-Prot'):
        out = mg.querymany(ensembl, scopes = 'ensembl.gene',
                           fields = 'uniprot', species = species,
                           as_dataframe = True)
        d = {}
        for k, r in out.iterrows():
            if isinstance(r['uniprot'], float):
                d[k] = None
            else:
                if not entry_type in r['uniprot']:
                    next = r['uniprot']['TrEMBL'] if not isinstance(r['uniprot']['TrEMBL'], list) else r['uniprot']['TrEMBL'][0]
                    d[k] = next
                else:
                    next = r['uniprot'][entry_type] if not isinstance(r['uniprot'][entry_type], list) else r['uniprot'][entry_type][0]
                    d[k] = next
        return d


    @staticmethod
    def get_homologs_mouse_from_human(ensg_human_name):
        # mouse gene names

        query_genes = set(ensg_human_name)
        orthologs_tsv = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/human_mm10_homologs.tsv'
        human_orthologs = pd.read_csv(orthologs_tsv,
                                      sep='\t', header=None)
        genes_mm10 = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/gene_names.txt'
        ensmusg_by_gene = {t[2].upper(): t[0] for t in
                           [r.split(",") for r in open(genes_mm10)]}
        human_orthologs.columns = ['human', 'transcript.id', 'mouse']
        ensg_mouse = human_orthologs[[k in query_genes for k in human_orthologs['human']]]['mouse']

        d = DataFrameAnalyzer.get_dict(human_orthologs, 'human', 'mouse')
        return {k : d[k] if k in d else None for k in query_genes}

    @staticmethod
    def get_homologs_human_from_mouse(ensg_mouse_name):
        '''
        Return a dictionary mapping each mouse gene to its respective human homolog, or None
        :param ensg_mouse_name:
        :return:
        '''

        query_genes = set(ensg_mouse_name)
        orthologs_tsv = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/human_mm10_homologs.tsv'
        human_homologs = pd.read_csv(orthologs_tsv , sep='\t')
        genes_mm10 = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/gene_names.txt'

        human_homologs.columns = ['human', 'transcript.id', 'mouse']
        d = DataFrameAnalyzer.get_dict(human_homologs, 'mouse', 'human')

        return {k: d[k] if k in d else None for k in query_genes}