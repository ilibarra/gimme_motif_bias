'''
Created on 

DESCRIPTION

@author: ignacio
'''

import multiprocessing
from .path_functions import *
from multiprocessing import Manager
from .DataFrameAnalyzer import DataFrameAnalyzer
from .func import *

class TabulaMuris:

    def __init__(self, method=None, tabula_muris_rootdir=''):
        if method not in {'FACS', 'droplet'}:
            print('please specific a method (FACS/droplet)')
            assert 1 > 2

        if len(tabula_muris_rootdir) == 0:
            print('Tabula Muris directory is given (parm tabula_muris_rootdir). Some functions might not work')
        self.method = method
        self.basedir = tabula_muris_rootdir
        self.cpm_dir = join(tabula_muris_rootdir, '00_data_ingest/05_robj_to_dataframes')
        self.annot = self.get_annotation()

    def get_annotation(self):
        p = join(self.basedir, '00_data_ingest/18_global_annotation_csv/annotations_%s.csv' % self.method.lower())
        if not exists(p):
            p = join("input", 'annotations_%s.csv.gz' % self.method.lower())

        print(p)
        df = DataFrameAnalyzer.read_tsv_gz(p, sep=',')
        return df

    @staticmethod
    def get_sg_scores():
        p = '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data/full_results_tabula_muris.tsv.gz'
        df = DataFrameAnalyzer.read_tsv_gz(p)
        df['sg'] = np.sqrt(((df['z.score']) ** 2) + (df['z.score.all'] ** 2))
        return df

    def get_cpm(self, tissue=None, nrows=None):
        print('loading from method: %s' % self.method)
        res = []
        for f in listdir(self.cpm_dir):
            p = join(self.cpm_dir, f)
            if tissue != None and not tissue in f:
                continue
            # print self.method.lower() in f
            if not self.method.lower() in f:
                continue
            print('next tissue: %s' % f)
            print('n rows', nrows)
            print(exists(p), p)
            df = DataFrameAnalyzer.read_tsv_gz(p, sep='\t', index_col=0, nrows=nrows, keep_default_na=False)
            res.append(df)
        res = pd.concat(res)
        return res

    def get_raw_counts(self, tissue=None):
        print('loading from method: %s' % self.method)
        res = []
        for f in listdir(join(self.basedir, self.method)):
            p = join(self.basedir, self.method, f)
            if tissue != None and not tissue in f:
                print('skip tissue', f, '...')
                continue
            print('next tissue: %s' % f)
            df = DataFrameAnalyzer.read_tsv(p, sep=',', index_col=0)
            print(df.head())
            res.append(df)
        res = pd.concat(res)
        return res

    def get_tissue_names(self):
        return set(self.annot['tissue'])

    def get_cell_ontology_class_names(self):
        return set(self.annot['cell_ontology_class'])

    @staticmethod
    def get_known_reprogramming_protocols():
        p = '/home/rio/data/tabula-muris/23_tf_analysis/reprogramming/known_reprogramming_protocols.csv'
        df = DataFrameAnalyzer.read_tsv(p, sep=',')
        trans_df = TabulaMuris.get_translate_celltype_cellOntology_table()
        ont_by_celltype = DataFrameAnalyzer.get_dict(trans_df, 'CellType', 'Cell_Ontology_Class')
        df['ont'] = [{ont_by_celltype[k] for k in ont_by_celltype if k in celltype} for celltype in df['CellType']]
        df['ont'] = [str(list(k)[0]).split(";")[0] if len(k) > 0 else None for k in df['ont']]
        df['Genes'] = df['Genes'].str.upper().str.split(";")
        return df

    @staticmethod
    def get_translate_celltype_cellOntology_table():
        p = '/home/rio/data/tabula-muris/23_tf_analysis/reprogramming/translate_celltype_cellOntology.csv'
        trans_df = DataFrameAnalyzer.read_tsv(p, sep=',')
        return trans_df


    def get_cpm_by_tissue_and_ont(self):
        from lib.ThreadingUtils import ThreadingUtils
        import multiprocessing
        from multiprocessing import Manager
        manager = multiprocessing.Manager()
        cpm_by_tissue = manager.dict()
        cpm_by_tissue_and_ont = manager.dict()

        tm = TabulaMuris(method=self.method)
        tissue_names = tm.get_tissue_names()

        def load_dataframes(tiss):
            # load dataframes in parallel
            print(tiss)
            cpm = tm.get_cpm(tissue=tiss)
            cpm_by_tissue[tiss] = cpm
            cpm['k'] = cpm.index.to_series().map(DataFrameAnalyzer.get_dict(tm.annot, 'cell', 'cell_ontology_class'))
            # print cpm[['k']].head()
            for cell_ont_class, grp in cpm.groupby('k'):
                k = tiss + ";" + cell_ont_class
                if not k in list(cpm_by_tissue_and_ont.keys()):
                    print('adding %s, %s, %i' % (tiss, cell_ont_class, grp.shape[0]))
                    cpm_by_tissue_and_ont[k] = grp
            print(len(list(cpm_by_tissue.keys())))

        ThreadingUtils.run(load_dataframes, [[t] for t in list(tissue_names)], min(len(tissue_names), 12))

        print('getting genes...')
        genes = set()
        for tiss in list(cpm_by_tissue.keys()):
            genes = genes.union(set(cpm_by_tissue[tiss].columns))

        cpm_by_tissue_and_ont = {k: cpm_by_tissue_and_ont[k] for k in list(cpm_by_tissue_and_ont.keys())}
        return genes, cpm_by_tissue_and_ont

    def group_genes_by_ont(self, genes, cpm_by_tissue_and_ont):
        manager = multiprocessing.Manager()
        res_by_core = manager.dict()

        def get_avr_by_genes(pi, genes):
            next_res = []
            for gi, g in enumerate(genes):
                if gi % 100 == 0:
                    print(pi, gi, g, len(genes))
                if g == 'k':
                    continue
                for next_ont in {k.split(";")[1] for k in cpm_by_tissue_and_ont}:
                    next = [cpm_by_tissue_and_ont[k] for k in cpm_by_tissue_and_ont if
                            k.split(";")[1] == next_ont]
                    n_cells = sum([n.shape[0] for n in next])
                    expr_mean = np.mean(pd.concat([n[g] for n in next]))
                    expr_median = np.nanmedian(pd.concat([n[g] for n in next]))
                    next_res.append([next_ont, g, n_cells, expr_mean, expr_median])
            next_res = pd.DataFrame(next_res, columns=['ont', 'gene', 'n', 'avr.expr', 'median.expr'])
            res_by_core[pi] = next_res

        ThreadingUtils.run(get_avr_by_genes, [[qi, q] for qi, q in enumerate(ThreadingUtils.chunks(list(genes), 1000))], 10)
        res = pd.concat([res_by_core[k] for k in list(res_by_core.keys())])
        return res

    def get_ont_by_ncell(self, n_cells_cutoff=100):
        if not hasattr(self, 'expr'):
            self.expr = self.get_avr_expr_by_gene()
        expr = self.expr
        return expr[expr['n'] >= n_cells_cutoff]

    def get_zscores(self, ont_subset=None, n_cells_cutoff=None, get_raw=False,
                    by='median.expr'):
        '''
        Given a set of ontologies and a subset for cell counts, return all possible zscores for each gene
        :param ont_subset:
        :return:
        '''
        if not hasattr(self, 'expr'):
            self.expr = self.get_avr_expr_by_gene()

        sel = self.expr
        if n_cells_cutoff is not None:
            sel = sel[sel['n'] >= n_cells_cutoff]
        if ont_subset is not None:
            sel = sel[sel['ont'].isin(ont_subset)]
        hm = sel.pivot('ont', 'gene', by)
        means = hm.mean(axis=0)
        stddev = hm.std(axis=0)
        if get_raw:
            return hm
        else:
            zscores = (hm - means) / stddev
            return zscores

    def get_ont_group_names(self):
        return {'shortlist1', 'shortlist2', 'all'}

    def get_ont_definitions(self):
        df = pd.read_excel('/g/scb2/zaugg/rio/data/tabula-muris/00_data_ingest/shortlist_cell_onts_counts.xlsx')
        return {label: set(grp['ont']) for label, grp in df.groupby('group')}

    def get_zscores_by_shortlist(self, shortlist_name, by='median.expr', **kwargs):
        ont = self.get_ont_definitions()
        if not shortlist_name in ont:
            print('shortlist not found. Returning ALL')
        zscores = self.get_zscores(ont_subset=ont[shortlist_name] if shortlist_name in ont else None, by=by, **kwargs)
        return zscores

    def get_avr_expr_by_gene(self, overwrite=True, query=None, add_external=True): # 'Heart_and_Aorta'):
        avr_expr_path = join("../../data/00_data_ingest", 'expr_by_gene_%s.tsv.gz' % (self.method))
        if not exists(avr_expr_path):
            genes, cpm_by_tissue_and_ont = self.get_cpm_by_tissue_and_ont()
            res = self.group_genes_by_ont(genes, cpm_by_tissue_and_ont)
            DataFrameAnalyzer.to_tsv_gz(res, avr_expr_path)

        print('reading...')
        print(avr_expr_path)

        res = DataFrameAnalyzer.read_tsv_gz(avr_expr_path)
        if add_external:
            import scanpy as sc
            res2 = []
            bkp_path_ext = avr_expr_path.replace(".tsv.gz", '_external.tsv.gz')
            if not exists(bkp_path_ext):
                for next_dir in ['../../data/dobnikar_et_al/normed_ln_cpm_p1',
                                 '../../data/MCA']:
                    for f in listdir(next_dir):
                        if not f.endswith('.h5'):
                            continue
                        print('reading external', next_dir, f)
                        p = join(next_dir, f)
                        df = sc.read_h5ad(p).to_df()
                        # df = DataFrameAnalyzer.read_tsv_gz(p)
                        df_long = df.mean(axis=0).reset_index()
                        df_long.columns = ['gene', 'avr.expr']
                        df_long['n'] = df.shape[0]
                        df_long['median.expr'] = list(df.median(axis=0))
                        df_long['ont'] = f.replace(".tsv.gz", '')
                        res2.append(df_long)
                res2 = pd.concat(res2)
                DataFrameAnalyzer.to_tsv_gz(res2, bkp_path_ext)
            res2 = DataFrameAnalyzer.read_tsv_gz(bkp_path_ext)
            res = pd.concat([res, res2])

        return res

    @staticmethod
    def get_topn_genes_by_zscores(n_genes, method, shortlist_label, by='median.expr', tfs_only=False):
        zscores = TabulaMuris.get_zscores_three_ontology_groups(by=by, tfs_only=tfs_only)
        zscores = zscores[(zscores['method'] == method)]
        genes_by_ont = {}
        # define top-NGENES genes by k
        label = 'z.score.' + shortlist_label
        for k, grp in zscores.groupby('ont'):
            sel = grp.sort_values(label, ascending=False).head(n_genes)
            sel = sel[sel[label] > 0]
            genes_by_ont[k] = set(sel['gene'])
        return genes_by_ont

    @staticmethod
    def get_zscores_three_ontology_groups(tfs_only=False, by='avr.expr', overwrite=False):
        '''
        Given three definitions for ontologies, return the respective Z-scores
        :param tfs_only:
        :return:

        '''

        bkp_path = "/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data_2018-07-09_00-00-00/data/zscores_tabula_muris_%s_TFs_%s.tsv.gz" %\
                   (by, tfs_only)
        if not exists(bkp_path) or overwrite:
            if tfs_only:
                from lib.HumanTFs import HumanTFs
                tfs = HumanTFs.get_tf_names()
                from lib.MyGeneAnalyzer import MyGeneAnalyzer
                ensgmus_by_ensghuman = MyGeneAnalyzer.get_homologs_mouse_from_human(set(tfs['Ensembl ID']))
                mouse_tf_symbols = set(MyGeneAnalyzer.get_symbol_by_ensembl(list(ensgmus_by_ensghuman.values()), species='mouse').values())

            all = []
            z_by_k = {}
            raw_by_method = {}
            for method in ['FACS', 'droplet']:
                print(method)
                tm = TabulaMuris(method=method)
                for shortlist_label in ['all', 'shortlist1', 'shortlist2']:
                    print(shortlist_label)
                    k = method + "." + shortlist_label
                    if not k in z_by_k:
                        z = tm.get_zscores_by_shortlist(shortlist_label, n_cells_cutoff=100, by=by)
                        if tfs_only:
                            z = z[[c for c in z if c in mouse_tf_symbols]]
                        z = z.transpose()
                        z_by_k[k] = z

                        if shortlist_label == 'all':
                            raw = tm.get_zscores_by_shortlist(shortlist_label, n_cells_cutoff=100, get_raw=True, by=by)
                            raw_by_method[method] = raw.transpose()
            res = []
            for method in ['FACS', 'droplet']:
                for ont in z_by_k[method + ".all"].columns:
                    print(method, ont)
                    for g in z_by_k[method + ".all"].index:
                        if g not in {'Gpr84', 'Aida'}:
                            continue
                        v = [method, ont, g, raw_by_method[method][ont][g]]
                        for shortlist_label in ['all', 'shortlist1', 'shortlist2']:
                            v += [z_by_k[method + "." + shortlist_label][ont][g] if ont in z_by_k[
                                method + "." + shortlist_label] else np.nan]
                        res.append(v)

            res = pd.DataFrame(res, columns=['method', 'ont', 'gene', 'cpm', 'z.score.all', 'z.score.shortlist1',
                                             'z.score.shortlist2'])

            DataFrameAnalyzer.to_tsv_gz(res, bkp_path)

        return DataFrameAnalyzer.read_tsv_gz(bkp_path)

    def get_zscores_by_gene(self, n_cells_cutoff=None, overwrite=False, add_external=False, inputdir='input'):
        zscores_path = join(inputdir, 'expr_by_gene_zscores%s_%s.tsv.gz' %
                            (("_" + str(n_cells_cutoff)) if n_cells_cutoff is not None else '',
                            'with_external_' + str(int(add_external))))

        # zscores_path = join("/home/rio/data/tabula-muris/00_data_ingest", 'expr_by_gene_%s.tsv.gz' %
        #                     (str(self.method)))
        print(exists(zscores_path), zscores_path)
        if not exists(zscores_path):
            expr = self.get_avr_expr_by_gene(add_external=add_external)
            expr['k'] = expr['ont']
            mu_n_sigma_by_gene = {}

            table = []
            n_genes = 0
            n_genes_tot = len(set(expr['gene']))
            print('# genes to scan', n_genes_tot)
            for g, grp in expr.groupby('gene'):
                if n_genes % 100 == 0:
                    print(g, n_genes, 'out of', n_genes_tot, ', # rows', len(table))
                n_genes += 1
                for ont, grp2 in grp.groupby('ont'):
                    avr_by_ont, n_cells = sum(grp2['avr.expr'] * grp2['n']) / sum(grp2['n']), sum(grp2['n'])
                    table.append([ont, g, n_cells, avr_by_ont])

            table = pd.DataFrame(table, columns=['ont', 'gene.name', 'n.cells', 'avr.expr'])
            for g, grp in table[table['n.cells'] >= n_cells_cutoff].groupby('gene.name'):
                mu = np.mean(grp['avr.expr'])
                sigma = np.std(grp['avr.expr'])
                if sigma == 0:
                    sigma = 1.0
                mu_n_sigma_by_gene[g] = [mu, sigma]

            mu = np.array([x[0] for x in table['gene.name'].map(mu_n_sigma_by_gene)])
            sigma = np.array([x[1] for x in table['gene.name'].map(mu_n_sigma_by_gene)])
            table['z.score'] = (table['avr.expr'] - mu) / sigma
            table = table[table['n.cells'] >= n_cells_cutoff]
            DataFrameAnalyzer.to_tsv_gz(table, zscores_path)
        return DataFrameAnalyzer.read_tsv_gz(zscores_path)

    def get_genes_by_ont(self, n_genes, n_cells_cutoff=100, field='gene.name', add_external=False):
        plot = False
        # print('loading zscores expression...')

        zscores = self.get_zscores_by_gene(n_cells_cutoff=n_cells_cutoff, add_external=add_external)

        genes = set(zscores['gene.name'])
        ens_by_gene_bkp = "input/ensembl_by_symbol_tabula_muris.bkp"
        if not exists(ens_by_gene_bkp):
            name_by_gene = MyGeneAnalyzer.get_ensembl_by_symbol(genes, species='mouse')
            DataFrameAnalyzer.to_pickle(name_by_gene, ens_by_gene_bkp)
        name_by_gene = DataFrameAnalyzer.read_pickle(ens_by_gene_bkp)
        zscores['ensembl'] = zscores['gene.name'].map(name_by_gene)

        query_tissues = None  # {'Brain', 'Pancreas', 'Limb_Muscle', 'Skin', 'Liver', 'Fat', 'Heart', 'Lung'}

        genes_by_ont = {}
        # define top-NGENES genes by k
        for k, grp in zscores.groupby('ont'):
            reject = True
            if query_tissues is None:
                reject = False
            else:
                for s in query_tissues:
                    if s in k:
                        reject = False
            if reject:
                continue
            sel = grp.sort_values('z.score', ascending=False).head(n_genes)
            sel = sel[sel['z.score'] > 0]
            genes_by_ont[k] = set([s for s in sel[field] if str(s) != 'nan'])

        return genes_by_ont