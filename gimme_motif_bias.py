'''
Created on 04/19/2020

Calculate motif biases for motif pairs, using standalone Python libraries
@author: ignacio
'''


import utilities
from utilities import *

# Main script function
def calculate_motif_bias(a, b, motif_id, **kwargs):
    tm = TabulaMuris(method='FACS')

    plot = False
    print('loading zscores expression...')
    zscores = tm.get_zscores_by_gene(n_cells_cutoff=10, add_external=True, overwrite=False,
                                     inputdir='input')
    genes = set(zscores['gene.name'])

    input_dir, output_dir = kwargs.get('inputdir'), kwargs.get('outputdir')

    ens_by_gene_bkp = join(input_dir, "ensembl_by_symbol_tabula_muris.bkp")
    if not exists(ens_by_gene_bkp):
        name_by_gene = MyGeneAnalyzer.get_ensembl_by_symbol(genes, species='mouse')
        DataFrameAnalyzer.to_pickle(name_by_gene, ens_by_gene_bkp)
    name_by_gene = DataFrameAnalyzer.read_pickle(ens_by_gene_bkp)
    zscores['ensembl'] = zscores['gene.name'].map(name_by_gene)

    query_tissues = None # {'Brain', 'Pancreas', 'Limb_Muscle', 'Skin', 'Liver', 'Fat', 'Heart', 'Lung'}
    N_GENES = 1000
    genes_by_ont = tm.get_genes_by_ont(N_GENES, add_external=True, n_cells_cutoff=10)

    # save these genes in each tissue, assuming files do not exist
    output_dir_genes = "input/genes_by_ont_tabula_muris"
    if not exists(output_dir_genes):
        mkdir(output_dir_genes)
    for ont in genes_by_ont:
        output_path = join(output_dir_genes, ont + ".txt")
        # print(exists(output_path), output_path)
        if not exists(output_path):
            DataFrameAnalyzer.write_list(genes_by_ont[ont], output_path)

    print('reading results motif hits')

    # cisbp motif hits =
    motifs_path_all = 'input/motif_hits_cisbp_build_1.94d_mm10.tsv.gz'

    print('reading all motifs paths...')
    print(motifs_path_all)
    all_motifs_df = DataFrameAnalyzer.read_tsv_gz(motifs_path_all)

    for m, grp in all_motifs_df.groupby('#pattern name'):
        print(m, grp.shape[0])
        DataFrameAnalyzer.to_tsv_gz(grp, join('input/motif_hits_cisbp_build_1.94d_mm10/%s' % m.replace(".txt", '.tsv.gz')))


    print('done...')
    motif_ids = set(all_motifs_df['filename'].str.replace(".tsv.gz", ""))
    print(all_motifs_df.head())
    print(len(motif_ids))
    print(list(motif_ids)[:10])

    tfs = HumanTFs.get_tf_motifs_cisbp(datadir="../../data/")
    ensembl_by_model = {m: set(grp['Ensembl ID']) for m, grp in tfs.groupby('CIS-BP ID')}
    tfname_by_ensembl = DataFrameAnalyzer.get_dict(tfs, 'Ensembl ID', 'HGNC symbol')

    human_orthologs = DataFrameAnalyzer.read_tsv("../../data/human_mm10_homologs.tsv", sep='\t')
    # print human_orthologs.head()
    engmus_by_enghuman = DataFrameAnalyzer.get_dict(human_orthologs, 'Gene stable ID', 'Mouse gene stable ID')

    # using new ENSEMBL identifiers to map between old and new db
    updated_gtex_identifiers = GTEXAnalyzer.get_updated_ensembl_identifiers(datadir='../../data')
    query = kwargs.get('query', None)
    n_group = kwargs.get('ngroup', 1)

    counter = -1

    counter += 1

    print(counter, 'next motif', m)

    grp = None
    if all_motifs_df is None:
        grp = DataFrameAnalyzer.read_tsv_gz(join(motifs_dir, m + ".tsv.gz"))
    else:
        grp = all_motifs_df[all_motifs_df['filename'] == m]

    grp.columns = ['motif.id', 'ensembl'] + list(grp.columns[2:])
    # put gene name into res column
    from lib.SequenceMethods import SequenceMethods
    tss = MyGeneAnalyzer.get_gene_tss('mouse', 2000, datadir="../../data")
    tss = SequenceMethods.parse_range2coordinate(tss)
    grp['gene.name'] = grp['ensembl'].map(DataFrameAnalyzer.get_dict(tss, 'range', 'SYMBOL'))

    # if m != 'NR2E3_HUMAN.H11MO.0.C':
    #     continue
    # print m, grp.shape[0]
    # print updated_gtex_identifiers.head()

    # some cisbp motifs are associated to more than one gene
    for ensg_human in ensembl_by_model[m]:
        code = m + "_" + ensg_human
        print(code)
        if query is not None:
            if not query in code:
                continue

        print('analyzing', ensg_human)

        # if code != 'M06660_1.94d_ENSG00000082641':
        #     continue
        pkl_path = join('../../data/enrichment_heatmaps_cisbp_mm10', 'tabula_muris',
                        m + "_" + ensg_human + ".pkl")
        if exists(pkl_path):
            print('pkl path exists for engs_human', ensg_human, 'skip...')
            continue
        output_dir_mm10 = "../../data/figures/enrichment_n_depletion_clustermaps_expr_atlas_CISBP_build_1.94d_mm10"
        output_dir_mm10_motif_biases = "../../data/figures/enrichment_n_depletion_clustermaps_tabula_muris_CISBP_build_1.94d_mm10"
        if not exists(output_dir_mm10):
            mkdir(output_dir_mm10)
        output_path = join(output_dir_mm10, m + "_" + ensg_human)

        if not exists(output_dir_mm10_motif_biases):
            mkdir(output_dir_mm10_motif_biases)
        output_path_motif_biases = join(output_dir_mm10_motif_biases, m + ".tsv.gz")

        # print updated_gtex_identifiers[updated_gtex_identifiers['Description'] == 'DUX4']
        sel = updated_gtex_identifiers[(updated_gtex_identifiers['ensembl.mygene'] == ensg_human) |
                                       (updated_gtex_identifiers['ensembl.gtex'] == ensg_human)]

        if sel.shape[0] < 1 and len(set(sel['ensembl.gtex'])) != 1:
            print(m, ensg_human)
            print(sel)
            print('skip...')
            continue
        else:
            print(sel)

        ensembl_gtex = list(sel['Name'])[0]
        print(ensembl_gtex)

        fig = plt.figure(figsize=(kwargs.get('w', 16), kwargs.get('h', 5.5)))
        plot_i = 0
        for zscores, label in zip([zscores], ['tabula-muris']):
            print('next check:', label)
            if not 'Mouse gene stable ID' in zscores:
                zscores['Mouse gene stable ID'] = [engmus_by_enghuman[idx]
                                                   if idx in engmus_by_enghuman else None
                                                   for idx in zscores.index]
            mouse_ensg_set = set(zscores['Mouse gene stable ID'])
            print(output_path)
            # print jobid, counter, m, ensembl_by_model[m]

            # print genes_df.keys()
            sel_keys = list(genes_by_ont.keys())
            # sel_keys = {'Pancreas_leukocyte', 'Liver_hepatocyte', 'Brain_Non-Myeloid_neuron', 'Limb_Muscle_skeletal muscle satellite cell',
            #             'Heart_cardiac muscle cell'}

            if not exists(output_path_motif_biases):
                pairwise_combinations = [[c1, c2] for c1, c2 in product(sel_keys, repeat=2)]
                # calculate enrichmments using a new method
                print('calculating pair wise enrichments...')
                try:
                    enrichments = EnrichmentAnalyzer.get_motif_enrichments_by_pairwise_grouping(genes_by_ont, grp,
                                                                                                label=m,
                                                                                                column_gene='gene.name')
                    enrichments['odds.ratio'] = np.where(enrichments['a'] == enrichments['b'], 1.0,
                                                         enrichments['odds.ratio'])
                    enrichments['log2FC'] = np.where(enrichments['a'] == enrichments['b'], 0.0,
                                                     enrichments['log2FC'])
                    # print list(genes_df.keys())
                    # print enrichments
                    print('enrichments calculated...')
                    if 'ensembl' in enrichments:
                        del enrichments['ensembl']

                    DataFrameAnalyzer.to_tsv_gz(enrichments, output_path_motif_biases)
                except Exception:
                    print('an error happened while calculating pairwise enrichments. Please check...')
                    continue

            enrichments = DataFrameAnalyzer.read_tsv_gz(output_path_motif_biases)
            # get the values from the actual Z-scores from expression
            # print m, m in ensg_by_motif_id
            # print 'here...', len(ensg_by_motif_id)
            zscore_by_tissue = {}
            reject = False

            ensg_human = ensg_human.replace(".", '')
            print(ensg_human)
            ensg_mouse = engmus_by_enghuman[ensg_human] if ensg_human in engmus_by_enghuman else None

            # print ensg_mouse

            print('label', ensembl_gtex)
            if label == 'tabula-muris':
                columns = []  # ['Brain_Non-Myeloid_neuron', 'Pancreas_endothelial cell']
                cell_types = []  # ['Brain', 'Pancreas']

            # print list(zscores.columns)
            n_mouse = zscores[zscores['ensembl'] == ensg_mouse].shape[0]
            n_human = zscores[zscores['ensembl'] == ensg_human].shape[0]

            print('# mouse', n_mouse)
            print('# human', n_human)
            # if label == 'GTEx' and n_mouse == 0 or label == 'Expression Atlas' and n_human == 0:
            #     reject = True
            # if reject:
            #     break

            print('label', label)
            finish = False

            for ri, r in zscores[zscores['ensembl'] == ensg_mouse].iterrows():
                zscore_by_tissue[r['ont']] = r['z.score']

            if finish:
                if label == 'GTEx':
                    continue
                else:
                    print('breaking...')
                    break

                    # print zscore_by_tissue
                    # if reject: # ensg was not found in requested table: skip
                    # continue
            enrichments['p.adj'] = RFacade.get_bh_pvalues(enrichments['p.val'])

            # print zscore_by_tissue
            enrichments['z.score.expr'] = [zscore_by_tissue[k] if k in zscore_by_tissue else np.nan
                                           for k in enrichments['a']]
            # print enrichments
            # print enrichments
            # print enrichments.sort_values('p.adj', ascending=True)
            # print '\nsignificant:'
            # print enrichments[enrichments['p.adj'] < 0.05]

            enrichments['p.adj.symbol'] = RFacade.get_pval_asterisks(enrichments['p.adj'])
            # print enrichments

            # print enrichments
            hm = DataFrameAnalyzer.dataframe_to_matrix(enrichments[['b', 'a', 'log2FC']])
            # print enrichments.sort_values('log2FC')

            enrichments['k'] = [str(a) + "\n" + ("%.2f" % b) for a, b in zip(enrichments['p.adj.symbol'],
                                                                             enrichments['z.score.expr'])]
            # print enrichments

            symbols = DataFrameAnalyzer.dataframe_to_matrix(enrichments[['b', 'a', 'z.score.expr']])
            # DataFrameAnalyzer.to_tsv_gz(hm, '../../data/playtest_sizes.tsv.gz', index=True)
            # DataFrameAnalyzer.to_tsv_gz(symbols, '../../data/playtest_color.tsv.gz', index=True)

            # symbols = DataFrameAnalyzer.dataframe_to_matrix(enrichments[['b', 'a', 'k']])


            old_labels = symbols.columns
            new_label_by_old = {}
            for k in old_labels:
                new_name = k
                for a, b in zip(columns, cell_types):
                    new_name = new_name.replace(a, b)
                new_label_by_old[k] = new_name

            symbols.columns = [new_label_by_old[k] if k in new_label_by_old else k for k in symbols.columns]
            hm.columns = [new_label_by_old[k] if k in new_label_by_old else k for k in hm.columns]
            symbols.index = [new_label_by_old[k] if k in new_label_by_old else k
                             for k in symbols.index]
            hm.index = [new_label_by_old[k] if k in new_label_by_old else k
                        for k in hm.index]

            symbols.sort_index(inplace=True)
            hm.sort_index(inplace=True)
            symbols = symbols[sorted(symbols.columns.astype(str))]
            hm = hm[sorted(hm.columns.astype(str))]

            # we save a heatmap with the obtained values and pvalues
            # print symbols
            print(hm.shape)

            if True or not exists(pkl_path):
                DataFrameAnalyzer.to_pickle([symbols, hm, enrichments], pkl_path)

            print('done...')
            # exit()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=int, default='hepatocyte', help='minimum primer length (def. 20)')
    parser.add_argument("--b", type=int, default='neuron', help='minimum primer length (def. 24)')
    parser.add_argument("--inputdir", type=int, default='input', help='input directory')
    parser.add_argument("--output", type=int, default='output', help='output directory')

    parser.add_argument("--motifid", type=float, default=None,
                        help='motif.id')

    opts = parser.parse_args()

    calculate_motif_bias(opts.a, opts.b, opts.motifid, inputdir=opts.inputdir, output=opts.output)