'''
Created on 04/19/2020

Calculate motif biases for motif pairs, using standalone Python libraries
@author: ignacio
'''


import utilities
from utilities import *
from itertools import product
import numpy as np
import pandas as pd

# Main script function
def calculate_motif_bias(a, b, motif_id, **kwargs):

    # save these genes in each tissue, assuming files do not exist
    input_dir, output_dir = kwargs.get('indir'), kwargs.get('outdir')

    ontdir = join(input_dir, "genes_by_ont")
    if kwargs.get('listont'):
        print('available cell types')
        for f in listdir(ontdir):
            d = join(ontdir, f)
            if isdir(d):
                print('ontgroup:%s' % f)
                for f2 in listdir(d):
                    print("\t" + f2.replace(".txt", ''))
            else:
                print(f.replace(".txt", ''))
        return
    if kwargs.get('listmotifs') is not None:
        tfs = HumanTFs.get_tf_motifs_cisbp(datadir="input")
        print(tfs[tfs['HGNC symbol'].str.lower().str.contains(kwargs.get('listmotifs').lower())])
        return

    is_group_a = isdir(join(ontdir, a))
    is_group_b = isdir(join(ontdir, b))
    names_a = {a} if not is_group_a else {f.replace(".txt", '') for f in listdir(join(ontdir, a))}
    names_b = {b} if not is_group_b else {f.replace(".txt", '') for f in listdir(join(ontdir, b))}


    print('Comparing %s vs %s' % (a, b))
    tm = TabulaMuris(method='FACS')

    plot = False
    print('loading zscores expression...')
    zscores = tm.get_zscores_by_gene(n_cells_cutoff=10, add_external=True, overwrite=False,
                                     inputdir='input')
    genes = set(zscores['gene.name'])


    ens_by_gene_bkp = join(input_dir, "ensembl_by_symbol_tabula_muris.bkp")
    if not exists(ens_by_gene_bkp):
        name_by_gene = MyGeneAnalyzer.get_ensembl_by_symbol(genes, species='mouse')
        DataFrameAnalyzer.to_pickle(name_by_gene, ens_by_gene_bkp)
    name_by_gene = DataFrameAnalyzer.read_pickle(ens_by_gene_bkp)
    zscores['ensembl'] = zscores['gene.name'].map(name_by_gene)

    query_tissues = None # {'Brain', 'Pancreas', 'Limb_Muscle', 'Skin', 'Liver', 'Fat', 'Heart', 'Lung'}
    N_GENES = kwargs.get('ngenes')

    genes_by_ont = tm.get_genes_by_ont(N_GENES, add_external=True, n_cells_cutoff=10)

    if not exists(ontdir):
        mkdir(ontdir)
    for ont in genes_by_ont:
        output_path = join(ontdir, ont + ".txt")
        # print(exists(output_path), output_path)
        if not exists(output_path):
            DataFrameAnalyzer.write_list(genes_by_ont[ont], output_path)

    print('Comparing %s versus %s' % (a, b))

    print('reading results motif hits')

    # cisbp motif hits =
    motifs_path_all = 'input/motif_hits_cisbp_build_1.94d_mm10.tsv.gz'

    print('reading all motifs paths...')
    print(motifs_path_all)


    print('done...')
    motif_ids = set()
    for f in listdir(join(input_dir, 'motif_hits_cisbp_build_1.94d_mm10')):
        motif_ids.add(f.replace(".tsv.gz", ''))

    if not motif_id in motif_ids:
        print('Your query motif ID is not in the mapped motifs directory. Cannot execute...')
        assert not motif_id in motif_ids
    else:
        print('%s query motif found in motifs directory' % motif_id)

    tfs = HumanTFs.get_tf_motifs_cisbp(datadir=input_dir)
    ensembl_by_model = {m: set(grp['Ensembl ID']) for m, grp in tfs.groupby('CIS-BP ID')}
    tfname_by_ensembl = DataFrameAnalyzer.get_dict(tfs, 'Ensembl ID', 'HGNC symbol')

    human_orthologs = DataFrameAnalyzer.read_tsv("input/human_mm10_homologs.tsv", sep='\t')
    # print human_orthologs.head()
    engmus_by_enghuman = DataFrameAnalyzer.get_dict(human_orthologs, 'Gene stable ID', 'Mouse gene stable ID')

    # using new ENSEMBL identifiers to map between old and new db
    updated_gtex_identifiers = GTEXAnalyzer.get_updated_ensembl_identifiers(datadir='input')

    motifs_path = join('input', 'motif_hits_cisbp_build_1.94d_mm10', motif_id + ".tsv.gz")
    print(exists(motifs_path)), motifs_path
    assert exists(motifs_path)
    grp = DataFrameAnalyzer.read_tsv_gz(join('input', 'motif_hits_cisbp_build_1.94d_mm10', motif_id + ".tsv.gz"))

    grp.columns = ['motif.id', 'ensembl'] + list(grp.columns[2:])
    # put gene name into res column
    tss = MyGeneAnalyzer.get_gene_tss('mouse', 2000, datadir=input_dir)
    tss = SequenceMethods.parse_range2coordinate(tss)
    grp['gene.name'] = grp['ensembl'].map(DataFrameAnalyzer.get_dict(tss, 'range', 'SYMBOL'))

    # if m != 'NR2E3_HUMAN.H11MO.0.C':
    #     continue
    # print m, grp.shape[0]
    # print updated_gtex_identifiers.head()

    # some cisbp motifs are associated to more than one gene
    for ensg_human in ensembl_by_model[motif_id]:
        code = motif_id + "_" + ensg_human
        print(code)
        print('analyzing', ensg_human)

        # if code != 'M06660_1.94d_ENSG00000082641':
        #     continue
        output_dir_enrichments = join(output_dir, 'enrichment_heatmaps_cisbp_mm10', 'tabula_muris')
        if not exists(output_dir_enrichments):
            makedirs(output_dir_enrichments)

        pkl_path = join(output_dir_enrichments, motif_id + "%s_%s_%s" % (ensg_human, a, b) + ".pkl")
        print(exists(pkl_path), pkl_path)
        if exists(pkl_path) and not kwargs.get('overwrite'):
            print('pkl path exists for engs_human', ensg_human, 'skip...')
            continue
        output_dir_mm10 = join("%s" % output_dir,
                               "enrichment_n_depletion_clustermaps_expr_atlas_CISBP_build_1.94d_mm10")
        output_dir_mm10_motif_biases = join("%s" % output_dir,
                                            "figures/enrichment_n_depletion_clustermaps_tabula_muris_CISBP_build_1.94d_mm10")
        if not exists(output_dir_mm10):
            makedirs(output_dir_mm10)
        output_path = join(output_dir_mm10, motif_id + "_" + ensg_human)

        if not exists(output_dir_mm10_motif_biases):
            makedirs(output_dir_mm10_motif_biases)
        output_path_motif_biases = join(output_dir_mm10_motif_biases, motif_id + ".tsv.gz")

        # print updated_gtex_identifiers[updated_gtex_identifiers['Description'] == 'DUX4']
        sel = updated_gtex_identifiers[(updated_gtex_identifiers['ensembl.mygene'] == ensg_human) |
                                       (updated_gtex_identifiers['ensembl.gtex'] == ensg_human)]

        if sel.shape[0] < 1 and len(set(sel['ensembl.gtex'])) != 1:
            print(motif_id, ensg_human)
            print(sel)
            print('skip...')
            continue
        else:
            print(sel)

        ensembl_gtex = list(sel['Name'])[0]
        print(ensembl_gtex)

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

            if not exists(output_path_motif_biases) or kwargs.get('overwrite'):
                pairwise_combinations = [[c1, c2] for c1, c2 in product(sel_keys, repeat=2)]
                # calculate enrichmments using a new method
                print('calculating pair wise enrichments...')
                try:
                    sub_genes_ont = {query: genes_by_ont[query] for query in names_a.union(names_b) if query in genes_by_ont}
                    print(len(sub_genes_ont), sub_genes_ont.keys())
                    enrichments = EnrichmentAnalyzer.get_motif_enrichments_by_pairwise_grouping(sub_genes_ont,
                                                                                                grp,
                                                                                                label=motif_id,
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

            if not exists(pkl_path) or kwargs.get('overwrite'):
                DataFrameAnalyzer.to_pickle([symbols, hm, enrichments], pkl_path)
                DataFrameAnalyzer.to_tsv_gz(enrichments, pkl_path.replace(".pkl", '.tsv.gz'))
                if kwargs.get('xlsx'):
                    enrichments.to_excel(pkl_path.replace(".pkl", '.xlsx'))





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--listont", action='store_true', help='List available ont and finish', default=False)
    parser.add_argument("--listmotifs", type=str, default=None, help='Get all motifs associated with TF and finish')
    parser.add_argument("-a", type=str, default='hepatocyte', help='Ontology label A')
    parser.add_argument("-b", type=str, default='neuron', help='Ontology label B or group lable (e.g. shortlist1)')
    parser.add_argument("--ngenes", type=int, help='set number of topN genes for comparison', default=1000)
    parser.add_argument("--indir", type=str, default='input', help='input directory')
    parser.add_argument("--outdir", type=str, default='output', help='output directory')
    parser.add_argument("--overwrite", action='store_true', help='Force writing')
    parser.add_argument("--xlsx", action='store_true', help='Save additional copy as Excel')
    parser.add_argument("--motifid", type=str, default=None, help='motif.id to be used (please run with listmotifs to see which ones are available)')

    opts = parser.parse_args()

    calculate_motif_bias(opts.a, opts.b, opts.motifid, indir=opts.indir, outdir=opts.outdir, listont=opts.listont,
                         listmotifs=opts.listmotifs, ngenes=opts.ngenes, overwrite=opts.overwrite, xlsx=opts.xlsx)