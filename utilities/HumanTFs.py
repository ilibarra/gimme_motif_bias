'''
Created on 6/24/2018

DESCRIPTION

@author: ignacio
'''

from .DataFrameAnalyzer import DataFrameAnalyzer

class HumanTFs:
    @staticmethod
    def get_tf_names(**kwargs):
        p = join(kwargs.get('datadir', '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data'), 'DatabaseExtract_v_1.01.csv')
        return DataFrameAnalyzer.read_tsv(p, sep=',')

    @staticmethod
    def get_tf_motifs_cisbp(**kwargs):
        p = join(kwargs.get('datadir', '/g/scb2/zaugg/rio/EclipseProjects/zaugglab/moritz_collaboration/data'), 'Human_TF_MotifList_v_1.01.csv')
        return DataFrameAnalyzer.read_tsv(p, sep=',')

    @staticmethod
    def get_tf_motifs_cisbp_best():
        df = HumanTFs.get_tf_motifs_cisbp()
        return df[df['Best Motif(s)? (Figure 2A)'] == True]

    @staticmethod
    def get_ensembl_by_dbd():
        dbd_df = HumanTFs.get_tf_names()
        return {dbd_id: set(grp['Ensembl ID']) for dbd_id, grp in dbd_df.groupby('DBD')}

    @staticmethod
    def get_dbd_by_ensembl(**kwargs):
        dbd_df = HumanTFs.get_tf_names(**kwargs)
        return DataFrameAnalyzer.get_dict(dbd_df, 'Ensembl ID', 'DBD')

    @staticmethod
    def get_ppm(cisbp_id):
        model_path = '/g/scb2/zaugg/zaugg_shared/annotations/motifs/CISBP/build_1.94d/PWMs/%s.txt' % cisbp_id
        if exists(model_path):
            ppm = DataFrameAnalyzer.read_tsv(model_path)
            ppm = ppm[ppm.columns[-4:]].transpose()
            return ppm
        else:
            assert exists(model_path)

    @staticmethod
    def plot_pwm_model(model_id, show_complementary=False, ax=None, title=None, ppm=None):
        if ppm is None:
            model_path = '/g/scb2/zaugg/zaugg_shared/annotations/motifs/CISBP/build_1.94d/PWMs/%s.txt' % model_id
            assert exists(model_path)
            ppm = DataFrameAnalyzer.read_tsv(model_path)
            ppm = ppm[ppm.columns[-4:]].transpose()

        from lib.motif_plotter import ConsensusMotifPlotter
        from lib.Motif.MotifAnalyzer import MotifAnalyzer
        motif_analyzer = MotifAnalyzer()
        if show_complementary:
            ppm = motif_analyzer.get_complementary_ppm(ppm)

        cbp = ConsensusMotifPlotter.from_ppm(ppm)
        ax = plt.subplot() if ax is None else ax
        cbp.plot(ax)
        remove_top_n_right_ticks(ax)
        if title is not None:
            plt.title(title, fontsize=8)

    @staticmethod
    def get_motif_scores(bed_coordinates, genome='hg19', stop_at=None, n_cores=1,
                         query_column=None, query_value=None):

        cisbp = HumanTFs.get_tf_motifs_cisbp_best().reset_index(drop=True)

        if query_column is not None and query_value is not None:
            cisbp = cisbp[cisbp[query_column] == query_value].reset_index(drop=True)

        from lib.FastaAnalyzer import FastaAnalyzer

        fa =  FastaAnalyzer.get_sequences_from_bed(bed_coordinates, genome=genome)
        fa_path = tempfile.mkstemp()[1]
        FastaAnalyzer.write_fasta_from_sequences(fa, fa_path)

        from lib.ThreadingUtils import ThreadingUtils
        from multiprocessing import Manager
        manager = Manager()
        out = manager.dict()

        def get_motif_hits(r, ri, n_cores=1):
            print((ri, 'out of', cisbp.shape[0]))
            # SCORE MOTIFS FOR BCL6 using CIS-BP
            cisbp_id = r['CIS-BP ID']
            cisbp_path = '/g/scb2/zaugg/zaugg_shared/annotations/motifs/CISBP/build_1.94d/PWMs/' + cisbp_id + ".txt"
            print((exists(cisbp_path), cisbp_path))

            # TRANSFAC motifs do not contain information: skip those
            if not exists(cisbp_path):
                return
            nlines = [r for r in open(cisbp_path)]
            if len(nlines) == 1:
                print('# of lines is len than 1 (no PWM info. Probably TRANSFAC)')

            motif_output_dir = '/tmp/motif_output_%i' % (ri)
            if not exists(motif_output_dir):
                mkdir(motif_output_dir)
            motif_code = basename(cisbp_path).replace(".txt", '')
            motif_output_path = join(motif_output_dir, motif_code + ".tsv.gz")
            print((exists(motif_output_path), motif_output_path))
            reject = False
            if exists(motif_output_path):
                # last test: the motif id has to be the same as the motif id
                df = DataFrameAnalyzer.read_tsv_gz(motif_output_path)
                if list(set(df['motif_id'].str.replace(".txt", '')))[0] != motif_code:
                    reject = True
                print('output exists. Skip...')
                print(('reject?', reject))
                print(motif_output_path)

            meme_tpm_path = join('/tmp/', motif_code + ".meme")
            from lib.Motif.MotifConverter import MotifConverter
            MotifConverter.convert_cisbp_to_meme(cisbp_path, meme_tpm_path)

            fimo_output_dir = '/tmp/fimo_output_%i'
            if not exists(fimo_output_dir):
                mkdir(fimo_output_dir)
            fimo_output_dir = join(fimo_output_dir, motif_code)

            fimo_cmd = ' '.join(['fimo', '-bgfile', '--uniform--', '--o', fimo_output_dir, meme_tpm_path, fa_path])
            print(fimo_cmd)
            system(fimo_cmd)

            # collect output
            res = DataFrameAnalyzer.read_tsv(join(fimo_output_dir, 'fimo.tsv'))
            DataFrameAnalyzer.to_tsv_gz(res, motif_output_path)
            remove(meme_tpm_path)
            system('rm -rf ' + fimo_output_dir)
            system('rm -rf ' + motif_output_dir)
            out[ri] = res

        ThreadingUtils.run(get_motif_hits, [[r, ri] for ri, r in cisbp.iterrows()][:stop_at if stop_at is not None else cisbp.shape[0]],
                           n_cores=n_cores)

        out = pd.concat([out[pi] for pi in list(out.keys())])
        return out[~out['motif_id'].str.startswith("#")].reset_index(drop=True)


