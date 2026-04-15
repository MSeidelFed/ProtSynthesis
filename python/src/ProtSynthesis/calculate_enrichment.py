import polars as pl
import scipy.stats

import matplotlib.pyplot as plt

def calculate_enrichment(gene_enrichment_groups, df_gene):

    """
    Function to calculate enrichment
    gene_enrichment_groups: dict (containing enriched and all genes for each experiment)
    df_gene: pl.DataFrame (containing codon and aa counts generated in count_cds.py)
    function returns a pl.DataFrame consisting of the [sample_id, codon, fold_change, pval_sf, pval_cdf]
    """

    enrichment = []

    for sample, genes in gene_enrichment_groups.items():
        enriched = genes['enriched']
        control = genes['all']


        enriched_counts = df_gene.filter(pl.col('gene_id').is_in(enriched)).select(pl.selectors.numeric()).sum().to_dicts()[0]
        control_counts = df_gene.filter(pl.col('gene_id').is_in(control)).select(pl.selectors.numeric()).sum().to_dicts()[0]

        enriched_codon_count = enriched_counts.pop('codon_count')
        control_codon_count = control_counts.pop('codon_count')
        

        for k in enriched_counts:
            e = enriched_counts[k]
            c = control_counts[k]


            fc = (e/enriched_codon_count) / (c/control_codon_count) if c != 0 else np.nan 
            pval_sf = scipy.stats.hypergeom.sf(e-1,control_codon_count, c, enriched_codon_count)
            pval_cdf = scipy.stats.hypergeom.cdf(e,control_codon_count, c, enriched_codon_count)

            enrichment.append([sample, k, fc, pval_sf, pval_cdf])

    df = pl.DataFrame(enrichment, schema = ['sample_id','codon','fold_change', 'pval_sf', 'pval_cdf'], orient='row').drop_nans()

    return df



def plot_fold_changes(df):

    """
    Generates bar plots of log2 fold change for codons and amino acids derived for the pl.DataFrame returned by calculate_enrichment()
    
    """

    sample_ids = df.select('sample_id').unique().to_numpy().T.tolist()[0]


    for i,s in enumerate(sample_ids[:-1]):

        fig_aa,ax_aa = plt.subplots(1)
        df_aa = df.sort('fold_change').filter((pl.col('codon').str.len_chars() == 1) & (pl.col('sample_id') == s))    
        ax_aa.bar(df_aa.select('codon').to_numpy().T.tolist()[0],np.log2(df_aa.select('fold_change').to_numpy().T[0]))
        ax_aa.spines['right'].set_visible(False)
        ax_aa.spines['top'].set_visible(False)
        fig_aa.suptitle(s)
        fig_aa.savefig(f'{s}_shoots_aas.pdf')


        fig_codon,ax_codon = plt.subplots(1)
        df_codon = df.sort('fold_change').filter((pl.col('codon').str.len_chars() == 3) & (pl.col('sample_id') == s))
        ax_codon.bar(df_codon.select('codon').to_numpy().T.tolist()[0],np.log2(df_codon.select('fold_change').to_numpy().T[0]))
        ax_codon.tick_params('x',rotation = 90)
        ax_codon.spines['right'].set_visible(False)
        ax_codon.spines['top'].set_visible(False)
        fig_codon.suptitle(s)

        fig_codon.savefig(f'{s}_shoots_codons.pdf')




def plot_fold_change_accross_samples(df):

    """
    Still working on these plots
    """


    fig,ax = plt.subplots(1)

    df_aa = df.sort('codon').filter((pl.col('codon').str.len_chars() == 1))
    scipy.stats.pearsonr(np.log2(df_aa.filter(pl.col('sample_id') == sample_ids[1]).select('fold_change').to_numpy().T[0]),np.log2(df_aa.filter(pl.col('sample_id') == sample_ids[2]).select('fold_change').to_numpy().T[0]))


