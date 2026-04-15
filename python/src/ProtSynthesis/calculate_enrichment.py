import polars as pl
import scipy.stats
import itertools

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
        depleted = genes['depleted']


        enriched_counts = df_gene.filter(pl.col('gene_id').is_in(enriched)).select(pl.selectors.numeric()).sum().to_dicts()[0]
        control_counts = df_gene.filter(pl.col('gene_id').is_in(control)).select(pl.selectors.numeric()).sum().to_dicts()[0]
        enriched_counts = df_gene.filter(pl.col('gene_id').is_in(enriched)).select(pl.selectors.numeric()).sum().to_dicts()[0]
        depleted_counts = df_gene.filter(pl.col('gene_id').is_in(depleted)).select(pl.selectors.numeric()).sum().to_dicts()[0]


        enriched_codon_count = enriched_counts.pop('codon_count')
        control_codon_count = control_counts.pop('codon_count')
        depleted_codon_count = depleted_counts.pop('codon_count')

        for k in enriched_counts:
            e = enriched_counts[k]
            d = depleted_counts[k]
            c = control_counts[k]


            fc_e = (e/enriched_codon_count) / (c/control_codon_count) if c != 0 else np.nan 
            pval_sf_e = scipy.stats.hypergeom.sf(e-1,control_codon_count, c, enriched_codon_count)
            pval_cdf_e = scipy.stats.hypergeom.cdf(e,control_codon_count, c, enriched_codon_count)

            fc_d = (d/depleted_codon_count) / (c/control_codon_count) if c != 0 else np.nan 
            pval_sf_d = scipy.stats.hypergeom.sf(d-1,control_codon_count, c, depleted_codon_count)
            pval_cdf_d = scipy.stats.hypergeom.cdf(d,control_codon_count, c, depleted_codon_count)




            enrichment.append([sample,'enriched', k, e/enriched_codon_count * 100, fc_e, pval_sf_e, pval_cdf_e])
            enrichment.append([sample,'depleted', k, d/depleted_codon_count * 100,fc_d, pval_sf_d, pval_cdf_d])




    df = pl.DataFrame(enrichment, schema = ['sample_id','fc_type','codon','norm_enrich_counts' ,'fold_change', 'pval_sf', 'pval_cdf'], orient='row').drop_nans()

    return df




def calculate_enrichment_accross_conditions(gene_enrichment_groups, df_gene):

    """
    Function to calculate enrichment
    gene_enrichment_groups: dict (containing enriched and all genes for each experiment)
    df_gene: pl.DataFrame (containing codon and aa counts generated in count_cds.py)
    function returns a pl.DataFrame consisting of the [sample_id, codon, fold_change, pval_sf, pval_cdf]
    """

    enrichment = []

    for e1,e2 in itertools.combinations(gene_enrichment_groups.keys(),2):
        e1_genes = gene_enrichment_groups[e1]['enriched']
        e2_genes = gene_enrichment_groups[e2]['enriched']

        

        enriched_counts = df_gene.filter(pl.col('gene_id').is_in(e1_genes)).select(pl.selectors.numeric()).sum().to_dicts()[0]
        control_counts = df_gene.filter(pl.col('gene_id').is_in(e2_genes)).select(pl.selectors.numeric()).sum().to_dicts()[0]

        enriched_codon_count = enriched_counts.pop('codon_count')
        control_codon_count = control_counts.pop('codon_count')
        


        for k in enriched_counts:
            e = enriched_counts[k]
            c = control_counts[k]



            fc = (e/enriched_codon_count) / (c/control_codon_count) if c != 0 else np.nan 
            pval_sf = scipy.stats.hypergeom.sf(e-1,control_codon_count, c, enriched_codon_count)
            pval_cdf = scipy.stats.hypergeom.cdf(e,control_codon_count, c, enriched_codon_count)

            enrichment.append([f'{e1}-{e2}', k, fc])

            print(enrichment)

    df = pl.DataFrame(enrichment, schema = ['sample_id','codon','fold_change'], orient='row').drop_nans()

    return df





def plot_fold_changes(df):

    """
    Generates bar plots of log2 fold change for codons and amino acids derived for the pl.DataFrame returned by calculate_enrichment()
    
    """
    
    samp, fct = df.select('sample_id','fc_type').unique().to_numpy().T.tolist()


    for i,(s,f) in enumerate(zip(samp,fct)):

        fig_aa,ax_aa = plt.subplots(1)
        df_aa = df.sort('fold_change').filter((pl.col('codon').str.len_chars() == 1) & (pl.col('sample_id') == s) & (pl.col('fc_type') == f))    
        ax_aa.bar(df_aa.select('codon').to_numpy().T.tolist()[0],np.log2(df_aa.select('fold_change').to_numpy().T[0]))
        ax_aa.spines['right'].set_visible(False)
        ax_aa.spines['top'].set_visible(False)
        ax_aa.set_ylabel(f'log2(AA-{f} / AA-all)')
        fig_aa.suptitle(f'{s}-{f}')
        fig_aa.savefig(f'{s.replace('/','-')}-{f}_roots_fc.pdf')

        fig_aa,ax_aa = plt.subplots(1)
        df_aa = df.sort('norm_enrich_counts').filter((pl.col('codon').str.len_chars() == 1) & (pl.col('sample_id') == s) & (pl.col('fc_type') == f))    
        ax_aa.bar(df_aa.select('codon').to_numpy().T.tolist()[0],(df_aa.select('norm_enrich_counts').to_numpy().T[0]))
        ax_aa.spines['right'].set_visible(False)
        ax_aa.spines['top'].set_visible(False)
        ax_aa.set_ylabel(f"{s}")
        fig_aa.suptitle(f'{s}-{f}')
        fig_aa.savefig(f'{s.replace('/','-')}-{f}_roots_aas.pdf')


        fig_codon,ax_codon = plt.subplots(1)
        df_codon = df.sort('fold_change').filter((pl.col('codon').str.len_chars() == 3) & (pl.col('sample_id') == s))
        ax_codon.bar(df_codon.select('codon').to_numpy().T.tolist()[0],np.log2(df_codon.select('fold_change').to_numpy().T[0]))
        ax_codon.tick_params('x',rotation = 90)
        ax_codon.spines['right'].set_visible(False)
        ax_codon.spines['top'].set_visible(False)
        fig_codon.suptitle(s)

        #fig_codon.savefig(f'{s}_roots_codons.pdf')




def plot_fold_change_accross_samples(df):

    """
    Still working on these plots
    """


    fig,ax = plt.subplots(1)

    df_aa = df.sort('codon').filter((pl.col('codon').str.len_chars() == 1))
    scipy.stats.pearsonr(np.log2(df_aa.filter(pl.col('sample_id') == sample_ids[1]).select('fold_change').to_numpy().T[0]),np.log2(df_aa.filter(pl.col('sample_id') == sample_ids[2]).select('fold_change').to_numpy().T[0]))


