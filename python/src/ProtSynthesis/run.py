import argparse
import json
import polars as pl
import toml 
import os

import src.ProtSynthesis.run_isocorr as run_isocorr
import src.ProtSynthesis.make_enrichment_groups as make_enrichment_groups
import src.ProtSynthesis.count_cds as count_cds
import src.ProtSynthesis.calculate_enrichment as calculate_enrichment




def main():

    parser = argparse.ArgumentParser(description= "Python wrapper around isocor")
    parser.add_argument("-p", type= str, help="Path to params.toml")
    args = parser.parse_args()

    config = toml.load(open(args.p,'r'))
    config = toml.load(open('misc/params.toml','r'))

    config_paths = config['paths']
    config_params = config['parameters']


    results_path = config_paths['results_path']
    output_path = config_paths['output_path']
    isotope_path = config_paths['isotope_path']
    schema_path = config_paths['schema_path']
    fasta_path = config_paths['fasta_path']
    idmap_path = config_paths['idmap_path']


    experiment_id = config_params.pop('experiment_id')

    figure_path = os.path.join(output_path,'figures/')
    data_path = os.path.join(output_path,'processed_data/')

    if not os.path.exists(figure_path):
        os.makedirs(figure_path)

    if not os.path.exists(data_path):
        os.mkdir(data_path)


    df = run_isocorr.run(results_path=results_path, isotopes_path=isotope_path, schema_path=schema_path, data_path=data_path, experimental_id = experiment_id)

    df = pl.read_csv('output/total_proteome/roots/processed_data/roots_2026_isocorrected_data.csv')


    df, enrichment_groups = make_enrichment_groups.run(df=df, idmap_path= idmap_path, figure_path=figure_path, params=config_params)
    enrichment_group_path = os.path.join(data_path,f'{experiment_id}_enrichment_groups.json')
    json.dump(enrichment_groups, open(enrichment_group_path,'w'))


    df_gene = count_cds.run(fasta_path=fasta_path)
    calculate_enrichment.calculate_enrichment(gene_enrichment_groups=enrichment_groups, df_gene= df_gene)





if __name__ == "__main__":
    main()