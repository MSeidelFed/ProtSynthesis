import argparse
import glob

import pandas as pd
import numpy as np

import isocor
import count_glyser



def build_isotope_lookup(isotopes_path):

    '''
    Function to generated isotope dictionary from isotope.dat file

    Dictionary structure 
    '''

    isotopes = {}

    with open(isotopes_path, 'r') as f:
        for i,l in enumerate(f.readlines()):
            l = l.strip('\n')
            if i == 0:
                headers = l.split(',')
            else:
                element, mass, abundance = l.split(",")
                if element not in isotopes:
                    isotopes[element] = {"mass":[], "abundance": []}
                isotopes[element]["mass"].append(float(mass))
                isotopes[element]["abundance"].append(float(abundance))

    return isotopes


def consolidate_data(file_list):

    '''
    Combines raw_outputs(.tsv) to one output where file directory is used for sample name
    '''

    
    for i,f in enumerate(file_list):
        sample_name = f.split('/')[-2]

        tmpdf = pd.read_csv(f, sep = '\t')
        tmpdf.loc[:, 'sample'] = sample_name

        if i == 0:
            df = tmpdf
        else:
            df = pd.concat([df, tmpdf])

    df.columns = [c.replace(' ','_').lower() for c in df.columns]
    df = df.fillna(0)

    return df


def make_inputs_to_isocorr(df, isotopes, tracer_mass = 15, tracer_element = "N", tracer_purity = [0.01,0.99]):
    '''
    Parses raw dataframe to run isocorr
    returns dataframe containing the following columns [sample, protein_key, protein_group, mean_enrichment]
    '''

    tracer = f'{tracer_mass}{tracer_element}'

    metabolties = {}
    measurements = []

    mean_enrichment = []

    intensity_columns = sorted([c for c in df.columns if 'intensity' in c],key=lambda a: int(a.split('_')[1]))

    df = df[df.loc[:,intensity_columns].sum(axis = 1) != 0]


    for i in range(df.shape[0]):
        row = df.iloc[i,:]
        peptide = row.loc['peptide_key']

        if peptide not in metabolties:
            metabolties[peptide] = isocor.LowResMetaboliteCorrector(formula=row.loc['formula'].replace(' ',''), tracer= tracer, label = peptide, tracer_purity = tracer_purity,  correct_NA_tracer = True, data_isotopes = isotopes, derivative_formula = None)

        number_of_peaks = metabolties[peptide].formula[tracer_element] + 1

        if number_of_peaks <= len(intensity_columns):
            peaks = row[intensity_columns[:number_of_peaks]]
        else:
            extension = [0] * (number_of_peaks-len(intensity_columns))
            peaks = row[intensity_columns[:number_of_peaks]].to_list() + extension

        _,_,_, menrich = metabolties[peptide].correct(peaks)
        mean_enrichment.append(menrich)

    df_isocorr = df.loc[:,['file','peptide_key','protein_gorup']]
    
    df_isocorr.loc[:,'mean_enrichment'] = mean_enrichment
    
    return df_isocorr



def main():

    parser = argparse.ArgumentParser(description= "Python wrapper around isocor")
    parser.add_argument("-f", type= str, help="Path to results.dat")
    parser.add_argument("-i", type = str, help='isotopes.dat path')
    parser.add_argument("-o", type = str, help="Output path")
    parser.add_argument("-s", type=str, help = "path to schema" )
    args = parser.parse_args()

    results_path = args.f
    raw_data = pd.read_csv(results_path, sep = '\t')
    raw_data.columns = [c.replace(' ','_').lower() for c in raw_data.columns]
    raw_data = raw_data.fillna(0)

    output_path = args.o
    isotopes_path = args.i
    schema_path = args.s

    schema_path = 'data/total_proteome_labeled/2026_schema_shoot.csv'
    schema_df = pd.read_csv(schema_path)

    sample_lookup = {file:f'{sample}_{file.split('_')[-1]}' for _, sample, file in schema_df.loc[:,['Sample','File']].itertuples()}
    

    isotopes = build_isotope_lookup(isotopes_path)

    isocorrected_data = make_inputs_to_isocorr(raw_data, isotopes= isotopes)
    isocorrected_data.loc[:,'sample'] = isocorrected_data.loc[:,'file'].map(sample_lookup)
    #isocorrected_data = isocorrected_data[isocorrected_data.loc[:,"mean_enrichment"] != 0]

    isocorrected_data.loc[:,'gly_ser_count'] = isocorrected_data.loc[:,'peptide_key'].apply(count_glyser.count_gly_ser)
    #isocorrected_data = isocorrected_data[isocorrected_data.loc[:,"gly_ser_count"] != 0]

    isocorrected_data.loc[:,'normalized_mean_enrichment'] = (isocorrected_data.loc[:,'mean_enrichment'] / isocorrected_data.loc[:,'gly_ser_count'].replace(0,np.nan)).replace(np.nan, 0)
    isocorrected_data.to_csv(output_path)


if __name__ == "__main__":
    main()
