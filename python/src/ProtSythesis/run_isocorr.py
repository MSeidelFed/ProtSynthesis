import argparse
import glob

import pandas as pd
import numpy as np

import isocor



def build_isotope_lookup():

    '''
    Function to generated isotope dictionary from isotope.dat file

    Dictionary structure 
    '''

    isotopes = {}

    with open("Isotopes.dat", 'r') as f:
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


def make_inputs_to_isocorr(df, tracer_mass = 15, tracer_element = "N", tracer_purity = [0.01,0.99]):
    '''
    Parses raw dataframe to run isocorr
    returns dataframe containing the following columns [sample, protein_key, protein_group, mean_enrichment]
    '''

    tracer = f'{tracer_mass}{tracer_element}'

    isotopes = build_isotope_lookup()
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

    df_isocorr = df.loc[:,['sample','peptide_key','protein_gorup']]
    
    df_isocorr.columns = ['sample','peptide_key', 'protein_group']

    df_isocorr.loc[:,'mean_enrichment'] = mean_enrichment
    
    return df_isocorr



def main():

    parser = argparse.ArgumentParser(description= "Python wrapper around isocor")
    parser.add_argument("-f", type= str, help="Path to files with wild card")
    parser.add_argument("-o", type = str, help="Output path")
    args = parser.parse_args()

    file_list = glob.glob(args.f)

    raw_data = consolidate_data(file_list = file_list)
    isocorrected_data = make_inputs_to_isocorr(raw_data)


if __name__ == "__main__":
    main()
