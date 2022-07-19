import utilities

import collections
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from BCBio.GFF import (GFFExaminer, GFFParser, DiscoGFFParser)
from matplotlib import pyplot as plt
from config import config


def gff_blast_parse():
    dict_name_n_location_rec = {}

    parser = GFFParser()
    a = parser.parse(config['gff_file'])
    rec_dict = SeqIO.to_dict(a)
    for name_rec, rec in rec_dict.items():
        rec_features = rec.features
        for rec in rec_features:
            if rec.type == 'gene':
                print(rec.location)
                if name_rec not in dict_name_n_location_rec.keys():
                    dict_name_n_location_rec[name_rec] = []
                dict_name_n_location_rec[name_rec].append([rec.location.nofuzzy_start, rec.location.nofuzzy_end])
    return dict_name_n_location_rec


def parse_blast_results():
    uncommented = 'C:\\Users\\Inbal\\PycharmProjects\\CrisprAnalyzeDataCow\\Data\\results_table.txt'
    blast_result = SearchIO.parse(uncommented, 'blast-tab')

    dict_amount_of_matches_lengthx_to_each_grna = {}
    dict_amount_for_each_match_for_each_grna = {}
    dict_amount_for_each_match_for_each_grna_hit_id_n_range = {}
    for all_gene_hits in blast_result:
        dict_amount_of_matches_lengthx_to_each_grna[all_gene_hits.id] = {}
        for gene_i in all_gene_hits.hsps:
            gene_i_hit = gene_i.aln_span - gene_i.mismatch_num
            if gene_i_hit not in dict_amount_for_each_match_for_each_grna:
                dict_amount_for_each_match_for_each_grna[gene_i_hit] = {}
            try:
                dict_amount_of_matches_lengthx_to_each_grna[all_gene_hits.id][gene_i_hit] += 1
            except:
                dict_amount_of_matches_lengthx_to_each_grna[all_gene_hits.id][gene_i_hit] = 1
            try:
                dict_amount_for_each_match_for_each_grna[gene_i_hit][all_gene_hits.id] += 1
            except:
                dict_amount_for_each_match_for_each_grna[gene_i_hit][all_gene_hits.id] = 1
            if gene_i.hit_id not in dict_amount_for_each_match_for_each_grna_hit_id_n_range.keys():
                dict_amount_for_each_match_for_each_grna_hit_id_n_range[gene_i.hit_id] = []
            dict_amount_for_each_match_for_each_grna_hit_id_n_range[gene_i.hit_id].append(
                [gene_i.hit_start, gene_i.hit_end])

    dict_amount_for_each_match_for_each_grna = collections.OrderedDict(
        sorted(dict_amount_for_each_match_for_each_grna.items()))

    df_amount_of_matches_lengthx_to_each_grna = utilities.dict_to_csv(dict_amount_of_matches_lengthx_to_each_grna, config[
        'blast_results'] + '/dict_amount_of_matches_lengthx_to_each_grna.csv')
    df_amount_for_each_match_for_each_grna = utilities.dict_to_csv(dict_amount_for_each_match_for_each_grna, config[
        'blast_results'] + '/dict_amount_for_each_match_for_each_grna.csv')

    df_amount_of_matches_lengthx_to_each_grna = df_amount_of_matches_lengthx_to_each_grna.fillna(0)
    df_amount_for_each_match_for_each_grna = df_amount_for_each_match_for_each_grna.fillna(0)
    return df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna, dict_amount_for_each_match_for_each_grna_hit_id_n_range


def amount_of_grna_for_each_match(df):
    df_sum_each_dist_each_grna = df.sum(axis=1)
    plt.bar(list(df_sum_each_dist_each_grna.keys()), df_sum_each_dist_each_grna.values)
    plt.xlabel('Similarity rate')
    plt.ylabel('# gRNA')
    plt.title('Amount of gRNA matches to each similarity rate in the cow genome')
    plt.savefig(config['blast_results'] + '/Amount of gRNA matches to each similarity rate in the cow genome')
    plt.close()


def sum_amount_for_each_match_for_each_grna_sum_from_20_to_i(df):
    results_dir_sum_name = utilities.create_dir(
        config['blast_results'] + config['sum_amount_for_each_match_for_each_grna_sum_from_20_to_i'])

    df_col_names = list(df.columns)
    df_col_names.sort(reverse=True)

    df_sum_each_dist_each_grna_concat = df[df_col_names[0]]
    for i_match_idx, i_match in enumerate(df_col_names):
        if len(df_col_names) == i_match_idx + 1:
            break
        df_sum_each_dist_each_grna_temp = df[df_col_names[i_match_idx + 1]]
        df_sum_each_dist_each_grna_concat = pd.concat(
            (df_sum_each_dist_each_grna_concat, df_sum_each_dist_each_grna_temp), axis=1).sum(axis=1)
        # plt.bar(list(df_sum_each_dist_each_grna.keys()), df_sum_each_dist_each_grna.values)
        plt.bar(range(df_sum_each_dist_each_grna_concat.keys().__len__()), df_sum_each_dist_each_grna_concat.values)
        plt.ylabel('amount of matches')
        plt.xlabel('gRNA')
        plt.title('Amount of gRNA for each match sum from(' + str(df_col_names[0]) + ')to(' + str(i_match - 1) + ')')
        plt.savefig(results_dir_sum_name + '/Sum amount for each match for each grna sum from(' + str(
            df_col_names[0]) + ')to(' + str(i_match - 1) + ').svg')
        plt.close()

        utilities.df_to_csv(df_sum_each_dist_each_grna_concat, results_dir_sum_name +
                  '/Sum amount for each match for each grna sum from(' + str(df_col_names[0]) + ')to(' + str(
            i_match - 1) + ').csv')


def for_each_match_for_each_grna(df):
    results_dir_name = utilities.create_dir(config['blast_results'] + config['for_each_match_for_each_grna'])

    df_col_names = list(df.columns)
    df_col_names.sort(reverse=True)
    for i_match in df_col_names:
        df_sum_each_dist_each_grna = df[i_match]
        # plt.bar(list(df_sum_each_dist_each_grna.keys()), df_sum_each_dist_each_grna.values)
        plt.figure(figsize=(15, 15))
        plt.bar(range(df_sum_each_dist_each_grna.keys().__len__()), df_sum_each_dist_each_grna.values)
        plt.ylabel('amount of matches')
        plt.xlabel('gRNA')
        plt.title('Amount of gRNA for each match ' + str(i_match))
        plt.savefig(results_dir_name + '/for_each_match_for_each_grna_' + str(i_match) + '.svg')
        plt.close()


def groups_of_amount_of_guides_per_match(df_amount_for_each_match_for_each_grna):
    results_dir = utilities.create_dir(config['blast_results'] + config['sum_amount_of_gRNA_with_match_in_bin_x_i'])

    df_amount_for_each_match_for_each_grna_columns = df_amount_for_each_match_for_each_grna.columns

    for columns_name in df_amount_for_each_match_for_each_grna_columns:
        dict_amount_for_each_match_for_each_grna = {}
        array_amount_for_each_match_for_each_grna_columns_name_values = df_amount_for_each_match_for_each_grna[
            columns_name].values
        for count in array_amount_for_each_match_for_each_grna_columns_name_values:
            try:
                dict_amount_for_each_match_for_each_grna[count] += 1
            except:
                dict_amount_for_each_match_for_each_grna[count] = 1

        dict_amount_for_each_match_for_each_grna = collections.OrderedDict(
            sorted(dict_amount_for_each_match_for_each_grna.items()))

        plt.bar(list(dict_amount_for_each_match_for_each_grna.keys()),
                list(dict_amount_for_each_match_for_each_grna.values()), width=0.3)
        plt.ylabel('Amount of gRNA')
        plt.xlabel('Sum of gRNA with ' + str(columns_name) + 'nuc match')
        plt.title('Amount of gRNA with (' + str(columns_name) + ') nuc match')

        plt.savefig(results_dir + '/Amount of gRNA with (' + str(columns_name) + ') nuc match.svg')
        plt.close()

        utilities.df_to_csv(df_amount_for_each_match_for_each_grna[columns_name], results_dir
                  + '/Amount of gRNA with (' + str(columns_name) + ') nuc match.csv')


def histograms(df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna):
    amount_of_grna_for_each_match(df_amount_of_matches_lengthx_to_each_grna)

    # Drop gRNA that has too many matches, so the histograms will be clearer
    df_amount_for_each_match_for_each_grna_after_drop_grna = df_amount_for_each_match_for_each_grna.drop(
        config['drop_grna_that_has_too_many_matches_in_cow_genome'])

    for_each_match_for_each_grna(df_amount_for_each_match_for_each_grna_after_drop_grna)
    sum_amount_for_each_match_for_each_grna_sum_from_20_to_i(df_amount_for_each_match_for_each_grna_after_drop_grna)
    groups_of_amount_of_guides_per_match(df_amount_for_each_match_for_each_grna_after_drop_grna)

    # Dropped gRNA that has too many matches, write them in csv
    df_amount_for_each_match_for_each_grna_too_many_matches = df_amount_for_each_match_for_each_grna.loc[
        config['drop_grna_that_has_too_many_matches_in_cow_genome']]
    utilities.df_to_csv(df_amount_for_each_match_for_each_grna_too_many_matches,
              config['blast_results'] + config['grna_that_has_too_many_matches_in_cow_genome'])



def compare_cow_genome_to_grna_locations(dict_gff_cow_genome_name_n_location_rec,
                                         dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range):
    print('compare_cow_genome_to_grna_locations')

    # Array with all the locations blast and gff are
    location_in_blast_and_not_in_gff_header = ['name']
    location_in_blast_and_not_in_gff_data = []
    file_path_location_in_blast_and_not_in_gff = config['blast_results'] + '/location_in_blast_and_not_in_gff.csv'

    location_in_blast_and_in_gff_header = ['name', 'location_blast', 'location_gff']
    location_in_blast_and_in_gff_data = []
    file_path_location_in_blast_and_in_gff = config['blast_results'] + '/location_in_blast_and_in_gff.csv'

    dict_grna_in_cow_genome = {}

    set_gff_cow_genome_name_n_location_rec = set(dict_gff_cow_genome_name_n_location_rec)
    set_blast_amount_for_each_match_for_each_grna_hit_id_n_range = set(
        dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range)

    intersection = set_gff_cow_genome_name_n_location_rec.intersection(
        set_blast_amount_for_each_match_for_each_grna_hit_id_n_range)
    for name_idx, name in enumerate(intersection):
        is_in_blast = False
        for location_blast in dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range[name]:
            for location_gff in dict_gff_cow_genome_name_n_location_rec[name]:
                if (location_blast[0] >= location_gff[0]) and (location_blast[1] <= location_gff[1]):
                    print(
                        f'{name}, ({location_blast[0]} >= {location_gff[0]}) and ({location_blast[1]} <= {location_gff[1]})')
                    dict_grna_in_cow_genome[name] = (location_blast, location_gff)
                    location_in_blast_and_in_gff_data.append([name, location_blast, location_gff])
                    is_in_blast = True
                    break
        if is_in_blast:
            location_in_blast_and_not_in_gff_data.append([name])

    utilities.write_rows_to_csv(file_path_location_in_blast_and_not_in_gff, location_in_blast_and_not_in_gff_header,
                      location_in_blast_and_not_in_gff_data)
    utilities.write_rows_to_csv(file_path_location_in_blast_and_in_gff, location_in_blast_and_in_gff_header,
                      location_in_blast_and_in_gff_data)

    return dict_grna_in_cow_genome


if __name__ == '__main__':
    print('parse_blast_results')
    df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna, dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range = parse_blast_results()
    print('histograms')
    histograms(df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna)
    print('gff_blast_parse')
    dict_gff_cow_genome_name_n_location_rec = gff_blast_parse()
    print('compare_cow_genome_to_grna_locations')
    compare_cow_genome_to_grna_locations(dict_gff_cow_genome_name_n_location_rec,
                                         dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range)

    print('Finished')
