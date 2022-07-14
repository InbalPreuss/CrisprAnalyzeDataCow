import collections
import os

import numpy as np
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from BCBio import GFF
from BCBio.GFF import (GFFExaminer, GFFParser, DiscoGFFParser)
from matplotlib import pyplot as plt
from config import config


def dict_to_csv(dict, file_name):
    df = pd.DataFrame(dict)
    df.to_csv(file_name)
    return df


def df_to_csv(df, file_name):
    df.to_csv(file_name)
    return df


def create_dir(dir_name):
    results_dir = dir_name
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    return dir_name


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

    df_amount_of_matches_lengthx_to_each_grna = dict_to_csv(dict_amount_of_matches_lengthx_to_each_grna, config[
        'blast_results'] + '/dict_amount_of_matches_lengthx_to_each_grna.csv')
    df_amount_for_each_match_for_each_grna = dict_to_csv(dict_amount_for_each_match_for_each_grna, config[
        'blast_results'] + '/dict_amount_for_each_match_for_each_grna.csv')
    return df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna, dict_amount_for_each_match_for_each_grna_hit_id_n_range


def amount_of_grna_for_each_match(df):
    df_sum_each_dist_each_grna = df.sum(axis=1)
    print(df_sum_each_dist_each_grna)
    plt.bar(list(df_sum_each_dist_each_grna.keys()), df_sum_each_dist_each_grna.values)
    plt.xlabel('match_num')
    plt.ylabel('# gRNA')
    plt.title('Amount of gRNA for each match')
    plt.savefig(config['blast_results'] + '/Amount of gRNA for each match')
    plt.close()


def sum_amount_for_each_match_for_each_grna_sum_from_20_to_i(df):
    results_dir_sum_name = create_dir(
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

        df_to_csv(df_sum_each_dist_each_grna_concat, results_dir_sum_name +
                  '/Sum amount for each match for each grna sum from(' + str(df_col_names[0]) + ')to(' + str(
            i_match - 1) + ').csv')


def for_each_match_for_each_grna(df):
    results_dir_name = create_dir(config['blast_results'] + config['for_each_match_for_each_grna'])

    df_col_names = list(df.columns)
    df_col_names.sort(reverse=True)
    df = df.fillna(0)
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
    results_dir = create_dir(config['blast_results'] + config['sum_amount_of_gRNA_with_match_in_bin_x_i'])

    df_amount_for_each_match_for_each_grna = df_amount_for_each_match_for_each_grna.fillna(0)
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
        # dict_amount_for_each_match_for_each_grna = {str(int(key)): value for key, value in dict_amount_for_each_match_for_each_grna.items()}

        plt.bar(list(dict_amount_for_each_match_for_each_grna.keys()),
                list(dict_amount_for_each_match_for_each_grna.values()), width=0.8)
        plt.ylabel('Amount of gRNA')
        plt.xlabel('Sum of gRNA with value ' + str(columns_name))

        # font_size = int(300 / np.power(dict_amount_for_each_match_for_each_grna.__len__(),(4/3)))
        # plt.xticks(fontsize=font_size)

        # plt.title('Amount of gRNA for each match sum from(' + str(df_col_names[0]) + ')to(' + str(i_match - 1) + ')')
        plt.savefig(results_dir + '/Amount of gRNA with match in bin x_i(' + str(columns_name) + ').svg')
        plt.close()

        df_to_csv(df_amount_for_each_match_for_each_grna[columns_name], results_dir
                  + '/Amount of gRNA with match in bin x_i(' + str(columns_name) + ').csv')


def histograms(df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna):
    amount_of_grna_for_each_match(df_amount_of_matches_lengthx_to_each_grna)
    for_each_match_for_each_grna(df_amount_for_each_match_for_each_grna)
    sum_amount_for_each_match_for_each_grna_sum_from_20_to_i(df_amount_for_each_match_for_each_grna)
    groups_of_amount_of_guides_per_match(df_amount_for_each_match_for_each_grna)


def compare_cow_genome_to_grna_locations(dict_gff_cow_genome_name_n_location_rec,
                                         dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range):
    print('compare_cow_genome_to_grna_locations')
    dict_grna_in_cow_genome = {}

    set_gff_cow_genome_name_n_location_rec = set(dict_gff_cow_genome_name_n_location_rec)
    set_blast_amount_for_each_match_for_each_grna_hit_id_n_range = set(
        dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range)

    for name in set_gff_cow_genome_name_n_location_rec.intersection(
            set_blast_amount_for_each_match_for_each_grna_hit_id_n_range):
        for location_blast in dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range[name]:
            for location_gff in dict_gff_cow_genome_name_n_location_rec[name]:
                if (location_blast[0] >= location_gff[0]) and (location_blast[0] <= location_gff[0]):
                    print(
                        f'{name}, ({location_blast[0]} >= {location_gff[0]}) and ({location_blast[0]} <= {location_gff[0]})')
                    dict_grna_in_cow_genome[name] = (location_blast, location_gff)

    return dict_grna_in_cow_genome


if __name__ == '__main__':
    df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna, dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range = parse_blast_results()
    histograms(df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna)
    print(5)
    dict_gff_cow_genome_name_n_location_rec = gff_blast_parse()
    compare_cow_genome_to_grna_locations(dict_gff_cow_genome_name_n_location_rec,
                                         dict_blast_amount_for_each_match_for_each_grna_hit_id_n_range)
