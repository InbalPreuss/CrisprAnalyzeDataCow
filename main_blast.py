import pandas as pd
from Bio import SearchIO
from Bio.Blast import NCBIXML
from matplotlib import pyplot as plt

from config import config


def parse_blast_results():
    uncommented = 'C:\\Users\\Inbal\\PycharmProjects\\CrisprAnalyzeDataCow\\Data\\results_table.txt'
    blast_result = SearchIO.parse(uncommented, 'blast-tab')

    dict_amount_of_matches_lengthx_to_each_grna = {}
    dict_amount_for_each_match_for_each_grna = {}
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

    df_amount_of_matches_lengthx_to_each_grna = dict_to_csv(dict_amount_of_matches_lengthx_to_each_grna, config[
        'blast_results'] + '/dict_amount_of_matches_lengthx_to_each_grna.csv')
    df_amount_for_each_match_for_each_grna = dict_to_csv(dict_amount_for_each_match_for_each_grna, config[
        'blast_results'] + '/dict_amount_for_each_match_for_each_grna.csv')
    return df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna


def dict_to_csv(dict, file_name):
    df = pd.DataFrame(dict)
    df.to_csv(file_name)
    return df


def df_to_csv(df, file_name):
    df.to_csv(file_name)
    return df


def amount_of_grna_for_each_match(df):
    df_sum_each_dist_each_grna = df.sum(axis=1)
    print(df_sum_each_dist_each_grna)
    plt.bar(list(df_sum_each_dist_each_grna.keys()), df_sum_each_dist_each_grna.values)
    plt.xlabel('match_num')
    plt.ylabel('# gRNA')
    plt.title('Amount of gRNA for each match')
    plt.savefig(config['blast_results'] + '/Amount of gRNA for each match')
    plt.close()


def for_each_match_for_each_grna(df):
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
        plt.savefig(config['blast_results'] + '/Amount for_each_match_for_each_grna_' + str(i_match) + '.svg')
        plt.close()
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
        plt.savefig(config['blast_results'] + '/Amount for each match for each grna sum from(' + str(
            df_col_names[0]) + ')to(' + str(i_match - 1) + ').svg')
        plt.close()

        df_to_csv(df_sum_each_dist_each_grna_concat, config[
        'blast_results'] + '/Amount for each match for each grna sum from(' + str(
            df_col_names[0]) + ')to(' + str(i_match - 1) + ').csv')


def histograms(df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna):
    amount_of_grna_for_each_match(df_amount_of_matches_lengthx_to_each_grna)
    for_each_match_for_each_grna(df_amount_for_each_match_for_each_grna)


if __name__ == '__main__':
    df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna = parse_blast_results()
    histograms(df_amount_of_matches_lengthx_to_each_grna, df_amount_for_each_match_for_each_grna)
