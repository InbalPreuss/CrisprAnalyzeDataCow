import math

from config import config
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Levenshtein import distance as lev


def open_folder():
    df_guides_3k_file_name = pd.read_excel(config['guides_3k_file_name'])
    return df_guides_3k_file_name


def calculate_dist_2_3(df_guides_3k_file_name, dist_calculation_name):
    dict_dist_2_3 = {}
    for dist in range(2, 4):
        dict_dist_2_3[dist] = {}

    for seq1_idx, seq1 in enumerate(df_guides_3k_file_name.seq):
        for seq2_idx, seq2 in enumerate(df_guides_3k_file_name.seq):
            if seq1_idx == seq2_idx:
                continue
            if 'hamming' == dist_calculation_name:
                distance_seq1_seq2 = hamming_distance(seq1, seq2)
            elif 'levenshtein' == dist_calculation_name:
                distance_seq1_seq2 = lev(seq1, seq2)
            if distance_seq1_seq2 == 2 or distance_seq1_seq2 == 3:
                seq1_name = seq1 + "," + df_guides_3k_file_name.Name[seq1_idx]
                seq2_name = seq2 + "," + df_guides_3k_file_name.Name[seq2_idx]
                if seq2_name in dict_dist_2_3[distance_seq1_seq2]:
                    continue
                if seq1_name in dict_dist_2_3[distance_seq1_seq2]:
                    dict_dist_2_3[distance_seq1_seq2][seq1_name].add(seq2_name)
                else:
                    dict_dist_2_3[distance_seq1_seq2][seq1_name] = {seq2_name}


    return dict_dist_2_3


def calculate_dist(df_guides_3k_file_name, dist_calculation_name):
    dict_dist = {}
    dict_dist_count = {}
    dict_dist_2_3 = {}
    for dist in range(0, 21):
        dict_dist[dist] = {}
        dict_dist_2_3[dist] = {}
        dict_dist_count[dist] = 0

    rows, cols = (df_guides_3k_file_name.shape[0], df_guides_3k_file_name.shape[0])
    hamming_distance_seq1_seq2_matrix_temp = [[0] * cols] * rows
    distance_seq1_seq2_matrix = np.copy(hamming_distance_seq1_seq2_matrix_temp)

    for seq1_idx, seq1 in enumerate(df_guides_3k_file_name.seq):
        for seq2_idx, seq2 in enumerate(df_guides_3k_file_name.seq):
            if seq1_idx == seq2_idx:
                continue
            if 'hamming' == dist_calculation_name:
                distance_seq1_seq2 = hamming_distance(seq1, seq2)
            elif 'levenshtein' == dist_calculation_name:
                distance_seq1_seq2 = lev(seq1, seq2)
            distance_seq1_seq2_matrix[seq1_idx][seq2_idx] = distance_seq1_seq2
            if seq1 in dict_dist[distance_seq1_seq2]:
                dict_dist[distance_seq1_seq2][seq1].add(seq2)
            else:
                dict_dist[distance_seq1_seq2][seq1] = {seq2}

            dict_dist_count[distance_seq1_seq2] = dict_dist_count[distance_seq1_seq2] + 1

    dict_dist_count_each_pair = {k: int(v / 2) for k, v in dict_dist_count.items()}
    print(f'dict_dist_count_each_pair: {dict_dist_count_each_pair}')
    return distance_seq1_seq2_matrix, dict_dist, dict_dist_count_each_pair


def hamming_distance(string1, string2):
    # Start with a distance of zero, and count up
    distance = 0
    # Loop over the indices of the string
    L = len(string1)
    for i in range(L):
        # Add 1 to the distance if these two characters are not equal
        if string1[i] != string2[i]:
            distance += 1
    # Return the final count of differences
    return distance


def matrix_to_csv(dist_matrix, df_guides_3k_file_name, dist_calculation_name):
    df = pd.DataFrame(dist_matrix)
    idx = 0
    new_col = df_guides_3k_file_name.seq.values
    df.insert(loc=idx, column='GuideName', value=new_col)
    if 'hamming' == dist_calculation_name:
        newfile = np.savetxt(config['hamming_dist_path'], df.to_numpy(), delimiter=',',
                             header='GuideName,' + ','.join(df_guides_3k_file_name.seq.values), fmt="%s")
    elif 'levenshtein' == dist_calculation_name:
        newfile = np.savetxt(config['levenshtein_dist_path'], df.to_numpy(), delimiter=',',
                         header='GuideName,' + ','.join(df_guides_3k_file_name.seq.values), fmt="%s")


def dict_dist_to_csv(dict_dist, dist_calculation_name, file_name):
    df = pd.DataFrame(dict_dist)
    df.to_csv(file_name)


def matrix_to_heatmap(dist_calculation_name):
    if 'hamming' == dist_calculation_name:
        df = pd.read_csv(config['hamming_dist_path'], index_col=0)
        plt.figure(figsize=(20, 15))
        plt.xlabel('gRNA', fontsize=30)
        plt.ylabel('gRNA', fontsize=30)
        plt.suptitle('3k gRNA - distance between each pair', fontsize=40)
        plt.imshow(df, cmap='hot', interpolation='nearest')
        cbar = plt.colorbar()
        cbar.set_label(label="red - low dist. yellow - high dist.", size=20)
        plt.savefig(config['heatmap_hamming_dist'])
        plt.close()
    elif 'levenshtein' == dist_calculation_name:
        df = pd.read_csv(config['levenshtein_dist_path'], index_col=0)
        plt.figure(figsize=(20, 15))
        plt.xlabel('gRNA', fontsize=30)
        plt.ylabel('gRNA', fontsize=30)
        plt.suptitle('3k gRNA - distance between each pair', fontsize=40)
        plt.imshow(df, cmap='hot', interpolation='nearest')
        cbar = plt.colorbar()
        cbar.set_label(label="red - low dist. yellow - high dist.", size=20)
        plt.savefig(config['heatmap_levenshtein_dist'])
        plt.close()


def hist_per_dist(dict_dist_count, dist_calculation_name):
    dict_count = {}

    plt.bar(dict_dist_count.keys(), dict_dist_count.values())
    plt.xlabel('distance')
    plt.ylabel('amount of sequences')
    if 'hamming' == dist_calculation_name:
        plt.title('Hamming')
        plt.suptitle('Distance between each pair 3k gRNA')
        plt.savefig(config['hist_amount_seq_per_dist_hamming'])
        plt.close()
    elif 'levenshtein' == dist_calculation_name:
        plt.title('levenshtein')
        plt.suptitle('Distance between each pair 3k gRNA')
        plt.savefig(config['hist_amount_seq_per_dist_levenshtein'])
        plt.close()


if __name__ == '__main__':
    # df_guides_3k_file_name = open_folder()
    # Hamming
    # hamming_dist_matrix, dict_per_dist_hamming, dict_dist_count_hamming = calculate_dist(df_guides_3k_file_name, 'hamming')
    # dict_dist_2_3_hamming = calculate_dist_2_3(df_guides_3k_file_name, 'hamming')
    # dict_dist_to_csv(dict_dist_2_3_hamming, 'hamming', config['hamming_dict_dist_2_3_path'])
    # matrix_to_csv(hamming_dist_matrix, df_guides_3k_file_name, 'hamming')
    # matrix_to_heatmap('hamming')
    # dict_dist_to_csv(dict_per_dist_hamming, 'hamming', config['hamming_per_dist_path'])
    # hist_per_dist(dict_dist_count_hamming, 'hamming')


    # Levenshtein
    # levenshtein_dist_matrix, dict_per_dist_levenshtein, dict_dist_count_levenshtein = calculate_dist(df_guides_3k_file_name, 'levenshtein')
    # dict_dist_2_3_levenshtein = calculate_dist_2_3(df_guides_3k_file_name, 'levenshtein')
    # dict_dist_to_csv(dict_dist_2_3_levenshtein, 'levenshtein', config['levenshtein_dict_dist_2_3_path'])
    # matrix_to_csv(levenshtein_dist_matrix, df_guides_3k_file_name, 'levenshtein')
    # matrix_to_heatmap('levenshtein')
    # dict_dist_to_csv(dict_per_dist_levenshtein, 'levenshtein', config['levenshtein_per_dist_path'])
    # hist_per_dist(dict_dist_count_levenshtein, 'levenshtein')
    
    print(4)
