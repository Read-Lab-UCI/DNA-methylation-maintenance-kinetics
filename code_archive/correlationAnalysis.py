# -*- coding:utf-8 -*-
import math,os,collections
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from localreg import *

sns.set(color_codes=True)
BASE_DIR = '..'
EACH_SUB_FIG_SIZE = 5
D_MAX = 1000
is_inter_with_other_cpg = True
n_category = 3
def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def read_bed_file_and_store_pos_to_a_struct(bed_tsv_file_path, K_RD, spearman=False):
    """
    Read tsv format bed file (No chr prefix at 1st colum), and store the chr, pos, k or methylation level info into dict
    :param bed_tsv_file_path: tsv format bed file input path
    :return: the dict with chr, pos, k or methylation level info
    """
    print("reading %s" % bed_tsv_file_path)
    df = pd.read_csv(bed_tsv_file_path, sep='\t', header=None).values
    chrs = df[:, 0]
    unique_chr, _ = np.unique(chrs, return_counts=True)
    dict_to_store = {chr_i: {} for chr_i in unique_chr}
    val_col = 4 if K_RD else 3
    if spearman:
        df_vals = df[:, val_col]
        bins = [np.percentile(df_vals, i * (100. / n_category)) for i in range(n_category + 1)]
        vals = np.digitize(df_vals, bins, right=True)
        print(df_vals)
        print(bins)
        print(vals)
    else:
        vals = df[:, val_col]
    for idx, item in enumerate(df):
        chr_i, pos, val = item[0], item[1], vals[idx]
        dict_to_store[chr_i][int(pos)] = val
    return [dict_to_store]


def filter_d_length_to_generate_CpG_pairs(Chr_CpG_pos_and_methy_dict, D_MAX):
    """
    Extract the CpG pairs with distance of d, and there can be CpGs between CpG pairs with d distance
    :param Chr_CpG_pos_and_methy_dict: The dict with chr, pos, k or methylation level info provided
    :param d: the distance that required for CpG pairs
    :return: An array that stores all CpG pairs with distance of d
    """
    CpG_Pairs_Dict = {}
    D_RANGE = range(2, D_MAX)
    Chr_CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[0]
    for d in D_RANGE:
        if d % 10 == 0:
            print("Calculating CpG pairs with %d distance" % d)
        for chr_i in Chr_CpG_pos_and_methy_dict.keys():
            CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[chr_i]
            for key in CpG_pos_and_methy_dict.keys():
                pos = key
                methy_level_1 = CpG_pos_and_methy_dict[key]
                pos_off_by_d = pos + d
                if pos_off_by_d in CpG_pos_and_methy_dict.keys():
                    methy_level_2 = CpG_pos_and_methy_dict[pos_off_by_d]
                    methy_levels_of_CpG_pairs =[methy_level_1, methy_level_2]
                    if d not in CpG_Pairs_Dict.keys():
                        CpG_Pairs_Dict[d] = []
                    CpG_Pairs_Dict[d].append(methy_levels_of_CpG_pairs)
    return [CpG_Pairs_Dict]

def filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(Chr_CpG_pos_and_methy_dict, D_MAX):
    """
    Extract the CpG pairs with distance of d, and required no CpGs between CpG pairs with d distance
    :param d: the distance that required for CpG pairs
    :param length_of_od:
    :param od_keys:
    :param od_vals:
    :return:
    """
    CpG_Pairs_list = {}
    Chr_CpG_pos_and_methy_dict=Chr_CpG_pos_and_methy_dict[0]
    for chr_i in Chr_CpG_pos_and_methy_dict.keys():
        CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[chr_i]
        sorted_struct = collections.OrderedDict(sorted(CpG_pos_and_methy_dict.items()))
        od_keys = list(sorted_struct.keys())
        od_vals = list(sorted_struct.values())
        for i in range(0, len(sorted_struct) - 1):
            # 后一个位点与前一个位点之间的距离为d
            pre_od_key = od_keys[i]
            post_od_key = od_keys[i + 1]
            dist = post_od_key - pre_od_key
            if 0 < dist <D_MAX:
                methy_levels_of_CpG_pairs = [od_vals[i], od_vals[i + 1]]
                if dist not in CpG_Pairs_list.keys():
                    CpG_Pairs_list[dist] = []
                CpG_Pairs_list[dist].append(methy_levels_of_CpG_pairs)
    return [CpG_Pairs_list]

def calc_C_d_by_correlation(CpG_Pairs_Dict, spearman = False):
    """
    calc the Pearson correlation of CpG pairs list
    :param CpG_pairs: The input CpG pairs list
    :return: The dict with {d: correlation}
    """
    RD_dict = {}
    CpG_Pairs_Dict = CpG_Pairs_Dict[0]
    for DISTANCE, CpG_pairs in CpG_Pairs_Dict.items():
        if len(CpG_pairs) > 5:
            CpG_arr = np.array(CpG_pairs)
            if spearman:
                r_d, _ = spearmanr(CpG_arr[:, 0], CpG_arr[:, 1])
            else:
                r_d = np.corrcoef(CpG_arr[:, 0], CpG_arr[:, 1])[0][1]
            RD_dict[DISTANCE] = r_d
    return [RD_dict] if RD_dict else []

def write_RD_into(RD_dict, out_fp):
    if RD_dict:
        RD_dict = RD_dict[0]
        sorted_struct = collections.OrderedDict(sorted(RD_dict.items()))
        od_keys = list(sorted_struct.keys())
        ltw = "\n".join(["%d\t%.3f" %(d, RD_dict[d]) for d in od_keys]) + "\n"
        with open(out_fp, 'w') as RD_file:
            RD_file.write(ltw)
def calc_correlation(bed_tsv_file_path, out_R_d_correlation_path,  d_max, is_inter_with_other_cpg, K_RD=1, spearman= False):
    [Chr_CpG_pos_and_methy_dict] = read_bed_file_and_store_pos_to_a_struct(bed_tsv_file_path, K_RD, spearman)
    if is_inter_with_other_cpg:
        CpG_pairs = filter_d_length_to_generate_CpG_pairs([Chr_CpG_pos_and_methy_dict], d_max)
    else:
        CpG_pairs = filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg([Chr_CpG_pos_and_methy_dict], d_max)
    RD_dict = calc_C_d_by_correlation(CpG_pairs, spearman)
    write_RD_into(RD_dict, out_R_d_correlation_path)

def plot_local_regression_and_RD_with_path(fig_dir, RD_fp, label, color, spearman = False):
    mkdir(fig_dir)
    fig, ax = plt.subplots(1, 1, figsize=(1 * EACH_SUB_FIG_SIZE, 1 * EACH_SUB_FIG_SIZE))
    RD_df = pd.read_csv(RD_fp, sep="\t", header=None).values
    x = RD_df[:, 0]
    y = RD_df[:, 1]
    y2 = localreg(x, y, degree=2, kernel=tricube, width=100)
    ax.scatter(x, y, s=8, color= color, label=label)
    ax.plot(x, y2, "k-", linewidth=2)
    ax.set_xticks(range(0, D_MAX + 1, 200))
    ax.set_xlim(0, D_MAX)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Genomic Distance(bp)")
    ax.set_ylabel("%s Correlation" % "Pearson" if not spearman else "Spearman")
    plt.savefig(os.path.join(fig_dir, "correlation_function.png"), dpi=300)


if __name__ == "__main__":
    bed_file_path = os.path.join(BASE_DIR, "data", "MyData", "METHY_AND_K.bed")
    out_Rd_corr_fp = os.path.join(BASE_DIR, "data", "MyData", "categorical_RD.bed")
    K_RD = 1
    spearman = True
    calc_correlation(bed_file_path, out_Rd_corr_fp, D_MAX, is_inter_with_other_cpg, K_RD=K_RD, spearman=spearman)
    fig_dir = os.path.join(BASE_DIR, "figures")
    plot_local_regression_and_RD_with_path(fig_dir, out_Rd_corr_fp, label="categorical R(d)", color="blue")
