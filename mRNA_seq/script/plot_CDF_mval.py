from typing import List
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style("ticks")


def flatten_list(input_list: List[list]):
    return [item for sublist in input_list for item in sublist]


if __name__ == "__main__":
    histone_types = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
    histone_col = sns.color_palette("Set2", 6)
    # histone_col_light = sns.color_palette("hls", 8)
    path = '/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_annotated/'
    # path = '/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/majiq_res/correl/temp/'

    fig = plt.figure()
    ax = fig.add_subplot(111)

    all_his_data = []
    for i, his in enumerate(histone_types):
        all_files = glob.glob(path + his + "/*.bed")
        his_data = []
        for file in all_files:
            print(file)
            data = pd.read_csv(file, header=None, sep=' ').iloc[:, 0].tolist()
            his_data.append(data)
            sns.ecdfplot(data=data, color=histone_col[i], alpha=0.008, linewidth=0.75)
        all_his_data.append(flatten_list(his_data))

    no_nonzero = [-1*(np.count_nonzero(his_data)) for his_data in all_his_data]
    new_idx = np.argsort(no_nonzero)
    sorted_all_his_data = [all_his_data[i] for i in new_idx]
    sorted_histone_types = [histone_types[i] for i in new_idx]

    for i, his_data in enumerate(sorted_all_his_data):
        sns.ecdfplot(data=his_data, color=histone_col[i], alpha=0.9, label=i, linewidth=2)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("M-value", fontsize=16)
    plt.ylabel("Cumulative proportion", fontsize=16)

    h, l = ax.get_legend_handles_labels()
    plt.legend(handles=h, labels=sorted_histone_types, loc='lower right', title="Histone type")
    plt.savefig('mval_cdf2.tiff', dpi=300)
