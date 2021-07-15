from typing import List
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from statannot import add_stat_annotation

sns.set_style("ticks")


def flatten_list(input_list: List[list]):
    return [item for sublist in input_list for item in sublist]


def add_blank_subplot(ax):
    ax.tick_params(axis='both', colors='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['left'].set_color('white')
    return ax

if __name__ == "__main__":
    histone_types = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
    histone_col = sns.color_palette("Set2", 6)
    # histone_col_light = sns.color_palette("hls", 8)
    path = '/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_annotated/'
    # path = '/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/majiq_res/correl/temp/'

    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot2grid((5, 4), (0, 0), colspan=2, rowspan=2)

    all_his_data = []
    for i, his in enumerate(histone_types):
        all_files = glob.glob(path + his + "/*.bed")
        his_data = []
        for file in all_files:
            print(file)
            data = pd.read_csv(file, header=None, sep=' ').iloc[:, 0].tolist()
            his_data.append(data)
            # ECDF plot
            sns.ecdfplot(data=data, color=histone_col[i], alpha=0.01, linewidth=0.75)
        all_his_data.append(flatten_list(his_data))

    no_nonzero = [1 * (np.count_nonzero(his_data)) for his_data in all_his_data]
    new_idx = np.argsort(no_nonzero)
    sorted_histone_types = [histone_types[i] for i in new_idx]
    sorted_all_his_data = [all_his_data[i] for i in new_idx]

    print('Drawing first plot...')
    for i, his_data in enumerate(sorted_all_his_data):
        sns.ecdfplot(data=his_data, color=histone_col[i], alpha=0.9, label=i, linewidth=2)

    ax1.set_title('A', fontweight='bold', loc='left', fontsize=18)
    ax1.set_xlabel("M-value", fontsize=16)
    ax1.set_ylabel("Cumulative proportion", fontsize=16)
    h, l = ax1.get_legend_handles_labels()

    # Remove zero
    sorted_all_his_data = [[i for i in his_data if i != 0] for his_data in sorted_all_his_data]
    sorted_all_his_data_his = [len(data) * [sorted_histone_types[i]] for i, data in enumerate(sorted_all_his_data)]
    sorted_all_his_data = pd.DataFrame(flatten_list(sorted_all_his_data))
    sorted_all_his_data_his = pd.DataFrame(flatten_list(sorted_all_his_data_his))
    sorted_all_his_df = pd.concat([sorted_all_his_data, sorted_all_his_data_his], axis=1, ignore_index=True)
    sorted_all_his_df.columns = ['M-val', 'Histone']

    # Add another plot
    print('Drawing second plot...')
    # boxpairs = [('H3K27ac', 'H3K27me3'), ('H3K27ac', 'H3K36me3')]
    boxpairs = [(sorted_histone_types[x], sorted_histone_types[y])
                for x in range(len(sorted_histone_types)) for y in range(len(sorted_histone_types)) if x < y]

    ax2 = plt.subplot2grid((5, 4), (0, 2), colspan=2, rowspan=4)
    sns.violinplot(data=sorted_all_his_df, x='Histone', y='M-val', hue='Histone', palette='Set2', cut=2,
                   linewidth=0.75, inner=None, dodge=False)
    ax2.legend([], [], frameon=False)
    ax2, test_results = add_stat_annotation(ax2, data=sorted_all_his_df, x='Histone', y='M-val', box_pairs=boxpairs,
                                           test='Mann-Whitney', text_format='star', loc='inside', linewidth=0.9)

    ax2.set_title('B', loc='left', fontweight='bold', fontsize=18)
    ax2.set_xlabel(None)
    ax2.set_ylabel("M-value", fontsize=16)
    ax2.set_xticklabels(labels=sorted_histone_types, rotation=45)

    ax3 = plt.subplot2grid((5, 4), (2, 0), colspan=2, rowspan=3)
    ax3 = add_blank_subplot(ax3)
    ax4 = plt.subplot2grid((5, 4), (4, 2), colspan=2, rowspan=1)
    ax4 = add_blank_subplot(ax4)

    plt.legend(handles=h, labels=sorted_histone_types, edgecolor='white', loc='lower right', bbox_to_anchor=(1, 0),
                          title="Histone type", ncol=2, framealpha=0)

    plt.tight_layout()
    plt.xticks(fontsize=12)
    plt.savefig('mval_whitney.tiff', dpi=300)
