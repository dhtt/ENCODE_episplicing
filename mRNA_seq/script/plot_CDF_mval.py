from typing import List
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("ticks")


def flatten_list(input_list: List[list]):
    return [item for sublist in input_list for item in sublist]


if __name__ == "__main__":
    histone_types = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
    histone_col = sns.color_palette("husl", 6)
    histone_col_light = sns.color_palette("hls", 8)
    path = '/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_annotated/'  # use your path

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
            sns.ecdfplot(data=data, color=histone_col[i], alpha=0.00075)
        all_his_data.append(flatten_list(his_data))

    for i, his_data in enumerate(all_his_data):
        sns.ecdfplot(data=his_data, color=histone_col_light[i], alpha=0.9, label=i)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("M-value", fontsize=16)
    plt.ylabel("Cumulative proportion", fontsize=16)

    h, l = ax.get_legend_handles_labels()
    plt.legend(handles=h, labels=histone_types, loc='lower right', title="Histone type")
    plt.savefig('mval_cdf2.tiff', dpi=300)
