import re
import numpy as np
from itertools import combinations
import json

if __name__ == "__main__":
    his = 'H3K36me3'
    res_dir = '/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/majiq_res/correl/'\
              + his + '.bed'
    with open(res_dir) as reader:
        lines = [line.rstrip().split('\t') for line in reader]

    gene_mval = []
    gene_dpsi = []
    result = {}
    current_gene = lines[0][4].split('"')[1]
    for k, line in enumerate(lines):
        print(k)
        chr, start, end, strand = line[0], line[1], line[2], line[3]
        flank_type, gene_id, exon_id, mval, dhm_pval = str, str, str, float, float
        lsv_type, lsv_pval, dpsi, a5ss, a3ss, es = '_', 1.0, 0.0, False, False, False
        lsv_start = line[5].split(',')
        for i, txt in enumerate(line):
            if i == 4:
                txt = re.split(';|"', txt)
                flank_type, gene_id, exon_id, mval, dhm_pval = txt[0], txt[2], txt[8], float(txt[10]), float(
                    txt[11])
            elif i == 7:
                txt_ = [re.split(';', t) for t in re.split(',', txt)]
                if len(lsv_start) >= 2:
                    # print('*' * 20, k, lsv_start, flank_type, len(lsv_start), '*' * 20)

                    all_pairs = list(combinations(range(0, len(lsv_start)), 2))
                    same_lsv = [[txt_[l[0]], txt_[l[1]]] for l in all_pairs if
                                lsv_start[l[0]] == lsv_start[l[1]]]
                    if len(same_lsv) != 0:
                        if flank_type == 'flank_end':
                            a = [l for l in same_lsv[0] if l[0] == 's'][0]
                        else:
                            a = [l for l in same_lsv[0] if l[0] == 't'][0]
                    else:
                        sources = [l for l in txt_ if l[0] == 's']
                        targets = [l for l in txt_ if l[0] == 't']
                        if flank_type == 'flank_end' and len(sources) != 0:
                            a = sources[np.nanargmax([abs(float(l[2])) for l in sources])]
                        elif flank_type == 'flank_start' and len(targets) != 0:
                            a = targets[np.nanargmax([abs(float(l[2])) for l in targets])]
                        else:
                            a = txt_[np.nanargmax([abs(float(l[2])) for l in txt_])]
                    # print('a: ', a)
                    lsv_type, lsv_p_val, dpsi, a5ss, a3ss, es = a[0], a[1], a[2], a[5], a[6], a[7]

        # print(chr, start, end, strand, flank_type, gene_id, exon_id, mval, dhm_pval, lsv_type, lsv_pval, dpsi, a5ss, a3ss, es)

        if current_gene != gene_id:
            result[current_gene] = (gene_mval, gene_dpsi)
            gene_mval, gene_dpsi  = [], []
            current_gene = gene_id
        else:
            gene_mval.append(mval)
            gene_dpsi.append(dpsi)

    with open(str('/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/majiq_res/correl/res/'
                  + his + '.json'), 'w') as fp:
        json.dump(result, fp)