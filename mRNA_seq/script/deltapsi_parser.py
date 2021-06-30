from collections import defaultdict
from typing import List
import pandas as pd
import glob


def print_dict(a_dict):
    for k_, v_ in a_dict.items():
        print(k_, ": ", v_)


def flatten_list(input_list: List[list]):
    return [item for sublist in input_list for item in sublist]


class LSV:
    def __init__(self, gene_id, LSV_info):
        self.gene_id = gene_id
        self.loc_LSV = LSV_info['LSV_ID']
        self.p_dPSI = LSV_info['P(|dPSI|>=0.20)_per_LSV_junction']
        self.e_dPSI = LSV_info['E(dPSI)_per_LSV_junction']
        self.no_junctions = [int(x) for x in LSV_info['Num_Junctions']]
        self.no_exons = [int(x) for x in LSV_info['Num_Exons']]
        self.A5SS = [LSV_info['A5SS'][i] for i in range(len(LSV_info['LSV_ID']))]
        self.A3SS = [LSV_info['A3SS'][i] for i in range(len(LSV_info['LSV_ID']))]
        self.ES = [LSV_info['ES'][i] for i in range(len(LSV_info['LSV_ID']))]

        self.sig_LSV = {'loc_LSV': [], 'chr': [], 'strand': [], 'gene': [], 'LSV_type': [], 'start': [], 'end': [],
                        'p_dPSI': [], 'e_dPSI': [], 'no_junctions': [], 'no_exons': [], 'A5SS': [], 'A3SS': [], 'ES': []
                        }

    def get_sig_LSV(self, mapping_info, PSI_threshold=0.2, PSI_probability=0.95):
        for i, LSV_ in enumerate(self.e_dPSI):
            sig_e_dPSI = abs(max(LSV_, key=abs))
            sig_idx = LSV_.index(max(LSV_, key=abs))
            if sig_e_dPSI > PSI_threshold and self.p_dPSI[i][sig_idx] > PSI_probability:
                loc_LSV = self.loc_LSV[i].split(':')
                loc_LSV_coord = loc_LSV[2].split('-')
                self.sig_LSV['loc_LSV'].append(self.loc_LSV[i])

                gene_chr = mapping_info.query('gene == "' + loc_LSV[0] + '"')
                if len(gene_chr) != 0:
                    self.sig_LSV['chr'].append(gene_chr['chr'].values[0])
                    self.sig_LSV['strand'].append(gene_chr['strand'].values[0])
                else:
                    self.sig_LSV['chr'].append("chr_")
                    self.sig_LSV['strand'].append("strand_")

                self.sig_LSV['gene'].append(loc_LSV[0])
                self.sig_LSV['LSV_type'].append(loc_LSV[1])
                self.sig_LSV['start'].append(loc_LSV_coord[0])
                self.sig_LSV['end'].append(loc_LSV_coord[1])
                self.sig_LSV['p_dPSI'].append(self.p_dPSI[i][sig_idx])
                self.sig_LSV['e_dPSI'].append(self.e_dPSI[i][sig_idx])
                self.sig_LSV['no_junctions'].append(self.no_junctions[i])
                self.sig_LSV['no_exons'].append(self.no_exons[i])
                self.sig_LSV['A5SS'].append(self.A5SS[i])
                self.sig_LSV['A3SS'].append(self.A3SS[i])
                self.sig_LSV['ES'].append(self.ES[i])

        return self.sig_LSV

    # def trim_sig_LSV(self):
    #     remove_idx = []
    #     if len(self.sig_LSV['loc_LSV']) != 0:
    #         all_coords = [x.split(':')[2] for x in self.sig_LSV['loc_LSV']]
    #         for i, x in enumerate(all_coords):
    #             for j, y in enumerate(all_coords):
    #                 if j > i and x == y:
    #                     remove_idx.append(i) if self.sig_LSV['e_dPSI'][i] < self.sig_LSV['e_dPSI'][j] \
    #                         else remove_idx.append(j)
    #     remove_idx.sort()
    #
    #     for k_, v_ in self.sig_LSV.items():
    #         self.sig_LSV[k_] = [x for i, x in enumerate(v_) if i not in remove_idx]
    #
    #     return self.sig_LSV

    def trim_sig_LSV(self):
        remove_idx = []
        if len(self.sig_LSV['loc_LSV']) != 0:
            all_coords = [x.split(':')[2] for x in self.sig_LSV['loc_LSV']]
            for i, x in enumerate(all_coords):
                for j, y in enumerate(all_coords):
                    if j > i and x == y:
                        remove_idx.append(i) if self.sig_LSV['e_dPSI'][i] < self.sig_LSV['e_dPSI'][j] \
                            else remove_idx.append(j)
        remove_idx.sort()

        for k_, v_ in self.sig_LSV.items():
            self.sig_LSV[k_] = [x for i, x in enumerate(v_) if i not in remove_idx]

        return self.sig_LSV


class LSVParser:
    def __init__(self,
                 delta_psi_path='/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/majiq'
                                '/aorta_adiposetissue/aorta_adiposetissue.deltapsi.tsv',
                 mapping_path='/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/refgen/gene_chr.txt'
                 ):
        print(delta_psi_path)
        self.mapping_info = pd.read_csv(mapping_path, '\t')

        with open(delta_psi_path) as reader:
            lines = [line.rstrip().split('\t') for line in reader]

        LSV_info = defaultdict()
        for i, line in enumerate(lines):
            col_name = lines[0]
            num_col = 15
            if i != 0:
                gene_name = line[0]
                LSV_info[gene_name] = defaultdict(list) if gene_name not in LSV_info.keys() else LSV_info[gene_name]

                if 'na' not in line[1] and 'na' not in line[2]:
                    for j in range(1, num_col - 1):
                        info = line[j]
                        if j in [3, 4, 5, 6, 7]:
                            info = [float(number) for number in info.split(';')]
                        LSV_info[gene_name][col_name[j]].append(info)

                    # If there are IR, store coordinates, else store "_"
                    ir_coords = line[num_col - 1] if len(line) == num_col else "_"
                    LSV_info[gene_name][col_name[num_col - 1]].append(ir_coords)
                    
        self.LSV = {gene: LSV(gene, info) for gene, info in LSV_info.items()}  # Each gene has an LSV info list
        self.all_genes = self.LSV.keys()
        self.no_genes = len(self.all_genes)

        self.sig_LSV = defaultdict()
        self.no_sig_genes = None
        self.sig_LSV_df = None

    def get_sig_LSV(self, trim=False):
        for gene, info in self.LSV.items():
            LSV_info = info.get_sig_LSV(mapping_info=self.mapping_info)
            if trim:
                LSV_info = info.trim_sig_LSV()
            self.sig_LSV[gene] = LSV_info if len(LSV_info['gene']) != 0 else None
        self.sig_LSV = {k: x for k, x in self.sig_LSV.items() if x is not None}
        self.no_sig_genes = len(self.sig_LSV)
        return self.sig_LSV

    def sig_LSV_to_DF(self):
        sig_LSV_list = [pd.DataFrame.from_dict(info, orient='columns') for gene, info in self.sig_LSV.items()]
        self.sig_LSV_df = pd.concat(sig_LSV_list)
        return self.sig_LSV_df

    def to_gtf(self, output_path):
        self.sig_LSV_df['none'] = '.'
        self.sig_LSV_df = self.sig_LSV_df[['chr', 'none', 'none', 'start', 'end', 'none', 'strand', 'none', 'gene',
                                           'LSV_type', 'p_dPSI', 'e_dPSI', 'no_junctions', 'no_exons', 'A5SS', 'A3SS',
                                           'ES']]
        self.sig_LSV_df.to_csv(output_path, sep='\t', header=False, index=False)


if __name__ == "__main__":
    # result_folder = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/majiq/temp_res/"
    # all_deltapsi_dirs = glob.glob("/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/
    # majiq/*/*.deltapsi.tsv")
    result_folder = '/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_to_be_mapped/'
    all_deltapsi_dirs = glob.glob('/home/dhthutrang/ENCODE/mRNA_seq/script/majiq_deltapsi_res/*/*.deltapsi.tsv')
    refgen_map = '/home/dhthutrang/ENCODE/refgen/gene_chr.txt'
    for input_dir in all_deltapsi_dirs:
        LSV1 = LSVParser(delta_psi_path=input_dir, mapping_path=refgen_map)

        LSV_list = LSV1.LSV
        LSV1.get_sig_LSV(trim=True)
        LSV1.sig_LSV_to_DF()

        output_dir = result_folder + input_dir.split('/')[-1]
        LSV1.to_gtf(output_dir)

    # # TESTING
    # LSV1 = LSVParser(delta_psi_path='/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/majiq/aorta_adiposetissue/temp.tsv')
    # LSV_list = LSV1.LSV
    #
    # # for key, val in LSV_list.items():
    # #     print('_' * 80)
    # #     print(key)
    # #     print(val.loc_LSV, val.no_exons, val.AS_types)
    #
    # # setd5 = LSV_list['SETD5']
    # # setd5.get_sig_LSV()
    # # setd5.trim_sig_LSV()
    # sig_LSV = LSV1.get_sig_LSV(trim=False)
    # temp = len(sig_LSV)
    # df = LSV1.sig_LSV_to_DF()
    # temp = len(df)
    # print(df.head(25))
    # LSV1.to_gtf('/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/majiq/temp_res/temp.txt')

