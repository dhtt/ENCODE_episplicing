import re
import numpy as np
from itertools import combinations
import json
import pandas as pd
from scipy.stats import linregress
import statsmodels.stats.multitest


def parse_annotated_result(res_dir, out_dir):
    his = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
    his_df = [pd.read_csv(str('temp/' + h + '.bed'), header=None, sep=' ').iloc[:, 0].to_list() for h in his]

    with open(res_dir) as reader:
        lines = [line.rstrip().split('\t') for line in reader]

    gene_dpsi, gene_as = [], []
    result = {}
    current_gene = lines[0][6].split('"')[1]
    for k, line in enumerate(lines):
        # print('*' * 20, k, line, '*' * 20)
        chrom, flank_type, start, end, strand, gene_id, exon_id = line[0], line[2], line[3], line[4], line[5], str, str
        lsv_type, lsv_pval, dpsi, a5ss, a3ss, es = '_', 1.0, 0.0, "False", "False", "False"
        lsv_start, lsv_end = line[7].split(','), line[8].split(',')
        # print(lsv_start, lsv_end)

        for i, txt in enumerate(line):
            # print(i, txt)
            if i == 6:
                txt = re.split(';|"', txt)
                gene_id, exon_id = txt[1], txt[7]
            elif i == 9:
                txt_ = [re.split(';', t) for t in re.split(',', txt)]
                # print(txt_)
                if len(lsv_start) >= 2:
                    # print('*' * 20, k, lsv_start, flank_type, len(lsv_start), '*' * 20)

                    all_pairs = list(combinations(range(0, len(lsv_start)), 2))
                    same_lsv = [[txt_[pair[0]], txt_[pair[1]]] for pair in all_pairs if
                                lsv_start[pair[0]] == lsv_start[pair[1]] and lsv_end[pair[0]] == lsv_end[pair[1]]]
                    if len(same_lsv) != 0:
                        if flank_type == 'flank_end':
                            a = [pair for pair in same_lsv[0] if pair[0] == 's'][0]
                        else:
                            a = [pair for pair in same_lsv[0] if pair[0] == 't'][0]
                    else:
                        sources = [fl for fl in txt_ if fl[0] == 's']
                        targets = [fl for fl in txt_ if fl[0] == 't']
                        if flank_type == 'flank_end' and len(sources) != 0:
                            a = sources[np.nanargmax([abs(float(fl[2])) for fl in sources])]
                        elif flank_type == 'flank_start' and len(targets) != 0:
                            a = targets[np.nanargmax([abs(float(fl[2])) for fl in targets])]
                        else:
                            a = txt_[np.nanargmax([abs(float(fl[2])) for fl in txt_])]
                    # print('a: ', a)
                    lsv_type, lsv_p_val, dpsi, a5ss, a3ss, es = a[0], a[1], a[2], a[5], a[6], a[7]
                elif len(lsv_start) == 1 and str(lsv_start[0]) != str(-1):
                    a = txt_[0]
                    lsv_type, lsv_p_val, dpsi, a5ss, a3ss, es = a[0], a[1], a[2], a[5], a[6], a[7]

        # print(chrom, start, end, strand, flank_type, gene_id, exon_id, lsv_type, lsv_pval, dpsi, a5ss, a3ss, es)
        # print(k)

        if k == len(lines) - 1:
            result[current_gene] = {'majiq': gene_dpsi,
                                    'histone': [his_[max(0, (k - len(gene_dpsi))): k + 1] for his_ in his_df],
                                    'as': gene_as, 'length': len(gene_dpsi)}

        if current_gene != gene_id:
            result[current_gene] = {'majiq': gene_dpsi,
                                    'histone': [his_[max(0, (k - len(gene_dpsi))): k] for his_ in his_df],
                                    'as': gene_as, 'length': len(gene_dpsi)}
            gene_dpsi = [float(dpsi)]
            gene_as = [[a5ss, a3ss, es]]
            current_gene = gene_id
        else:
            gene_dpsi.append(float(dpsi))
            gene_as.append([a5ss, a3ss, es])

    with open(out_dir, 'w') as outfile:
        json.dump(result, outfile)


def get_mask(rval, pval):
    all_r_bool = [np.where(np.asarray(abs(rval.loc[:, his])) >= 0.5, True, False) for his in histone_types]
    all_padj_bool = [statsmodels.stats.multitest.multipletests(pval.loc[:, his], method='fdr_bh', alpha=0.05)[0]
                     for his in histone_types]
    mask_bool = np.logical_and(all_r_bool, all_padj_bool)
    return mask_bool


if __name__ == "__main__":
    histone_types = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
    # ==================================
    # PART 1: PARSING
    # ==================================
    # parse_annotated_result(
    #     res_dir='/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/majiq_res/correl/temp/majiq.bed',
    #     out_dir='/Users/dhthutrang/Documents''/BIOINFO/Episplicing/ENCODE_episplicing/mRNA_seq/script/majiq_res/correl/res2/majiq.js')

    # ==================================
    # PART 2: CORRELATION
    # ==================================
    # with open('res2/majiq.js') as json_file:
    #     data = json.load(json_file)
    #
    # all_cor = {}
    # for gene, info in data.items():
    #     # print(gene)
    #     majiq_res = info['majiq']
    #     his_res = info['histone']
    #     no_exon = info['length']
    #     as_event = info['as']
    #     # print(majiq_res, len(majiq_res))
    #     # print(his_res[0], len(his_res[0]))
    #     # print(len(majiq_res), len(his_res[0]))
    #     if no_exon > 3:
    #         # |M-value| >= 1
    #         for i, res in enumerate(his_res):
    #             his_res[i] = [abs(np.where(abs(np.asarray(k)) >= 1, k, 0)) for k in res]
    #         # |dPSI| >= 0.2
    #         for i, res in enumerate(majiq_res):
    #             majiq_res[i] = abs(np.where(abs(np.asarray(res)) >= 0.2, res, 0))
    #         gene_his_cor = [linregress(majiq_res, i)[2:4] for i in his_res]
    #     else:
    #         gene_his_cor = [[0.0, 1.0] for i in his_res]
    #     all_cor[gene] = gene_his_cor
    #
    # with open('res2/majiq_cor_abs.js', 'w') as fp:
    #     json.dump(all_cor, fp)

    # ==================================
    # PART 3: RESULT
    # ==================================
    with open('res2/majiq_cor_abs.js') as json_file:
        data = json.load(json_file)
    all_r, all_p = [], []
    for all_cor in data.values():
        all_r.append([cor_his[0] for cor_his in all_cor])
        all_p.append([cor_his[1] for cor_his in all_cor])
    all_r = pd.DataFrame(all_r, index=data.keys(), columns=histone_types)
    all_p = pd.DataFrame(all_p, index=data.keys(), columns=histone_types)
    mask_gene = get_mask(rval=all_r, pval=all_p)

    all_sig_genes = {his: [] for his in histone_types}
    for x, mask in enumerate(mask_gene):
        his = histone_types[x]
        true_idx = np.where(mask)[0]
        all_sig_genes[his] = all_r.iloc[true_idx, x].index.values.tolist()

    with open('res2/majiq_sig_genes_abs.js', 'w') as fp:
        json.dump(all_sig_genes, fp)

    # ==================================
    # PART 4: Compare DEXSEQ
    # ==================================
    with open('res2/majiq_sig_genes_abs.js') as fp:
        majiq_res = json.load(fp)

    # dexseq_res = {histone_types[0]: ['ACOT11','ACP5','ACYP2','AEBP1','AFG3L2','AIF1','AIFM3','ANGPT1','ANK1','ANKS3','AOX1','AP4M1','ARHGAP6','ARMC9','ARNT2','ARSA','ATF3','ATP5PD','BICRA','BNIPL','C11orf91+CD59','C1orf162+ATP5PB','CC2D2B','CCDC114','CCDC13','CCDC192','CCDC33','CCDC73','CCDC78','CDC14B','CDCA5','CDT1','CELF4','CEP112','CFAP157','CHMP7','CKMT2','CLIC1','COMMD4','COQ5','CR1','CRACR2A','CRB1','CREB3L2','CRHR1+LINC02210-CRHR1','CRTAP','CSNK2A2','CUEDC2','CUL4A','DHTKD1','DIAPH3','DNAAF1','DNPH1','DPF1','DRAM1','DSCAML1','ECHDC2','EIF4EBP3+ANKHD1-EIF4EBP3+ANKHD1','EML1','EPHA4','ERAP1','ESPN','EVI5L','EYA1','FAM107A','FAM111A','FAM178B','FAM72D','FBXO41','FEZF1','FNDC3A','FYB1','GGT5','GIGYF1','GPLD1','GPR176','GPR19','GPRIN3','GTF2I','HCK','HCRTR1','HNRNPUL1','HPSE2','IFNGR2','IGSF22','IKZF1','IL12RB1','IL23A','IL32','IQCC','JOSD2','KCTD5','KIAA1522','KLHL18','KLHL2','KLHL26','LAMC1','LETM1','LETMD1','LIMK1','LOC102723750','MAD1L1','MAD2L2','MAFK','MAML1','MAML3','MBOAT2','ME1','MED6','MELTF','MFAP3','MFRP+C1QTNF5','MGAT4B','MIF4GD','MMP14','MMRN2','MRC2','MVB12A','MYLK4','MYO19','MYOM2','NDOR1','NECAB3','NEU3','NISCH','NMT2','NOXO1','NSFL1C','NT5C3A','NUP62+IL4I1','OGFOD3','P4HA1','PBK','PCDH1','PDZD2','PHF23','PHF7','PIF1','PIK3CG','PKD2L1','PODXL','POFUT2','POLR2H','PPP1R21','PPP3R1','PPP5D1','PPP6R1','PRAG1','PRKACB','PRKG1','PRPF4','PRPF8','PSMB10','PTPN9','RAB11A','RASA2','RASA3','RBM7','RBPMS2','RNF32','ROBO3','RPL18','RPL22','RPS27','RPS4X','RTF1','SAXO1','SCAI','SERPINB9','SHANK1','SHCBP1','SHF','SIAH1','SIMC1','SLC15A2','SLC16A11','SLC25A21','SLC4A8','SLC52A2','SLC6A8','SMO','SOBP','SPTAN1','SRGAP2C','SSBP4','ST3GAL6','ST7','STAB1','STPG2','STX12','SURF6','SYT5','TBC1D10B','TCN2','TCTE3','TENM2','TEX264','TJP2','TLK1','TLN1','TM9SF1','TMEM116','TMEM45A','TMEM8B','TNFRSF10D','TNS3','TOMM20','TRIM2','TXK','UBE2O','USP36','USP53','VIPR1','VPS35L','WNT5A','XPNPEP1','YPEL1','YTHDF1','ZNF316','ZNF385C','ZNF554','ZNF627','ZNF660+ZNF197+ZNF660-ZNF197','ZNF780B','ZNF827','ZNF829'],
    #  histone_types[1]: ['ACOXL','AFAP1L2','ARHGAP6','ARMC9','ASTN2','ATG9B','ATP6V1B1','ATP8A2','AXL','BAHD1','BCL6B','BLK','BNIP5','C2CD4C','C3orf86+ZNF445','CACNA1G','CARHSP1','CARMIL3','CCDC136','CCDC149','CCDC57','CD207','CDC14B','CDK5RAP1','CEP170B','CFAP157','CHMP7','CNTN4','COL1A2','CR1','CRACR2A','CRAMP1','DIDO1','DNM1','ECHDC2','EML1','ERBB2','ESYT3','FGFRL1','FN3K','FOXP2','FRMD4B','FTCD','GPSM1','GRAMD1C','GRIK1','HCRTR1','HIP1','HIPK1','HLA-DPA1','HPCAL1','IDUA','IGF2BP1','IL12RB1','IL17RC','IRAK4','ITIH4','KAZN','KCNIP1','KRT12','LHX6','LIN7B','LMO3','LRRC2','LRRC36','LYPD6B','MAD2L2','MBD2','MCF2L2','MCTP1','MKLN1','MMP15','MTMR1','MYH3+SCO1','NBAS','NDUFS2','NINJ2','NRP2','NRXN1','NTRK1','NUP62+IL4I1','OGFOD1','PCED1B','PIWIL1','PLCL2','PLEKHM2','PODXL','PRKG1','PRPF6','PXDNL','RACGAP1','RBM44','RIPOR2','S100A2','SCYL1','SEC24D','SEMA6B','SEPTIN12','SEPTIN5','SERINC5','SH3RF2','SIRPA','SLC4A8','STIL','STOX2','TCHP','TEX15','TEX30','TLR5','TRIM65','TRPM6','TSC22D3','TUBA4A','USH2A','WNT5A','XYLB','ZNF333','ZNF740'],
    #  histone_types[2]: ['AAGAB','ACOT11','ACOX2','ACP5','ACP6','ACTR3C','ADCY7','ADRA1A','AGBL2','ANKRD11','ANKRD34B','ANO7','AP1G2','AP4M1','APH1A','APLP2','ARHGEF17','ARMC8','ARNT2','ASCC2','ASMTL','AXIN2','B3GALNT1','BEND3','C16orf46','C1orf127','C3orf20','CAPN12','CATSPERG','CCDC39','CCDC65','CCDC73','CD37','CDT1','CEP170','CEP44','CFAP44','CHKA','CHST12','CIT','CLDN10','CLMN','CNOT8','CNTD1','CNTNAP2','COL6A2','CORO6','CTSD','CTSK','CYTH2','DHX38','DHX8','DLG3','DNAAF1','DNAH9','DNPEP','DNPH1','E2F1','ECT2L','EEF1AKMT4-ECE2+ECE2+EEF1AKMT4','EFCAB7','ELAC2','EML6','ENPP5','FAM149B1','FAM166A','FAM184B','FBXL17','FBXL5','FBXO42','FCGR2A','FRMD4A','GALT','GAS2L3','GDI2','GIGYF1','GLB1L3','GRAMD1B','GSKIP','HCRTR1','HEATR1','HELZ2','HK1','ICA1L','IGF2BP1','IKZF1','IL16','IL7','IPO11','ITGB1BP1','JOSD2','KIAA0319','KIAA0895','KLHL14','KLHL32','KLRG2','LCORL','LHFPL2','LNX1','LOC100289561','LPAR1','LRBA','LRRC27','LRRC37B','LRRC39','LRRC6','LRRC75A','MAMSTR','MASP2','MCMDC2','MED14','MED15','MEF2A','MEGF11','MIF4GD','MLLT10','MMP14','MOCS1','MROH2A','MUTYH','MYT1','NAALADL1','NADK','NDUFA9','NEMP1','NLGN4X','NRG4','NUAK1','NUDT17','OAS3','OCLN','OPRL1','P3R3URF+PIK3R3+P3R3URF-PIK3R3','PAIP2B','PBK','PCED1A','PDGFD','PEX11G','PLEKHG3','PM20D2','PPIP5K2','PRKAA2','PXDNL','RAB21','RAB8A+HSH2D','RAC2','RANBP6','RASSF8','RBFOX3','RBP1','RDH11','RELB','REPS2','RFNG','RHEBL1','RHPN2','RNF150','RNF24','RNF32','RPGRIP1','RPL13A','RPL23','RPL26L1','RPS11','RPS4X','RPS8','RUNX1','SACS','SETD9','SHCBP1','SHISA5','SHISA8','SLC26A5','SLC2A12','SLC35A1','SLC39A5','SLC4A11','SLFN11','SMAD5','SMO','SORD','SP140','SPAG17','SPOCK1','SPPL2B','SPTAN1','SRPK2','STK25','STK32C','STX1B','SULF2','SURF6','SYN1','SYN3','SYNDIG1L','SYNPO','TAF6L','TBC1D1','TBKBP1','TDRKH','THNSL2','TKFC','TLDC2','TMEM168','TMEM175','TMEM225B+ZNF655','TMEM45B','TPSG1','TTC16','TTC23L','TUT7','TXNDC15','USP44','VOPP1','XBP1','ZACN','ZDHHC1','ZDHHC6','ZEB1','ZNF112','ZNF141','ZNF227'],
    #  histone_types[3]: ['ABCC3','ABLIM2','ACSL5','ADAM23','ADAMTS13','AGBL3','ANGPT1','AOX1','AP1G2','AP3D1','ARL8A','ATP8A2','AVPR2','BEST1','BICRA','BZW2','C12orf75','C16orf46','C16orf89','C2CD4C','C2orf81','CABP5+PLA2G4C','CALHM2','CAPG','CARHSP1','CBX3','CC2D1B','CCDC107','CCDC130','CCDC73','CCDC82','CCNK','CD79B','CDC7','CDCA2','CENPT','CEP76','CGN','CGREF1','CHD5','CHN1','CHURC1+FNTB+CHURC1-FNTB','CKMT2','CLIC1','CNOT8','COBLL1','COL12A1','COL1A2','COPZ1','CPA5','CPXM1','CR1','CSF2RB','CSGALNACT2','CTCFL','CUL5','CUX1','CYTH2','DCAF1','DCC','DDX39B','DENND10','DENND5B','DGKZ','DHX38','DIDO1','DLK2','DNAAF1','DNAJC17','DOK5','DPM2','EEF2KMT','EIF2B3','EIF4E','EIF4EBP3+ANKHD1-EIF4EBP3+ANKHD1','EIF4ENIF1','ELANE','EPB41L5','EPHA7','ERICH6','EVC','EXOC4','FABP6','FAM178B','FAM193A','FAM98A','FBXO41','FDFT1','FGD1','FHL3','FKBP6','FLNA','FLNC','FOXM1','GALT','GATD1','GK5','GLB1L3','GPAT2','GPSM1','GRIA2','GRK1','HAUS3+POLN','HCFC1','HCK','HLA-DMA','HTR3A','IGSF22','IL12RB1','ILK','IRF8','JOSD2','KCTD16','KDM5D','KHDRBS2','KIAA1549','KMT5B','KRI1','LARP4B','LCN12','LIG1','LYPD5','MAOA','MAST1','MCMDC2','MED13L','MGST1','MORC4','MPHOSPH10','MPHOSPH9','MROH2A','MYO19','MYO1D','MYO5A','NAIP','NAT1','NCAPH2','NDUFV3','NFATC3','NPM3','NPR3','NT5C3A','NUMB','NUP62+IL4I1','OFD1','OS9','OTUB1','OXTR','PACSIN1','PARD3','PBK','PDE11A','PHGDH','PIGBOS1','PIP4K2B','PLAAT3','PLCL1','POGLUT1','PPP2R3A','PPP5D1','PRPF8','PTPRE','PUS7','RABGAP1L','RALBP1','RANBP2','RASL10B','RNF212','RPL17-C18orf32+RPL17+C18orf32','RPS11','RPS29','RTN4','RWDD3+TLCD4-RWDD3+TLCD4','SAMD3','SCN2A','SEC14L2','SEMA5B','SEPTIN3','SERINC5','SFI1','SGK2','SHANK1','SHCBP1','SHISA4','SIRPA','SLC23A3','SLC25A21','SLC26A4','SLC36A4','SLC44A5','SLC5A10','SLC5A6','SLC7A6','SLC7A9','SLFNL1','SLMAP','SMG7','SMO','SNPH','SNX30','SRPRA','SRSF10','SSB','ST7','STAC3','STK32A','STPG2','STT3B','STYXL1','SUPT3H','TADA2A','TAF4B','TAF6L','TATDN1','TBC1D19','TBC1D9B','TCTE1','TGDS','THRAP3','TLN1','TMEM123','TMEM225B+ZNF655','TMEM255B','TMEM74B','TOB2','TRAPPC4','TRIM26','TRPT1','TSPOAP1','TTBK2','TTC16','TUT7','UBE2D4','UGT3A2','UNG','UPF2','USP42','VCL','WDR37','XRN2','XYLB','YLPM1','ZBTB41','ZC3H14','ZDHHC6','ZNF333','ZNF44','ZNF580','ZNHIT6','ZSCAN32'],
    #  histone_types[4]: ['AATK','ABCC10','ACYP2','ADIRF+AGAP11','ADK','ADPRHL1','ADSS1','AFG3L2','AIF1','AKAP8L','AKNAD1','ANGPT1','ANKLE2','ANKRD11','ANKRD46','AP1G2','AP4M1','APLP2','ARHGAP6','ARHGEF2','ARL11','ARSA','ASB2','ATP5PD','BACH2','BAHCC1','BICRA','BNIPL','C12orf75','CACNA1A','CALHM2','CC2D2B','CCDC120','CCDC130','CCDC22','CD163L1','CDC14B','CDCA5','CELSR2','CEP112','CFAP157','CFLAR','CGN','CHMP7','CHST12','CKMT2','COL1A2','COL6A1','CORIN','CR1','CRHR1+LINC02210-CRHR1','CROCC','CRYBB2','CSF3R','CTNNBIP1','CTSB','CYP26B1','DAZAP2','DENND2C','DERL2','DHRS2','DHTKD1','DHX38','DNAH5','DNAJC6','DPP6','EBF1','ECHDC3','EIF2AK1','EIF3G','EMC9','EML1','ENPP3','EPS15','ERAP1','ERBB2','EYA1','FABP6','FAM107A','FAM166A','FAM72D','FBXW7','FN3K','FNDC3A','FXYD7','FZR1','GALNT2','GALT','GANC','GFAP','GIGYF1','GNG12','GPR176','GPR19','GPRIN3','GPSM1','GRAMD1C','GRHL3','GRIK1','GTF2I','HDAC9','HGSNAT','IGF2BP1','IGSF22','IL12RB1','IL32','IZUMO1','JOSD2','KAZN','KCNIP1','KCNMA1','KCTD5','KIAA0895L','KLHL2','LAMC1','LGR4','LGR6','LIMK1','LOC102723750','LRP8','LRRC36','LTBP4','LYPD6B','MAD1L1','MAD2L2','MAML1','MAML3','MAP4K1','MARK2','MARK3','MAU2','MBOAT2','MCRIP1','ME1','ME3','MED6','MELTF','MFAP3','MFRP+C1QTNF5','MMP14','MMP15','MMP24OS+MMP24-AS1-EDEM2+EDEM2','MPHOSPH9','MRC2','MROH7','MTUS1','MYLK4','NBAS','NFE2L1','NICN1+AMT','NINJ2','NIPSNAP2','NMT2','NSD3','NTRK1','OAS3','OCIAD1','P4HA1','PAIP2B','PARP3','PBK','PCDHA9+PCDHA6+PCDHAC1+PCDHA2+PCDHA7+PCDHA4+PCDHA3+PCDHA8+PCDHA1+PCDHAC2+PCDHA5+PCDHA13+PCDHA12+PCDHA10+PCDHA11','PKD2L1','PKLR','PKN1','PKP2','PLAUR','PLCH1','PLEC','PLEKHB2','PLS3','PNKD','POLR3H','PPP3R1','PPP5D1','PRKACB','PRPF4','PRSS36','PTCH2','PTDSS1','PTPN9','RAB11A','RABGAP1L','RAC2','RALGDS','RASA4B','RBPMS2','RGPD2','RGS14','RHCG','RIMS2','RNLS','ROR2','RTF1','RTN1','SAXO1','SCN2A','SEPTIN9','SERPINB8','SF3A2','SGCA','SH2B3','SH3GL1','SHANK1','SHCBP1','SHF','SIAH1','SIMC1','SLC25A21','SLC4A8','SMIM12','SMO','SMYD2','SRC','SRGAP2C','SSBP4','ST7','STPG2','STX12','SYTL3','TAF6L','TANK','TCN2','TCTE3','TDRKH','TENM2','TLK1','TLL2','TM9SF1','TMEM234','TMEM241','TMEM45A','TNFRSF10D','TNPO3','TNS3','TNS4','TOMM20','TRIM11','TRIM2','TRIM24','TRPM3','TSPAN33','TTLL6','UBXN11','UQCC2','USF1','USP53','VIPR1','VLDLR','WDR46','WNT10B','XPNPEP1','YPEL1','YTHDF1','ZFYVE28','ZNF254','ZNF354C','ZNF362','ZNF540','ZNF554','ZNF575','ZNF780B','ZNF827','ZNF829','ZNF83'],
    #  histone_types[5]: ['ABCB11','ABCC2','ABCD4','ACAP3','ACTR3C','ADAMTS13','ADM2+MIOX','AFF3','AK8','ANKRD13B','AOAH','APLP2','APOBEC3B','ARHGAP23','ARHGEF15','ARHGEF6','ARL13B','ASH1L','ASPM','ATAD3A','ATG9A','ATP10D','BCAR1','BCAT1','BEST1','BICDL1','BRCA2','BRSK2','BTBD7','C19orf57','CADM4','CCDC186','CD200','CDC5L','CDHR1','CDHR3','CECR2','CELSR2','CERS3','CFAP100','CLSPN','CMTM1+CKLF-CMTM1+CKLF','CNOT2','CNTFR','COLEC12','COMMD5','CPSF6','CROT','CRYBB2','CTIF','CUL1','CYREN','DAAM2','DAB2IP','DBX1','DCAF1','DCAF12','DCDC1','DENND5B','DHCR24','DPP6','EEF1AKMT4-ECE2+ECE2+EEF1AKMT4','EIF4E3','ESRRA','EVL','EXT1','EYA1','F11R','F3','FAM81A','FANK1','FGF17','FKBP9','FNDC1','FOXP1','GCNA','GOSR1','GPR176','GRIK1','GSTM3','HGSNAT','HNRNPAB','HOOK2','HSPB11','IFT81','IFTAP','IGF2BP1','ITPKC','KCNA6','KCNH2','KCNT2','KDM8','KIF1C','KIF21B','LGSN','LIN7B','LOC101930307','LOC107986035+SEMA3F','LRRC27','MAD2L2','MAML3','MANEA','MAP1LC3C','MAP3K20','MAPKAPK5','MBD2','MCF2L','MCM2','MGAM2','MIB2','MMP25','MRC2','MRPL3','MSANTD4','MTRR','MUTYH','NAF1','NLE1','NSUN7','NUDC','ORC1','OSBPL7','OTULIN','PARP8','PIP5KL1','PKD1','PLBD2','PLXDC2','POMT2','PRDM10','PRKN','PRPF19','PSMD12','PSMD5','PSMG4','PTCH2','PTP4A1','QPCTL','RAB3GAP1','RACGAP1','RARRES1','RASSF2','RBM39','RNH1','ROBO2','RTN4','RYK','SARM1','SDR42E2','SEC31B','SELENOW','SLC25A21','SLC38A9','SLC49A3','SLC7A1','SPARCL1','SPOCK1','ST6GALNAC1','STIL','STK32C','STX6','TANC1','TARBP1','TBC1D8','TCF7L2','TLN1','TMCC2','TMEM163','TMEM169','TMEM219','TMEM266','TMEM87A','TNFRSF10D','TOX2','TRAF7','TRIP10','TRIP12','TRPA1','TUBGCP2','UPF2','USP5','UTY','VWA8','WDR37','WDR97','WLS','YARS1','ZBTB41','ZBTB8A','ZC3H12B','ZC3HC1','ZDHHC11','ZDHHC13','ZNF148','ZNF225','ZNF230','ZNF333','ZNF337','ZNF33B','ZNF426','ZNF521','ZNF619','ZNF667','ZNF888+ZNF816+ZNF816-ZNF321P']}
    #
    # overlap_res = {}
    # for his in histone_types:
    #     overlap_res[his] = [gene for gene in majiq_res[his] if gene in dexseq_res[his]]
    #
    # with open('res2/gene_list_abs.txt', 'w') as filehandle:
    #     for his, res in majiq_res.items():
    #         filehandle.write('%s\n' % str('=====' + his + '====='))
    #         for i in res:
    #             filehandle.write('%s, ' % i)
    #         filehandle.write('\n')

# TODO 1: filter by AS events, flank end/start
# TODO 2: for all pairs
# TODO 3: find threshold for padj and val
# TODO 4: filter first exon
