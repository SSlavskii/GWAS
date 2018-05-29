import subprocess
import pandas as pd
import argparse

plink_path = '/home/ubuntu/gwas/old_gwas/tools/plink/plink'
ref_path = '/home/ubuntu/gwas/old_gwas/data/ref/eur_genotype_data/by_chr/'

RS_ID = 'SNP'
CHR = 'Chr'

snp_table_path = './ng.3097-S2.csv'


def get_snp_list(chr_num, table_path, bim_path):
    df = pd.read_table(table_path)
    df = df[df[CHR] == chr_num]
    bim = pd.read_table(bim_path, header=None)  # rs_id column has index 1
    snp_list = list(df[df[RS_ID].isin(bim[1])][RS_ID])
    return snp_list


def get_matrix_for_chromosome(chr_num, table_path):
    bfile_path = f"{ref_path}ALL.chr{chr_num}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_nodup"
    snp_list = get_snp_list(chr_num, table_path, bfile_path+'.bim')
    snps = ' '.join(snp_list)
    chr_num = str(chr_num)



    if len(snp_list) < 2:
        print(f"less than two SNPs on chromosome {chr_num}")
        return -1

    query = f"{plink_path} --bfile {bfile_path} --r square -out {chr_num} --snps {snps}"

    subprocess.call(query, shell=True)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-chr', type=int, dest="chr_num")
    return parser


if __name__ == "__main__":
    args, _ = create_parser().parse_known_args()
    get_matrix_for_chromosome(args.chr_num, snp_table_path)
