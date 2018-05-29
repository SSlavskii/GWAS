import subprocess
import pandas as pd
import argparse
import os

plink_path = '/home/ubuntu/gwas/old_gwas/tools/plink/plink'
ref_path = '/home/ubuntu/gwas/old_gwas/data/ref/eur_genotype_data/by_chr/'

RS_ID = 'ref_snp_id'
CHR = 'ref_chrom'
table_path = '/home/ubuntu/gwas/old_gwas/sega/ld_real/core_snps.csv'
R_threshold = 0.2

def get_snp_list(chr_num, table_path, bim_path):
    df = pd.read_table(table_path)
    df = df[df[CHR] == chr_num]
    bim = pd.read_table(bim_path, header=None)  # rs_id column has index 1
    snp_list = list(df[df[RS_ID].isin(bim[1])][RS_ID])
    return snp_list


def get_table_for_chr(chr_num):
    bfile_path = f"{ref_path}ALL.chr{chr_num}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_nodup"
    snp_list = get_snp_list(chr_num, table_path, bfile_path+'.bim')
    snps = ' '.join(snp_list)
    chr_num = str(chr_num)

    query = f"{plink_path} --bfile {bfile_path} --r -out {chr_num} --ld-snps {snps}" \
            f" --ld-window-kb 150 --ld-window 10000000"

    subprocess.call(query, shell=True)
    df = pd.read_csv(f"{chr_num}.ld", sep='\s+')
    df.dropna(inplace=True)
    df.drop(df[abs(df["R"]) < R_threshold].index, inplace=True)

    df.to_csv(f"{chr_num}.ld", index=False, sep='\t')
    os.remove(f"{chr_num}.log")
    os.remove(f"{chr_num}.nosex")


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-chr', type=int, default=0, dest="chr_num")
    return parser


if __name__ == "__main__":
    args, _ = create_parser().parse_known_args()
    if args.chr_num > 0:
        print(f"Calculation for chromosome {args.chr_num}")
        get_table_for_chr(args.chr_num)
    else:
        print('Calculation for all chromosomes')
        for i in range(1,23):
            get_table_for_chr(i)
        frames = [pd.read_csv(f"{i}.ld", sep='\t') for i in range(1, 23)]
        result = pd.concat(frames)
        result.to_csv("all_chr.ld", index=False, sep='\t')
        for i in range(1,23):
            os.remove(f"{i}.ld")