import pandahouse as ph
import subprocess
import os
import pandas as pd
import argparse
from scipy import stats

CLICKHOUSE_CONNECTION_PARAMS = {'host': 'http://localhost:8123', 'database': 'gwas'}
plink_path = '/home/ubuntu/gwas/old_gwas/tools/plink/plink'
ref_path  = '/home/ubuntu/gwas/old_gwas/data/ref/eur_genotype_data/by_chr/'

if os.environ.get("TEST"):
    CLICK_DB = "test"
    CLICK_TABLE = "gwas_snp"
else:
    CLICK_DB = "gwas"
    CLICK_TABLE = "test_gwas_snp"


class SNP:
    def __init__(self, rs_id):
        df = self.get_snp_info(rs_id)
        self.chr = df['chrom'].item()
        self.id  = df['snp_id'].item()
        self.ea = df['ea'].item()
        self.ra = df['ra'].item()
        self.bp = df['bp'].item()

    def get_snp_info(self, sr_id):
        query = "select chrom, snp_id, ea, ra, bp from {0}.{1} where snp_num=={2}  and  outlier=0  limit 1".\
            format(CLICK_DB, CLICK_TABLE, sr_id)
        print(query)
        row = ph.read_clickhouse(query, connection=CLICKHOUSE_CONNECTION_PARAMS)
        return row


def get_z_correlation(rs_id1, rs_id2):
    query = "SELECT z_1, z_2  \
       FROM (\
            SELECT gwas_1.gwas_id, gwas_1.z AS z_1\
            FROM {0}.{1} as gwas_1\
            WHERE gwas_1.snp_num = {2} and  outlier=0 \
        ) ALL INNER JOIN (\
            SELECT gwas_2.gwas_id, gwas_2.z AS z_2\
            FROM {0}.{1} as gwas_2\
            WHERE gwas_2.snp_num = {3} and  outlier=0 \
        ) USING gwas_id".format(CLICK_DB, CLICK_TABLE, rs_id1, rs_id2)

    z_df = ph.read_clickhouse(query, connection=CLICKHOUSE_CONNECTION_PARAMS)
    z_df_filtered = z_df.loc[(abs(z_df['z_1']) < 2) & (abs(z_df['z_2']) < 2)]

    z_corr_spearman = stats.spearmanr(z_df['z_1'], z_df['z_2']).correlation
    z_corr_less_than_2 = stats.pearsonr(z_df_filtered['z_1'], z_df_filtered['z_2'])[0]
    z_corr_ordinary = stats.pearsonr(z_df['z_1'], z_df['z_2'])[0]

    return [z_corr_ordinary, z_corr_less_than_2, z_corr_spearman]


def get_plink_output(snp1, snp2, chr):
    out_file_prefix = str(snp1) + '_' + str(snp2)
    try:
        os.remove(out_file_prefix + '.ld')
        os.remove(out_file_prefix + '.frq')
    except FileNotFoundError:
        pass

    subprocess.call([plink_path, '--bfile', ref_path +
                     'ALL.chr' + str(chr) + '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_nodup',
                     '--snps', 'rs' + str(snp1), 'rs' + str(snp2), '-r', 'in-phase',
                     '--freq', '-out', out_file_prefix])

    ld_table = pd.read_table(out_file_prefix + '.ld', sep='\s+')
    freq_table = pd.read_table(out_file_prefix + '.frq', sep='\s+')

    ld_table['EA_A'] = freq_table[freq_table['SNP'] == ld_table['SNP_A'].item()]['A1'].item()
    ld_table['RA_A'] = freq_table[freq_table['SNP'] == ld_table['SNP_A'].item()]['A2'].item()
    ld_table['EA_B'] = freq_table[freq_table['SNP'] == ld_table['SNP_B'].item()]['A1'].item()
    ld_table['RA_B'] = freq_table[freq_table['SNP'] == ld_table['SNP_B'].item()]['A2'].item()

    ld_table['EAF_A'] = freq_table[freq_table['SNP'] == ld_table['SNP_A'].item()]['MAF'].item()
    ld_table['EAF_B'] = freq_table[freq_table['SNP'] == ld_table['SNP_B'].item()]['MAF'].item()

    os.remove(out_file_prefix + '.ld')
    os.remove(out_file_prefix + '.frq')
    os.remove(out_file_prefix + '.log')
    os.remove(out_file_prefix + '.nosex')
    return ld_table


def generate_df_for_snp_pair(snp1, snp2):
    first_snp = SNP(snp1)
    second_snp = SNP(snp2)
    list_of_z_corr = get_z_correlation(snp1, snp2)
    dist = abs(first_snp.bp - second_snp.bp)
    plink_info = get_plink_output(snp1, snp2, first_snp.chr)  # this is dataframe
    plink_info['dist'] = dist

    # restore coordination between first and second snp and A and B in plink_info
    if first_snp.id == plink_info['SNP_B'].item():
        first_snp, second_snp = second_snp, first_snp

    # add correlation column to dataframe
    plink_info['z_corr_raw'] = list_of_z_corr[0]
    plink_info['z_corr_less_2'] = list_of_z_corr[1]
    plink_info['z_corr_spearman'] = list_of_z_corr[2]


    factor = 1
    # check if alleles coincide with clickhouse
    if plink_info['EA_A'].item() != first_snp.ea:
        print('EA_A != clickhouse')
        plink_info['EA_A'], plink_info['RA_A'] = plink_info['RA_A'].item(), plink_info['EA_A'].item()
        plink_info['EAF_A'] = 1 - plink_info['EAF_A'].item()
        factor = factor*-1
    if plink_info['EA_B'].item() != second_snp.ea:
        print('EA_B != clickhouse')
        plink_info['EA_B'], plink_info['RA_B'] = plink_info['RA_B'].item(), plink_info['EA_B'].item()
        plink_info['EAF_B'] = 1 - plink_info['EAF_B'].item()
        factor = factor*-1

    plink_info['R'] = plink_info['R'].item()*factor

    return plink_info


def get_query_table():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    results = parser.parse_args()
    return results.input


if __name__ == "__main__":
    input_path = get_query_table()
    query_table = pd.read_table(input_path, header=None, sep='\t')
    rs_pairs = [data.tolist() for index, data in query_table.iterrows()]
    list_of_dfs = [generate_df_for_snp_pair(pair[0], pair[1]) for pair in rs_pairs]
    df = pd.concat(list_of_dfs)
    df.to_csv('ld_table.csv', sep='\t', index=False)

