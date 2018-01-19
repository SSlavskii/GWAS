import pandahouse as ph
import subprocess
import os

CLICKHOUSE_CONNECTION_PARAMS = {'host': 'http://localhost:8123', 'database': 'gwas'}

if os.environ.get("TEST"):
    CLICK_DB = "test"
    CLICK_TABLE = "gwas_snp"
else:
    CLICK_DB = "gwas"
    CLICK_TABLE = "test_gwas_snp"


snp1 = 12740374
snp2 = 602633
chr = '1'


def get_z_correlation(sr_id1, sr_id2):
    query = "SELECT corr(z_1, z_2)\
       FROM (\
            SELECT gwas_1.gwas_id, gwas_1.z AS z_1\
            FROM {0}.{1} as gwas_1\
            WHERE gwas_1.snp_num = {2}\
        ) ALL INNER JOIN (\
            SELECT gwas_2.gwas_id, gwas_2.z AS z_2\
            FROM {0}.{1} as gwas_2\
            WHERE gwas_2.snp_num = {3}\
        ) USING gwas_id".format(CLICK_DB, CLICK_TABLE, sr_id1, sr_id2)
    z_corr = ph.read_clickhouse(query, connection=CLICKHOUSE_CONNECTION_PARAMS)
    return z_corr


def get_snp_alleles(sr_id):
    query = "select snp_id, ea, ra from {0}.{1} where snp_num=={2}  limit 1".format(CLICK_DB, CLICK_TABLE, sr_id)
    print(query)
    eff_allele = ph.read_clickhouse(query, connection=CLICKHOUSE_CONNECTION_PARAMS)
    return eff_allele


def get_plink_ld(snp1, snp2):
    subprocess.check_output(['plink', '--bfile',
                            '/home/ubuntu/storage/data/ref/eur_genotype_data/by_chr/ALL.chr' +
                             str(chr) + '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_nodup',
                             '--snps', 'rs' + str(snp1), 'rs' + str(snp2), '-r'])
    return


get_plink_ld(snp1, snp2)
