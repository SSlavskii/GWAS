import pandas as pd
import argparse
from table_creation import get_z_correlation

def get_query_table():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    results = parser.parse_args()
    return results.input





if __name__=="__main__":
    input_path = get_query_table()
    query_table = pd.read_table(input_path, header=None, sep='\t')
    sr_pairs = [data.tolist() for index, data in query_table.iterrows()]
    list_of_dfs = [get_z_correlation(pair[0], pair[1]) for pair in sr_pairs]
    print(list_of_dfs)
