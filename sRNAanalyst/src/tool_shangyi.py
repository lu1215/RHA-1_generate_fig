import argparse
import pandas as pd 
from collections import defaultdict
import json


bwt_dict = {
    0: 'read_id',   # read id and read count
    1: 'reverse_complement',    # +: means no reverse complement, -: means reverse complement
    2: 'ref_id',     # References sequence NAME
    3: 'target_pos', # 0-based leftmost mapping POSition
    4: 'read_seq',   # segment SEQuence
    5: 'ASCII-encoded read qualities', 
    6: 'read mutation times',   # times of same read and same mutation condition
    7: 'mutation_str'
}

def bwt_to_df(data):
    df = pd.read_csv(data, sep='\t', header=None, low_memory=False)
    # df = df[[0, 1, 2, 3, 4, 7]]
    df = df.rename(columns=bwt_dict)
    df[['read_id', 'read_count']] = df['read_id'].str.split('|', expand=True)
    df['read_count'] = pd.to_numeric(df['read_count'], errors='coerce')
    df['reverse_complement'] = df['reverse_complement'].apply(lambda x: True if x == '-' else False)
    df['read_len'] = df["read_seq"].apply(len)
    df['end_pos'] = df['target_pos'] + df['read_len'] - 1
    df['init_pos'] = df['target_pos']
    return df

def bwt_to_df_yiting_file(data):
    df_data = pd.read_csv(data, sep='\t', header=None, low_memory=False)
    # df = df[[0, 1, 2, 3, 4, 7]]
    df_data = df_data.rename(columns=bwt_dict)
    print(data)
    # csv = data.name.replace(".bwt", "_22g.csv")
    print(data.name.split('_'))
    if 'mis' in data.name.split('/')[-1]:
        replace_str = '_' + data.name.split('_')[-1]
    else:
        replace_str = ".bwt"
    print(replace_str)
    csv = data.name.replace(replace_str, "_22g.csv")
    print(csv)
    
    df_csv = pd.read_csv(csv)
    # merge readcount column
    df_data = pd.merge(df_data, df_csv, left_on='read_seq', right_on='input_seq', how='inner')
    df_data.drop(columns=["input_seq"], inplace=True)
    # distribute readcount
    df_data['seq_count'] = df_data.groupby('read_seq')['read_seq'].transform('count')
    df_data["read_count"] = df_data['read_count'] / df_data['seq_count']
    df_data['read_len'] = df_data["read_seq"].apply(len)
    # +1 because 0-indexed in .bwt file
    df_data['end_pos'] = df_data['target_pos'] + df_data['read_len'] - 1 + 1
    df_data['init_pos'] = df_data['target_pos'] + 1
    return df_data

def filter_from_list(target_list: str, df: pd.DataFrame) -> pd.DataFrame:
    with open(target_list, 'r') as f:
        target_dict = { line.strip(): True for line in f}
    df_filtered = df[df['ref_id'].isin(target_dict)]
    return df_filtered

def filter_len(df: pd.DataFrame, length: int) -> pd.DataFrame:
    if 'read_len' in df.columns:
        return df[df["read_len"] == length]
    else:
        df['read_len'] = df["read_seq"].apply(len)
        return df[df["read_len"] == length]
        
def calculating_mutation_distirbute(df):
    mutation_count = defaultdict(float)
    ## drop nan, prevent error occur ##
    # df = df[df['mutation_str'] != 'nan']
    total_read_count = df["read_count"].sum()
    df = df.dropna()
    for _, row in df.iterrows():
        mutations = row['mutation_str'].split(',')
        readcount = row['read_count']
        
        for mutation in mutations:
            position = int(mutation.split(':')[0])
            mutation_count[position] += readcount
    result_dict = dict(mutation_count)
    df_ret = pd.DataFrame(
        [(k, v, v / total_read_count) for k, v in result_dict.items()],
        columns=["location", "mutation_readcount", "ratio"]
    )
    # df_ret = pd.DataFrame.from_dict(result_dict, orient='index')
    return df_ret, total_read_count

def calculating_mutation_cond(df) -> dict:
    mutation_count = defaultdict(float)
    ## drop nan, prevent error occur ##
    # df = df[df['mutation_str'] != 'nan']
    df = df.dropna()
    result = {}
    for _, row in df.iterrows():
        mutations = row['mutation_str'].split(',')
        readcount = row['read_count']
        
        for mutation in mutations:
            position, change = mutation.split(':')
            position = int(position)
            
            if position not in result:
                result[position] = {}
            
            if change not in result[position]:
                result[position][change] = 0
            
            # 累加 readcount
            result[position][change] += readcount
    return result

def calculating_read_mutation_cond(df) -> dict:
    df = df.dropna()
    result = {}
    ################################################################################
    ## ignoring reverse complement, We choosing nucleotides for the original read ##
    ################################################################################
    complement_dict = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    for _, row in df.iterrows():
        mutations = row['mutation_str'].split(',')
        readcount = row['read_count']
        
        for mutation in mutations:
            position, change = mutation.split(':')
            changed_nuc = change.split(">")[1]
            changed_nuc = complement_dict[changed_nuc] if row["reverse_complement"] else changed_nuc
            position = int(position)
            
            if position not in result:
                result[position] = {}
            
            if changed_nuc not in result[position]:
                result[position][changed_nuc] = 0
            
            # 累加 readcount
            result[position][changed_nuc] += readcount
    return result

#############################################
# Distribute Read-Count by Duplicated Reads #
#############################################
def distribute(df):

    cmt = "# distribute\n"

    df_tmp = pd.DataFrame()
    df_tmp['dup'] = df.pivot_table(columns=['read_id'], aggfunc='size')
    df_tmp = df_tmp.reset_index()
    df = pd.merge(df, df_tmp, how='left', on='read_id')
    df['read_count'] = df['read_count']/df['dup']
    del df['dup']

    return df, cmt

############################
# Read Count Normalization #
############################
def normalize(df, factor):
    df['read_count'] = df['read_count']*float(factor)
    cmt = "# normalize={}\n".format(factor)
    return df, cmt

def main():
    parser = argparse.ArgumentParser(description="Process some files.")
    
    # 定義全域參數
    parser.add_argument('--flen', action='store_true', help='Enable filter specific len')
    parser.add_argument('--calmut', action='store_true', help='Enable calculating mutation distribute')
    parser.add_argument('--calmutcond', action='store_true', help='Enable calculating mutation condition')
    parser.add_argument('--calrmc', action='store_true', help='Enable calculating mutation condition(choosing original read nuc)')
    parser.add_argument('--len',type=int, help='len of filter')
    parser.add_argument('--fl', action='store_true', help='Enable filter from list')
    parser.add_argument('-l', '--list', type=str, help='Path to the target list file')
    parser.add_argument("-i", "--input",
                        nargs='?',
                        type=argparse.FileType('r'),
                        help="input file in 'csv', 'tsv', 'fasta', 'fastq' or 'sam' format")
    parser.add_argument("-o", "--output",
                        nargs='?',
                        type=argparse.FileType('w'),
                        help="output file in 'csv' or 'fasta' format")
    parser.add_argument('--normalize')
    parser.add_argument('--norm_factor')
    parser.add_argument('--bwt')
    parser.add_argument('--process_yiting_file')
    # 解析參數
    args = parser.parse_args()
    data = args.input if args.input else None
    cmt = ""
    if args.input and not args.bwt:
        for line in data:
            if line[0]=='#': cmt += line   
            else: break
        header = line[:-1].split(',')
        df_in = pd.read_csv(data, comment='#', names=header)

    # 如果 --fl 被指定，則檢查 -l, -i, -o 是否存在並執行函式
    if args.fl:
        if not (args.list and args.input and args.output):
            parser.error("--fl requires -l, -i, and -o to be specified")
        ret_df = filter_from_list(df=df_in, target_list = args.list)
        args.output.write(cmt)
        ret_df.to_csv(args.output, index=False)
    
    if args.flen:
        if not (args.input and args.output and args.len):
            parser.error("--flen requires -i, and -o to be specified")
        ret_df = filter_len(df=df_in, length=args.len)
        args.output.write(cmt)
        ret_df.to_csv(args.output, index=False)
    
    if args.calmut:
        if not (args.input and args.output):
            parser.error("--calmut requires -i, and -o to be specified")
        ret_df, total_read_count = calculating_mutation_distirbute(df=df_in)
        cmt += '# total_read_count={}\n'.format(total_read_count)
        args.output.write(cmt)
        ret_df.to_csv(args.output, index=False)
        # ret_dict = calculating_mutation_distirbute(df=df_in)
        # json.dump(ret_dict, args.output, indent=4)
        # args.output.close()

    if args.calmutcond:
        if not (args.input and args.output):
            parser.error("--calmutcond requires -i, and -o to be specified")
        ret_dict = calculating_mutation_cond(df=df_in)
        json.dump(ret_dict, args.output, indent=4)
        args.output.close()
    
    if args.calrmc:
        if not (args.input and args.output):
            parser.error("--calrmc requires -i, and -o to be specified")
        ret_dict = calculating_read_mutation_cond(df=df_in)
        json.dump(ret_dict, args.output, indent=4)
        args.output.close()

    if args.normalize and args.norm_factor and args.input and args.output and args.bwt:
        print(args.input)
        df = bwt_to_df(args.input)
        cmt = ""
        df, tmp = distribute(df)
        cmt += tmp
        df, tmp = normalize(df, args.norm_factor)
        cmt += tmp
        args.output.write(cmt)
        df.to_csv(args.output, index=False)
    
    if args.process_yiting_file and args.input and args.output and args.norm_factor:
        ## 2.7822835313855494
        cmt = "# normalize={}\n".format(args.norm_factor)
        df = bwt_to_df_yiting_file(args.input)
        args.output.write(cmt)
        df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
