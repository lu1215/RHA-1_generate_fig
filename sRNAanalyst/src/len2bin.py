import argparse
import multiprocessing
import pandas as pd

def generate_bin(df, index, bin_size):
    region = []
    for num in range(bin_size):
        init_pos = int(df.at[index,'bin_len']*num) + 1
        end_pos = int(df.at[index,'bin_len']*(num+1))
        region.append(str(init_pos)+'-'+str(end_pos))
    return region

def change_format(df, bin_size):
    # cpu limitation
    cpu = multiprocessing.cpu_count()
    if cpu > 8:
        cpu = int(cpu/2) 
    elif cpu > 1:
        cpu = cpu - 1

    # multiprocess every bin
    df['bin_len'] = df['length']/bin_size
    args = [(df,i,bin_size) for i in range(len(df))]
    pool = multiprocessing.Pool(cpu)
    bin_lst = pool.starmap(generate_bin,args)
    pool.close()
    pool.join()
    
    meta_df = pd.DataFrame(bin_lst)
    meta_df.columns = meta_df.columns + 1
    meta_df.columns = meta_df.columns.astype(str)
    meta_df['ref_id'] = df['ref_id']
    
    return meta_df

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="This is a tool to change file format from length to regions.",
    )
    parser.add_argument("-i", "--input", help="input file in with 'ref_id' and 'length'")
    parser.add_argument("-o", "--output", help="output file in with 'ref_id' and 'N regions'")
    parser.add_argument("-b", "--bin", default=100, help="set bin size (region size)")
    args = parser.parse_args()
    
    df = pd.read_csv(args.input)
    df = change_format(df, args.bin)
    df.to_csv(args.output, index=False)

    print("Program complete with success.")
