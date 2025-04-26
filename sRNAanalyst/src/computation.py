##################################
# Copyright (C) 2023 Ryan Chung  #
#                                #
# History:                       #
# 2023/04/10 Ryan Chung          #
##################################
 
import itertools
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import pandas as pd
import seaborn as sns
import sys
from tqdm import tqdm
if sys.version_info >= (3,6):
    from statannotations.Annotator import Annotator
else:
    from statannot import add_stat_annotation

#########################
# Analyze Single Region #
#########################
def analyze_single_region(
    read_df=None,      # dataframe of reads
    ref_df=None,       # dataframe of regions
    column=None,       # column name which want to analyze
    set_region=False,  # output with "_region" column
    set_count=False,   # output with "_count" column
    set_length=False,  # output with "_length" column
    set_density=False, # output with "_density" column
):       
    
    # inner function to find length of overlap-region
    # x is a dataframe with four columns: 
    # [init_pos_x, end_pos_x, init_pos_y, end_pos_y] 
    def _overlap(x):
        sorted_x = sorted(x)
        return sorted_x[2] - sorted_x[1] + 1
        
    # split positions from column of ref_df
    position = ref_df[column].str.split('-', expand=True)
    df_ref = pd.DataFrame()
    if set_region:
        df_ref[column] = ref_df[column].copy()
    df_ref['init_pos'] = position[0].fillna(0).astype('uint32')
    df_ref['end_pos'] = position[1].fillna(0).astype('uint32')
    df_ref['ref_id'] = ref_df['ref_id'].copy()
    df_ref.set_index('ref_id', inplace=True)
        
    # copy columns from read_df
    df = read_df[['ref_id','read_count','init_pos','end_pos']]
        
    # merge two dataframe
    # read_df: init_pos_x, end_pos_x
    # ref_df: init_pos_y, end_pos_y
    df = pd.merge(df, df_ref, on='ref_id')
    df['len'] = df['end_pos_x'] - df['init_pos_x'] + 1
    df['len_overlap'] = 0
        
    # find reads that overlap with region
    df_tmp = df.loc[
        (df['init_pos_x']<=df['end_pos_y']) & (df['init_pos_y']<=df['end_pos_x']),
        ['init_pos_x','end_pos_x','init_pos_y','end_pos_y']
    ]
    len_overlap = df_tmp.apply(lambda x: _overlap(x), axis=1)
    df.loc[len_overlap.index, 'len_overlap'] = len_overlap
        
    # calculate the ratio of read-count
    df[column+'_count'] = df['read_count'] * df['len_overlap'] / df['len']
        
    # merge read-count to ref_df
    df = df.groupby('ref_id')[column+'_count'].sum()
    target_df = pd.merge(df_ref, df, on='ref_id', how='left')
        
    # add columns of 'region', read-count', 'length', 'density'
    target_df[column+'_count'].fillna(0, inplace=True)
    if set_length:
        target_df.loc[ target_df['init_pos']>0, column+'_len'] = target_df['end_pos'] - target_df['init_pos'] + 1
    if set_density:
        target_df.loc[ target_df['init_pos']>0, column+'_den'] = target_df[column+'_count'] / (target_df['end_pos'] - target_df['init_pos'] + 1)
    if set_count:
        target_df.loc[ target_df['init_pos']==0, column+'_count'] = np.nan
    else:
        target_df.drop(columns=[column+'_count'], inplace=True)
        
    del target_df['init_pos']
    del target_df['end_pos']

    return target_df


############################
# Analyze Multiple Regions #
############################
def analyze_multiple_region(
    read_df=None,      # dataframe of reads
    ref_df=None,       # dataframe of regions
    columns=None,      # column names which want to analyze
    set_region=False,  # output with "_region" column
    set_count=False,   # output with "_count" column
    set_length=False,  # output with "_length" column
    set_density=False, # output with "_density" column
):
    # initailize
    if columns==None:
        columns = ref_df.columns.tolist()
        columns.remove('ref_id')

    # multiprocessing
    cpu = multiprocessing.cpu_count()
    if cpu > 8:
        cpu = int(cpu/2) 
    elif cpu > 1:
        cpu = cpu - 1
    pool = multiprocessing.Pool(cpu)
    args = [(read_df,ref_df,col,set_region,set_count,set_length,set_density) for col in columns]
    df_list = pool.starmap(analyze_single_region, tqdm(args))
    pool.close()
    pool.join()

    return df_list


###########################
# Analyze Single Position #
###########################
def analyze_single_position(
    read_df=None, # dataframe of reads
    ref_df=None,  # dataframe of regions
    column=None,  # column name which want to analyze
    limit=None,   # left-limit and right-limit in list
):

    # inner function to find positioins of overlap-region
    # x is a dataframe with four columns: 
    # [init_pos_x, end_pos_x, init_pos_y, end_pos_y]
    def _overlap(x):
        sorted_x = sorted(x)
        return sorted_x[1], sorted_x[2]

    # initailize
    left_limit,right_limit = limit
    df_ref = ref_df.copy()
    df_ref['index'] = df_ref.index
    try:
        df_ref['target'] = df_ref[column]
        df_ref['init_pos'] = df_ref['target'] + left_limit
        df_ref['end_pos'] = df_ref['target'] + right_limit
    except:
        print('[Error] Column "{}" is not found.'.format(column))
        sys.exit(1)
        
    # check window is in mRNA boundary (head,tail), 1-based
    check_boundary = True 
    if 'head' not in df_ref.columns:
        df_ref['head'] = 1
    if 'tail' not in df_ref.columns:
        if 'length' not in df_ref.columns:
            print('[Warning] Columns "tail" or "length" are not found, won\'t check mRNA boundary.')
            check_boundary = False
        else:
            df_ref['tail'] = df_ref['length']
    if check_boundary:
        df_ref.loc[ df_ref['init_pos']<df_ref['head'], 'init_pos'] = df_ref['head']
        df_ref.loc[ df_ref['end_pos']>df_ref['tail'], 'end_pos'] = df_ref['tail']
        
    # remain necessary columns
    df_ref = df_ref[['ref_id','index','init_pos','end_pos','head','tail','target']]
    df = read_df[['ref_id','read_count','init_pos','end_pos']]
        
    # merge two dataframe
    # read_df: init_pos_x, end_pos_x
    # ref_df: init_pos_y, end_pos_y
    df = pd.merge(df, df_ref, on='ref_id')
        
    # find reads that overlap with region
    df = df[ (df['init_pos_x']<=df['end_pos_y']) & (df['init_pos_y']<=df['end_pos_x']) ]
        
    # shift index as 0-based array
    if check_boundary:
        for col in ['head','tail']:
            df_ref[col] = df_ref[col] - df_ref['target'] - left_limit
    for col in ['init_pos_x','init_pos_y','end_pos_x','end_pos_y']:
        df[col] = df[col] - df['target'] - left_limit
        
    # find overlap positoins
    df_tmp = df[['init_pos_x','end_pos_x','init_pos_y','end_pos_y']].copy()
    overlap = df_tmp.apply(lambda x: _overlap(x), axis=1)
    df_tmp2 = pd.DataFrame(overlap.to_list(), index=df_tmp.index)
    df_tmp['init_pos'] = df_tmp2[0]   # 0-based initial position
    df_tmp['end_pos'] = df_tmp2[1]+1  # terminal of end position

    # remain necessary columns
    # [ref_id, index, read-count, init_pos, end_pos]
    df = pd.concat([df, df_tmp], axis=1).reset_index(drop=True)
    df.drop(columns=['init_pos_x','end_pos_x','init_pos_y','end_pos_y'], inplace=True)

    # build array for boundary plot
    arr = np.zeros((len(df_ref),(right_limit-left_limit+1)))
    for i in range(len(df)):
        arr[df.at[i,'index'], df.at[i,'init_pos']:df.at[i,'end_pos']] += df.at[i,'read_count']

    # build dataframe
    # set values outside the boundaries to NaN
    target_df = pd.DataFrame(arr, index=df_ref['ref_id'], 
                            columns=[*range(right_limit-left_limit+1)]).round(3)
    if check_boundary:
        mask = ((target_df.columns.values < df_ref['head'].values[:, None]) |
                (target_df.columns.values > df_ref['tail'].values[:, None]))
        target_df[mask] = np.nan
    target_df.columns = [*range(left_limit,right_limit+1)]
    target_df = target_df.reset_index()

    return target_df


##############################
# Analyze Multiple Positions #
##############################
def analyze_multiple_position(
    read_df=None, # dataframe of reads
    ref_df=None,  # dataframe of regions
    columns=None, # column names which want to analyze
    limit=None,   # left-limit and right-limit in list
):
    # multiprocessing
    cpu = multiprocessing.cpu_count()
    if cpu > 8:
        cpu = int(cpu/2) 
    elif cpu > 1:
        cpu = cpu - 1
    pool = multiprocessing.Pool(cpu)
    args = [(read_df,ref_df,col,lim) for col,lim in zip(columns,limit)]
    df_list = pool.starmap(analyze_single_position, tqdm(args))
    pool.close()
    pool.join()

    return df_list


############
# Box Plot #
############
def box_plot(
    data=None,            # dataframe of reads
    filter=None,          # reference names which want to analyze
    columns=None,         # column names which want to analyze
    title=None,           # title of figure
    detail=False,         # show mean and median at x-ticks
    scale='normal',       # scale of y-axis
    test=None,            # hypothesis testing
    test_format='simple', # format of testing
    x_axis='data',        # field of x-axis in condition 4
    y_label=None,         # label of y-axis
    style='darkgrid',     # backgroud style of figure
    color='coolwarm',     # color palette of figure
):      
    # determine the condition
    # 1. One-One:     data(1) and filter(1)
    # 2. One-Multi:   data(1) and filter(N)
    # 3. Multi-One:   data(N) and filter(1)
    # 4. Multi-Multi: data(N) and filter(N)
    data_num = len(data) if data else 0
    filter_num = len(filter) if filter else 0
    if data_num<=1 and filter_num<=1:
        condition = 1
        x = 'region'
        hue = None
    elif data_num<=1:
        condition = 2
        x = 'filter'
        hue = None
    elif filter_num<=1:
        condition = 3
        x = 'data'
        hue = None
    else:
        condition = 4
        if x_axis=='data':
            x = 'data'
            hue = 'filter'
        elif x_axis=='filter':
            x = 'filter'
            hue = 'data'
        else:
            print("[Error]")
            print("Wrong value (x_axis={}) of Box_Plot in stylesheet.yml".format(x_axis))
            sys.exit(1)

    # initialize necessary parameters
    # in order to prevent misuse of the configuration file
    detail = False if detail is None else detail
    scale = 'normal' if scale is None else scale
    test_format = 'simple' if test_format is None else test_format
    x_axis = 'data' if x_axis is None else x_axis
    style = 'darkgrid' if style is None else style
    color = 'coolwarm' if color is None else color
    if columns is None:
        # all columns in dataframe
        if condition <= 2:
            columns = data[0]['dataframe'].columns.tolist()
            columns.remove('ref_id')
        # intersection columns in dataframes
        else:
            columns = []
            for d in data:
                if columns:
                    columns = list(set(columns) & set(d['dataframe'].columns.tolist()))
                else:
                    columns = d['dataframe'].columns.tolist()
            columns.remove('ref_id')
    orignal_columns = [col.split('_')[0] for col in columns]

    # merge different dataframes frome wide-format to long-format
    tmp_df = pd.DataFrame()
    for d in data:
        df = d['dataframe']
        df = df.melt(id_vars=['ref_id'], value_vars=columns, var_name='region')
        df = df.dropna().reset_index(drop=True)
        df['region'] = df['region'].str.split('_',expand=True)[0]
        df['data'] = d['name']
        tmp_df = pd.concat([tmp_df,df])
    tmp_df = tmp_df.reset_index(drop=True)
    columns = orignal_columns

    # merge dataframe with different filters
    plot_df = pd.DataFrame()
    if filter:
        for f in filter:
            df = tmp_df[ tmp_df['ref_id'].isin(f['id']) ].reset_index(drop=True)
            df['filter'] = f['name']
            plot_df = pd.concat([plot_df,df])
    else:
        plot_df = tmp_df
        plot_df['filter'] = 'all_mRNA'
    plot_df = plot_df.reset_index(drop=True)

    if scale=='log2':
        plot_df['value'] = plot_df['value'].apply(lambda x: np.log2(x) if x!=0 else x)
    elif scale=='log10':
        plot_df['value'] = plot_df['value'].apply(lambda x: np.log10(x*pow(10,6)) if x!=0 else x)

    sns.set(style=style)
    ax_num = 1 if condition==1 else len(columns)
    fig, axes = plt.subplots(1, ax_num, figsize=(1+5*ax_num, 5), dpi=200, sharey=True)

    # set limit for y-ticks
    max_yticks, min_yticks = 0,0
    for i in range(ax_num):
        sub_df = plot_df if condition==1 else plot_df[ plot_df['region']==columns[i] ]
        x_lst = sub_df[x].drop_duplicates().to_list()
        for j in range(len(x_lst)):
            tmp = sub_df[ sub_df[x]==x_lst[j] ]
            max_val = tmp['value'].max()
            min_val = tmp['value'].min()
            q1 = tmp['value'].quantile(0.25)
            q3 = tmp['value'].quantile(0.75)
            iqr = q3 - q1
            upper_limit = q3 + 1.5 * iqr
            upper_limit = max_val if max_val < upper_limit else upper_limit
            lower_limit = q1 - 1.5 * iqr
            lower_limit = min_val if min_val > lower_limit else lower_limit
            max_yticks = upper_limit if upper_limit > max_yticks else max_yticks
            min_yticks = lower_limit if lower_limit < min_yticks else min_yticks
    ymax = max_yticks + 0.05*(max_yticks-min_yticks)
    ymin = min_yticks - 0.05*(max_yticks-min_yticks)
    plt.ylim(ymin, ymax)

    # plot figure
    for i in range(ax_num):
        ax = axes if ax_num<=1 else axes[i]
        sub_df = plot_df if condition==1 else plot_df[ plot_df['region']==columns[i] ]
        sns.boxplot(ax=ax, data=sub_df, x=x, y='value', hue=hue,
                    width=0.3, showfliers=False, showmeans=False,
                    medianprops=dict(color='orange'), palette=color)
        ax.set_xlabel('')
        ax.set_ylabel('')
        if i==0:
            ax.set_ylabel(y_label, fontsize=14)
        if condition==4:
            if i==ax_num-1:
                ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            else:
                ax.get_legend().remove()
        if ax_num==1:
            ax.set_title('')
        else:
            ax.set_title(orignal_columns[i])

        # add test to axes
        add_test = False
        pairs = []
        xticks = [xtick.get_text() for xtick in ax.get_xticklabels()]
        if test!=None and test!='None' and (condition>1 or len(columns)>1):
            if test=='U-test':
                add_test = True
                test_method = 'Mann-Whitney'
            elif test=='T-test':
                add_test = True
                test_method = 't-test_welch'
            else:
                print("Wrong Testing Method")
            if add_test:
                try:
                    pairs = list(itertools.combinations(xticks, 2))
                    if condition==4:
                        hue_labels = sub_df[hue].drop_duplicates().to_list()
                        #pairs = [tuple((x, hue) for x in pair) for pair in pairs for hue in hue_labels]
                        hue_labels = list(itertools.combinations(hue_labels, 2))
                        pairs = [tuple((xtick, h) for h in hue) for hue in hue_labels for xtick in xticks]
                        ax_title = orignal_columns[i] + '\n'*(len(pairs)-((len(xticks)-1)*len(hue_labels)))*2
                    else:
                        ax_title = orignal_columns[i] + '\n'*len(pairs)*2
                    if ax_num!=1:
                        ax.set_title(ax_title)
                    
                    if sys.version_info >= (3,6):    
                        annot = Annotator(ax, pairs, data=sub_df, x=x, y='value', hue=hue)
                        annot.configure(test=test_method, text_format=test_format, loc='outside', verbose=0)
                        annot.apply_and_annotate()
                    else:
                        add_stat_annotation(ax, data=sub_df, x=x, y='value', hue=hue,
                                            box_pairs=pairs, loc='outside', verbose=0,
                                            test=test_method, text_format=test_format)

                except ValueError as error:
                    print("[Error] {}".format(error))
    
        # add mean and median to x-ticks
        if detail:
            num_lst = [len(sub_df[sub_df[x]==xtick]) for xtick in xticks]
            mean_list = [round(sub_df.loc[sub_df[x]==xtick,'value'].mean(),3) for xtick in xticks]
            median_list = [round(sub_df.loc[sub_df[x]==xtick,'value'].median(),3) for xtick in xticks]
            new_xticks = ["{}\nnumber:{}\navg:{}\nmedian:{}".format(xtick,num,mean,median) for (xtick,num,mean,median) in zip(xticks,num_lst,mean_list,median_list)]
            xticks_dict = dict(zip(xticks, new_xticks))
            ax.set_xticklabels([ xticks_dict[xtick.get_text()] for xtick in ax.get_xticklabels() ], fontsize=10)
    
    if title:
        if detail:
            fig.suptitle(title, y=-0.07, fontsize=16)
        else:
            fig.suptitle(title, y=0.01, fontsize=16)
    
    # change gap size between subplots
    plt.subplots_adjust(wspace=0.1)

    return fig, plot_df


#############
# Line Plot #
#############
def line_plot(
    data=None,          # dataframe of reads
    filter=None,        # reference names which want to analyze
    columns=None,       # column names which want to analyze
    hue=None,           # hue of figure  in condition 4
    title=None,         # title of figure
    xlabel=None,        # xlabel of figure
    vertical_line=None, # positions of vertical line in list
    style='darkgrid',   # backgroud style of figure
    color='deep',       # color palette of figure
):
    # count average and standard-error of a dataframe
    def _count_avg(data):
        data.columns = data.columns.astype(int)
        df = pd.DataFrame()
        df['ratio'] = data[data!=0].count().sort_index() / data.count().sort_index()
        df['avg'] = data.mean(axis = 0).sort_index()
        df['avg_dis'] = df['avg'] / sum(df['avg'])
        df['ste'] = np.std(data, ddof=1) / np.sqrt(len(data))
        df['avg_plus_ste'] = df['avg'] + df['ste']
        df['avg_minus_ste'] = df['avg'] - df['ste']
        df['avg_plus_ste_dis'] = df['avg_plus_ste'] / sum(df['avg'])
        df['avg_minus_ste_dis'] = df['avg_minus_ste'] / sum(df['avg'])
        df = df.reset_index()
        df['index'] = df['index']
        return df
    
    # determine the condition
    # 1. One-One:     data(1) and filter(1)
    # 2. One-Multi:   data(1) and filter(N)
    # 3. Multi-One:   data(N) and filter(1)
    # 4. Multi-Multi: data(N) and filter(N)
    data_num = len(data) if data else 0
    filter_num = len(filter) if filter else 0
    if data_num<=1 and filter_num<=1:
        condition = 1
        hue = 'data'
        hue2 = None
    elif data_num<=1:
        condition = 2
        hue = 'filter'
        hue2 = None
    elif filter_num<=1:
        condition = 3
        hue = 'data'
        hue2 = None
    else:
        condition = 4
        if hue==None:
            hue = 'data'
            hue2 = 'filter'
        elif hue=='data':
            hue2 = 'filter'
        elif hue=='filter':
            hue2 = 'data'
        else:
            print("[Error]")
            print("Wrong value (hue={}) of Line_Plot in stylesheet.yml".format(hue))
            sys.exit(1)

    # initialize necessary parameters
    if columns is None:
        # all columns in dataframe
        if condition <= 2:
            columns = data[0]['dataframe'].columns.tolist()
            columns.remove('ref_id')
        # intersection columns in dataframes
        else:
            columns = []
            for d in data:
                if columns:
                    columns = list(set(columns) & set(d['dataframe'].columns.tolist()))
                else:
                    columns = d['dataframe'].columns.tolist()
            columns.remove('ref_id')
    columns_num = len(columns)

    # merge different dataframes and filters
    plot_df = pd.DataFrame()
    for d in data:
        tmp_df = d['dataframe']
        if filter:
            for f in filter:
                df = tmp_df.loc[tmp_df['ref_id'].isin(f['id']), columns].reset_index(drop=True)
                df = _count_avg(df)
                df['data'] = d['name']
                df['filter'] = f['name']
                plot_df = pd.concat([plot_df,df])
        else:
            df = tmp_df[columns]
            df = _count_avg(df)
            df['data'] = d['name']
            df['filter'] = 'all_mRNA'
            plot_df = pd.concat([plot_df,df])
    plot_df = plot_df.reset_index(drop=True)

    name = plot_df[hue].drop_duplicates().to_list()
    if len(name) == 2:
        # only two data
        palette = ['blue', 'red']
    else:
        # more than two data using sns.color_palette
        palette = sns.color_palette(color, len(name)).as_hex()
    # plot figure
    sns.set(style=style)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(7,12), dpi=200)
    ax1 = sns.lineplot(ax=ax1, data=plot_df, x='index', y='ratio', hue=hue, style=hue2, palette=palette, ci=None)
    ax2 = sns.lineplot(ax=ax2, data=plot_df, x='index', y='avg', hue=hue, style=hue2, palette=palette, ci=None)
    ax3 = sns.lineplot(ax=ax3, data=plot_df, x='index', y='avg_dis', hue=hue, style=hue2, palette=palette, ci=None)
    
    # fill colors between (+/- 1STE)
    
    # palette = sns.color_palette(color,len(name)).as_hex()
    for i in range(int(len(plot_df)/columns_num)):
        start = i*columns_num
        end = (i+1)*columns_num
        tmp = plot_df[start:end]
        c = palette[name.index(tmp[hue].iloc[0])]
        ax2.fill_between(x=tmp['index'], y1=tmp['avg_plus_ste'], y2=tmp['avg_minus_ste'], alpha=0.3, color=c)
        ax3.fill_between(x=tmp['index'], y1=tmp['avg_plus_ste_dis'], y2=tmp['avg_minus_ste_dis'], alpha=0.3, color=c)
        
    ax1.tick_params(axis='y', labelsize=12)
    ax1.set_ylabel('binding ratio', fontsize=14)
    ax1.set_title("mRNA site ratio (site/position)", fontsize=14)

    ax2.tick_params(axis='y', labelsize=12)
    ax2.set_ylabel('read counts', fontsize=14)
    ax2.set_title("Average (+/- 1STE)", fontsize=14)

    ax3.tick_params(axis='x', labelsize=12)
    ax3.tick_params(axis='y', labelsize=12)
    ax3.set_xlabel(xlabel, fontsize=14)
    ax3.set_ylabel('read counts distribution', fontsize=14)
    ax3.set_title("Average (+/- 1STE)", fontsize=14)

    if condition==1:
        ax1.get_legend().remove()
        ax2.get_legend().remove()
        ax3.get_legend().remove() 
    else:
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        ax2.get_legend().remove()
        ax3.get_legend().remove()   

    if vertical_line:
        for pos in vertical_line:
            ax1.axvline(x=pos,c='k',linestyle='dashed')
            ax2.axvline(x=pos,c='k',linestyle='dashed')
            ax3.axvline(x=pos,c='k',linestyle='dashed')
    
    if title:
        fig.suptitle(title, y=0.05, fontsize=16)

    return fig, plot_df


################
# Scatter Plot #
################
def scatter_plot(
    data=None,        # dataframe of reads
    filter=None,      # reference names which want to analyze
    columns=None,     # column names which want to analyze
    hue=None,         # hue of figure  in condition 4
    title=None,       # title of figure
    scale='log2',     # scale of y-axis
    x_label=None,     # label of x-axis
    y_label=None,     # label of y-axis
    show_others=False,# show other values which not in filters
    style='darkgrid', # backgroud style of figure
    color='deep',     # color palette of figure
    BCV=None,         # biological coefficient of variation 
):
    # determine the condition
    # 1. One-One:     data(1) and filter(1)
    # 2. One-Multi:   data(1) and filter(N)
    # 3. Multi-One:   data(N) and filter(1)
    # 4. Multi-Multi: data(N) and filter(N)
    data_num = len(data) if data else 0
    filter_num = len(filter) if filter else 0
    if data_num<=1:
        condition = 2
        hue = 'filter'
        hue2 = None
    else:
        condition = 4
        if hue==None:
            hue = 'data'
            hue2 = 'filter'
        elif hue=='data':
            hue2 = 'filter'
        elif hue=='filter':
            hue2 = 'data'
        else:
            print("[Error]")
            print("Wrong value (hue={}) of Line_Plot in stylesheet.yml".format(hue))
            sys.exit(1)

    # initialize necessary parameters
    # in order to prevent misuse of the configuration file
    if columns is None:
        # all columns in dataframe
        if condition <= 2:
            columns = data[0]['dataframe'].columns.tolist()
            columns.remove('ref_id')
        # intersection columns in dataframes
        else:
            columns = []
            for d in data:
                if columns:
                    columns = list(set(columns) & set(d['dataframe'].columns.tolist()))
                else:
                    columns = d['dataframe'].columns.tolist()
            columns.remove('ref_id')
    orignal_columns = [col.split('_')[0] for col in columns]
    
    # merge different dataframes frome wide-format to long-format
    tmp_df = pd.DataFrame()
    for d in data:
        df = d['dataframe']
        tmp_df2 = pd.DataFrame()
        for col in orignal_columns:
            new_columns = {
                col+'_count': 'x',
                col+'_count_y': 'y',
                col+'_filter': 'filter'
            }
            df2 = df[['ref_id',col+'_count',col+'_count_y',col+'_filter']]
            df2 = df2.rename(columns=new_columns)
            df2['region'] = col
            df2 = df2.dropna().reset_index(drop=True)
            tmp_df2 = pd.concat([tmp_df2,df2])
        tmp_df2['data'] = d['name']
        tmp_df = pd.concat([tmp_df,tmp_df2])
    tmp_df = tmp_df.reset_index(drop=True)
    
    # merge dataframe with different filters
    plot_df = pd.DataFrame()
    valid = 'p < 0.05' if (BCV and BCV!='None') else 'valid'
    invalid_df = tmp_df[tmp_df['filter']!=valid].reset_index(drop=True)
    valid_df = tmp_df[tmp_df['filter']==valid].reset_index(drop=True)
    valid_df['others'] = True
    if filter:
        for f in filter:
            valid_df.loc[ valid_df['ref_id'].isin(f['id']), 'others'] = False
            df = valid_df[ valid_df['ref_id'].isin(f['id']) ].drop(columns='others').reset_index(drop=True)
            df['filter'] = f['name']
            plot_df = pd.concat([plot_df,df])
        df = valid_df[ valid_df['others']==True ].drop(columns='others').reset_index(drop=True)
        if show_others:
            df['filter'] = 'others'
            plot_df = pd.concat([df,plot_df])
    else:
        plot_df = valid_df
    plot_df = pd.concat([invalid_df,plot_df]).reset_index(drop=True) 

    if scale=='log2':
        for col in ['x','y']:
            plot_df[col] = plot_df[col].apply(lambda x: np.log2(x) if x!=0 else x)
    elif scale=='log10':
        for col in ['x','y']:
            plot_df[col] = plot_df[col].apply(lambda x: np.log10(x) if x!=0 else x)
    
    # plot figure
    sns.set(style=style)
    ax_num = 1 if condition==1 else len(columns)
    # print("ax_num: ", ax_num)
    # # setting colors
    # if ax_num == 2:
    #     # only two data
    #     palette = ['blue', 'black']
    # else:
    #     # more than two data using sns.color_palette
    #     palette = color
    palette = ['gray', 'black']
    fig = plt.figure(figsize=(1+5*ax_num, 5), dpi=200)
    for i in range(ax_num):
        sub_df = plot_df if condition==1 else plot_df[ plot_df['region']==orignal_columns[i] ]
        ax = fig.add_subplot(1,ax_num,i+1)
        ax = sns.scatterplot(data=sub_df, x='x', y='y', hue=hue, style=hue2, palette=palette, linewidth=0.3, s=10)

        # add base line
        diag_x = list(ax.get_xlim())
        diag_y = list(ax.get_ylim())
        ax.plot(diag_x, diag_y, c='black')
        ax.plot([diag_x[0]+1, diag_x[1]], [diag_y[0], diag_y[1]-1], c='black' ,linestyle="--")
        ax.plot([diag_x[0], diag_x[1]-1], [diag_y[0]+1, diag_y[1]], c='black' ,linestyle="--")
        
        ax.set_xlabel('')
        ax.set_ylabel('')
        if i==0:
            ax.set_ylabel(y_label, fontsize=14)
        if condition>1:
            if i==ax_num-1:
                ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            else:
                ax.get_legend().remove()
            ax.set_title(orignal_columns[i])
    fig.text(0.5, -0.01, x_label, ha='center',size=14)
    if title:
        fig.suptitle(title, y=-0.05, fontsize=16)

    return fig, plot_df
    

