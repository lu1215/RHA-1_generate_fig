##################################
# Copyright (C) 2023 Ryan Chung  #
#                                #
# History:                       #
# 2023/04/10 Ryan Chung          #
##################################

import computation as cp
import itertools
import os
import numpy as np
import pandas as pd
import sys

## Main Tools

################
# Density Plot #
################
def density_plot(run_config=None, plot_config=None):
    
    # analyze
    if run_config and run_config['run']:
    
        config_list = recombine_dict(run_config, iter=['data'])

        for config in config_list:
            data_df, cmt = read_data(config['data'], config['BASE_DIR'])
            reference_df, _ = read_data(config['reference'], config['BASE_DIR'])
            config['comment'] = cmt

            df_list = cp.analyze_multiple_region(
                read_df=data_df,            # dataframe of reads
                ref_df=reference_df,        # dataframe of regions
                columns=config['columns'],  # column names which want to analyze
                set_region=True,            # output with "_region" column
                set_count=True,             # output with "_count" column
                set_density=True,           # output with "_density" column
            )

            column_name = '-'.join(config['columns'])
            metadata = pd.concat(df_list, axis=1).reset_index()
            metapath = os.path.join(run_config['csv_path'],
                                    "{}_region_{}.csv".format(config['data']['name'],column_name))
            if plot_config and plot_config['data']:
                append_config(plot_config['data'], config['data']['name'], path=metapath, dataframe=metadata)
            with open(os.path.join(config['BASE_DIR'],metapath), "w") as f:
                config['comment'] += "# analyze_regions={}\n".format(config['columns'])
                f.write(config['comment'])
                metadata.to_csv(f, index=False)
    
    # plot
    if plot_config and plot_config['plot']:

        # read datas and filters into the configuration file
        for i, data in enumerate(plot_config['data']):
            plot_config['data'][i]['dataframe'], _ = read_data(data, plot_config['BASE_DIR'])
        read_filter(plot_config['filter'], plot_config['BASE_DIR'])
        
        # special setting
        plot_config['columns'] = [ x+'_den' for x in plot_config['columns'] ]
        if plot_config['scale']=='log2': 
            y_label = 'log2( counts/nt )'
        elif plot_config['scale']=='log10': 
            y_label = 'log10( counts*M/nt )'
        else:
            y_label = 'counts/nt'

        seg = plot_config['segmentation']
        if plot_config['merge'] or plot_config['data']==None or plot_config['filter']==None:
            config_list = [ plot_config ]
        elif plot_config['data'] and len(plot_config['data'])<=1:
            config_list = [ plot_config ]
        elif plot_config['filter'] and len(plot_config['filter'])<=1:
            config_list = [ plot_config ]
        else: 
            config_list = recombine_dict(plot_config, iter=[seg])
            for i, config in enumerate(config_list):
                if config['title']:
                    config_list[i]['title'] = config[seg]['name'] + ' ' + config['title']
                else:
                    config_list[i]['title'] = config[seg]['name']
                config_list[i][seg] = [ config[seg] ]

        for i, config in enumerate(config_list):

            fig, plot_df = cp.box_plot(
                data=config['data'],                # dataframe of reads
                filter=config['filter'],            # reference names which want to analyze
                columns=config['columns'],          # column names which want to analyze
                title=config['title'],              # title of figure
                scale=config['scale'],              # scale of y-axis
                test=config['test'],                # hypothesis testing
                test_format=config['test_format'],  # format of testing
                detail=config['detail'],            # show mean and median at x-ticks
                x_axis=config['segmentation'],      # field of x-axis in condition 4
                y_label=y_label,                    # label of y-axis
                style=config['style'],              # backgroud style of figure
                color=config['color'],              # color palette of figure
            )
            fig.savefig(os.path.join(config['BASE_DIR'],"{}_{}.{}".format(config['fig_path'],i,config['fig_format'])), bbox_inches='tight')


#################
# Metagene Plot #
################# 
def metagene_plot(run_config=None, plot_config=None):

    # analyze
    if run_config and run_config['run']:
    
        config_list = recombine_dict(run_config, iter=['data'])

        for config in config_list:
            data_df, cmt = read_data(config['data'], config['BASE_DIR'])
            reference_df, _ = read_data(config['reference'], config['BASE_DIR'])
            config['comment'] = cmt

            df_list = cp.analyze_multiple_region(
                read_df=data_df,           # dataframe of reads
                ref_df=reference_df,       # dataframe of regions
                set_count=True,            # output with "_count" column
            )
            metadata = pd.concat(df_list, axis=1).reset_index()

            # remove "_count" at the end of column's names
            metadata.columns = [ col.split('_')[0] if col!='ref_id' else col for col in metadata.columns ]
            metapath = os.path.join(run_config['csv_path'],
                                    "{}_metagene.csv".format(config['data']['name']))
            if plot_config and plot_config['data']:
                append_config(plot_config['data'], config['data']['name'], path=metapath, dataframe=metadata)
            with open(os.path.join(config['BASE_DIR'],metapath), "w") as f:
                config['comment'] += "# metagene_bin=100\n"
                f.write(config['comment'])
                metadata.to_csv(f, index=False)

    # plot
    if plot_config and plot_config['plot']:
        
        # read datas and filters into the configuration file
        for i, data in enumerate(plot_config['data']):
            plot_config['data'][i]['dataframe'], _ = read_data(data, plot_config['BASE_DIR'])
        read_filter(plot_config['filter'], plot_config['BASE_DIR'])

        seg = plot_config['segmentation']
        if plot_config['merge'] or len(plot_config['data'])<=1 or len(plot_config['filter'])<=1:
            config_list = [ plot_config ]
        else: 
            config_list = recombine_dict(plot_config, iter=[seg])
            for i, config in enumerate(config_list):
                if config['title']:
                    config_list[i]['title'] = config[seg]['name'] + ' ' + config['title']
                else:
                    config_list[i]['title'] = config[seg]['name']
                config_list[i][seg] = [ config[seg] ]
            
        for i, config in enumerate(config_list):
            fig, plot_df = cp.line_plot(
                data=config['data'],         # dataframe of reads
                filter=config['filter'],     # reference names which want to analyze
                hue=config['segmentation'],  # hue of figure in condition 4
                title=config['title'],       # title of figure
                style=config['style'],       # backgroud style of figure
                color=config['color'],       # color palette of figure
                xlabel='region(%)',          # xlabel of figure
            )
            fig.savefig(os.path.join(config['BASE_DIR'],"{}_{}.{}".format(config['fig_path'],i,config['fig_format'])), bbox_inches='tight')


#################
# Position Plot #
#################
def position_plot(run_config=None, plot_config=None):

    # analyze
    if run_config and run_config['run']:
    
        config_list = recombine_dict(run_config, iter=['data'])

        for config in config_list:
            data_df, cmt = read_data(config['data'], config['BASE_DIR'])
            reference_df, _ = read_data(config['reference'], config['BASE_DIR'])
            config['comment'] = cmt

            df_list = cp.analyze_multiple_position(
                read_df=data_df,           # dataframe of reads
                ref_df=reference_df,       # dataframe of regions
                columns=config['columns'], # column name (position) which want to analyze
                limit=config['limit'],     # left-limit and right-limit in list
            )
            for df,col,limit in zip(df_list,config['columns'],config['limit']):
                metapath = os.path.join(config['BASE_DIR'],run_config['csv_path'],
                                        "{}_position_{}({}-{}).csv".format(config['data']['name'],col,limit[0],limit[1]))
                if plot_config and plot_config['data']:
                    append_config(plot_config['data'], config['data']['name'], path=metapath, dataframe=df)
                with open(metapath, "w") as f:
                    f.write(config['comment'] + "# position={}({}-{})\n".format(col,limit[0],limit[1]))
                    df.to_csv(f, index=False)

    # plot
    if plot_config and plot_config['plot']:
        
        # read datas and filters into the configuration file
        for i, data in enumerate(plot_config['data']):
            plot_config['data'][i]['dataframe'], _ = read_data(data, plot_config['BASE_DIR'])
        read_filter(plot_config['filter'], plot_config['BASE_DIR'])

        seg = plot_config['segmentation']
        if plot_config['merge'] or len(plot_config['data'])<=1 or len(plot_config['filter'])<=1:
            config_list = [ plot_config ]
        else: 
            config_list = recombine_dict(plot_config, iter=[seg])
            for i, config in enumerate(config_list):
                if config['title']:
                    config_list[i]['title'] = config[seg]['name'] + ' ' + config['title']
                else:
                    config_list[i]['title'] = config[seg]['name']
                config_list[i][seg] = [ config[seg] ]
            
        for i, config in enumerate(config_list):
            fig, plot_df = cp.line_plot(
                data=config['data'],         # dataframe of reads
                filter=config['filter'],     # reference names which want to analyze
                hue=config['segmentation'],  # hue of figure in condition 4
                title=config['title'],       # title of figure
                style=config['style'],       # backgroud style of figure
                color=config['color'],       # color palette of figure
                xlabel='nt',                 # xlabel of figure
                vertical_line=[0],           # positions of vertical line in list
            )
            fig.savefig(os.path.join(config['BASE_DIR'],"{}_{}.{}".format(config['fig_path'],i,config['fig_format'])), bbox_inches='tight')


####################
# Fold Change Plot #
####################
def fold_change_plot(run_config=None, plot_config=None):
    
    # analyze
    if run_config and run_config['run']:
        
        config_list = recombine_dict(run_config, iter=['data'])

        for config in config_list:
            data_df, cmt = read_data(config['data'], config['BASE_DIR'])
            reference_df, _ = read_data(config['reference'], config['BASE_DIR'])
            config['comment'] = cmt

            df_list = cp.analyze_multiple_region(
                read_df=data_df,            # dataframe of reads
                ref_df=reference_df,        # dataframe of regions
                columns=config['columns'],  # column names which want to analyze
                set_region=True,            # output with "_region" column
                set_count=True,             # output with "_count" column
                set_density=True,           # output with "_density" column
            )
            column_name = '-'.join(config['columns'])
            metadata = pd.concat(df_list, axis=1).reset_index()
            metapath = os.path.join(config['BASE_DIR'],run_config['csv_path'],
                                    "{}_region_{}.csv".format(config['data']['name'],column_name))
            if plot_config and plot_config['data']:
                append_config(plot_config['data'], config['data']['name'], path=metapath, dataframe=metadata)
            with open(metapath, "w") as f:
                config['comment'] += "# analyze_regions={}\n".format(config['columns'])
                f.write(config['comment'])
                metadata.to_csv(f, index=False)

    # plot
    if plot_config and plot_config['plot']:

        # read datas and filters into the configuration file
        # find alpha
        alpha = np.inf
        plot_config['columns'] = [ x+'_count' for x in plot_config['columns'] ]
        for i, data in enumerate(plot_config['data']):
            df, _ = read_data(data, plot_config['BASE_DIR'])
            df = df.set_index('ref_id')
            df = df[ plot_config['columns'] ]
            if plot_config['delete_zero']:
                    df = df.replace(0.0, np.nan)
            plot_config['data'][i]['dataframe'] = df
            alpha = min(alpha, df.replace(0.0, np.nan).min().min())
        read_filter(plot_config['filter'], plot_config['BASE_DIR'])

        # special setting
        if plot_config['scale']=='log2': 
            y_label = 'log2( mutant+α / control+α )\nα={:.3f}'.format(alpha)
        elif plot_config['scale']=='log10': 
            y_label = 'log10( mutant+α / control+α )\nα={:.3f}'.format(alpha)
        else:
            y_label = 'mutant+α / control+α\nα={:.3f}'.format(alpha)

        # split datas into "control group" and "experimental group"
        new_config = []
        control = True
        tmp_df = pd.DataFrame()
        for data in plot_config['data']:
            df = data['dataframe']
            name = data['name']
            if control:
                tmp_name = name
                tmp_df = df
                control = False
            else:
                # add smallest value as alpha
                df = (df+alpha) / (tmp_df+alpha)
                df = df.replace([np.inf, -np.inf], np.nan)
                data['dataframe'] = df.reset_index()
                data['name'] = name +'\n'+ tmp_name
                new_config.append(data)
                control = True
        plot_config['data'] = new_config

        seg = plot_config['segmentation']
        if plot_config['merge'] or len(plot_config['data'])<=1 or len(plot_config['filter'])<=1:
            config_list = [ plot_config ]
        else: 
            config_list = recombine_dict(plot_config, iter=[seg])
            for i, config in enumerate(config_list):
                if config['title']:
                    config_list[i]['title'] = config[seg]['name'] + ' ' + config['title']
                else:
                    config_list[i]['title'] = config[seg]['name']
                config_list[i][seg] = [ config[seg] ]
            
        for i, config in enumerate(config_list):
            fig, plot_df = cp.box_plot(
                data=config['data'],                # dataframe of reads
                filter=config['filter'],            # reference names which want to analyze
                columns=config['columns'],          # column names which want to analyze
                title=config['title'],              # title of figure
                scale=config['scale'],              # scale of y-axis
                test=config['test'],                # hypothesis testing
                test_format=config['test_format'],  # format of testing
                detail=config['detail'],            # show mean and median at x-ticks
                x_axis=config['segmentation'],      # field of x-axis in condition 4
                y_label=y_label,                    # label of y-axis
                style=config['style'],              # backgroud style of figure
                color=config['color'],              # color palette of figure
            )
            fig.savefig(os.path.join(config['BASE_DIR'],"{}_{}.{}".format(config['fig_path'],i,config['fig_format'])), bbox_inches='tight')
    

################
# Scatter Plot #
################
def scatter_plot(run_config=None, plot_config=None):

    # analyze
    if run_config and run_config['run']:
    
        config_list = recombine_dict(run_config, iter=['data'])

        for config in config_list:
            data_df, cmt = read_data(config['data'], config['BASE_DIR'])
            reference_df, _ = read_data(config['reference'], config['BASE_DIR'])
            config['comment'] = cmt

            df_list = cp.analyze_multiple_region(
                read_df=data_df,            # dataframe of reads
                ref_df=reference_df,        # dataframe of regions
                columns=config['columns'],  # column names which want to analyze
                set_region=True,            # output with "_region" column
                set_count=True,             # output with "_count" column
                set_density=True,           # output with "_density" column
            )
            column_name = '-'.join(config['columns'])
            metadata = pd.concat(df_list, axis=1).reset_index()
            metapath = os.path.join(config['BASE_DIR'],run_config['csv_path'],
                                    "{}_region_{}.csv".format(config['data']['name'],column_name))
            if plot_config and plot_config['data']:
                append_config(plot_config['data'], config['data']['name'], path=metapath, dataframe=metadata)
            with open(metapath, "w") as f:
                config['comment'] += "# analyze_regions={}\n".format(config['columns'])
                f.write(config['comment'])
                metadata.to_csv(f, index=False)

    # plot
    if plot_config and plot_config['plot']:

        # read datas and filters into the configuration file
        # find alpha
        alpha = np.inf
        plot_config['columns'] = [ x+'_count' for x in plot_config['columns'] ]
        for i, data in enumerate(plot_config['data']):
            df, _ = read_data(data, plot_config['BASE_DIR'])
            df = df.set_index('ref_id')
            df = df[ plot_config['columns'] ]
            plot_config['data'][i]['dataframe'] = df
            alpha = min(alpha, df.replace(0.0, np.nan).min().min())
        read_filter(plot_config['filter'], plot_config['BASE_DIR'])

        # special setting
        if plot_config['scale']=='log2': 
            x_label = 'log2( control+α ), α={:.3f}'.format(alpha)
            y_label = 'log2( mutant+α )\nα={:.3f}'.format(alpha)
        elif plot_config['scale']=='log10': 
            x_label = 'log10( control+α ), α={:.3f}'.format(alpha)
            y_label = 'log10( mutant+α )\nα={:.3f}'.format(alpha)
        else:
            x_label = 'control+α , α={:.3f}'.format(alpha)
            y_label = 'mutant+α\nα={:.3f}'.format(alpha)

        # split datas into "control group" and "experimental group"
        new_config = []
        control = True
        tmp_df = pd.DataFrame()
        for data in plot_config['data']:
            df = data['dataframe']
            name = data['name']
            if control:
                tmp_df = df
                tmp_name = name
                tmp_col = df.columns.to_list()
                control = False
            else:
                # add smallest value as alpha
                df.columns = [col+'_y' for col in df.columns.to_list()]
                tmp_col += df.columns.to_list()
                df = pd.concat([tmp_df, df], axis=1)
                for col in plot_config['columns']:
                    original_col = col.split('_')[0]
                    df[original_col+'_filter'] = 'valid'
                    df.loc[(df[col]==0) | (df[col].isnull()), original_col+'_filter'] = 'invalid'
                    df.loc[(df[col+'_y']==0) | (df[col+'_y'].isnull()), original_col+'_filter'] = 'invalid'
                    
                    # add p-value
                    if plot_config['BCV'] and plot_config['BCV']!='None':
                        compare = [col,col+'_y']
                        pvalue_df = df[ df[original_col+'_filter']=='valid' ]
                        pvalue_df = pvalue_df[compare]
                        pvalue_df = add_pvalue(pvalue_df, compare, float(plot_config['BCV']))[['PValue']]
                        df = pd.merge(df, pvalue_df, on='ref_id', how='left')
                        df.loc[df[original_col+'_filter']=='valid', original_col+'_filter'] = 'p < 0.05'
                        df.loc[df['PValue']>=0.05, original_col+'_filter'] = 'p ≥ 0.05'
                        del df['PValue']
                
                for col in tmp_col:
                    df[col] = df[col] + alpha
                data['dataframe'] = df.reset_index()
                data['name'] = name +'\n'+ tmp_name
                new_config.append(data)
                control = True
        plot_config['data'] = new_config

        seg = plot_config['segmentation']
        if plot_config['merge'] or len(plot_config['data'])<=1 or len(plot_config['filter'])<=1:
            config_list = [ plot_config ]
        else: 
            config_list = recombine_dict(plot_config, iter=[seg])
            for i, config in enumerate(config_list):
                if config['title']:
                    config_list[i]['title'] = config[seg]['name'] + ' ' + config['title']
                else:
                    config_list[i]['title'] = config[seg]['name']
                config_list[i][seg] = [ config[seg] ]
            
        for i, config in enumerate(config_list):
            fig, plot_df = cp.scatter_plot(
                data=config['data'],                # dataframe of reads
                filter=config['filter'],            # reference names which want to analyze
                columns=config['columns'],          # column names which want to analyze
                hue=config['segmentation'],         # hue of figure in condition 4
                title=config['title'],              # title of figure
                scale=config['scale'],              # scale of y-axis
                x_label=x_label,                    # label of x-axis
                y_label=y_label,                    # label of y-axis
                show_others=config['show_others'],  # show other values which not in filters
                style=config['style'],              # backgroud style of figure
                color=config['color'],              # color palette of figure
                BCV=config['BCV'],                  # biological coefficient of variation
            )
            fig.savefig(os.path.join(config['BASE_DIR'],"{}_{}.{}".format(config['fig_path'],i,config['fig_format'])), bbox_inches='tight')
        

## Other Functions

##################################
# Recombine Dictionary to a List #
##################################
def recombine_dict(dict_, iter=[]):

    if iter==[]:
        return dict_
    
    # Split Static and Iterable Fields
    static_dict = {}
    iterable_dict = {}
    for key, value in dict_.items():
        if key in iter:
            iterable_dict[key] = value
        else:
            static_dict[key] = value

    # Iter Dictionary
    # eg.
    #   input: {'a':[1, 2], 'b':[5, 6]}
    #   output: [{'a':1, 'b':5},
    #            {'a':1, 'b':6},
    #            {'a':2, 'b':5},
    #            {'a':2, 'b':6}]
    values_list = []
    for value in iterable_dict.values():
        if isinstance(value, list):
            values_list.append(value)
        else:
            values_list.append([value])
    combinations = itertools.product(*values_list)
    new_dict = [dict(zip(iterable_dict.keys(), com)) for com in combinations]
    
    # Add Static Fields Back
    for dic in new_dict:
        dic.update(static_dict)

    return new_dict


#############################
# Auto-detect and Read Data #
#############################
def read_data(data, BASE_DIR=None):

    cmt,path = "",""
    if isinstance(data, str):
        path = os.path.join(BASE_DIR,data)
    elif isinstance(data, dict):
        if 'dataframe' in data:
            df = data['dataframe']
            return df, cmt
        else:
            try:
                path = os.path.join(BASE_DIR,data['path'])
            except:
                print("[Error] Unknown path of file: {}".format(data['name']))
                sys.exit(1)
    else: 
        print("[Error] Unknown type of {}, which can only be 'str' or 'dict'.".format(data))
        sys.exit(1)
    
    _, file_extension = os.path.splitext(path)
    file_extension = file_extension[1:]

    if file_extension=='csv':
        with open(path, 'r') as f:
            for line in f:
                if line[0]=='#': cmt += line   
                else: break
            header = line[:-1].split(',')
            df = pd.read_csv(f, comment='#', names=header)
    else:
        print("[Error] Unknown file extension {}, which can only be 'csv'.".format(file_extension))
        sys.exit(1)
    
    return df, cmt


########################
# Read ID into Filters #
########################
def read_filter(filter, BASE_DIR):
    if filter:
        for ft in filter:
            if 'id' not in ft:
                if 'path' in ft:
                    with open(os.path.join(BASE_DIR,ft['path']), 'r') as f:
                        ft['id'] = f.read().splitlines()


####################################
# Append Run-Config to Plot-Config #
####################################
def append_config(config, name, path=None, dataframe=None):
    for data in config:
        if data['name']==name:
            if path is not None:
                data['path'] = path
            if dataframe is not None:
                data['dataframe'] = dataframe
            break


####################################
# Calculate DEG P-value with edgeR #
####################################
def add_pvalue(data, compare, bcv):

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    import rpy2.robjects as robjects

    edgeR = importr('edgeR')

    # From Pandas to R
    with localconverter(robjects.default_converter + pandas2ri.converter):
        matrix = robjects.conversion.py2rpy(data)
    
    # DEG exact-test with BCV
    group = robjects.StrVector(compare)
    dge_list = edgeR.DGEList(matrix, group=group)
    bcv = bcv ** 2
    et = edgeR.exactTest(dge_list, dispersion=bcv)
    
    # From R to Pandas
    with localconverter(robjects.default_converter + pandas2ri.converter):
        df = robjects.conversion.rpy2py(et.rx2('table'))
    df.index.name = 'ref_id'
    
    return df