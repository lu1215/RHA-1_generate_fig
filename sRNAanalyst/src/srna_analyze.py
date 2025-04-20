##################################
# Copyright (C) 2023 Ryan Chung  #
#                                #
# History:                       #
# 2023/04/10 Ryan Chung          #
##################################

import argparse
import multiprocessing
import os
import time
import yaml
from utility import density_plot, metagene_plot, position_plot, fold_change_plot, scatter_plot

__version__ = "version 1.0"

###############################
# Read and Merge Config Files #
###############################
def read_config(files=None):

    config = {}
    if files != None:
        for file in files:
            # read config file
            with open(file, 'r') as f:
                try:
                    new_config = yaml.safe_load(f)
                except yaml.YAMLError as exc:
                    print(exc)
            
            # check if keys are duplicated
            for key in new_config.keys():
                if key in config:
                    config[key].update(new_config[key])
                else:
                    config[key] = new_config[key]
    
    return config


##############################
# Merge two Configs on a Key #
##############################
def merge_config(target, source, on):
    if target:
        for config_t in target:
            for config_s in source:
                if config_t[on]==config_s[on]:
                    config_t.update(config_s)
                    break
    return target


####################################
# Append Run-Config to Plot-Config #
####################################
def append_path(target, source, column, tool=None):
    for i in range(len(target['data'])):
        if (tool=='boundary') or (tool=='codon'):
            annot = '_position_{}({}-{})'.format(source['columns'][column],source['limit'][column][0],source['limit'][column][1])
        else:
            annot = ''
        target['data'][i]['path'] = os.path.join(source['csv_path'], target['data'][i]['name'] + annot + '.csv')


#################################
# Set Config Files for Analysis #
#################################
def set_config(config, tool, set_stylesheet=True):
    
    if tool=='density':
        config['run_config']['Density']['data'] = merge_config(
            target = config['run_config']['Density']['data'], 
            source = config['samples']['data'],
            on = 'name')
        config['plot_config']['Density']['filter'] = merge_config(
            target = config['plot_config']['Density']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['run_config']['Density']['BASE_DIR'] = config['run_config']['BASE_DIR']
        config['plot_config']['Density']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Density']['fig_path'] += 'Density'
        if set_stylesheet:
            config['plot_config']['Density'].update(config['stylesheet']['General'])
            config['plot_config']['Density'].update(config['stylesheet']['Box_Plot'])
    
    elif tool=='metagene':
        config['run_config']['Metagene']['data'] = merge_config(
            target = config['run_config']['Metagene']['data'], 
            source = config['samples']['data'],
            on = 'name')
        config['plot_config']['Metagene']['filter'] = merge_config(
            target = config['plot_config']['Metagene']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['run_config']['Metagene']['BASE_DIR'] = config['run_config']['BASE_DIR']
        config['plot_config']['Metagene']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Metagene']['fig_path'] += 'Metagene'
        if set_stylesheet:
            config['plot_config']['Metagene'].update(config['stylesheet']['General'])
            config['plot_config']['Metagene'].update(config['stylesheet']['Line_Plot'])

    elif tool=='boundary':
        config['run_config']['Boundary']['data'] = merge_config(
            target = config['run_config']['Boundary']['data'], 
            source = config['samples']['data'],
            on = 'name')
        config['plot_config']['Head']['filter'] = merge_config(
            target = config['plot_config']['Head']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['plot_config']['Tail']['filter'] = merge_config(
            target = config['plot_config']['Tail']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['run_config']['Boundary']['BASE_DIR'] = config['run_config']['BASE_DIR']
        config['plot_config']['Head']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Tail']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Head']['fig_path'] += 'Head'
        config['plot_config']['Tail']['fig_path'] += 'Tail'
        if set_stylesheet:
            config['plot_config']['Head'].update(config['stylesheet']['General'])
            config['plot_config']['Head'].update(config['stylesheet']['Line_Plot'])
            config['plot_config']['Tail'].update(config['stylesheet']['General'])
            config['plot_config']['Tail'].update(config['stylesheet']['Line_Plot'])
    
    elif tool=='codon':
        config['run_config']['Codon']['data'] = merge_config(
            target = config['run_config']['Codon']['data'], 
            source = config['samples']['data'],
            on = 'name')
        config['plot_config']['Start_Codon']['filter'] = merge_config(
            target = config['plot_config']['Start_Codon']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['plot_config']['Stop_Codon']['filter'] = merge_config(
            target = config['plot_config']['Stop_Codon']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['run_config']['Codon']['BASE_DIR'] = config['run_config']['BASE_DIR']
        config['plot_config']['Start_Codon']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Stop_Codon']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Start_Codon']['fig_path'] += 'Start_Codon'
        config['plot_config']['Stop_Codon']['fig_path'] += 'Stop_Codon'
        if set_stylesheet:
            config['plot_config']['Start_Codon'].update(config['stylesheet']['General'])
            config['plot_config']['Start_Codon'].update(config['stylesheet']['Line_Plot'])
            config['plot_config']['Stop_Codon'].update(config['stylesheet']['General'])
            config['plot_config']['Stop_Codon'].update(config['stylesheet']['Line_Plot'])

    elif tool=='fold':
        config['run_config']['Fold_Change']['data'] = merge_config(
            target = config['run_config']['Fold_Change']['data'], 
            source = config['samples']['data'],
            on = 'name')
        config['plot_config']['Fold_Change']['filter'] = merge_config(
            target = config['plot_config']['Fold_Change']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['run_config']['Fold_Change']['BASE_DIR'] = config['run_config']['BASE_DIR']
        config['plot_config']['Fold_Change']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Fold_Change']['fig_path'] += 'Fold_Change'
        if set_stylesheet:
            config['plot_config']['Fold_Change'].update(config['stylesheet']['General'])
            config['plot_config']['Fold_Change'].update(config['stylesheet']['Box_Plot'])

    elif tool=='scatter':
        config['run_config']['Scatter']['data'] = merge_config(
            target = config['run_config']['Scatter']['data'], 
            source = config['samples']['data'],
            on = 'name')
        config['plot_config']['Scatter']['filter'] = merge_config(
            target = config['plot_config']['Scatter']['filter'], 
            source = config['samples']['filter'],
            on = 'name')
        config['run_config']['Scatter']['BASE_DIR'] = config['run_config']['BASE_DIR']
        config['plot_config']['Scatter']['BASE_DIR'] = config['plot_config']['BASE_DIR']
        config['plot_config']['Scatter']['fig_path'] += 'Scatter'
        if set_stylesheet:
            config['plot_config']['Scatter'].update(config['stylesheet']['General'])
            config['plot_config']['Scatter'].update(config['stylesheet']['Scatter_Plot'])
    
    return config


################
# Run Analysis #
################
def run(config, tool):

    cpu = multiprocessing.cpu_count()
    cpu = 2 if cpu>2 else cpu

    if tool=='density':
        density_plot(config['run_config']['Density'], config['plot_config']['Density'])
    
    elif tool=='metagene':
        metagene_plot(config['run_config']['Metagene'], config['plot_config']['Metagene'])

    elif tool=='boundary':
        position_plot(run_config=config['run_config']['Boundary'])
        if config['run_config']['Boundary']['run']:
            for i, col in enumerate(['Head','Tail']):
                append_path(config['plot_config'][col], config['run_config']['Boundary'], i, tool)
        arg = [(None,config['plot_config']['Head']),
                (None,config['plot_config']['Tail'])]
        pool = multiprocessing.Pool(cpu)
        pool.starmap(position_plot,arg)
        pool.close()
        pool.join()
    
    elif tool=='codon':
        position_plot(run_config=config['run_config']['Codon'])
        if config['run_config']['Codon']['run']:
            for i, col in enumerate(['Start_Codon','Stop_Codon']):
                append_path(config['plot_config'][col], config['run_config']['Codon'], i, tool)
        arg = [(None,config['plot_config']['Start_Codon']),
                (None,config['plot_config']['Stop_Codon'])]
        pool = multiprocessing.Pool(cpu)
        pool.starmap(position_plot,arg)
        pool.close()
        pool.join()

    elif tool=='fold':
        fold_change_plot(config['run_config']['Fold_Change'], config['plot_config']['Fold_Change'])

    elif tool=='scatter':
        scatter_plot(config['run_config']['Scatter'], config['plot_config']['Scatter'])
    

################
# Main Program #
################
if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser(
        description="This program is for analyzing NGS-Seq data.",
    )
    parser.add_argument("-v", "--version",
                        action="version",
                        version="%(prog)s " + __version__)
    parser.add_argument("--config",
                        help="path for storing configuration files (YAML)")
    parser.add_argument("--density",
                        action="store_true",
                        help="analyze read-count into mutiple regions")
    parser.add_argument("--metagene",
                        action="store_true",
                        help="analyze read-count into metagene")
    parser.add_argument("--boundary",
                        action="store_true",
                        help="analyze read-count near boundary")
    parser.add_argument("--codon",
                        action="store_true",
                        help="analyze read-count near start and stop codon")
    parser.add_argument("--fold",
                        action="store_true",
                        help="analyze read-count fold-change in two conditions")
    parser.add_argument("--scatter",
                        action="store_true",
                        help="analyze read-count scatter plot in two conditions")
    args = parser.parse_args()

    # main program
    T = time.time()

    # initialize config
    if args.config!=None:
        base_path = args.config
    else:
        base_path = os.path.abspath(__file__+'/../config')
    config = {
        'run_config': read_config([base_path+'/run_config.yml']),
        'plot_config': read_config([base_path+'/plot_config.yml']),
        'samples': read_config([base_path+'/samples.yml']),
        'stylesheet': read_config([base_path+'/stylesheet.yml']),
    }
    
    if args.density: 
        config = set_config(config, 'density')
        print("Start running Density...")
        run(config, 'density')

    if args.metagene:
        config = set_config(config, 'metagene')
        print("Start running Metagene...")
        run(config, 'metagene')
    
    if args.boundary:
        config = set_config(config, 'boundary')
        print("Start running Boundary...")
        run(config, 'boundary')
        
    if args.codon:
        config = set_config(config, 'codon')
        print("Start running Codon...")
        run(config, 'codon')
    
    if args.fold:
        config = set_config(config, 'fold')
        print("Start running Fold-change...")
        run(config, 'fold')
    
    if args.scatter:
        config = set_config(config, 'scatter')
        print("Start running Scatter...")
        run(config, 'scatter')

    print("Program complete.")
    print("Time:{:.3f}s".format(time.time()-T))
    