import argparse
import os
import sys
import json
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from os import walk

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def load_env():
    fits = os.getenv("FITS_FOLDER") + "/"
    def_file = os.getenv("DEFAULT_DIR")
    tables_loc = os.getenv("TABLES_LOC")
    minard_loc = os.getenv("PLOTS_TABLES")
    runtime_loc = os.getenv("RUNTIME_LOC")
    mpl.rcParams['font.size'] = 18
    save_plots = 1
    return fits, def_file, tables_loc, minard_loc, runtime_loc, save_plots

def list_of_tables():
    f = []
    table_files = []
    for (dirpath, dirnames, filenames) in walk(tables_loc):
        f.extend(filenames)
        break
    for table in f:
        if ".ratdb" in table:
            ### this here removes old (very first) datasets
            if int(table[:6]) > 111000:
                table_files.append( table )
    return table_files

def parse_data():
    data = [None] * len(table_files)
    i = 0
    table_files.sort()
    for table in table_files:
        with open(tables_loc+table) as f:
            table_data = json.load(f)
            data[i] = table_data
            i += 1
    return data, i

def make_ipw():
    plt.figure(0, figsize=(16,8))
    plt.title('IPW')
    plt.ylabel('IPW', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fibre index', horizontalalignment='right', x=1.0, weight='bold')
    for set in range(0,sets):
        int_map = map(int, data[set]['fibre_index'])
        int_list = list(int_map)
        plt.plot(int_list, data[set]['pulse_width'], 'o', alpha=0.4, label=min(data[set]['run_number_used']))
    plt.grid(True)
    plt.xticks(np.arange(0, 97, step=5))
    plt.legend()
    if save_plots == 1:
        plt.savefig('all_ipw.png')
    return

def make_pca():
    plt.figure(1, figsize=(16,8))
    plt.title('PCA offset')
    plt.ylabel('PCA_offset [ns]', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fibre index', horizontalalignment='right', x=1.0, weight='bold')
    for set in range(0,sets):
        int_map = map(int, data[set]['fibre_index'])
        int_list = list(int_map)
        plt.errorbar(int_list, data[set]['pca_offset'], yerr=data[set]['pca_offset_width'], fmt='o', alpha=0.4, label=min(data[set]['run_number_used']))
    plt.grid(True)
    plt.xticks(np.arange(0, 97, step=5))
    plt.legend()
    if save_plots == 1:
        plt.savefig('all_pca.png')
    return

def make_angb():
    plt.figure(2, figsize=(16,8))
    plt.title('Angular B parameter (slope)')
    plt.ylabel('ang b', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fibre index', horizontalalignment='right', x=1.0, weight='bold')
    for set in range(0,sets):
        int_map = map(int, data[set]['fibre_index'])
        int_list = list(int_map)
        plt.errorbar(int_list, data[set]['ang_b'], yerr=data[set]['ang_b_err'], fmt='o', alpha=0.4, label=min(data[set]['run_number_used']))
    plt.grid(True)
    plt.xticks(np.arange(0, 97, step=5))
    plt.legend()
    if save_plots == 1:
        plt.savefig('all_angb.png')
    return

def pca_get_relative():
    pca_rels = [None] * len(table_files)
    for set in range(0,sets):
        pca_rel = []
        pmin = min(data[set]['pca_offset'])
        for i in range(0, len(data[set]['fibre_index'])):
            pca_rel.append( data[set]['pca_offset'][i]-pmin )
        pca_rels[set] = pca_rel
    return pca_rels

def make_pca_relative(pca_rels):
    plt.figure(3, figsize=(16,8))
    plt.title('PCA offset - relative')
    plt.ylabel('PCA_offset [ns]', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fibre index', horizontalalignment='right', x=1.0, weight='bold')
    for set in range(0,sets):
        int_map = map(int, data[set]['fibre_index'])
        int_list = list(int_map)
        plt.plot(int_list, pca_rels[set], 'o', alpha=0.4, label=min(data[set]['run_number_used']))
    plt.grid(True)
    plt.xticks(np.arange(0, 97, step=5))
    plt.legend()
    if save_plots == 1:
        plt.savefig('all_pca_rel.png')
    return

def make_plots():
    make_ipw()
    make_pca()
    make_angb()
    pca_rels = pca_get_relative()
    make_pca_relative(pca_rels)
    return

def move_plots():
    print "Moving plots:"
    cmd = 'mv ' + runtime_loc + 'all*.png ' + minard_loc + "all"
    print cmd
    os.system( cmd ) #using os.system here due to wildcard*
    return

if __name__=="__main__":

    # parser set-up
    parser = argparse.ArgumentParser(description='compare ALL TELLIE PCA ratdb tables')

    # load env
    fits, def_file, tables_loc, minard_loc, runtime_loc, save_plots = load_env()

    # get list of all tables
    table_files = list_of_tables()

    # parse data from tables
    data, sets = parse_data()

    # make plots
    make_plots()

    # move plots
    move_plots()
