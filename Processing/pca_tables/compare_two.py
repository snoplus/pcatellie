import argparse
import os
import sys
import json
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
    return fits, def_file, tables_loc, minard_loc, runtime_loc

def load_default():
    def_fibres = []
    def_u = []
    def_v = []
    def_w = []
    with open(fits+def_file) as f:
        def_data = json.load(f)

    def_fibres = def_data["fibre"]
    def_u = def_data["u"]
    def_v = def_data["v"]
    def_w = def_data["w"]
    return def_fibres, def_u, def_v, def_w

def load_tables(file1, file2):
    with open(tables_loc+file1) as data_file1:
        data1 = json.load(data_file1)
    with open(tables_loc+file2) as data_file2:
        data2 = json.load(data_file2)
    return data1, data2

def parse_table(table):
    index = table['fibre_index']
    ipw = table['pulse_width']
    pca_offset = table['pca_offset']
    width = table['pca_offset_width']
    angb = table['ang_b']
    angerr = table['ang_b_err']
    u = table['u']
    v = table['v']
    w = table['w']
    int_map = map(int, index)
    int_list = list(int_map)
    index_int = int_list
    return index_int, ipw, pca_offset, width, angb, angerr, u, v, w

def set_up_legend():
    l1 = mpatches.Patch(color='blue',  label=name1)
    l2 = mpatches.Patch(color='red',   label=name2)
    return l1, l2

def make_ipw(l1, l2):
    plt.figure(0, figsize=(16,8))
    plt.title('IPW')
    plt.ylabel('IPW', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fibre index', horizontalalignment='right', x=1.0, weight='bold')
    plt.plot(index, ipw, 'bx', alpha=0.4)
    plt.plot(index2, ipw2, 'rx', alpha=0.4)

    plt.legend(handles=[l1, l2])
    plt.grid(True)
    plt.xticks(np.arange(0, 97, step=5))
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_ipw.png')
    return

def make_ipw_hist(l1, l2):
    pca_diff = []
    for i in range (0, len(index)):
    	for j in range (0, len(index2)):
    		if index[i] == index2[j]:
    			pca_diff.append( ipw[i] - ipw2[j] )
    plt.figure(1, figsize=(12,10))
    plt.title('IPW difference: ' + name1 + " vs " + name2)
    plt.ylabel('count', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('IPW difference', horizontalalignment='right', x=1.0, weight='bold')
    plt.hist(pca_diff, bins = 14)
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_ipw_H.png')
    return

def make_pca(l1, l2):
    plt.figure(2, figsize=(16,8))
    plt.title('PCA offset')
    plt.ylabel('PCA_offset [ns]', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fibre index', horizontalalignment='right', x=1.0, weight='bold')
    plt.errorbar(index, pca_offset, yerr=width, fmt='o', alpha=0.4, mfc='blue')
    plt.errorbar(index2, pca_offset2, yerr=width2, fmt='o', alpha=0.4, mfc='red')

    plt.legend(handles=[l1, l2])
    plt.grid(True)
    plt.xticks(np.arange(0, 97, step=5))
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_pca_offset.png')
    return

def make_pca_hist(l1, l2):
    pca_diff = []
    for i in range (0, len(index)):
    	for j in range (0, len(index2)):
    		if index[i] == index2[j]:
    			pca_diff.append( pca_offset[i] - pca_offset2[j] )

    plt.figure(3, figsize=(12,10))
    plt.title('PCA offset difference: ' + name1 + " vs " + name2)
    plt.ylabel('count', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('pca offset difference [ns]', horizontalalignment='right', x=1.0, weight='bold')
    plt.hist(pca_diff, bins = 14)
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_pca_offset_H.png')
    return

def make_angb(l1, l2):
    plt.figure(4, figsize=(16,8))
    plt.title('Angular B parameter (slope)')
    plt.ylabel('ang b', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fibre index', horizontalalignment='right', x=1.0, weight='bold')
    plt.errorbar(index, angb, yerr=angerr, fmt='o', alpha=0.4, mfc='blue')
    plt.errorbar(index2, angb2, yerr=angerr2, fmt='o', alpha=0.4, mfc='red')
    plt.legend(handles=[l1, l2])
    plt.grid(True)
    plt.xticks(np.arange(0, 97, step=5))
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_angB.png')
    return

def make_angb_hist(l1, l2):
    pca_diff = []
    for i in range (0, len(index)):
    	for j in range (0, len(index2)):
    		if index[i] == index2[j]:
    			pca_diff.append( angb[i] - angb2[j] )

    plt.figure(5, figsize=(12,10))
    plt.title('Angular B parameter difference: ' + name1 + " vs " + name2)
    plt.ylabel('count', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('ang b difference', horizontalalignment='right', x=1.0, weight='bold')
    plt.hist(pca_diff, bins = 14)
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_angB_H.png')
    return

def calc_dir_diffs():
    twoMINUSone = []
    oneMINUSdb = []
    twoMINUSdb = []

    for i in range (0, len(index)):
    	for j in range (0, len(index2)):
    		if index[i] == index2[j]:
    		    dir = np.array([u[i],v[i],w[i]])
    		    dir2 = np.array([u2[j],v2[j],w2[j]])
    		    dirDB = np.array([def_u[i], def_v[i], def_w[i]])
    		    ang = angle(dir,dir2)
    		    ang2 = angle(dir,dirDB)
    		    ang3 = angle(dir2,dirDB)
    		    if ((math.degrees(ang) > 10) or (math.degrees(ang2) > 10) or (math.degrees(ang3) > 10)):
    			    print "???!"
    		    twoMINUSone.append( math.degrees(ang) )
    		    oneMINUSdb.append( math.degrees(ang2) )
    		    twoMINUSdb.append( math.degrees(ang3) )
    return twoMINUSone, oneMINUSdb, twoMINUSdb

def dir_hist1(twoMINUSone):
    plt.figure(6, figsize=(12,10))
    plt.title('Fitted fibre direction: ' + name1 + " vs " + name2)
    plt.ylabel('count', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fitted direction (u,v,w) - Difference [deg]', horizontalalignment='right', x=1.0, weight='bold')
    plt.hist(twoMINUSone, bins = 14)
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_dir1vs2.png')
    return

def dir_hist2(oneMINUSdb):
    plt.figure(7, figsize=(12,10))
    plt.title('Fitted fibre direction: ' + name1 + " vs DB")
    plt.ylabel('count', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fitted direction (u,v,w) - Difference [deg]', horizontalalignment='right', x=1.0, weight='bold')
    plt.hist(oneMINUSdb, bins = 14)
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_dir1vsDB.png')
    return

def dir_hist3(twoMINUSdb):
    plt.figure(8, figsize=(12,10))
    plt.title('Fitted fibre direction: ' + name2 + " vs DB")
    plt.ylabel('count', horizontalalignment='right', y=1.0, weight='bold')
    plt.xlabel('fitted direction (u,v,w) - Difference [deg]', horizontalalignment='right', x=1.0, weight='bold')
    plt.hist(twoMINUSdb, bins = 14)
    #plt.show()
    if save_plots == 1:
        plt.savefig('tab_dir2vsDB.png')
    return

def make_plots():
    l1, l2 = set_up_legend()
    make_ipw(l1, l2)
    make_ipw_hist(l1, l2)
    make_pca(l1, l2)
    make_pca_hist(l1, l2)
    make_angb(l1, l2)
    make_angb_hist(l1, l2)
    twoMINUSone, oneMINUSdb, twoMINUSdb = calc_dir_diffs()
    dir_hist1(twoMINUSone)
    dir_hist2(oneMINUSdb)
    dir_hist3(twoMINUSdb)
    return

def move_plots():
    print "Moving plots:"
    create_minard_dir()
    cmd = 'mv ' + runtime_loc + 'tab_*.png ' + minard_loc + name1
    print cmd
    os.system( cmd ) #using os.system here due to wildcard*
    return

def create_minard_dir():
    print "Creating minard dir:"
    cmd = 'mkdir ' +  minard_loc + name1
    print cmd
    os.system( cmd ) #using os.system here due to wildcard*
    return

if __name__=="__main__":

    # need two ratdb tables to compare
    if len(sys.argv) != 3:
        print "provide two tables as arguments"
        exit()
    else:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        name1 = sys.argv[1][0:6]
        name2 = sys.argv[2][0:6]
        print "Comparing: " + name1 + " and " + name2
        save_plots = 1

    # parser set-up
    parser = argparse.ArgumentParser(description='compare two TELLIE PCA ratdb tables, passed as arguments')

    # load env
    fits, def_file, tables_loc, minard_loc, runtime_loc = load_env()

    # load default dirs
    def_fibres, def_u, def_v, def_w = load_default()

    # load the tables to compare
    data1, data2 = load_tables(file1, file2)

    # parse table files
    index, ipw, pca_offset, width, angb, angerr, u, v, w = parse_table(data1)
    index2, ipw2, pca_offset2, width2, angb2, angerr2, u2, v2, w2 = parse_table(data2)

    make_plots()
    #plt.show()

    move_plots()
