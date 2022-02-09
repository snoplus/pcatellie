import os
import sys

def load_env():
    pca_cons_loc = os.getenv("PCA_CONS")
    runtime_loc = os.getenv("RUNTIME_LOC")
    return pca_cons_loc, runtime_loc

def get_args():
    if len(sys.argv) != 3:
        print "Please provide PCA_TW and PCA_GF files!"
        exit()
    else:
        pca_tw = sys.argv[1]
        pca_gf = sys.argv[2]
        return pca_tw, pca_gf

def create_content():
    macro = '/rat/physics_list/OmitAll true\n'
    macro += '/rat/db/set DETECTOR geo_file "geo/snoplusnative.geo"\n\n'
    macro += '/rat/db/load ' + pca_cons_loc + pca_tw+'\n'
    macro += '/rat/db/load ' + pca_cons_loc + pca_gf+'\n\n'
    macro += '/rat/inzdab/load ' + pca_cons_loc + 'LB/lb3_0.zdab\n'
    macro += '/rat/inzdab/load ' + pca_cons_loc + 'LB/lb3_1.zdab\n'
    macro += '/rat/inzdab/load ' + pca_cons_loc + 'LB/lb3_2.zdab\n\n'
    macro += '/run/initialize\n\n'
    macro += '/rat/proc calibratePMT\n'
    macro += '/rat/procset eca 1\n' # Apply ECA and crosstalk flag
    macro += '/rat/procset pca 1\n\n' # Apply PCA
    macro += '/rat/proc count\n'
    macro += '/rat/procset update 10000\n\n'
    macro += '/rat/proc user\n'
    macro += '/rat/procset file "' + run_n + '.root"\n\n'
    macro += '/rat/inzdab/read\n\n'
    macro += 'exit'
    return macro

def save_file(content):
    file = open(runtime_loc + run_n + '_apply.mac', 'w')
    file.write(content)
    file.close()
    return

def get_run_n():
    run_n = pca_tw.split("_")[1]
    return run_n

if __name__=="__main__":
    pca_cons_loc, runtime_loc = load_env()
    pca_tw, pca_gf = get_args()
    run_n = get_run_n()
    content = create_content()
    save_file(content)
    print "DONE creating bench apply mac."
