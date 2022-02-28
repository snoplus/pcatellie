import os
import sys

def load_env():
    runtime_loc = os.getenv("RUNTIME_LOC")
    data_loc = os.getenv("DATA_LOC")
    return runtime_loc, data_loc

def get_args():
    if len(sys.argv) < 3:
        print "Please provide the runlist"
        exit()
    else:
        runs = []
        for run in sys.argv[1:]:
            runs.append( run )
        return runs

def get_run_n(runs):
    runs.sort()
    run_n = int(runs[0].split("_")[1])
    return run_n, runs

def create_content(runs, data_loc):
    macro = '/rat/physics_list/OmitAll true\n\n'
    macro += '/rat/db/set DETECTOR geo_file "geo/snoplusnative.geo"\n\n'
    macro += '/rat/db/set PCA_GENERATION pca_source 1\n'
    macro += '/rat/db/set PCA_GENERATION low_occ_lim 300\n'
    macro += '/rat/db/set PCA_GENERATION pca_verbosity 2\n'
    macro += '/rat/db/set PCA_GENERATION pca_source_mode 0\n'
    macro += '/rat/db/set PCA_GENERATION applyAngSys 1\n\n'

    for run in runs:
        macro += '/rat/inzdab/load ' + data_loc + run + '\n'

    macro += '\n/run/initialize\n\n'
    macro += '/rat/proc/if trigTypeSelector\n'
    macro += '\t/rat/procset trigType "N100Low"\n'
    macro += '\t/rat/procset trigType "N100Med"\n'
    macro += '\t/rat/procset trigType "N100High"\n'
    macro += '\t/rat/procset trigType "N20"\n'
    macro += '\t/rat/procset trigType "N20LB"\n'
    macro += '\t/rat/procset trigType "Pedestal"\n'
    macro += '\t/rat/procset trigType "EXT8PulseAsy"\n'
    macro += '/rat/proc/else\n'
    macro += '\t/rat/proc/if trigTypeSelector\n'
    macro += '\t\t/rat/procset trigType "EXTASY"\n\n'
    macro += '\t\t/rat/proc/if nhitCut\n'
    macro += '\t\t\t/rat/procset nhit 700\n\n'
    macro += '\t\t/rat/proc/else\n'
    macro += '\t\t\t/rat/proc calibratePMT\n'
    macro += '\t\t\t/rat/procset eca 1\n'
    macro += '\t\t\t/rat/procset pca 0\n\n'
    macro += '\t\t\t/rat/proc count\n'
    macro += '\t\t\t/rat/procset update 10000\n\n'
    macro += '\t\t\t/rat/proc genPCA\n\n'
    macro += '\t\t/rat/proc/endif\n'
    macro += '\t/rat/proc/endif\n'
    macro += '/rat/proc/endif\n\n'
    macro += '/rat/inzdab/read\n\n'
    macro += 'exit'
    return macro

def save_file(content):
    file = open(runtime_loc + str(run_n) + '_pca_proc.mac', 'w')
    file.write(content)
    file.close()
    return

if __name__=="__main__":
    runtime_loc, data_loc = load_env()
    runs = get_args()
    run_n, runs = get_run_n(runs)
    content = create_content(runs, data_loc)
    save_file(content)
    print "DONE creating pca proc mac."
