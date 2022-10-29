import subprocess
import shlex
import time
import argparse
import os
import sys
import re
import couchdb
import json
import datetime as dt
from os import listdir
from os.path import isfile, join
from datetime import datetime
from rat.ratdb import RATDBConnector

def load_env():
    scripts_loc = os.getenv("SCRIPTS_LOC")
    upl_env = os.getenv("UPL_ENV")
    upl_val1 = os.getenv("UPL_VAL1")
    upl_val2 = os.getenv("UPL_VAL2")
    upl_fits = os.getenv("UPL_FITS")
    upl_ratdb = os.getenv("UPL_RATDB")
    make_table = os.getenv("MAKE_TABLE")
    fits_folder = os.getenv("FITS_FOLDER")
    dir_fit_file = os.getenv("DIR_FIT_FILE")
    ang_fit_file = os.getenv("ANG_FIT_FILE")
    offset_fit_file = os.getenv("OFFSET_FIT_FILE")
    default_dir = os.getenv("DEFAULT_DIR")
    tables_loc = os.getenv("TABLES_LOC")
    tables_scripts = os.getenv("TABLES_SCRIPTS_LOC")
    val1_log = os.getenv("VAL1_LOG")
    val2_log = os.getenv("VAL2_LOG")
    dir_log = os.getenv("DIR_LOG")
    as_log = os.getenv("AS_LOG")
    pca_log = os.getenv("PCA_LOG")
    bench_log = os.getenv("BENCH_LOG")
    data_loc = os.getenv("DATA_LOC")
    logs_loc = os.getenv("LOGS")
    plots = os.getenv("PLOTS")
    runtime_loc = os.getenv("RUNTIME_LOC")
    compare_table = os.getenv("COMPARE_TABLE")
    compare_all = os.getenv("COMPARE_ALL")
    cdb_add = os.getenv("COUCH_ADD")
    cdb_db = os.getenv("COUCH_DB")
    cdb_user = os.getenv("COUCHDB_USER")
    cdb_pw = os.getenv("COUCHDB_PW")
    create_bench_apply = os.getenv("CREATE_BENCH_APPLY")
    bench_cd = os.getenv("BENCH_CD")
    bench_tw = os.getenv("BENCH_TW")
    bench_peak = os.getenv("BENCH_PEAK")
    pca_cons = os.getenv("PCA_CONS")
    bench_root = os.getenv("BENCH_ROOT")
    upl_bench = os.getenv("UPL_BENCH")
    create_pca_proc = os.getenv("CREATE_PCA_PROC")
    checkpca_loc = os.getenv("CHECKPCA_LOC")
    return scripts_loc, upl_env, upl_val1, upl_val2, upl_fits, upl_ratdb, make_table, fits_folder, dir_fit_file, ang_fit_file, offset_fit_file, default_dir, tables_loc, tables_scripts, val1_log, val2_log, dir_log, as_log, pca_log, bench_log, data_loc, logs_loc, plots, runtime_loc, compare_table, compare_all, cdb_add, cdb_db, cdb_user, cdb_pw, create_bench_apply, create_pca_proc, bench_cd, bench_tw, bench_peak, pca_cons, bench_root, upl_bench, checkpca_loc

def parse_arguments():
    ### Parse arguments
    parser = argparse.ArgumentParser(description='The master script for TELLIE automation - processing.')
    parser.add_argument("-f", dest="run_list_file", help="File containing the run list (full path)", type=str, required=True)
    parser.add_argument("-c", dest="cores", help="How many cores to use (ie consecutive jobs) ?", default=4, type=int)
    parser.add_argument("-e", dest="upload_env", help="Should reupload env ?", default=0, type=int)
    parser.add_argument("-val1", dest="val1", help="Should run validate (1) ?", default=1, type=int)
    parser.add_argument("-val2", dest="val2", help="Should run validate (2) ?", default=1, type=int)
    parser.add_argument("-fit1", dest="fit1", help="Should run position fit ?", default=1, type=int)
    parser.add_argument("-fit2", dest="fit2", help="Should run ang sys fit ?", default=1, type=int)
    parser.add_argument("-fit3", dest="fit3", help="Should run pca offset fit ?", default=1, type=int)
    parser.add_argument("-pca_tab", dest="pca_tab", help="Should create pca table ?", default=1, type=int)
    parser.add_argument("-pca_tab_comp", dest="pca_tab_comp", help="Should compare pca table ?", default=1, type=int)
    parser.add_argument("-pca_tab_upl", dest="pca_tab_upl", help="Should upload pca table ?", default=0, type=int)
    args = parser.parse_args()
    return args

def call_command(cmd, use_shell=False):
    args = shlex.split(cmd)
    process = subprocess.Popen( args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=use_shell )
    return process

def check_status(process):
    return process.poll()

def check_jobs(jobs, cmds_copy=[]):
    running = 0
    success = 0
    fail = 0
    for i in range(0, len(jobs)):
        if check_status(jobs[i]) == None:
            running += 1
        elif check_status(jobs[i]) == 0:
            success += 1
        else:
            fail += 1
            if len(cmds_copy) > 0:
                failed_jobs.append( cmds_copy[i] )
    return running, success, fail

def jobs_running(cmds, retry=False):
    cmds_copy = list(cmds)
    while len(cmds) > 0:                    # keep submitting jobs if there are any in q
        print "Jobs to submit: ", len(cmds)
        if check_jobs(jobs)[0] < cores:     # if we have empty slot(s), submit
            print "Submitting: ", cmds[0]
            process = call_command( cmds[0] )
            jobs.append( process )
            time.sleep(15)
            del cmds[0]
        # wait for slots                    
        print check_jobs(jobs)                
        time.sleep(10)
    print "All submitted, now waiting..."
    while check_jobs(jobs)[0] != 0:
        print check_jobs(jobs)
        time.sleep(10)
        insert_line()
    else: 
        print "DONE: ALL JOBS"
        print check_jobs(jobs)
        insert_line()
        update_all_jobs1(jobs)
        if retry == True:
            print check_jobs(jobs, cmds_copy)
            retry_failed_jobs(failed_jobs, jobs)
            update_all_jobs2(jobs)
        return

def retry_failed_jobs(failed_jobs, jobs):
    print "Failed jobs:", len(failed_jobs)
    print "Will attempt to retry..."
    jobs *= 0
    jobs_running(failed_jobs)
    return

def reset_job_counter(jobs):
    jobs *= 0
    jobs = []
    return

def update_all_jobs1(jobs):
    all_jobs[0] += check_jobs(jobs)[0]
    all_jobs[1] += check_jobs(jobs)[1]
    print "ALL jobs current:", all_jobs
    return

def update_all_jobs2(jobs):
    all_jobs[0] += check_jobs(jobs)[0]
    all_jobs[1] += check_jobs(jobs)[1]
    all_jobs[2] += check_jobs(jobs)[2]
    print "ALL jobs current:", all_jobs
    return

def reupload_env():
    # submit reupload
    cmd = scripts_loc + upl_env
    cmd_full = "python " + cmd
    process_reupload = call_command( cmd_full )

    counter = 0
    while check_status(process_reupload) != 0:
        time.sleep(0.5)
        counter += 1
        if counter > 10:
            print "Error uploading environment"
            sys.exit()
    else:
        print "DONE Uploading environment!"
        insert_line()
        return

def parse_run_list(runfile):
    with open(runfile) as f:
        data = f.readlines()
    runlist = [int(x.strip()) for x in data]
    return runlist

def check_data_exists(runlist, data_loc):
    supplied_runs_c = 0
    good_runs_c = 0
    good_runs = []
    for run in runlist:
        supplied_runs_c += 1
        filename = "SNOP_0000" + str(run) + "_000.zdab"
        if os.path.exists(data_loc + filename) == True:
            good_runs.append( filename )
            good_runs_c += 1
    return good_runs, supplied_runs_c, good_runs_c

def set_job_limit(cores_arg):
    cores = cores_arg
    return cores

def call_validate1(good_runs, plots):
    print "Calling validate1 script..."
    cmds = []
    for run in good_runs:
        print "Run: ", run
        cmd = val1_log + "myrat " + val1_log + "template.mac -i " + data_loc + run
        cmds.append( cmd )
        insert_line()

    print "Jobs in queue: ", len(cmds)
    jobs_running( cmds, True )
    reset_job_counter(jobs)
    return

def call_position_fit(good_runs):
    print "Calling position fit script..."
    cmds = []
    for run in good_runs:
        print "Run: ", run
        cmd = dir_log + "myrat " + dir_log + "template.mac -i " + data_loc + run
        cmds.append( cmd )
        insert_line()

    print "Jobs in queue: ", len(cmds)
    jobs_running( cmds, True )
    reset_job_counter(jobs)
    return

def call_angsys_fit(good_runs):
    print "Calling ang sys fit script..."
    cmds = []
    for run in good_runs:
        print "Run: ", run
        cmd = as_log + "myrat " + as_log + "template.mac -i " + data_loc + run
        cmds.append( cmd )
        insert_line()

    print "Jobs in queue: ", len(cmds)
    jobs_running( cmds, True )
    reset_job_counter(jobs)
    return

def call_offset_fit(good_runs):
    print "Calling pca offset script..."
    cmds = []
    for run in good_runs:
        print "Run: ", run
        cmd = pca_log + "myrat " + pca_log + "template.mac -i " + data_loc + run
        cmds.append( cmd )
        insert_line()

    print "Jobs in queue: ", len(cmds)
    jobs_running( cmds, True )
    reset_job_counter(jobs)
    return

def call_validate2(good_runs):
    print "Calling validate2 script..."
    cmds = []
    for run in good_runs:
        print "Run: ", run
        cmd = val2_log + "myrat " + val2_log + "template.mac -i " + data_loc + run
        cmds.append( cmd )
        insert_line()

    print "Jobs in queue: ", len(cmds)
    jobs_running( cmds, True )
    reset_job_counter(jobs)
    return

def get_final_job_count(all_jobs):
    print "Final count: (" + str(all_jobs[0]) + ", " + str(all_jobs[1]) + ", " + str(all_jobs[2]) + ")"
    return

def create_run_folder(plots, good_runs):
    print "Creating folder for plots:"
    for run in good_runs:
        runn = str(extrac_run_number(run))
        cmd = "mkdir " + plots + str(runn)
        print cmd
        call_command( cmd )
    insert_line()
    return

def move_plots(plots, runtime_loc, good_runs):
    print "Moving plots:"
    for run in good_runs:
        runn = str(extrac_run_number(run))
        cmd = "mv " + runtime_loc + runn + "*.png " + plots + runn
        print cmd
        os.system( cmd ) #using os.system here due to wildcard*
    insert_line()
    return

def move_fits(fits_folder, runtime_loc):
    print "Moving fits:"
    sorted_runlist = sorted(runlist)
    run_name = str(sorted_runlist[0])
    cmd = 'mv ' + runtime_loc + dir_fit_file + " " + fits_folder + run_name + "_" + dir_fit_file
    cmd2 = 'mv ' + runtime_loc + ang_fit_file + " " + fits_folder + run_name + "_" + ang_fit_file
    cmd3 = 'mv ' + runtime_loc + offset_fit_file + " " + fits_folder + run_name + "_" + offset_fit_file
    print cmd
    os.system( cmd ) #using os.system here due to wildcard*
    print cmd2
    os.system( cmd2 )
    print cmd3
    os.system( cmd3 )
    insert_line()
    return

def extrac_run_number(run):
    runn = int(run.split("_")[1])
    return runn

def upload_val1(upl_val1, scripts_loc, good_runs):
    print "Uploading val1:"
    for run in good_runs:
        runn = str(extrac_run_number(run))
        cmd = "python " + scripts_loc + upl_val1 + " " + runn + "_val1.log"
        print cmd
        #call_command( cmd )
    insert_line()
    return

def upload_val2(upl_val2, scripts_loc, good_runs):
    print "Uploading val2:"
    for run in good_runs:
        runn = str(extrac_run_number(run))
        cmd = "python " + scripts_loc + upl_val2 + " " + runn + "_val2.log"
        print cmd
        #call_command( cmd )
    insert_line()
    return

def upload_fits(upl_fits, scripts_loc, good_runs):
    print "Uploading fits:"
    for run in good_runs:
        runn = str(extrac_run_number(run))
        cmd = "python " + scripts_loc + upl_fits + " " + runn + "_pos.log " + runn + "_as.log " + runn + "_pca.log "
        print cmd
        #call_command( cmd )
    insert_line()
    return

def make_pca_table(make_table, scripts_loc, tables_loc, good_runs):
    print "Creating RATDB table:"
    cmd = "python " + scripts_loc + make_table
    print cmd
    job = call_command( cmd )
    wait_for_job(job, 60)

    runn = str(extrac_run_number(good_runs[0]))
    cmd2 = "mv new_table.ratdb " + tables_loc + runn + ".ratdb"
    print cmd2
    call_command( cmd2 )
    insert_line()
    new_table = runn + ".ratdb"
    return new_table

def compare_tables_two(tables_loc, tables_scripts, compare_table, good_runs):
    print "Comparing two tables:"
    runn = str(extrac_run_number(good_runs[0]))
    new_table = runn + ".ratdb"
    print "new table: ", new_table
    old_table = get_previous_table(tables_loc)
    print "old table: ", old_table
    cmd = "python " + tables_scripts + compare_table + " " + new_table + " " + old_table
    print cmd
    job = call_command( cmd )
    wait_for_job(job, 36)
    insert_line()
    return

def compare_tables_all(tables_scripts, compare_all):
    print "Comparing all tables:"
    cmd = "python " + tables_scripts + compare_all
    print cmd
    job = call_command( cmd )
    wait_for_job(job, 36)
    insert_line()
    return

def get_previous_table(tables_loc):
    onlyfiles = [f for f in listdir(tables_loc) if isfile(join(tables_loc, f))]
    tabs = []
    for tab in onlyfiles:
        tabs.append( int(re.search(r'\d+', tab).group()))
    tabs.sort()
    old_table = str(tabs[-2]) + ".ratdb"
    return old_table

def move_table_plots(tables_loc, plots):
    ### this is part of table scipts, no need to do anymore
    print "Moving table plots:"
    cmd = "mv " + tables_loc + "*.png " + plots + "tables/"
    print cmd
    os.system( cmd ) #using os.system here due to wildcard*
    insert_line()
    return

def upload_table(new_table, scripts_loc, upl_ratdb):
    print "Uploading new ratdb table:"
    cmd = "python " + scripts_loc + upl_ratdb + " " + new_table
    print cmd
    #job = call_command( cmd )
    #wait_for_job(job, 36)
    insert_line()
    return

def cleanup(runtime_loc):
    print "Cleanup:"
    new_run_n = str(new_table.split(".")[0])
    cmd = "mkdir " + logs_loc + new_run_n
    cmd2 = "mv " + new_run_n[0] + "*.log " + logs_loc + new_run_n
    print cmd
    print cmd2
    os.system( cmd ) #using os.system here due to wildcard*
    os.system( cmd2 )
    insert_line()
    cmd = "rm " + runtime_loc + "*.mac "
    print cmd
    os.system( cmd )
    cmd2 = "rm " + runtime_loc + "*.ratdb "
    print cmd2
    os.system( cmd2 )
    insert_line()
    return

def log_cleanup(runtime_loc):
    print "Log leanup:"
    cmd = "rm " + runtime_loc + "rat.*.log "
    print cmd
    os.system( cmd ) #using os.system here due to wildcard*
    insert_line()
    return

def wait_for_job(job, limit):
    counter = 0
    print "waiting for job to finish, limit: ", limit
    while check_status(job) is None:
        print counter
        time.sleep(5)
        counter += 1
        if counter > limit:
            print "Exceeded limit!"
            sys.exit()
    else:
        print "DONE, status: ", check_status(job)
        if check_status(job) == 1: sys.exit()
        insert_line()
        return

def insert_line():
    print ""
    return

def create_dataset_doc(runlist, cdb_add, cdb_db, cdb_user, cdb_pw):
    json_data = prepare_json(runlist, cdb_add, cdb_db, cdb_user, cdb_pw)
    upload_dataset_doc(json_data, cdb_add, cdb_db, cdb_user, cdb_pw)
    return

def prepare_json(runs, cdb_add, cdb_db, cdb_user, cdb_pw):
    runs.sort()
    fibres = get_fibres(runs, cdb_add, cdb_db, cdb_user, cdb_pw)
    runrange = [ runs[0], runs[-1] ]
    day_delta = get_run_date(runs[0])
    day = get_day(day_delta)
    name = get_name(day)
    json_data = {
      "type": "runlist",
      "runlist": runs,
      "fibres": fibres,
      "runrange": runrange,
      "name": name
    }
    return json_data

def get_fibres(runs, cdb_add, cdb_db, cdb_user, cdb_pw):
    couchserver = couchdb.Server("http://"+cdb_user+":BiPo214intheneck@couch.snopl.us")
    db = couchserver['telliedb']
    runRows = db.view("_design/runs/_view/run_by_number", startkey=runs[0], endkey=runs[-1])
    fibres = []
    for row in runRows:
        if row.key in runs:
            fibres.append( db[row.id]['sub_run_info'][0]['fibre'] )
    return fibres

def get_name(day):
    temp = str(day).split("-")
    year = temp[0]
    month_n = temp[1]
    datetime_object = datetime.strptime(month_n, "%m")
    month_name = datetime_object.strftime("%b")
    name = month_name + " " + year
    return name

def get_day(day_delta):
    time_delta = dt.timedelta(days=day_delta)
    day0 = datetime(2010, 1, 1, 0, 0, 0)
    run_day = day0 + time_delta
    return run_day

def get_run_date(run):
    db = RATDBConnector(server="postgres://snoplus@pgsql.snopl.us:5400/ratdb")
    try:
        result = db.fetch(obj_type="RUN", run=run)
        if len(result) == 0:
            # Nothing was found. Just exit
            print("Didn't find a run table for run {0}".format(run))
            exit()
    except Exception:
        e = exc_info()[1]
        print("==> Database Exception: {0}".format(e))
    finally:
        if db is not None:
            db.close_ratdb_connection()
    return result[0]['data']['start_day']

def upload_dataset_doc(runlist_doc, cdb_add, cdb_db, cdb_user, cdb_pw):
    couchserver = couchdb.Server("http://"+cdb_user+":"+cdb_pw+"@"+cdb_add)
    db = couchserver[cdb_db]
    doc = runlist_doc
    try:
        doc = db.save(doc)
        print "Uploaded dataset document."
        insert_line()
        return
    except:
        print "Couldn't upload the dataset document :/"
        exit()
    return

def prepare_timestamp():
    timestamp_temp = datetime.now()
    timestamp = timestamp_temp.strftime("%d-%b-%Y %H:%M:%S.%f")
    return timestamp

def create_bench_apply_mac(tw_table, gf_table) :
    print "Creating bench apply macro:"
    bench_apply_macro  = tw_table.split("_")[1] + '_apply.mac'
    cmd = "python " + scripts_loc + create_bench_apply + " " + tw_table + " " + gf_table
    print cmd
    job = call_command( cmd )
    wait_for_job(job, 2)
    insert_line()
    return bench_apply_macro

def call_bench_apply(bench_apply_macro):
    print "Calling bench apply:"
    cmd = bench_log + "myrat " + runtime_loc + bench_apply_macro
    print cmd
    #job = call_command( cmd )
    #wait_for_job(job, 10000)
    insert_line()
    return

def get_previous_tw(tables_loc):
    onlyfiles = [f for f in listdir(tables_loc) if isfile(join(tables_loc, f))]
    tabs = []
    for tab in onlyfiles:
        tabs.append( int(re.search(r'\d+', tab).group()))
    tabs.sort()
    old_table = "PCATW_" + str(tabs[-2]) + "_0.ratdb"
    return old_table

def extract_run_numbers(new_table, prew_tw):
    new_run_n = str(int(new_table.split("_")[1]) + 1)
    old_run_n = str(int(prew_tw.split("_")[1]) + 1)
    return new_run_n, old_run_n

def extract_run_numbers2(new_table, prew_tw):
    new_run_n = str(int(new_table.split(".")[0]) +1)
    old_run_n = str(int(prew_tw.split(".")[0]) +1)
    return new_run_n, old_run_n

def get_previous_peak(tables_loc):
    onlyfiles = [f for f in listdir(tables_loc) if isfile(join(tables_loc, f))]
    tabs = []
    for tab in onlyfiles:
        tabs.append( int(re.search(r'\d+', tab).group()))
    tabs.sort()
    old_root = str(tabs[-1]) + ".root"
    return old_root

def call_cd_compare(bench_cd, new_table):
    print "Running bench cd compare macro:"
    new_run_n = str(new_table.split(".")[0])
    old_run_n = str(get_previous_set(new_run_n))
    new_tw = "PCATW_" + new_run_n + "_0.ratdb"
    old_tw = "PCATW_" + old_run_n + "_0.ratdb"
    new_run = str(int(new_run_n)+1)
    old_run = str(int(old_run_n)+1)
    cmd = bench_cd + "Compare " + pca_cons + new_tw + " " + pca_cons + old_tw + " " + new_run + " " + old_run
    print cmd
    #job = call_command( cmd )
    #wait_for_job(job, 20)
    insert_line()
    log_name = new_run + "-" + old_run + "_cd_compare.log"
    return log_name

def call_tw_compare(bench_tw, new_table):
    print "Running bench tw compare macro:"
    new_run_n = str(new_table.split(".")[0])
    old_run_n = str(get_previous_set(new_run_n))
    new_tw = "PCATW_" + new_run_n + "_0.ratdb"
    old_tw = "PCATW_" + old_run_n + "_0.ratdb"
    new_run = str(int(new_run_n)+1)
    old_run = str(int(old_run_n)+1)
    cmd = bench_tw + "Compare " + pca_cons + new_tw + " " + pca_cons + old_tw + " " + new_run + " " + old_run
    print cmd
    #job = call_command( cmd )
    #wait_for_job(job, 8)
    insert_line()
    log_name = new_run + "-" + old_run + "_tw_compare.log"
    return log_name

def call_peak_compare(bench_peak, new_table):
    print "Running bench peak compare macro:"
    new_run_n = str(new_table.split(".")[0])
    old_run_n = str(get_previous_set(new_run_n))
    new_root = new_run_n + ".root"
    old_root = old_run_n + ".root"
    new_run = str(int(new_run_n)+1)
    old_run = str(int(old_run_n)+1)
    cmd = bench_peak + "Compare " + runtime_loc + new_root + " " + bench_root + old_root + " " + new_run + " " + old_run
    print cmd
    #job = call_command( cmd )
    #wait_for_job(job, 3)
    insert_line()
    log_name = new_run_n + "-" + old_run_n + "_peak_compare.log"
    return log_name

def upload_bench(cd_log, tw_log, peak_log, scripts_loc, upl_bench):
    print "Calling upload bench script:"
    cmd = "python " + scripts_loc + upl_bench + " " + cd_log + " " + tw_log + " " + peak_log
    print cmd
    job = call_command( cmd )
    wait_for_job(job, 3)
    insert_line()
    return

def move_bench_plots(plots, runlist):
    print "Creating move bench scripts:"
    sorted_runlist = sorted(runlist)
    run_name = sorted_runlist[0]
    # create dir
    cmd = "mkdir " + plots + "bench/" + str(run_name)
    print cmd
    call_command( cmd )
    insert_line()
    # move plots
    cmd = "mv *_comp_*.png " + plots + "bench/" + str(run_name)
    print cmd
    os.system( cmd )
    insert_line()
    return

def move_bench_root():
    print "Moving bench root file:"
    # move bench ROOT
    cmd = "mv *.root " + bench_root
    print cmd
    os.system( cmd )
    insert_line()

def set_new_names(new_table):
    new_set_name = new_table.split(".")[0]
    tw_table = "PCATW_" + new_set_name + "_0.ratdb" # make sure this works
    gf_table = "PCAGF_" + new_set_name + "_0.ratdb"
    pca_root = "PCA_" + new_set_name + "_0.root"
    pca_log_file = "PCA_log_" + new_set_name + "_0.ratdb"
    bench_root_file = new_set_name + ".root"
    return tw_table, gf_table, pca_root, pca_log_file, bench_root_file

def create_pca_proc_mac() :
    print "Creating pca proc macro:"
    pca_proc_macro  = new_table.split(".")[0] + '_pca_proc.mac'
    cmd = "python " + scripts_loc + create_pca_proc + " " + new_table + " " + " ".join([str(item) for item in good_runs])
    print cmd
    job = call_command( cmd )
    wait_for_job(job, 2)
    insert_line()
    return pca_proc_macro

def call_pca_proc(pca_proc_macro):
    print "Calling pca proc:"
    cmd = "rat " + runtime_loc + pca_proc_macro
    print cmd
    job = call_command( cmd )
    wait_for_job(job, 10800) # this is 15h, need to test
    insert_line()
    return

def get_global_offset():
    # get the global offset from log file
    log_files = []
    for file in os.listdir("."):
        if file.endswith(".log"):
            log_files.append( file )
    latest_file = max(log_files, key=os.path.getctime)

    with open(latest_file) as f:
        log_data = f.readlines()

    for line in log_data:
        if "Global time offset" in line:
            global_offset = float(line.split(" ")[4])
    print "Global time offset extracted:", global_offset
    return global_offset

def call_checkPCA(pca_root, tw_table, gf_table, global_offset):
    print "Calling checkPCA script:"
    run_n = pca_root.split("_")[1]
    run = str(int(run_n)+1)
    old_run_n = str(get_previous_set(run_n))
    old_run = str(int(old_run_n)+1)
    cmd = checkpca_loc + "CheckPCALaser " + pca_cons + pca_root + " " + pca_cons + tw_table + " " + pca_cons + gf_table + " 200 400 " + run + " " + str(global_offset) + " " + pca_cons + "PCATW_" + old_run_n + "_0.ratdb " + pca_cons + "PCAGF_" + old_run_n + "_0.ratdb " + old_run
    print cmd
    #job = call_command( cmd )
    #wait_for_job(job, 360)
    insert_line()
    return

def get_previous_set(run_n):
    onlyfiles = [f for f in listdir(pca_cons) if isfile(join(pca_cons, f))]
    tabs = []
    for tab in onlyfiles:
        tabs.append( int(re.search(r'\d+', tab).group()))
    sets = set(tabs)
    sets_list = list(sets)
    sets_list.sort()
    for i in range(0, len(sets_list)):
        if int(sets_list[i]) == int(run_n):
            index = i
    old_run = sets_list[index-1]
    return old_run

def call_compareTW(tw_table):
    print "Calling compareTW script:"
    run_n = tw_table.split("_")[1]
    old_run = get_previous_set(run_n)
    old_table = "PCATW_" + str(old_run) + "_0.ratdb"
    run = str(int(run_n)+1)
    old_run = str(int(old_run)+1)
    cmd = checkpca_loc + "CompareTW " + pca_cons + tw_table + " " + pca_cons + old_table + " 200 400 " + run + " " + old_run
    print cmd
    #job = call_command( cmd )
    #wait_for_job(job, 180)
    insert_line()
    return

def setup_tw_tables():
    run_n = tw_table.split("_")[1]
    tw_file = "PCATW_" + run_n + "_0.ratdb"
    gf_file = "PCAGF_" + run_n + "_0.ratdb"
    cmd = "cp " + pca_cons + tw_file + " ."
    os.system( cmd )
    cmd2 = "cp " + pca_cons + gf_file + " ."
    os.system( cmd2 )
    mod_local_tables(tw_file, gf_file)
    return

def mod_local_tables(tw_file, gf_file):
    print "Modding local tables:"
    # TW table
    filename = tw_file
    with open(filename, 'r') as f:
        data = json.load(f)
        data['run_range'] = [117566, 117568]

    os.remove(filename)
    with open(filename, 'w') as f:
        json.dump(data, f)
    # GF table
    filename = gf_file
    with open(filename, 'r') as f:
        data = json.load(f)
        data['run_range'] = [117566, 117568]

    os.remove(filename)
    with open(filename, 'w') as f:
        json.dump(data, f)
    return    

def move_pca_const():
    print "Moving PCA constants:"
    # move constants
    cmd = "mv *.root " + pca_cons
    print cmd
    os.system( cmd )
    insert_line()
    cmd = "mv *.ratdb " + pca_cons
    print cmd
    os.system( cmd )
    insert_line()
    # move log
    cmd = "mv rat*.log " + pca_cons
    print cmd
    os.system( cmd )
    insert_line()

def move_pca_plots(plots, runlist):
    print "Moving PCA plots:"
    # make dataset folder
    sorted_runlist = sorted(runlist)
    run_name = "SET_" + str(sorted_runlist[0])
    cmd = "mkdir " + plots + str(run_name)
    print cmd
    os.system( cmd )
    insert_line()
    # move PMT plots
    cmd = "mkdir PMTs"
    print cmd
    os.system( cmd )
    insert_line()
    cmd2 = "mv PMT-* PMTs"
    print cmd2
    os.system( cmd2 )
    insert_line()
    cmd3 = "mv TW_PMT-* PMTs"
    print cmd3
    os.system( cmd3 )
    insert_line()
    cmd4 = "mv PMTs " + plots + str(run_name)
    print cmd4
    os.system( cmd4 )
    insert_line()
    # move other plots
    cmd = "mv *.png " + plots + str(run_name)
    print cmd
    os.system( cmd )
    insert_line()

if __name__=="__main__":

    ### load environment
    scripts_loc, upl_env, upl_val1, upl_val2, upl_fits, upl_ratdb, make_table, fits_folder, dir_fit_file, ang_fit_file, offset_fit_file, default_dir, tables_loc, tables_scripts, val1_log, val2_log, dir_log, as_log, pca_log, bench_log, data_loc, logs_loc, plots, runtime_loc, compare_table, compare_all, cdb_add, cdb_db, cdb_user, cdb_pw, create_bench_apply, create_pca_proc, bench_cd, bench_tw, bench_peak, pca_cons, bench_root, upl_bench, checkpca_loc = load_env()

    ### parse arguments
    args = parse_arguments()

    ### check arguments
    # reupload env
    if args.upload_env == 1:
        reupload_env()

    ### for everything else, we need a runlist to process
    runlist = parse_run_list(args.run_list_file)
    #create_dataset_doc(runlist, cdb_add, cdb_db, cdb_user, cdb_pw)

    ### check runs exist
    good_runs, supplied_runs_c, good_runs_c = check_data_exists(runlist, data_loc)
    print "Supplied runs: ", supplied_runs_c
    print "Good runs: ", good_runs_c
    print ""

    ### set job limit here
    job_counter = 0
    cores = set_job_limit(args.cores)
    print "Will use " + str(cores) + " cores :)"
    print ""

    ### start processing here
    jobs = []
    failed_jobs = []
    all_jobs = [0, 0, 0]

    test_run = good_runs

    ### call processing scripts here
    # validation 1
    #if args.val1 == 1:
        #call_validate1(test_run, plots)
        #upload_val1(upl_val1, scripts_loc, test_run)

    ### call fits scripts
    #if args.fit1 == 1:
        #call_position_fit(test_run)

    #if args.fit2 == 1:
        #call_angsys_fit(test_run)

    #if args.fit3 == 1:
        #call_offset_fit(test_run)
        #upload_fits(upl_fits, scripts_loc, good_runs)

    # validation 2
    #if args.val2 == 1:
        #call_validate2(test_run)
        #upload_val2(upl_val2, scripts_loc, test_run)

    # rat log cleanup
    #log_cleanup(runtime_loc)

    # create run folder
    #create_run_folder(plots, test_run)

    # move plots
    #move_plots(plots, runtime_loc, test_run)

    # make pca table
    #if args.pca_tab == 1:
        #new_table = make_pca_table(make_table, scripts_loc, tables_loc, test_run)

    #move_fits(fits_folder, runtime_loc)

    # compare table & move plots
    #if args.pca_tab_comp == 1:
        #compare_tables_two(tables_loc, tables_scripts, compare_table, good_runs)
        #compare_tables_all(tables_scripts, compare_all)

    # upload table
    #if args.pca_tab_upl == 1:
        #upload_table(new_table, scripts_loc, upl_ratdb)

    new_table = "114670.ratdb"
    ### PCA Processor
    # create macro
    tw_table, gf_table, pca_root, pca_log_file, bench_root_file = set_new_names(new_table)
    #pca_proc_macro = create_pca_proc_mac()
    #call_pca_proc(pca_proc_macro)
    #global_offset = get_global_offset()
    #move_pca_const()
    #call_checkPCA(pca_root, tw_table, gf_table, global_offset)
    #call_compareTW(tw_table)
    #move_pca_plots(plots, runlist)

    ### Benchmarking scripts
    # benchmarking 1: apply
    #bench_apply_macro = create_bench_apply_mac(tw_table, gf_table)

    # benchmarking 2: process
    #call_bench_apply(bench_apply_macro)

    # benchmarking 3: compare scripts
    #setup_tw_tables()
    #cd_log = call_cd_compare(bench_cd, new_table)
    #tw_log = call_tw_compare(bench_tw, new_table)
    #peak_log = call_peak_compare(bench_peak, bench_root_file)
    #upload_bench(cd_log, tw_log, peak_log, scripts_loc, upl_bench)
    #move_bench_plots(plots, runlist)

    # move pca files
    #move_bench_root()

    # cleanup
    #cleanup(runtime_loc)
    #log_cleanup(runtime_loc)

    # check final job count
    #get_final_job_count(all_jobs)
    
