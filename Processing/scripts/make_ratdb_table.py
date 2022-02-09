import argparse
import sys
import couchdb
import json
import os
from datetime import datetime

def load_env():
    fits = os.getenv("FITS_FOLDER") + "/"
    dir_file = os.getenv("DIR_FIT_FILE")
    ang_file = os.getenv("ANG_FIT_FILE")
    pcaoffset_file = os.getenv("OFFSET_FIT_FILE")
    couch_add = os.getenv("COUCH_ADD")
    couch_user = os.getenv("COUCHDB_USER")
    couch_pw = os.getenv("COUCHDB_PW")
    couch_db = os.getenv("COUCH_DB_TELLIE")
    return dir_file, ang_file, pcaoffset_file, fits, couch_add, couch_user, couch_pw, couch_db

def load_pos_data():
    ### load pos data
    with open(fits+dir_fit) as f:
        pos_fit_data = f.readlines()
    pos_fit_data = [x.strip() for x in pos_fit_data]

    pos_fit_data_fibre = []
    pos_fit_data_x = []
    pos_fit_data_y = []
    pos_fit_data_z = []
    for line in pos_fit_data:
        split = line.split("\t")
        pos_fit_data_fibre.append( split[0] )
        pos_fit_data_x.append( float(split[4]) )
        pos_fit_data_y.append( float(split[5]) )
        pos_fit_data_z.append( float(split[6]) )

    print "Loaded " + str(len(pos_fit_data_fibre)) + " position fits"
    return pos_fit_data_fibre, pos_fit_data_x, pos_fit_data_y, pos_fit_data_z

def load_dir_data():
    ### load dir data
    with open(fits+dir_fit) as f:
        dir_fit_data = f.readlines()
    dir_fit_data = [x.strip() for x in dir_fit_data]

    dir_fit_data_fibre = []
    dir_fit_data_x = []
    dir_fit_data_y = []
    dir_fit_data_z = []
    for line in dir_fit_data:
        split = line.split("\t")
        dir_fit_data_fibre.append( split[0] )
        dir_fit_data_x.append( float(split[1]) )
        dir_fit_data_y.append( float(split[2]) )
        dir_fit_data_z.append( float(split[3]) )

    print "Loaded " + str(len(dir_fit_data_fibre)) + " direction fits"
    return dir_fit_data_fibre, dir_fit_data_x, dir_fit_data_y, dir_fit_data_z

def load_ang_data():
    ### load angular data
    with open(fits+ang_fit) as f:
        ang_fit_data = f.readlines()
    ang_fit_data = [x.strip() for x in ang_fit_data]

    ang_fit_data_fibre = []
    ang_fit_data_anga = []
    ang_fit_data_angaerr = []
    ang_fit_data_angb = []
    ang_fit_data_angberr = []
    for line in ang_fit_data:
        split = line.split()
        ang_fit_data_fibre.append( split[0] )
        ang_fit_data_anga.append( float(split[1]) )
        ang_fit_data_angaerr.append( float(split[2]) )
        ang_fit_data_angb.append( float(split[3]) )
        ang_fit_data_angberr.append( float(split[4]) )

    print "Loaded " + str(len(ang_fit_data_fibre)) + " angular fits"
    return ang_fit_data_fibre, ang_fit_data_anga, ang_fit_data_angaerr, ang_fit_data_angb, ang_fit_data_angberr

def load_pca_data():
    ### load pca data
    with open(fits+pca_fit) as f:
        pca_fit_data = f.readlines()
    pca_fit_data = [x.strip() for x in pca_fit_data]

    pca_fit_data_fibre = []
    pca_fit_data_run = []
    pca_fit_data_pca = []
    pca_fit_data_width = []
    for line in pca_fit_data:
        split = line.split("\t")
        pca_fit_data_fibre.append( split[0] )
        pca_fit_data_run.append( int(split[1]) )
        pca_fit_data_pca.append( float(split[2]) )
        pca_fit_data_width.append( float(split[3]) )

    print "Loaded " + str(len(pca_fit_data_fibre)) + " pca fits"
    return pca_fit_data_fibre, pca_fit_data_run, pca_fit_data_pca, pca_fit_data_width

def load_couchdb_data(runs):
    couch = couchdb.Server('https://' + couch_user + ':' + couch_pw + '@' + couch_add)
    db = couch[couch_db]

    fibres = []
    IPWs = []
    trig_dels = []
    fibre_dels = []

    print "Loading COUCHDB data..."
    for run in runs:
        print run
        for item in db.view('_design/runs/_view/run_by_number', key=run):
            doc = db[item.id]
            fibre = doc["sub_run_info"][0]['fibre']
            ipw = doc["sub_run_info"][0]['pulse_width']
            trig_del = doc["sub_run_info"][0]['trigger_delay']
            fibre_del = doc["sub_run_info"][0]['fibre_delay']
            fibres.append( str(fibre) )
            IPWs.append( int(ipw) )
            trig_dels.append( int(trig_del) )
            fibre_dels.append( float(fibre_del) )

    print "Loaded " + str(len(fibres)) + " couchdb entries"
    return fibres, IPWs, trig_dels, fibre_dels

def print_pos_data():
    print pos_fit_data_fibre
    print pos_fit_data_x
    print pos_fit_data_y
    print pos_fit_data_z

def print_dir_data():
    print dir_fit_data_fibre
    print dir_fit_data_x
    print dir_fit_data_y
    print dir_fit_data_z

def print_ang_data():
    print ang_fit_data_fibre
    print ang_fit_data_run
    print ang_fit_data_anga
    print ang_fit_data_angaerr
    print ang_fit_data_angb
    print ang_fit_data_angberr

def print_pca_data():
    print pca_fit_data_fibre
    print pca_fit_data_pca
    print pca_fit_data_width

def print_files():
    print dir_fit, ang_fit, pca_fit, fits


if __name__=="__main__":

    dir_fit = ""
    ang_fit = ""
    pca_fit = ""

    # parser set-up
    parser = argparse.ArgumentParser(description='create TELLIE PCA ratdb table')

    # load env
    dir_fit, ang_fit, pca_fit, fits, couch_add, couch_user, couch_pw, couch_db = load_env()

    print_files()

    # read in data from fit files
    pos_fit_data_fibre, pos_fit_data_x, pos_fit_data_y, pos_fit_data_z = load_pos_data()
    dir_fit_data_fibre, dir_fit_data_x, dir_fit_data_y, dir_fit_data_z = load_dir_data()
    ang_fit_data_fibre, ang_fit_data_anga, ang_fit_data_angaerr, ang_fit_data_angb, ang_fit_data_angberr = load_ang_data()
    pca_fit_data_fibre, pca_fit_data_run, pca_fit_data_pca, pca_fit_data_width = load_pca_data()

    #print_pos_data()
    #print_dir_data()
    #print_ang_data()
    #print_pca_data()

    #load couchdb data
    couch_fibres, IPWs, trig_dels, fibre_dels = load_couchdb_data(pca_fit_data_run)

    # Make final arrays
    final_fibre = []
    final_index = []
    final_run = []
    final_ipw = []
    final_trig = []
    final_del = []
    final_pca = []
    final_width = []
    final_anga = []
    final_angaerr = []
    final_angb = []
    final_angberr = []
    final_x = []
    final_y = []
    final_z = []
    final_u = []
    final_v = []
    final_w = []

    final_fibre = ang_fit_data_fibre
    final_run = pca_fit_data_run
    final_anga = ang_fit_data_anga
    final_angaerr = ang_fit_data_angaerr
    final_angb = ang_fit_data_angb
    final_angberr = ang_fit_data_angberr

    # Match from other vectors
    for fibre in final_fibre:
        for i in range (0, len(couch_fibres)):
            if fibre == couch_fibres[i]:
                final_ipw.append( IPWs[i] )
                final_trig.append( trig_dels[i] )
                final_del.append( fibre_dels[i] )

        for j in range (0, len(pca_fit_data_fibre)):
            if fibre == pca_fit_data_fibre[j]:
                final_pca.append( pca_fit_data_pca[j] )
                final_width.append( pca_fit_data_width[j] )

        for k in range (0, len(pos_fit_data_fibre)):
            if fibre == pos_fit_data_fibre[k]:
                final_x.append( pos_fit_data_x[k] )
                final_y.append( pos_fit_data_y[k] )
                final_z.append( pos_fit_data_z[k] )

        for l in range (0, len(dir_fit_data_fibre)):
            if fibre == dir_fit_data_fibre[l]:
                final_u.append( dir_fit_data_x[l] )
                final_v.append( dir_fit_data_y[l] )
                final_w.append( dir_fit_data_z[l] )

        final_index.append( str(fibre[2:5]) )

    # Check lengths
    print len(final_fibre), len(final_index), len(final_run), len(final_ipw), len(final_trig), len(final_del), len(final_pca), len(final_width), len(final_anga), len(final_angaerr), len(final_angb), len(final_angberr), len(final_x), len(final_y), len(final_z), len(final_u), len(final_v), len(final_w)

    # Create json data
    run_min = min(final_run)
    run_max = max(final_run)
    timestamp_date = datetime.today().strftime('%Y-%m-%d')
    timestamp_time = datetime.today().strftime('%H-%M-%S')
    timestamp = timestamp_date + "T" + timestamp_time

    json_data = {
    "fibre_index": final_index,
    "fibre": final_fibre,
    "run_number_used": final_run,
    "pulse_width": final_ipw,
    "trigger_delay": final_trig,
    "fibre_delay": final_del,
    "ang_a": final_anga,
    "ang_a_err": final_angaerr,
    "ang_b": final_angb,
    "ang_b_err": final_angberr,
    "pca_offset": final_pca,
    "pca_offset_width": final_width,
    "u": final_u,
    "v": final_v,
    "w": final_w,
    "run_range": [run_min, run_max],
    "comment": "",
    "pass": 1,
    "version": 4,
    "type": "TELLIE_PCA_FIBRE_OFFSETS",
    "index": "SLAVE",
    "timestamp": timestamp
    }

    print json_data

    # Save to ratdb file
    with open("new_table.ratdb", "w") as f:
        json.dump(json_data, f)
