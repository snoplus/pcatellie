import couchdb
import os
import sys
import re
import json
from datetime import datetime

def get_args():
    if len(sys.argv) != 2:
        print "Please provide log file name..."
        exit()
    else:
        filename = sys.argv[1]
        return filename

def get_runnumber(filename):
    run = int(re.search(r'\d+', filename).group())
    return run

def load_env():
    c_add = os.getenv("COUCH_ADD")
    c_db = os.getenv("COUCH_DB")
    c_user = os.getenv("COUCHDB_USER")
    c_pw = os.getenv("COUCHDB_PW")
    logs_folder = os.getenv("RUNTIME_LOC")
    return c_add, c_db, c_user, c_pw, logs_folder

def load_log(vars, filename):
    with open(vars[4] + filename) as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

def prepare_timestamp():
    timestamp_temp = datetime.now()
    timestamp = timestamp_temp.strftime("%d-%b-%Y %H:%M:%S.%f")
    return timestamp

def parse_content(content, timestamp, run):
    for line in content:
        if "TOF mean:" in line:
            TOF_mean = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "TOF RMS:" in line:
            TOF_rms = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Bucket mean:" in line:
            bucket_mean = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Bucket RMS:" in line:
            bucket_rms = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Ang sys mean:" in line:
            as_mean = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Ang sys RMS:" in line:
            as_rms = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "TOF min:" in line:
            TOF_min = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "TOF max:" in line:
            TOF_max = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Bucket min:" in line:
            bucket_min = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Bucket max:" in line:
            bucket_max = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Ang sys min:" in line:
            as_min = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Ang sys max:" in line:
            as_max = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Direct peaks:" in line:
            d_peaks = int(re.search(r'\d+', line).group())
        if "Direct peaktime:" in line:
            d_time = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Residual peaks:" in line:
            r_peaks = int(re.search(r'\d+', line).group())
        if "Residual peaktime:" in line:
            r_time = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Dir-Resid fit:" in line:
            dir_res = line.split(":")[1].split("\t")
        if "Direct vs Angle fit:" in line:
            dir_ang = line.split(":")[1].split("\t")
        if "Residual vs Angle fit:" in line:
            res_ang = line.split(":")[1].split("\t")
        if "TOF flag:" in line:
            FTOF = int(re.search(r'\d+', line).group())
        if "Bucket flag:" in line:
            FBuc = int(re.search(r'\d+', line).group())
        if "Ang sys flag:" in line:
            FAS = int(re.search(r'\d+', line).group())
        if "TOF min-max flag:" in line:
            FTOFmm = int(re.search(r'\d+', line).group())
        if "Bucket min-max flag:" in line:
            FBucmm = int(re.search(r'\d+', line).group())
        if "Ang sys min-max flag:" in line:
            FASmm = int(re.search(r'\d+', line).group())
        if "TOF progression flag:" in line:
            FTOFp = int(re.search(r'\d+', line).group())
        if "Ang sys progression flag:" in line:
            FASp = int(re.search(r'\d+', line).group())
        if "Resid peak flag:" in line:
            Fres = int(re.search(r'\d+', line).group())
        if "Dir-Resid fit flag:" in line:
            Fdr = int(re.search(r'\d+', line).group())
        if "Direct vs Angle fit flag:" in line:
            Fda = int(re.search(r'\d+', line).group())
        if "Residual vs Angle fit flag:" in line:
            Fra = int(re.search(r'\d+', line).group())
    bitword = str(FTOF)+str(FBuc)+str(FAS)+str(FTOFmm)+str(FBucmm)+str(FASmm)+str(FTOFp)+str(FASp)+str(Fres)+str(Fdr)+str(Fda)+str(Fra)
    json_data = {
    "run": run,
    "TOF mean": TOF_mean,
    "TOF RMS": TOF_rms,
    "Bucket mean": bucket_mean,
    "Bucket RMS": bucket_rms,
    "Ang sys mean": as_mean,
    "Ang sys RMS": as_rms,
    "TOF min": TOF_min,
    "TOF max": TOF_max,
    "Bucket min": bucket_min,
    "Bucket max": bucket_max,
    "Ang sys min": as_min,
    "Ang sys max": as_max,
    "Direct peaks": d_peaks,
    "Direct peaktime": d_time,
    "Residual peaks": r_peaks,
    "Residual peaktime": r_time,
    "Dir-Resid fit": dir_res,
    "Direct vs Angle fit": dir_ang,
    "Residual vs Angle fit": res_ang,
    "TOF flag": FTOF,
    "Bucket flag": FBuc,
    "Ang sys flag": FAS,
    "TOF min-max flag": FTOFmm,
    "Bucket min-max flag": FBucmm,
    "Ang sys min-max flag": FASmm,
    "TOF progression flag": FTOFp,
    "Ang sys progression flag": FASp,
    "Resid peak flag": Fres,
    "Dir-Resid fit flag": Fdr,
    "Direct vs Angle fit flag": Fda,
    "Residual vs Angle fit flag": Fra,
    "bitword": bitword,
    "type": 'val2',
    "timestamp": timestamp
    }
    return json_data

def upload_doc(vars, json_data):
    couchserver = couchdb.Server("http://"+vars[2]+":"+vars[3]+"@"+vars[0])
    db = couchserver[vars[1]]
    try:
        # save the new doc
        doc = db.save(json_data)
        return
    except:
        print "Couldn't create the doc :/"
        exit()

if __name__=="__main__":
    filename = get_args()
    run = get_runnumber(filename)
    vars = load_env()
    timestamp = prepare_timestamp()
    content = load_log(vars, filename)
    json_data = parse_content(content, timestamp, run)
    upload_doc(vars, json_data)
    print "DONE"
