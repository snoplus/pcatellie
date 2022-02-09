import couchdb
import os
import sys
import re
import json
from datetime import datetime

def get_args():
    if len(sys.argv) != 4:
        print "Please provide log file names (3), order is important.."
        exit()
    else:
        filename_dir = sys.argv[1]
        filename_as = sys.argv[2]
        filename_pca = sys.argv[3]
        return filename_dir, filename_as, filename_pca

def load_env():
    c_add = os.getenv("COUCH_ADD")
    c_db = os.getenv("COUCH_DB")
    c_user = os.getenv("COUCHDB_USER")
    c_pw = os.getenv("COUCHDB_PW")
    logs_folder1 = os.getenv("RUNTIME_LOC") # this allows for different location of the fit files (same now)
    logs_folder2 = os.getenv("RUNTIME_LOC")
    logs_folder3 = os.getenv("RUNTIME_LOC")
    return c_add, c_db, c_user, c_pw, logs_folder1, logs_folder2, logs_folder3

def load_log(vars, filename, z):
    with open(vars[z] + filename) as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

def prepare_timestamp():
    timestamp_temp = datetime.now()
    timestamp = timestamp_temp.strftime("%d-%b-%Y %H:%M:%S.%f")
    return timestamp

def parse_content(content, timestamp, type):
    if type == "dir":
        for line in content:
            if "Run:" in line:
                run = int(re.search(r'\d+', line).group())
            if "Fibre:" in line:
                fibre = line.split(" ")[1]
            if "Is affected by belly plate?" in line:
                belly = int(re.search(r'\d+', line).group())
            if "RATDB pos:" in line:
                ratdb_pos = line.split(":")[1].split("\t")
            if "RATDB dir:" in line:
                ratdb_dir = line.split(":")[1].split("\t")
            if "Correct fibre?" in line:
                correct = int(re.search(r'\d+', line).group())
            if "Final DIR fit:" in line:
                final_dir = line.split(":")[1].split("\t")
            if "DIR fit angular deviation from expected:" in line:
                DIR_dev = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
            if "Final REF fit:" in line:
                final_ref = line.split(":")[1].split("\t")
            if "REF fit angular deviation from expected:" in line:
                REF_dev = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
            if "New fitted fibre direction:" in line:
                new_fit = line.split(":")[1].split("\t")
            if "Angular difference:" in line:
                ang_dev = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        json_data = {
        "Run": run,
        "Fibre": fibre,
        "Is affected by belly plate?": belly,
        "RATDB pos": ratdb_pos,
        "RATDB dir": ratdb_dir,
        "Correct fibre?": correct,
        "Final DIR fit": final_dir,
        "DIR fit angular deviation from expected": DIR_dev,
        "Final REF fit": final_ref,
        "REF fit angular deviation from expected": REF_dev,
        "New fitted fibre direction": new_fit,
        "Angular difference": ang_dev,
        "type": "fits",
        "timestamp": timestamp
        }
        return json_data
    elif type == "as":
        for line in content:
            if "Mean hit time and RMS:" in line:
                mean_hit = line.split(":")[1].split("\t")
            if "Ang a:" in line:
                ang_a = line.split(":")[1].split("\t")
            if "Ang b:" in line:
                ang_b = line.split(":")[1].split("\t")
        json_data2 = {
        "Mean hit time and RMS": mean_hit,
        "Ang a": ang_a,
        "Ang b": ang_b
        }
        return json_data2
    elif type == "pca":
        for line in content:
            if "Injection time:" in line:
                inj = line.split(":")[1].split("\t")
        json_data3 = {
        "Injection time": inj
        }
        return json_data3

def merge_contents(json_data1, json_data2, json_data3):
    json_data1.update(json_data2)
    json_data1.update(json_data3)
    return json_data1

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
    filename_dir, filename_as, filename_pca = get_args()
    vars = load_env()
    timestamp = prepare_timestamp()
    content1 = load_log(vars, filename_dir, 4)
    content2 = load_log(vars, filename_as, 5)
    content3 = load_log(vars, filename_pca, 6)
    json_data1 = parse_content(content1, timestamp, "dir")
    json_data2 = parse_content(content2, timestamp, "as")
    json_data3 = parse_content(content3, timestamp, "pca")
    final_json = merge_contents(json_data1, json_data2, json_data3)
    upload_doc(vars, final_json)
    print "DONE"
