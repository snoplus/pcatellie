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
        filename_cd = sys.argv[1]
        filename_tw = sys.argv[2]
        filename_peak = sys.argv[3]
        return filename_cd, filename_tw, filename_peak

def load_env():
    c_add = os.getenv("COUCH_ADD")
    c_db = os.getenv("COUCH_DB")
    c_user = os.getenv("COUCHDB_USER")
    c_pw = os.getenv("COUCHDB_PW")
    logs_folder1 = os.getenv("RUNTIME_LOC")
    logs_folder2 = os.getenv("RUNTIME_LOC")
    logs_folder3 = os.getenv("RUNTIME_LOC")
    return c_add, c_db, c_user, c_pw, logs_folder1, logs_folder2, logs_folder3

def prepare_timestamp():
    timestamp_temp = datetime.now()
    timestamp = timestamp_temp.strftime("%d-%b-%Y %H:%M:%S.%f")
    return timestamp

def load_log(vars, filename, z):
    with open(vars[z] + filename) as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

def parse_content(content, timestamp, type):
    if type == "cd":
        out_ID = []
        out_diff = []
        nowGood = []
        nowBad = []
        for line in content:
            if "Run1:" in line:
                Run1 = int(re.search(r'\d+', line.split(":")[1]).group())
            if "Run2:" in line:
                Run2 = int(re.search(r'\d+', line.split(":")[1]).group())
            if "PMTOff:" in line:
                PMTOff = int(re.search(r'\d+', line).group())
            if "PMTZeroOc:" in line:
                PMTZeroOc = int(re.search(r'\d+', line).group())
            if "PMTLowOc:" in line:
                PMTLowOc = int(re.search(r'\d+', line).group())
            if "PMTGoodCal:" in line:
                PMTGoodCal = int(re.search(r'\d+', line).group())
            if "minDiff:" in line:
                minDiff = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
            if "maxDiff:" in line:
                maxDiff = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
            if "out_ID:" in line:
                for num in line.split(","):
                    out_ID.append( int(re.search(r'\d+', num).group()) )
            if "out_diff:" in line:
                for bit in line.split("[")[1:]:
                    temp_vec = []
                    q = bit.split("]")[0].split(",")
                    temp_vec.append(float(q[0]))
                    temp_vec.append(float(q[1]))
                    out_diff.append( temp_vec )
            if "nowGood:" in line:
                for num in line.split(","):
                    nowGood.append( int(re.search(r'\d+', num).group()) )
            if "nowBad:" in line:
                for num in line.split(","):
                    nowBad.append( int(re.search(r'\d+', num).group()) )
        json_data = {
        "Run1": Run1,
        "Run2": Run2,
        "PMTOff": PMTOff,
        "PMTZeroOc": PMTZeroOc,
        "PMTLowOc": PMTLowOc,
        "PMTGoodCal": PMTGoodCal,
        "minDiff": minDiff,
        "maxDiff": maxDiff,
        "out_ID": out_ID,
        "out_diff": out_diff,
        "nowGood": nowGood,
        "nowBad": nowBad,
        "type": "benchmark",
        "timestamp": timestamp
        }
        return json_data
    elif type == "tw":
        for line in content:
            if "IDs1:" in line:
                IDs1 = int(re.search(r'\d+', line.split(":")[1]).group())
            if "IDs2:" in line:
                IDs2 = int(re.search(r'\d+', line.split(":")[1]).group())
            if "off1:" in line:
                off1 = int(re.search(r'\d+', line.split(":")[1]).group())
            if "off2:" in line:
                off2 = int(re.search(r'\d+', line.split(":")[1]).group())
            if "goodCount:" in line:
                goodCount = int(re.search(r'\d+', line).group())
            if "badCount:" in line:
                badCount = int(re.search(r'\d+', line).group())
            if "Run1:" in line:
                Run1TW = int(re.search(r'\d+', line.split(":")[1]).group())
            if "Run2:" in line:
                Run2TW = int(re.search(r'\d+', line.split(":")[1]).group())
        json_data2 = {
        "Run1TW": Run1TW,
        "Run2TW": Run2TW,
        "IDs1": IDs1,
        "IDs2": IDs2,
        "off1": off1,
        "off2": off2,
        "goodCount": goodCount,
        "badCount": badCount
        }
        return json_data2
    elif type == "peak":
        for line in content:
            if "Run1:" in line:
                Run1Peak = int(re.search(r'\d+', line.split(":")[1]).group())
            if "Run2:" in line:
                Run2Peak = int(re.search(r'\d+', line.split(":")[1]).group())
            if "root1:" in line:
                root1 = line.split(":")[1][1:]
            if "root2:" in line:
                root2 = line.split(":")[1][1:]
        json_data3 = {
        "Run1Peak": Run1Peak,
        "Run2Peak": Run2Peak,
        "root1": root1,
        "root2": root2,
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
    filename_cd, filename_tw, filename_peak = get_args()
    vars = load_env()
    timestamp = prepare_timestamp()
    content1 = load_log(vars, filename_cd, 4)
    content2 = load_log(vars, filename_tw, 5)
    content3 = load_log(vars, filename_peak, 6)
    json_data1 = parse_content(content1, timestamp, "cd")
    json_data2 = parse_content(content2, timestamp, "tw")
    json_data3 = parse_content(content3, timestamp, "peak")
    final_json = merge_contents(json_data1, json_data2, json_data3)
    upload_doc(vars, final_json)
    print "DONE"
