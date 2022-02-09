import couchdb
import os
import json
from datetime import datetime

def load_os_vars():
    c_add = os.environ['COUCH_ADD']
    c_db = os.environ['COUCH_DB']
    c_user = os.environ['COUCHDB_USER']
    c_pw = os.environ['COUCHDB_PW']

    return c_add, c_db, c_user, c_pw

def load_couchdb_data():

    couch = couchdb.Server('https://' + c_user + ':' + c_pw + '@' + c_add)
    db = couch[c_db]

    belly_fibres = []

    print "Loading COUCHDB data..."

    for item in db.view('_design/belly/_view/belly'):
        doc = db[item.id]

    belly_fibres = doc["list"]

    return belly_fibres

def get_index(list):
    indexes = []
    for each in list:
        indexes.append( each[2:5] )

    return indexes

if __name__=="__main__":
    c_add, c_db, c_user, c_pw = load_os_vars()
    fibres = load_couchdb_data()
    indexes = get_index(fibres)

    ### Store to ratdb file
    # Generate current timestamp
    timestamp_date = datetime.today().strftime('%Y-%m-%d')
    timestamp_time = datetime.today().strftime('%H:%M:%S')
    timestamp = timestamp_date + "T" + timestamp_time

    json_data = {
    "list": fibres,
    "fibre_index": indexes,
    "run_range": [0, 2147483647],
    "comment": "list of TELLIE fibres affected by belly plates",
    "pass": 1,
    "version": 1,
    "type": "TELLIE_BELLY_PLATES",
    "index": "",
    "timestamp": timestamp
    }

    # Save to ratdb file
    with open("belly_fibres.ratdb", "w") as f:
        json.dump(json_data, f)
        print "stored to ratdb file successfully"
