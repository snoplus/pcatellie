import couchdb
import os
from datetime import datetime

def load_env():
    DIR_LIGHT_ANG = os.getenv("DIR_LIGHT_ANG")
    REF_LIGHT_ANG = os.getenv("REF_LIGHT_ANG")
    LOCALITY = os.getenv("LOCALITY")
    ANG_SYS_ANG = os.getenv("ANG_SYS_ANG")
    MIN_DIST = os.getenv("MIN_DIST")
    TOT_EVS = os.getenv("TOT_EVS")
    TOT_EVS_DEV = os.getenv("TOT_EVS_DEV")
    TOT_HITS = os.getenv("TOT_HITS")
    TH_EXTA = os.getenv("TH_EXTA")
    TH_BADCHAN = os.getenv("TH_BADCHAN")
    TH_BADECA = os.getenv("TH_BADECA")
    TH_BADPCA = os.getenv("TH_BADPCA")
    TH_XTALK = os.getenv("TH_XTALK")
    TH_NOTEN = os.getenv("TH_NOTEN")
    TH_OFFPMT = os.getenv("TH_OFFPMT")
    TH_OFFCHAN = os.getenv("TH_OFFCHAN")
    TH_DAQEN = os.getenv("TH_DAQEN")
    TH_NNORM = os.getenv("TH_NNORM")
    TH_BADPOS = os.getenv("TH_BADPOS")
    TH_LPCTIR = os.getenv("TH_LPCTIR")
    TH_LPCIP = os.getenv("TH_LPCIP")
    TH_LPCLOC = os.getenv("TH_LPCLOC")
    TH_LPCPATH = os.getenv("TH_LPCPATH")
    TH_ANG = os.getenv("TH_ANG")
    TH_NEARREF = os.getenv("TH_NEARREF")
    NHIT_DIV = os.getenv("NHIT_DIV")
    NHIT_MEAN = os.getenv("NHIT_MEAN")
    NHIT_MEAN_DEV = os.getenv("NHIT_MEAN_DEV")
    NHIT_RMS = os.getenv("NHIT_RMS")
    SUB_NHIT = os.getenv("SUB_NHIT")
    SUB_DEV = os.getenv("SUB_DEV")
    HIT_PEAK = os.getenv("HIT_PEAK")
    HIT_PEAK_DEV = os.getenv("HIT_PEAK_DEV")
    TIME_DIV = os.getenv("TIME_DIV")
    PMTS_BS = os.getenv("PMTS_BS")
    PMTS_BSP = os.getenv("PMTS_BSP")
    PMTS_BS_GO = os.getenv("PMTS_BS_GO")
    PMTS_BS_GOP = os.getenv("PMTS_BS_GOP")
    PMTS_BS_RAT = os.getenv("PMTS_BS_RAT")
    RUNTIME = os.getenv("RUNTIME")
    FREQ = os.getenv("FREQ")
    FREQ_DEV = os.getenv("FREQ_DEV")
    EV_SUB = os.getenv("EV_SUB")
    EV_SUB_DEV = os.getenv("EV_SUB_DEV")
    N_SUBS = os.getenv("N_SUBS")
    INTEG_ALL = os.getenv("INTEG_ALL")
    INTEG_BS = os.getenv("INTEG_BS")
    BS_PMTS = os.getenv("BS_PMTS")
    PIN_RMS = os.getenv("PIN_RMS")
    COV = os.getenv("COV")
    CORF = os.getenv("CORF")
    TOF_MEAN = os.getenv("TOF_MEAN")
    TOF_RMS = os.getenv("TOF_RMS")
    BUCK_MEAN = os.getenv("BUCK_MEAN")
    BUCK_RMS = os.getenv("BUCK_RMS")
    AS_MEAN_MIN = os.getenv("AS_MEAN_MIN")
    AS_MEAN_MAX = os.getenv("AS_MEAN_MAX")
    AS_RMS = os.getenv("AS_RMS")
    COR_DEV = os.getenv("COR_DEV")
    TOF_MIN = os.getenv("TOF_MIN")
    TOF_MAX = os.getenv("TOF_MAX")
    BUC_MIN = os.getenv("BUC_MIN")
    BUC_MAX = os.getenv("BUC_MAX")
    AS_MIN = os.getenv("AS_MIN")
    AS_MAX = os.getenv("AS_MAX")
    RESID_PEAK = os.getenv("RESID_PEAK")
    RESID_PEAK_DEV = os.getenv("RESID_PEAK_DEV")
    DIR_RESID_0 = os.getenv("DIR_RESID_0")
    DIR_RESID_1 = os.getenv("DIR_RESID_1")
    FIT_SLOPE = os.getenv("FIT_SLOPE")
    FIT_SLOPE_DEV = os.getenv("FIT_SLOPE_DEV")
    return DIR_LIGHT_ANG, REF_LIGHT_ANG, LOCALITY, ANG_SYS_ANG, MIN_DIST, TOT_EVS, TOT_EVS_DEV, TOT_HITS, TH_EXTA, TH_BADCHAN, TH_BADECA, TH_BADPCA, TH_XTALK, TH_NOTEN, TH_OFFPMT, TH_OFFCHAN, TH_DAQEN, TH_NNORM, TH_BADPOS, TH_LPCTIR, TH_LPCIP, TH_LPCLOC, TH_LPCPATH, TH_ANG, TH_NEARREF, NHIT_DIV, NHIT_MEAN, NHIT_MEAN_DEV, NHIT_RMS, SUB_NHIT, SUB_DEV, HIT_PEAK, HIT_PEAK_DEV, TIME_DIV, PMTS_BS, PMTS_BSP, PMTS_BS_GO, PMTS_BS_GOP, PMTS_BS_RAT, RUNTIME, FREQ, FREQ_DEV, EV_SUB, EV_SUB_DEV, N_SUBS, INTEG_ALL, INTEG_BS, BS_PMTS, PIN_RMS, COV, CORF, TOF_MEAN, TOF_RMS, BUCK_MEAN, BUCK_RMS, AS_MEAN_MIN, AS_MEAN_MAX, AS_RMS, COR_DEV, TOF_MIN, TOF_MAX, BUC_MIN, BUC_MAX, AS_MIN, AS_MAX, RESID_PEAK, RESID_PEAK_DEV, DIR_RESID_0, DIR_RESID_1, FIT_SLOPE, FIT_SLOPE_DEV

def prepare_timestamp():
    timestamp_temp = datetime.now()
    timestamp = timestamp_temp.strftime("%d-%b-%Y %H:%M:%S.%f")
    return timestamp

def prepare_json(vars, timestamp):
    json_data = {
    "DIR_LIGHT_ANG": vars[0],
    "REF_LIGHT_ANG": vars[1],
    "LOCALITY": vars[2],
    "ANG_SYS_ANG": vars[3],
    "MIN_DIST": vars[4],
    "TOT_EVS": vars[5],
    "TOT_EVS_DEV": vars[6],
    "TOT_HITS": vars[7],
    "TH_EXTA": vars[8],
    "TH_BADCHAN": vars[9],
    "TH_BADECA": vars[10],
    "TH_BADPCA": vars[11],
    "TH_XTALK": vars[12],
    "TH_NOTEN": vars[13],
    "TH_OFFPMT": vars[14],
    "TH_OFFCHAN": vars[15],
    "TH_DAQEN": vars[16],
    "TH_NNORM": vars[17],
    "TH_BADPOS": vars[18],
    "TH_LPCTIR": vars[19],
    "TH_LPCIP": vars[20],
    "TH_LPCLOC": vars[21],
    "TH_LPCPATH": vars[22],
    "TH_ANG": vars[23],
    "TH_NEARREF": vars[24],
    "NHIT_DIV": vars[25],
    "NHIT_MEAN": vars[26],
    "NHIT_MEAN_DEV": vars[27],
    "NHIT_RMS": vars[28],
    "SUB_NHIT": vars[29],
    "SUB_DEV": vars[30],
    "HIT_PEAK": vars[31],
    "HIT_PEAK_DEV": vars[32],
    "TIME_DIV": vars[33],
    "PMTS_BS": vars[34],
    "PMTS_BSP": vars[35],
    "PMTS_BS_GO": vars[36],
    "PMTS_BS_GOP": vars[37],
    "PMTS_BS_RAT": vars[38],
    "RUNTIME": vars[39],
    "FREQ": vars[40],
    "FREQ_DEV": vars[41],
    "EV_SUB": vars[42],
    "EV_SUB_DEV": vars[43],
    "N_SUBS": vars[44],
    "INTEG_ALL": vars[45],
    "INTEG_BS": vars[46],
    "BS_PMTS": vars[47],
    "PIN_RMS": vars[48],
    "COV": vars[49],
    "CORF": vars[50],
    "TOF_MEAN": vars[51],
    "TOF_RMS": vars[52],
    "BUCK_MEAN": vars[53],
    "BUCK_RMS": vars[54],
    "AS_MEAN_MIN": vars[55],
    "AS_MEAN_MAX": vars[56],
    "AS_RMS": vars[57],
    "COR_DEV": vars[58],
    "TOF_MIN": vars[59],
    "TOF_MAX": vars[60],
    "BUC_MIN": vars[61],
    "BUC_MAX": vars[62],
    "AS_MIN": vars[63],
    "AS_MAX": vars[64],
    "RESID_PEAK": vars[65],
    "RESID_PEAK_DEV": vars[66],
    "DIR_RESID_0": vars[67],
    "DIR_RESID_1": vars[68],
    "FIT_SLOPE": vars[69],
    "FIT_SLOPE_DEV": vars[70],
    "timestamp": timestamp,
    "type": "env"
    }
    return json_data

def load_chdb_vars():
    c_add = os.getenv("COUCH_ADD")
    c_db = os.getenv("COUCH_DB")
    c_user = os.getenv("COUCHDB_USER")
    c_pw = os.getenv("COUCHDB_PW")

    return c_add, c_db, c_user, c_pw

def get_existing_doc(c_add, c_db, c_user, c_pw):
    couchserver = couchdb.Server("http://"+c_user+":"+c_pw+"@"+c_add)
    db = couchserver[c_db]
    runRows = db.view("_design/env/_view/env", include_docs=True)
    IDs = []
    for runRow in runRows:
        id = runRow.id
        IDs.append(id)
    if len(IDs) == 1:
        doc_id = IDs[0]
        return doc_id
    else:
        "There are conflicting env documents!"
        exit()

def upload_to_chdb(c_add, c_db, c_user, c_pw, json_data, doc_id):
    couchserver = couchdb.Server("http://"+c_user+":"+c_pw+"@"+c_add)
    db = couchserver[c_db]
    doc = db[doc_id]
    # iterate through the json, update each element
    for key, value in json_data.items():
        doc[key] = value
    try:
        # save the updated doc
        doc = db.save(doc)
        return
    except:
        print "Couldn't update the doc :/"
        exit()

if __name__=="__main__":
    vars = load_env()
    timestamp = prepare_timestamp()
    json_data = prepare_json(vars, timestamp)
    c_add, c_db, c_user, c_pw = load_chdb_vars()
    doc_id = get_existing_doc(c_add, c_db, c_user, c_pw)
    upload_to_chdb(c_add, c_db, c_user, c_pw, json_data, doc_id)

    print "DONE"
