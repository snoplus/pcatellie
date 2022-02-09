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

def parse_content(content, timestamp):
    for line in content:
        if "Run:" in line:
            run = int(re.search(r'\d+', line).group())
        if "Channel:" in line:
            channel = int(re.search(r'\d+', line).group())
        if "Fibre:" in line:
            fibre = line.split(" ")[1]
        if "Mode:" in line:
            mode = line.split(" ")[1]
        if "Subs:" in line:
            subs = int(re.search(r'\d+', line).group())
        if "Total events:" in line:
            tot_evs = int(re.search(r'\d+', line).group())
        if "EXTA events:" in line:
            exta_evs = int(re.search(r'\d+', line).group())
        if "Passed events:" in line:
            pass_evs = int(re.search(r'\d+', line).group())
        if "CouchDB events:" in line:
            couch_evs = int(re.search(r'\d+', line).group())
        if "Total hits:" in line:
            tot_hits = int(re.search(r'\d+', line).group())
        if "Passed hits:" in line:
            pass_hits = int(re.search(r'\d+', line).group())
        if "Subruns:" in line:
            subruns = int(re.search(r'\d+', line).group())
        if "Subruns (loop):" in line:
            subs_loop = int(re.search(r'\d+', line).group())
        if "Passed evs:" in line:
            passed_evs = int(re.search(r'\d+', line).group())
        if "Events DB:" in line:
            evs_db = int(re.search(r'\d+', line).group())
        if "All hits:" in line:
            all_hits = int(re.search(r'\d+', line).group())
        if "CNEXTA:" in line:
            CNEXTA = line.split(" ")[1].split("\t")
        if "CCHS:" in line:
            CCHS = line.split(" ")[1].split("\t")
        if "CECA:" in line:
            CECA = line.split(" ")[1].split("\t")
        if "CPCA:" in line:
            CPCA = line.split(" ")[1].split("\t")
        if "CXT:" in line:
            CXT = line.split(" ")[1].split("\t")
        if "CE:" in line:
            CE = line.split(" ")[1].split("\t")
        if "COFF:" in line:
            COFF = line.split(" ")[1].split("\t")
        if "CCO:" in line:
            CCO = line.split(" ")[1].split("\t")
        if "CDAQ:" in line:
            CDAQ = line.split(" ")[1].split("\t")
        if "CNN:" in line:
            CNN = line.split(" ")[1].split("\t")
        if "CMAG:" in line:
            CMAG = line.split(" ")[1].split("\t")
        if "CTIR:" in line:
            CTIR = line.split(" ")[1].split("\t")
        if "CPV:" in line:
            CPV = line.split(" ")[1].split("\t")
        if "CRH:" in line:
            CRH = line.split(" ")[1].split("\t")
        if "CDAV:" in line:
            CDAV = line.split(" ")[1].split("\t")
        if "CANG:" in line:
            CANG = line.split(" ")[1].split("\t")
        if "CDIST:" in line:
            CDIST = line.split(" ")[1].split("\t")
        if "Rolling NHit mean:" in line:
            rol_n_mean = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Rolling NHit min:" in line:
            rol_n_min = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Rolling NHit max:" in line:
            rol_n_max = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "NHit distribution, mean:" in line:
            nD_mean = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "NHit distribution, rms:" in line:
            nD_rms = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Hit peak:" in line:
            hit_peak = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Mean of the time over event distribution:" in line:
            ToE_mean = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Direct hits, peaks:" in line:
            dh_peaks = int(re.search(r'\d+', line).group())
        if "Direct hits, peak time:" in line:
            dh_peak = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "PMTs in BS:" in line:
            PMT_BS = line.split(":")[1].split("\t")
        if "PMTs in BS, good occup:" in line:
            PMTS_BS_good = line.split(":")[1].split("\t")
        if "Good PMTs:" in line:
            good_PMT = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "EXTAlength:" in line:
            EXTAlength = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Frequency from length:" in line:
            f_from_length = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Occup 1-5 to all:" in line:
            occup_15 = re.findall(r"[-+]?\d*\.\d+|\d+", line.split(":")[1])[0]
        if "Occup BM 1-5 to all:" in line:
            occup_BS = re.findall(r"[-+]?\d*\.\d+|\d+", line.split(":")[1])[0]
        if "Good PMTs BS 1-5 occ:" in line:
            good_PMTS_BS = int(re.search(r'\d+', line.split(":")[1]).group())
        if "PIN RMS:" in line:
            pin_rms = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Covariance:" in line:
            cov = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Correlation factor:" in line:
            cor = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "Subruns flag:" in line:
            Fsubs = int(re.search(r'\d+', line).group())
        if "Event flag:" in line:
            Fevs = int(re.search(r'\d+', line).group())
        if "Hit flag:" in line:
            Fhit = int(re.search(r'\d+', line).group())
        if "EXTA flag:" in line:
            Fexta = int(re.search(r'\d+', line).group())
        if "Bad channel flag:" in line:
            FbadC = int(re.search(r'\d+', line).group())
        if "Bad ECA flag:" in line:
            FECA = int(re.search(r'\d+', line).group())
        if "BAD PCA flag:" in line:
            FPCA = int(re.search(r'\d+', line).group())
        if "X-talk flag:" in line:
            FX = int(re.search(r'\d+', line).group())
        if "Not enabled flag:" in line:
            FnEn = int(re.search(r'\d+', line).group())
        if "Offline PMT flag:" in line:
            FOffPMT = int(re.search(r'\d+', line).group())
        if "Offline channel flag:" in line:
            FOffC = int(re.search(r'\d+', line).group())
        if "Not DAQ enabled flag:" in line:
            FnDAQ = int(re.search(r'\d+', line).group())
        if "Not normal PMT flag:" in line:
            FnN = int(re.search(r'\d+', line).group())
        if "Bad PMT position flag:" in line:
            FbadPos = int(re.search(r'\d+', line).group())
        if "LPC TIR flag:" in line:
            FTIR = int(re.search(r'\d+', line).group())
        if "LPC invalid path flag:" in line:
            FLPCIP = int(re.search(r'\d+', line).group())
        if "LPC locality flag:" in line:
            FLPCloc = int(re.search(r'\d+', line).group())
        if "LPC weird path (not AV) flag:" in line:
            FLPCW = int(re.search(r'\d+', line).group())
        if "Angular cut flag:" in line:
            Fang = int(re.search(r'\d+', line).group())
        if "Near reflection flag:" in line:
            Fnear = int(re.search(r'\d+', line).group())
        if "Stable NHit flag:" in line:
            FstabN = int(re.search(r'\d+', line).group())
        if "NHit distribution flag:" in line:
            FNDist = int(re.search(r'\d+', line).group())
        if "NHit distribution over subruns flag:" in line:
            FNDistSubs = int(re.search(r'\d+', line).group())
        if "Trigger delay dev flag:" in line:
            Ftrig = int(re.search(r'\d+', line).group())
        if "Fibre delay dev flag:" in line:
            Ffib = int(re.search(r'\d+', line).group())
        if "Hit peak flag:" in line:
            FhitPeak = int(re.search(r'\d+', line).group())
        if "Time over event flag:" in line:
            FToE = int(re.search(r'\d+', line).group())
        if "Peak check flag:" in line:
            Fpeak = int(re.search(r'\d+', line).group())
        if "Check on PMTs (beamspot / occupancy) flag:" in line:
            FBSOcc = int(re.search(r'\d+', line).group())
        if "Runtime flag:" in line:
            Fruntime = int(re.search(r'\d+', line).group())
        if "Check on evens in subruns:" in line:
            FEvsSubs = int(re.search(r'\d+', line).group())
        if "Frequency flag:" in line:
            FFreq = int(re.search(r'\d+', line).group())
        if "Occupancy all flag:" in line:
            FOcc = int(re.search(r'\d+', line).group())
        if "Occupancy beamspot flag:" in line:
            FOccBS = int(re.search(r'\d+', line).group())
        if "Beamspot PMTs flag:" in line:
            FBSPMT = int(re.search(r'\d+', line).group())
        if "PIN RMS flag:" in line:
            FPINRMS = int(re.search(r'\d+', line).group())
        if "PIN-NHit cov and corr flag:" in line:
            FPINNHIT = int(re.search(r'\d+', line).group())

    bitword = str(Fsubs)+str(Fevs)+str(Fhit)+str(Fexta)+str(FbadC)+str(FECA)+str(FPCA)+str(FX)+str(FnEn)+str(FOffPMT)+str(FOffC)+str(FnDAQ)+str(FnN)+str(FbadPos)+str(FTIR)+str(FLPCIP)+str(FLPCloc)+str(FLPCW)+str(Fang)+str(Fnear)+str(FstabN)+str(FNDist)+str(FNDistSubs)+str(Ftrig)+str(Ffib)+str(FhitPeak)+str(FToE)+str(Fpeak)+str(FBSOcc)+str(Fruntime)+str(FEvsSubs)+str(FFreq)+str(FOcc)+str(FOccBS)+str(FBSPMT)+str(FPINRMS)+str(FPINNHIT)
    json_data = {
    "run": run,
    "channel": channel,
    "fibre": fibre,
    "mode": mode,
    "subs": subs,
    "Total events": tot_evs,
    "EXTA events": exta_evs,
    "Passed events": pass_evs,
    "CouchDB events": couch_evs,
    "Total hits": tot_hits,
    "Passed hits": pass_hits,
    "Subruns": subruns,
    "Subruns (loop)": subs_loop,
    "Passed evs": passed_evs,
    "Events DB": evs_db,
    "All hits": all_hits,
    "CNEXTA": CNEXTA,
    "CCHS": CCHS,
    "CECA": CECA,
    "CPCA": CPCA,
    "CXT": CXT,
    "CE": CE,
    "COFF": COFF,
    "CCO": CCO,
    "CDAQ": CDAQ,
    "CNN": CNN,
    "CMAG": CMAG,
    "CTIR": CTIR,
    "CPV": CPV,
    "CRH": CRH,
    "CDAV": CDAV,
    "CANG": CANG,
    "CDIST": CDIST,
    "Rolling NHit mean": rol_n_mean,
    "Rolling NHit min": rol_n_min,
    "Rolling NHit max": rol_n_max,
    "NHit distribution, mean": nD_mean,
    "NHit distribution, rms": nD_rms,
    "Hit peak": hit_peak,
    "Mean of the time over event distribution": ToE_mean,
    "Direct hits, peaks": dh_peaks,
    "Direct hits, peak time": dh_peak,
    "PMTs in BS": PMT_BS,
    "PMTs in BS, good occup": PMTS_BS_good,
    "Good PMTs": good_PMT,
    "EXTAlength": EXTAlength,
    "Frequency from length": f_from_length,
    "Occup 1-5 to all": occup_15,
    "Occup BM 1-5 to all": occup_BS,
    "Good PMTs BS 1-5 occ": good_PMTS_BS,
    "PIN RMS": pin_rms,
    "Covariance": cov,
    "Correlation factor": cor,
    "Subruns flag": Fsubs,
    "Event flag": Fevs,
    "Hit flag": Fhit,
    "EXTA flag": Fexta,
    "Bad channel flag": FbadC,
    "Bad ECA flag": FECA,
    "BAD PCA flag": FPCA,
    "X-talk flag": FX,
    "Not enabled flag": FnEn,
    "Offline PMT flag": FOffPMT,
    "Offline channel flag": FOffC,
    "Not DAQ enabled flag": FnDAQ,
    "Not normal PMT flag": FnN,
    "Bad PMT position flag": FbadPos,
    "LPC TIR flag": FTIR,
    "LPC invalid path flag": FLPCIP,
    "LPC locality flag": FLPCloc,
    "LPC weird path (not AV) flag": FLPCW,
    "Angular cut flag": Fang,
    "Near reflection flag": Fnear,
    "Stable NHit flag": FstabN,
    "NHit distribution flag": FNDist,
    "NHit distribution over subruns flag": FNDistSubs,
    "Trigger delay dev flag": Ftrig,
    "Fibre delay dev flag": Ffib,
    "Hit peak flag": FhitPeak,
    "Time over event flag": FToE,
    "Peak check flag": Fpeak,
    "Check on PMTs (beamspot / occupancy) flag": FBSOcc,
    "Runtime flag": Fruntime,
    "Check on evens in subruns": FEvsSubs,
    "Frequency flag": FFreq,
    "Occupancy all flag": FOcc,
    "Occupancy beamspot flag": FOccBS,
    "Beamspot PMTs flag": FBSPMT,
    "PIN RMS flag": FPINRMS,
    "PIN-NHit cov and corr flag": FPINNHIT,
    "bitword": bitword,
    "type": 'val1',
    "timestamp": timestamp
    }
    return json_data

def prepare_timestamp():
    timestamp_temp = datetime.now()
    timestamp = timestamp_temp.strftime("%d-%b-%Y %H:%M:%S.%f")
    return timestamp

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
    vars = load_env()
    timestamp = prepare_timestamp()
    content = load_log(vars, filename)
    json_data = parse_content(content, timestamp)
    upload_doc(vars, json_data)
    print "DONE"
