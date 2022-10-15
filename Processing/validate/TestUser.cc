///////////////////////////////////////////////////////////////////////////////
/// \file TestUser.cc
/// \brief An example User Processor
///
/// This example demonstrates how to write a User Processor. MyUserProc
/// (declared and defined below) just creates a 1-D histogram of the number of
/// photoelectrons in each event, filled it as events are generated and writes
/// it to a file called "numpe.root".
///
/// (Of course, you could easily do this just by writing the events to
/// disk using the outroot processor and making the histogram offline.
/// But then you wouldn't get to see how user processors work!)
///
///////////////////////////////////////////////////////////////////////////////

#include <RAT/Processor.hh>
#include <RAT/Log.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/ChanHWStatus.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/PhysicsUtil.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DB.hh>
#include <RAT/DBLink.hh>
#include <RAT/DBTable.hh>
#include <RAT/json.hh>
#include <RAT/DataQualityProc.hh>

#include <TFitResult.h>
#include <TFile.h>
#include <TMultiGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TVectorD.h>
#include <TGraph2D.h>
#include <TH2D.h>
#include <TView3D.h>
#include <TF2.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGraph2DErrors.h>
#include <TLegend.h>
#include <THStack.h>
#include <TSpectrum.h>

#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>
#include <numeric>
#include <stdint.h>

#include <G4PhysicalConstants.hh>

using namespace std;
using namespace RAT;
using CLHEP::nm;

// Class declaration
class MyUserProc : public Processor {
public:

  MyUserProc();

  virtual ~MyUserProc();

  virtual void BeginOfRun( DS::Run& run );

  virtual void EndOfRun( DS::Run& run );

  virtual Processor::Result DSEvent(DS::Run& run, DS::Entry& ds);

  virtual string printVector(const TVector3& v);

  virtual void Smooth(TGraph* g, TGraph* h, int width);

protected:

  // LPC
  double energy;
  double fLEDWavelength;
  double LOCALITY;

  // RAT stuff
  RAT::DU::GroupVelocity gv;
  RAT::DU::LightPathCalculator lpc;
  int NPMTS;
  int run_id;
  TVector2 flat;
  TVector3 pmtPosAll;
  RAT::DataQualityProc *dq;

  // Validation
  int val_evs;
  double val_evsDev;
  double val_hits;
  double val_exta;
  double val_badChan;
  double val_badECA;
  double val_badPCA;
  double val_Xtalk;
  double val_NotEn;
  double val_offPMT;
  double val_offChan;
  double val_daqEn;
  double val_Nnorm;
  double val_badPos;
  double val_LPCTIR;
  double val_LPCIP;
  double val_LPCLoc;
  double val_LPCPath;
  double val_Ang;
  double val_NearRef;
  double val_stableNHit;
  int val_NHitMean;
  int val_NHitMeanDev;
  double val_NHitRMS;
  double val_subNHit;
  double val_subNHitDev;
  double val_hitPeak;
  double val_hitPeakDev;
  double val_timeDev;
  int pmtsBS;
  int pmtsBSgoodOccup;
  double pmtsBSP;
  double pmtsBSgoodOccupP;
  int val_pmtsBS;
  int val_pmtsBSP;
  int val_pmtsBSGO;
  int val_pmtsBSGOP;
  double val_pmtsBSGO_ratio;
  ULong64_t firstEXTA;
  ULong64_t lastEXTA;
  ULong64_t prevEXTA10;
  ULong64_t prevEXTA50;
  ULong64_t EXTAlength;
  int val_runtime;
  double val_freq;
  double val_freqDev;
  int val_evsSub;
  double val_evsSubDev;
  int val_nsubs;
  int val_integ_all;
  int val_integ_bs;
  int val_bs_pmts;
  double val_pin_rms;
  double val_cov;
  double val_corf;

  int evsFlag;
  int hitFlag;
  int extaFlag;
  int badChanFlag;
  int badECAFlag;
  int badPCAFlag;
  int XtalkFlag;
  int NotEnFlag;
  int offPMTFlag;
  int offChanFlag;
  int daqEnFlag;
  int NnormFlag;
  int badPosFlag;
  int LPCTIRFlag;
  int LPCIPFlag;
  int LPCLocFlag;
  int LPCPathFlag;
  int angFlag;
  int nearRefFlag;
  int stableNHitFlag;
  int NHitDistFlag;
  int NHitSubsFlag;
  int TrigDelFlag;
  int FibDelFlag;
  int hitPeakFlag;
  int timeDevFlag;
  int timeRollFlag;
  int peakFlag;
  int PMTsFlag;
  int runtimeFlag;
  int freqRuntimeFlag;
  int evsSubFlag;
  int nsubsFlag;
  int integAllFlag;
  int integBSFlag;
  int BSpmtsFlag;
  int PINFlag;
  int PINNHitFlag;
  int flag;

  // Database
  RAT::DBLinkPtr tellie_run_data;
  std::string fibre_db;
  std::string run_mode_db;
  json::Value sub_run_info;
  RAT::DBLinkPtr fibre_data;
  TVector3 fibrepos, fibredir;
  int n_subruns;
  int channel;
  int temp_events_db;
  int total_events_db;
  std::vector<double> subruns;
  json::Value current_subrun;
  int PIN;
  float PIN_RMS;
  float fibre_delay;
  int trigger_delay;
  std::vector<double> PINs;
  std::vector<float> PIN_RMSes;
  std::vector<double> PIN_errs;
  std::vector<double> fibre_delays;
  std::vector<int> trigger_delays;

  // Event loop
  int count;
  int allEvs;
  int passedEvs;
  int trig;
  std::vector<double> pmtOccup;
  std::vector< std::vector < std::vector <double> > > fPMTs;
  int NHit;
  int NHit_cleaned;
  std::vector<int> NHits;
  std::vector<int> NHits_cleaned;
  int pmt_hits;
  int pmt_hits_passed;
  int pmt_hits_passed_total;
  double time_min;
  double time_max;
  std::vector< std::vector<double> > EXTAtimes10;
  std::vector< std::vector<double> > EXTAtimes50;
  int prevSubrun;
  int firstEXTAinSubrun;

  // Subruns
  std::vector< std::vector< std::vector<double> > > avgNHit;
  std::vector<double> avgNHits;
  std::vector<double> avgNHitsE;

  // Stats
  double hits_passed;
  int CNEXTA;          // not EXTA evs
  int CCHS;            // bad channel status
  int CECA;            // bad ECA bit
  int CPCA;            // bad PCA bit
  int CXT;             // x-talk
  int COFF;            // offline PMT
  int CE;              // not enabled PMT
  int CCO;             // offline channel
  int CDAQ;            // not DAQ enabled
  int CNN;             // not normal PMT
  int CMAG;            // bad PMT position
  int CTIR;            // lpc: total internal reflection
  int CPV;             // lpc: valid path
  int CRH;             // lpc: within locality
  int CANG;            // ang cut
  int CDIST;           // near reflection
  int CDAV;            // wierd path (not completely through AV)

  // ROOT
  TFile *RootFile;
  stringstream FileName;
  std::string outname;
  TView3D *view;

  TCanvas *fOccupC;
  TCanvas *fOccupBeamSpotC;
  TCanvas *fOccupCutC;
  TCanvas *fDirectHitsC;
  TCanvas *fNHitHistC;
  TCanvas *fNHitvsEventC;
  TCanvas *fNHitvsEventZOOMC;
  TCanvas *fTimevsAngleC;
  TCanvas *fTimevsEventC;
  TCanvas *fPINHistC;
  TCanvas *fPINC;
  TCanvas *fFibreDelaysC;
  TCanvas *fTriggerDelaysC;
  TCanvas *favgNHitHistC;
  TCanvas *favgNHitC;
  TCanvas *fPINvsNHitC;
  TCanvas *fXYThetaC;
  TCanvas *fXYThetaAtPMTC;
  TCanvas *fXYOccupAllC;
  TCanvas *fXYOccupBeamSpotC;
  TCanvas *fEvsSubrunC;
  TCanvas *fEXTADeltaC;

  TH1D *fOccup;
  TH1D *fOccupBeamSpot;
  TH1D *fOccupCut;
  TH1D *fDirectHits;
  TH1D *fNHitHist;
  TH1D *fCNHitHist;
  THStack *fSNHitHists;
  TH2D *fNHitvsEvent;
  TGraph *fNHitvsEventZOOM;
  TH2D *fTimevsAngle;
  TH2D *fTimevsEvent;
  TGraph *fTimevsEventFit;
  TH1D *fPINHist;
  TGraphErrors *fPIN;
  TGraph *fFibreDelays;
  TGraph *fTriggerDelays;
  TH1D *favgNHitHist;
  TGraphErrors *favgNHit;
  TGraphErrors *fPINvsNHit;
  TGraph2D *fXYTheta;
  TGraph2D *fXYThetaAtPMT;
  TGraph2D *fXYOccupAll;
  TGraph2D *fXYOccupBeamSpot;
  TH1D *fEvsSubrun;
  THStack *fEXTADelta;
  TH1D *fEXTADelta10;
  TH1D *fEXTADelta50;
  TH1D *fEXTAStolen;

  // File for results
  stringstream logFile_namess;
  string logFile_name;
  FILE *logFile;
  stringstream logFile_Flags;

};


namespace RAT {
  // Override user processor utility function
  Processor *construct_user_proc(std::string /*userProcName*/) {
    return new MyUserProc;
  }
} // namespace RAT

  // Class definition
  MyUserProc::MyUserProc() : Processor("user") {

    LOCALITY = atoi(getenv("LOCALITY"));
    fLEDWavelength = 506.0*1e-6;
    energy = util::WavelengthToEnergy(fLEDWavelength);
    view = new TView3D();

    // Validation
    val_evs = atoi(getenv("TOT_EVS"));
    val_evsDev = atof(getenv("TOT_EVS_DEV"));
    val_hits = atof(getenv("TOT_HITS"));
    val_exta = atof(getenv("TH_EXTA"));
    val_badChan = atof(getenv("TH_BADCHAN"));
    val_badECA = atof(getenv("TH_BADECA"));
    val_badPCA = atof(getenv("TH_BADPCA"));
    val_Xtalk = atof(getenv("TH_XTALK"));
    val_NotEn = atof(getenv("TH_NOTEN"));
    val_offPMT = atof(getenv("TH_OFFPMT"));
    val_offChan = atof(getenv("TH_OFFCHAN"));
    val_daqEn = atof(getenv("TH_DAQEN"));
    val_Nnorm = atof(getenv("TH_NNORM"));
    val_badPos = atof(getenv("TH_BADPOS"));
    val_LPCTIR = atof(getenv("TH_LPCTIR"));
    val_LPCIP = atof(getenv("TH_LPCIP"));
    val_LPCLoc = atof(getenv("TH_LPCLOC"));
    val_LPCPath = atof(getenv("TH_LPCPATH"));
    val_Ang = atof(getenv("TH_ANG"));
    val_NearRef = atof(getenv("TH_NEARREF"));
    val_stableNHit = atof(getenv("NHIT_DIV"));
    val_NHitMean = atoi(getenv("NHIT_MEAN"));
    val_NHitMeanDev = atoi(getenv("NHIT_MEAN_DEV"));
    val_NHitRMS = atof(getenv("NHIT_RMS"));
    val_subNHit = atof(getenv("SUB_NHIT"));
    val_subNHitDev = atof(getenv("SUB_DEV"));
    val_hitPeak = atof(getenv("HIT_PEAK"));
    val_hitPeakDev = atof(getenv("HIT_PEAK_DEV"));
    val_timeDev = atof(getenv("TIME_DIV"));
    val_pmtsBS = atoi(getenv("PMTS_BS"));
    val_pmtsBSP = atoi(getenv("PMTS_BSP"));
    val_pmtsBSGO = atoi(getenv("PMTS_BS_GO"));
    val_pmtsBSGOP = atoi(getenv("PMTS_BS_GOP"));
    val_pmtsBSGO_ratio = atof(getenv("PMTS_BS_RAT"));
    val_runtime = atoi(getenv("RUNTIME"));
    val_freq = atof(getenv("FREQ"));
    val_freqDev = atof(getenv("FREQ_DEV"))/100;
    val_evsSub = atoi(getenv("EV_SUB"));
    val_evsSubDev = atof(getenv("EV_SUB_DEV"));
    val_nsubs = atoi(getenv("N_SUBS"));
    val_integ_all = atoi(getenv("INTEG_ALL"));
    val_integ_bs = atoi(getenv("INTEG_BS"));
    val_bs_pmts = atoi(getenv("BS_PMTS"));
    val_pin_rms = atof(getenv("PIN_RMS"));
    val_cov = atof(getenv("COV"));
    val_corf = atof(getenv("CORF"));

  }

  MyUserProc::~MyUserProc() {
  }

  void MyUserProc::BeginOfRun( DS::Run& run ) {

    // Initialise counters
    count = 0;
    allEvs = 0;
    passedEvs = 0;
    run_id = run.GetRunID();
    pmt_hits = 0;
    pmt_hits_passed_total = 0;
    time_min = 9999;
    time_max = 0;
    firstEXTAinSubrun = 0;
    prevSubrun = -1;
    CNEXTA = 0;
    CCHS = 0;
    CECA = 0;
    CPCA = 0;
    CXT = 0;
    COFF = 0;
    CE = 0;
    CCO = 0;
    CDAQ = 0;
    CNN = 0;
    CMAG = 0;
    CTIR = 0;
    CPV = 0;
    CRH = 0;
    CANG = 0;
    CDIST = 0;
    CDAV = 0;

    // Validation
    evsFlag = 0;
    hitFlag = 0;
    extaFlag = 0;
    badChanFlag = 0;
    badECAFlag = 0;
    badPCAFlag = 0;
    XtalkFlag = 0;
    NotEnFlag = 0;
    offPMTFlag = 0;
    offChanFlag = 0;
    daqEnFlag = 0;
    NnormFlag = 0;
    badPosFlag = 0;
    LPCTIRFlag = 0;
    LPCIPFlag = 0;
    LPCLocFlag = 0;
    LPCPathFlag = 0;
    angFlag = 0;
    nearRefFlag = 0;
    NHitSubsFlag = 0;
    TrigDelFlag = 0;
    FibDelFlag = 0;
    hitPeakFlag = 0;
    timeDevFlag = 0;
    timeRollFlag = 0;
    peakFlag = 0;
    pmtsBS = 0;
    pmtsBSgoodOccup = 0;
    pmtsBSP = 0;
    pmtsBSgoodOccupP = 0;
    PMTsFlag = 0;
    runtimeFlag = 0;
    freqRuntimeFlag = 0;
    evsSubFlag = 0;
    nsubsFlag = 0;
    integAllFlag = 0;
    integBSFlag = 0;
    BSpmtsFlag = 0;
    PINFlag = 0;
    PINNHitFlag = 0;
    n_subruns = 0;

    // Initialise RAT vars
    gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    gv.BeginOfRun();
    lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lpc.BeginOfRun();
    lpc.SetELLIEEvent(true);
    dq = new RAT::DataQualityProc("");
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
    NPMTS = pmtinfo.GetCount();
    cout << "Initialised RAT. Number of PMTs is " << NPMTS << "." << endl;

    // Reset vectors
    pmtOccup.clear();
    pmtOccup.resize(NPMTS);
    fPMTs.clear();
    fPMTs.resize(NPMTS);

    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      pmtPosAll = pmtinfo.GetPosition(iPMT);
      int face;
      flat = dq->IcosProject(pmtPosAll, face);
      fPMTs[iPMT].resize(9);
      for (int g=0; g<9; g++){
        fPMTs[iPMT][g].resize(1);
      }
      fPMTs[iPMT][0][0] = iPMT;         // ID
      fPMTs[iPMT][1][0] = flat.X();     // x pos
      fPMTs[iPMT][2][0] = flat.Y();     // y pos
                                        // direct hits
      fPMTs[iPMT][4][0] = 0;            // in beam spot
      fPMTs[iPMT][5][0] = -5;           // angle at fibre
      fPMTs[iPMT][6][0] = -5;           // angle at bucket
      fPMTs[iPMT][7][0] = 0;            // overall occup (all hits)
      // this here marks offline (or otherwiser bad PMTs)
      if (chs.IsTubeOnline(iPMT) && chs.IsEnabled() && chs.IsChannelOnline(iPMT) && chs.IsDAQEnabled(iPMT)){
        fPMTs[iPMT][8][0] = 1;             // mark online
      } else {
        fPMTs[iPMT][8][0] = -1;            // mark as offline / bad
      }
    }

    // Load run fibre details
    try {
      tellie_run_data = DB::Get()->GetLink("TELLIE_RUN");
      sub_run_info = tellie_run_data->GetJSON("sub_run_info");
      n_subruns = sub_run_info.getArraySize();
      cout << "n_subs from couch: " << n_subruns << endl;

      // Access first subtun to get stuff
      run_mode_db = sub_run_info[0]["run_mode"].getString();
      fibre_db = sub_run_info[0]["fibre"].getString();
      channel = sub_run_info[0]["channel"].getInteger();

      total_events_db = 0;
      // Loop through subruns, store data
      for ( int i = 0; i < n_subruns; i++ ) {
        subruns.push_back((double)i);
        current_subrun = sub_run_info[i];
        temp_events_db = current_subrun["number_of_shots"].getInteger();
        PIN = current_subrun["pin_value"].getInteger();
        try {
          PIN_RMS = current_subrun["pin_rms"].getReal();
        } catch (...) {
          PIN_RMS = current_subrun["pin_rms"].getInteger();
        }
        try {
          fibre_delay = current_subrun["fibre_delay"].getReal();
        } catch (...) {
          fibre_delay = current_subrun["fibre_delay"].getInteger();
        }
        trigger_delay = current_subrun["trigger_delay"].getInteger();
        PINs.push_back((double)PIN);
        PIN_RMSes.push_back(PIN_RMS);
        PIN_errs.push_back((double)(PIN_RMS/sqrt(temp_events_db)));
        fibre_delays.push_back((double)fibre_delay);
        trigger_delays.push_back(trigger_delay);
        total_events_db += temp_events_db;
      }
    }
    catch (DBNotFoundError& e) {
      std::cout << "Proc(): Couldn't load TELLIE Run details from ratdb!" << std::endl;
    }
    cout << "DB: " << endl;
    cout << "FIBRE: " << fibre_db << endl;
    cout << "Channel: " << channel << endl;
    cout << "mode: " << run_mode_db << endl << endl;

    // Load fibre position and direction
    try {
      fibre_data = DB::Get()->GetLink("FIBRE", fibre_db);
      fibrepos.SetXYZ(fibre_data->GetD("x"), fibre_data->GetD("y"), fibre_data->GetD("z"));  // position
      fibredir.SetXYZ(fibre_data->GetD("u"), fibre_data->GetD("v"), fibre_data->GetD("w"));  // direction
    }
    catch (DBNotFoundError& e) {
      std::cout << "Proc(): Couldn't load FIBRE details from ratdb!" << std::endl;
    }
    cout << "RATDB: fibre " << fibre_db << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;

    // Subrun structure
    avgNHits.resize(n_subruns);
    avgNHitsE.resize(n_subruns);
    avgNHit.resize(n_subruns);
    EXTAtimes10.resize(n_subruns);
    EXTAtimes50.resize(n_subruns);
    for (size_t i=0; i<n_subruns; i++) {
      avgNHit[i].resize(3);             // for subruns: 0=events, 1=hits, 2=nhit by event
      for (size_t ii=0; ii<2; ii++) {
        avgNHit[i][ii].resize(1);
        avgNHit[i][ii][0] = 0;
      }
    }

    // Init event loop graphs
    fNHitvsEvent = new TH2D ("fNHitvsEvent", "", 100, 0, 200000, 80, 0, 80);
    fNHitvsEventZOOM = new TGraph();
    fTimevsEvent = new TH2D ("fTimevsEvent", "", 500, 0, 200000, 500, 270, 330);

    // Output file for stats
    logFile_namess.str("");
    logFile_Flags.str("");
    logFile_namess << run_id << "_val1.log";
    logFile_name = logFile_namess.str();
    logFile = fopen(logFile_name.c_str(), "w");

  }

  void MyUserProc::EndOfRun( DS::Run& run ) {

    cout << time_min << endl;
    cout << time_max << endl;

    // Stats
    double sumTime = 0;
    int countTime = 0;
    double meanTime = 0;
    hits_passed = (double)pmt_hits_passed_total/(double)pmt_hits*100;
    // Turn array of PMTs into percentage (occupancy)
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      pmtOccup[iPMT] = (float)fPMTs[iPMT][7][0]/(float)allEvs;
      if (fPMTs[iPMT][3].size() > 1){
        fPMTs[iPMT][3].erase(fPMTs[iPMT][3].begin());
        sumTime += std::accumulate(fPMTs[iPMT][3].begin(), fPMTs[iPMT][3].end(), 0.0);
        countTime += fPMTs[iPMT][3].size();
      }
    }
    meanTime = sumTime/(double)countTime;

    // Summary
    cout << "Number of events was " << count << endl;
    cout << "Number of EXTA events was " << allEvs << endl;
    cout << "Number of passed EXTA events was " << passedEvs << endl;
    cout << "Total events from CouchDB: " << total_events_db << endl;
    cout << "Total NHits: " << pmt_hits << endl;
    cout << "Total passed NHits: " << pmt_hits_passed_total << " " << hits_passed << "%" << endl << endl;
    // Stats summary
    cout << "Not EXTA: " << CNEXTA << " " << double(CNEXTA)/double(count)*100 << endl;
    cout << "Bad channel status: " << CCHS << " " << double(CCHS)/double(pmt_hits)*100 << endl;
    cout << "Bad ECA: " << CECA << " " << double(CECA)/double(pmt_hits)*100 << endl;
    cout << "Bad PCA: " << CPCA << " " << double(CPCA)/double(pmt_hits)*100 << endl;
    cout << "X-talk: " << CXT << " " << double(CXT)/double(pmt_hits)*100 << endl;
    cout << "Offline PMT: " << COFF << " " << double(COFF)/double(pmt_hits)*100 << endl;
    cout << "Not enabled: " << CE << " " << double(CE)/double(pmt_hits)*100 << endl;
    cout << "Offline channel: " << CCO << " " << double(CCO)/double(pmt_hits)*100 << endl;
    cout << "Not DAQ-enabled: " << CDAQ << " " << double(CDAQ)/double(pmt_hits)*100 << endl;
    cout << "Not normal PMT: " << CNN << " " << double(CNN)/double(pmt_hits)*100 << endl;
    cout << "Bad PMT position: " << CMAG << " " << double(CMAG)/double(pmt_hits)*100 << endl;
    cout << "LPC-TIR: " << CTIR << " " << double(CTIR)/double(pmt_hits)*100 << endl;
    cout << "LPC-invalid path: " << CPV << " " << double(CPV)/double(pmt_hits)*100 << endl;
    cout << "LPC-not within locality: " << CRH << " " << double(CRH)/double(pmt_hits)*100 << endl;
    cout << "Angular cut: " << CANG << " " << double(CANG)/double(pmt_hits)*100 << endl;
    cout << "Near reflection: " << CDIST << " " << double(CDIST)/double(pmt_hits)*100 << endl;
    cout << "Weird path (not through AV): " << CDAV << " " << double(CDAV)/double(pmt_hits)*100 << endl;
    // Print to log File
    // Run stats
    fprintf(logFile, "Run: %i\n", run_id);
    fprintf(logFile, "Channel: %i\n", channel);
    fprintf(logFile, "Fibre: %s\n", fibre_db.c_str());
    fprintf(logFile, "Mode: %s\n", run_mode_db.c_str());
    fprintf(logFile, "Subs: %i\n", n_subruns);
    fprintf(logFile, "\n");
    // Events stats
    fprintf(logFile, "Total events: %i\n", count);
    fprintf(logFile, "EXTA events: %i\n", allEvs);
    fprintf(logFile, "Passed events: %i\n", passedEvs);
    fprintf(logFile, "CouchDB events: %i\n", total_events_db);
    fprintf(logFile, "Total hits: %i\n", pmt_hits);
    fprintf(logFile, "Passed hits: %i\n", pmt_hits_passed_total);
    fprintf(logFile, "\n");

    // ROOT stuff
    fOccupC = new TCanvas ("fOccup", "fOccup", 1024, 768);
    fOccupBeamSpotC = new TCanvas ("fOccupBeamSpotC", "fOccupBeamSpotC", 1024, 768);
    fOccupCutC = new TCanvas ("fOccupCutC", "fOccupCutC", 1024, 768);
    fDirectHitsC = new TCanvas ("fDirectHits", "fDirectHits", 1024, 768);
    fNHitHistC = new TCanvas ("fNHitHist", "fNHitHist", 1024, 768);
    fNHitvsEventC = new TCanvas ("fNHitvsEvent", "fNHitvsEvent", 1024, 768);
    fNHitvsEventZOOMC = new TCanvas ("fNHitvsEventZOOM", "fNHitvsEventZOOM", 1024, 768);
    fTimevsAngleC = new TCanvas ("fTimevsAngleC", "fTimevsAngleC", 1024, 768);
    fTimevsEventC = new TCanvas ("fTimevsEventC", "fTimevsEventC", 1024, 768);
    fPINHistC = new TCanvas ("fPINHistC", "fPINHistC", 1024, 768);
    fPINC = new TCanvas ("fPIN", "fPIN", 1024, 768);
    fFibreDelaysC = new TCanvas ("fFibreDelays", "fFibreDelays", 1024, 768);
    fTriggerDelaysC = new TCanvas ("fTriggerDelays", "fTriggerDelays", 1024, 768);
    favgNHitHistC = new TCanvas ("favgNHitHist", "favgNHitHist", 1024, 768);
    favgNHitC = new TCanvas ("favgNHit", "favgNHit", 1024, 768);
    fPINvsNHitC = new TCanvas ("fPINvsNHit", "fPINvsNHit", 1024, 768);
    fXYThetaC = new TCanvas ("fXYTheta", "fXYTheta", 1024, 768);
    fXYThetaAtPMTC = new TCanvas ("fXYThetaAtPMT", "fXYThetaAtPMT", 1024, 768);
    fXYOccupAllC = new TCanvas ("fXYOccupAll", "fXYOccupAll", 1024, 768);
    fXYOccupBeamSpotC = new TCanvas ("fXYOccupBeamSpot", "fXYOccupBeamSpot", 1024, 768);
    fEvsSubrunC = new TCanvas ("fEvsSubrun", "fEvsSubrun", 1024, 768);
    fEXTADeltaC = new TCanvas ("fEXTADeltaC", "fEXTADeltaC", 1024, 768);

    fOccup = new TH1D ("fOccup", "", 100, 0, *max_element(pmtOccup.begin(), pmtOccup.end())*100+0.01 );
    fOccupBeamSpot = new TH1D ("fOccupBeamSpot", "", 100, 0, *max_element(pmtOccup.begin(), pmtOccup.end())*100+0.01 );
    fOccupCut = new TH1D ("fOccupCut", "", 100, 0.99, 5.01);
    fDirectHits = new TH1D ("fDirectHits", "", 1000, time_min-1, time_max+1);
    fNHitHist = new TH1D("fNHitHist", "", 100, *min_element(NHits.begin(), NHits.end()), *max_element(NHits.begin(), NHits.end()));
    fCNHitHist = new TH1D("fCNHitHist", "", 100, *min_element(NHits_cleaned.begin(), NHits_cleaned.end()), *max_element(NHits_cleaned.begin(), NHits_cleaned.end()));
    fSNHitHists = new THStack("fSNHitHists", "fSNHitHists");
    fTimevsAngle = new TH2D ("fTimevsAngle", "", 100, 0, 12, 500, meanTime-10, meanTime+10);
    fPINHist = new TH1D("fPINHist", "", 10, *min_element(PINs.begin(), PINs.end())+1, *max_element(PINs.begin(), PINs.end())+1);
    fFibreDelays = new TGraph();
    fTriggerDelays = new TGraph();
    fXYTheta = new TGraph2D();
    fXYThetaAtPMT = new TGraph2D();
    fXYOccupAll = new TGraph2D();
    fXYOccupBeamSpot = new TGraph2D();
    fTimevsEventFit = new TGraph();
    fEvsSubrun = new TH1D ("fEvsSubrun", "", 50, 4980, 5005);
    fEXTADelta = new THStack("fEXTADelta", "fEXTADelta");
    fEXTADelta10 = new TH1D ("fEXTADelta10", "", 200, 0.9, 3.1);
    fEXTADelta50 = new TH1D ("fEXTADelta50", "", 200, 0.9, 3.1);
    fEXTAStolen = new TH1D ("fEXTAStolen", "", 200, 0.9, 3.1);

    // Fill in plots
    // PMT loop plots
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      if (pmtOccup[iPMT] > 0){                                          // PMT with hits
        if (fPMTs[iPMT][8][0] != -1){                                   // this removes offline (bad) PMTs
          fOccup->Fill(pmtOccup[iPMT]*100);
        }
        if (fPMTs[iPMT][4][0] == 1){                                    // PMT in beamspot (within 12 deg)
          if (fPMTs[iPMT][8][0] != -1){
            fOccupBeamSpot->Fill(pmtOccup[iPMT]*100);
            pmtsBS++;
          }
          if ( pmtOccup[iPMT] >= 0.01 && pmtOccup[iPMT] <= 0.05){       // PMT in required occup range
            if (fPMTs[iPMT][8][0] != -1){
              fOccupCut->Fill(pmtOccup[iPMT]*100);
              pmtsBSgoodOccup++;
            }
            for (size_t i=0; i<fPMTs[iPMT][3].size(); i++) {
              if (fPMTs[iPMT][8][0] != -1){
                fDirectHits->Fill(fPMTs[iPMT][3][i]);
                fTimevsAngle->Fill(fPMTs[iPMT][5][0], fPMTs[iPMT][3][i]);
              }
            } // times loop
          } // occup cut range
        } // beamspot
      } // used PMT
      // XY plots
      if (fPMTs[iPMT][8][0] != -1){                                    // this marks offline (bad) PMT
        fXYTheta->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], fPMTs[iPMT][5][0]);
        fXYThetaAtPMT->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], fPMTs[iPMT][6][0]);
        fXYOccupAll->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], (double)fPMTs[iPMT][7][0]/(double)allEvs);
      } else {
        fXYTheta->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], -10);
        fXYThetaAtPMT->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], -10);
        fXYOccupAll->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], -0.01);
      }
      if (fPMTs[iPMT][4][0] == 1){
        if (fPMTs[iPMT][8][0] != -1){
          fXYOccupBeamSpot->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], (double)fPMTs[iPMT][3].size()/(double)allEvs);
        } else {
          fXYOccupBeamSpot->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], -0.01);
        }
      } else {
        if (fPMTs[iPMT][8][0] != -1){
          fXYOccupBeamSpot->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], 0);
        } else {
          fXYOccupBeamSpot->SetPoint(iPMT, fPMTs[iPMT][1][0], fPMTs[iPMT][2][0], -0.01);
        }
      }
    } // PMT loop

    // NHit plots
    for (size_t i=0; i<NHits.size(); i++) {
      fNHitHist->Fill(NHits[i]);
    }
    for (size_t i=0; i<NHits_cleaned.size(); i++) {
      fCNHitHist->Fill(NHits_cleaned[i]);
    }

    // CouchDB plots
    for (int sub = 0; sub < n_subruns; sub++){
      fFibreDelays->SetPoint(sub, sub, fibre_delays[sub]);
      fTriggerDelays->SetPoint(sub, sub, trigger_delays[sub]);
      fPINHist->Fill(PINs[sub]);
    }

    // Subruns
    for (size_t i=0; i<n_subruns; i++) {
      avgNHits[i] = (double)avgNHit[i][1][0]/(double)avgNHit[i][0][0];
      fEvsSubrun->Fill(EXTAtimes10[i].size()+1);
    }
    favgNHitHist = new TH1D("favgNHitHist", "", 20, *min_element(avgNHits.begin(), avgNHits.end())-1,*max_element(avgNHits.begin(), avgNHits.end())+1);
    for (size_t i=0; i<avgNHits.size(); i++) {
      favgNHitHist->Fill(avgNHits[i]);
    }
    for (size_t i=0; i<n_subruns; i++) {
      long double loopSum = 0;
      for (size_t ii=0; ii<avgNHit[i][2].size(); ii++) {
        loopSum += pow(avgNHit[i][2][ii] - avgNHits[i], 2);
      }
      avgNHitsE[i] = sqrt(loopSum/(double)(avgNHit[i][0][0]-1));
    }

    // PIN and NHit plots
    std::vector<double> x_errors (n_subruns, 0);
    double *subruns_array = &subruns[0];
    double *x_errors_array = &x_errors[0];
    double *PINs_array = &PINs[0];
    double *PIN_errs_array = &PIN_errs[0];
    fPIN = new TGraphErrors( n_subruns, subruns_array, PINs_array, x_errors_array, PIN_errs_array);
    double *avg_nhit_array = &avgNHits[0];
    double *errors_array = &avgNHitsE[0];
    favgNHit = new TGraphErrors( n_subruns, subruns_array, avg_nhit_array, x_errors_array, errors_array);
    fPINvsNHit = new TGraphErrors( n_subruns, avg_nhit_array, PINs_array, errors_array, PIN_errs_array);

    // Validation checks
    fprintf(logFile, "### VALIDATION ###\n");
    logFile_Flags << "### FLAGS ###" << endl;
    // 0: number of subruns (prevSubrun starts from -1, therefore the +1)
    if ( (n_subruns) == val_nsubs & (prevSubrun+1) == val_nsubs ){ nsubsFlag = 1; }
    logFile_Flags << "Subruns flag: " << nsubsFlag << endl;
    fprintf(logFile, "Subruns: %i\n", n_subruns);
    fprintf(logFile, "Subruns (loop): %i\n", prevSubrun+1);

    // 1: events
    if ( (passedEvs > val_evs-val_evs*val_evsDev) & (passedEvs < val_evs+val_evs*val_evsDev) & (passedEvs > total_events_db-val_evsDev*total_events_db) & (passedEvs < total_events_db+total_events_db*val_evsDev) ) {evsFlag = 1;}
    logFile_Flags << "Event flag: " << evsFlag << endl;
    fprintf(logFile, "Passed evs: %i\n", passedEvs);
    fprintf(logFile, "Events DB: %i\n", total_events_db);

    // 2: passed hits
    if (hits_passed > val_hits) {hitFlag = 1;}
    logFile_Flags << "Hit flag: " << hitFlag << endl;
    fprintf(logFile, "All hits: %i\n", pmt_hits);
    fprintf(logFile, "Passed hits: %i\t%f\n", pmt_hits_passed_total, hits_passed);

    // 3: cuts and thresholds
    if ( double(CNEXTA)/double(count)*100 < val_exta ) {extaFlag = 1;}
    if ( double(CCHS)/double(pmt_hits)*100 < val_badChan ) {badChanFlag = 1;}
    if ( double(CECA)/double(pmt_hits)*100 < val_badECA ) {badECAFlag = 1;}
    if ( double(CPCA)/double(pmt_hits)*100 < val_badPCA ) {badPCAFlag = 1;}
    if ( double(CXT)/double(pmt_hits)*100 < val_Xtalk ) {XtalkFlag = 1;}
    if ( double(CE)/double(pmt_hits)*100 < val_NotEn ) {NotEnFlag = 1;}
    if ( double(COFF)/double(pmt_hits)*100 < val_offPMT ) {offPMTFlag = 1;}
    if ( double(CCO)/double(pmt_hits)*100 < val_offChan ) {offChanFlag = 1;}
    if ( double(CDAQ)/double(pmt_hits)*100 < val_daqEn ) {daqEnFlag = 1;}
    if ( double(CNN)/double(pmt_hits)*100 < val_Nnorm ) {NnormFlag = 1;}
    if ( double(CMAG)/double(pmt_hits)*100 < val_badPos ) {badPosFlag = 1;}
    if ( double(CTIR)/double(pmt_hits)*100 < val_LPCTIR ) {LPCTIRFlag = 1;}
    if ( double(CPV)/double(pmt_hits)*100 < val_LPCIP ) {LPCIPFlag = 1;}
    if ( double(CRH)/double(pmt_hits)*100 < val_LPCLoc ) {LPCLocFlag = 1;}
    if ( double(CDAV)/double(pmt_hits)*100 < val_LPCPath ) {LPCPathFlag = 1;}
    if ( double(CANG)/double(pmt_hits)*100 < val_Ang ) {angFlag = 1;}
    if ( double(CDIST)/double(pmt_hits)*100 < val_NearRef ) {nearRefFlag = 1;}
    logFile_Flags << "EXTA flag: " << extaFlag << endl;
    logFile_Flags << "Bad channel flag: " << badChanFlag << endl;
    logFile_Flags << "Bad ECA flag: " << badECAFlag << endl;
    logFile_Flags << "BAD PCA flag: " << badPCAFlag << endl;
    logFile_Flags << "X-talk flag: " << XtalkFlag << endl;
    logFile_Flags << "Not enabled flag: " << NotEnFlag << endl;
    logFile_Flags << "Offline PMT flag: " << offPMTFlag << endl;
    logFile_Flags << "Offline channel flag: " << offChanFlag << endl;
    logFile_Flags << "Not DAQ enabled flag: " << daqEnFlag << endl;
    logFile_Flags << "Not normal PMT flag: " << NnormFlag << endl;
    logFile_Flags << "Bad PMT position flag: " << badPosFlag << endl;
    logFile_Flags << "LPC TIR flag: " << LPCTIRFlag << endl;
    logFile_Flags << "LPC invalid path flag: " << LPCIPFlag << endl;
    logFile_Flags << "LPC locality flag: " << LPCLocFlag << endl;
    logFile_Flags << "LPC weird path (not AV) flag: " << LPCPathFlag << endl;
    logFile_Flags << "Angular cut flag: " << angFlag << endl;
    logFile_Flags << "Near reflection flag: " << nearRefFlag << endl;
    fprintf(logFile, "CNEXTA: %i\t%f\n", CNEXTA, double(CNEXTA)/double(count)*100);
    fprintf(logFile, "CCHS: %i\t%f\n", CCHS, double(CCHS)/double(pmt_hits)*100);
    fprintf(logFile, "CECA: %i\t%f\n", CECA, double(CECA)/double(pmt_hits)*100);
    fprintf(logFile, "CPCA: %i\t%f\n", CPCA, double(CPCA)/double(pmt_hits)*100);
    fprintf(logFile, "CXT: %i\t%f\n", CXT, double(CXT)/double(pmt_hits)*100);
    fprintf(logFile, "CE: %i\t%f\n", CE, double(CE)/double(pmt_hits)*100);
    fprintf(logFile, "COFF: %i\t%f\n", COFF, double(COFF)/double(pmt_hits)*100);
    fprintf(logFile, "CCO: %i\t%f\n", CCO, double(CCO)/double(pmt_hits)*100);
    fprintf(logFile, "CDAQ: %i\t%f\n", CDAQ, double(CDAQ)/double(pmt_hits)*100);
    fprintf(logFile, "CNN: %i\t%f\n", CNN, double(CNN)/double(pmt_hits)*100);
    fprintf(logFile, "CMAG: %i\t%f\n", CMAG, double(CMAG)/double(pmt_hits)*100);
    fprintf(logFile, "CTIR: %i\t%f\n", CTIR, double(CTIR)/double(pmt_hits)*100);
    fprintf(logFile, "CPV: %i\t%f\n", CPV, double(CPV)/double(pmt_hits)*100);
    fprintf(logFile, "CRH: %i\t%f\n", CRH, double(CRH)/double(pmt_hits)*100);
    fprintf(logFile, "CDAV: %i\t%f\n", CDAV, double(CDAV)/double(pmt_hits)*100);
    fprintf(logFile, "CANG: %i\t%f\n", CANG, double(CANG)/double(pmt_hits)*100);
    fprintf(logFile, "CDIST: %i\t%f\n", CDIST, double(CDIST)/double(pmt_hits)*100);

    // 4: NHit over time - rolling average deviation
    size_t numPoints = NHits_cleaned.size();
    double *xarray = new double[numPoints];
    double *yarray = new double[numPoints];
    for (size_t i=0; i<numPoints; i++) {
      xarray[i] = i;
      yarray[i] = NHits_cleaned[i];
    }
    TGraph *nhitGraph = new TGraph(numPoints, xarray, yarray);
    delete[] xarray;
    delete[] yarray;
    TGraph* smooth = new TGraph();
    Smooth(smooth,nhitGraph,1000); // floating mean
    double meanS = smooth->GetMean(2);
    double *valS = smooth->GetY();
    double minS = 1e6;
    double maxS = -1e6;
    for (int i=0; i<smooth->GetN(); i++) {
      if (valS[i] > maxS) maxS = valS[i];
      if (valS[i] < minS) minS = valS[i];
    }
    double TOL = val_stableNHit;      // tolerance for fluctuations
    if (maxS/meanS > 1.+TOL || minS/meanS < 1.-TOL) {
      stableNHitFlag = 0;
    } else { stableNHitFlag = 1; }
    logFile_Flags << "Stable NHit flag: " << stableNHitFlag << endl;
    fprintf(logFile, "Rolling NHit mean: %f\n", meanS);
    fprintf(logFile, "Rolling NHit min: %f\n", minS);
    fprintf(logFile, "Rolling NHit max: %f\n", maxS);

    // 5: Cleaned NHit distribution: mean & rms
    double NHit_mean = fCNHitHist->GetMean();
    double NHit_rms = fCNHitHist->GetRMS();
    if ( NHit_mean > (val_NHitMean-val_NHitMeanDev) & NHit_mean < (val_NHitMean+val_NHitMeanDev) & NHit_rms < val_NHitRMS) {NHitDistFlag = 1;}
    logFile_Flags << "NHit distribution flag: " << NHitDistFlag << endl;
    fprintf(logFile, "NHit distribution, mean: %f\n", NHit_mean);
    fprintf(logFile, "NHit distribution, rms: %f\n", NHit_rms);

    // 6: NHit over subruns
    for (size_t i=0; i<avgNHits.size(); i++) {
      if ( avgNHits[i] > (val_subNHit-val_subNHitDev) & avgNHits[i] < (val_subNHit+val_subNHitDev) ) { NHitSubsFlag += 1; }
    }
    if ( NHitSubsFlag == val_nsubs ){ NHitSubsFlag = 1; }
    else { NHitSubsFlag = 0; }
    logFile_Flags << "NHit distribution over subruns flag: " << NHitSubsFlag << endl;

    // 7: software delays
    double sum = std::accumulate( fibre_delays.begin(), fibre_delays.end(), 0.0 );
    double mean = sum / fibre_delays.size();
    double dev = 0;
    for (size_t i=0; i<fibre_delays.size(); i++) {
      dev += abs(fibre_delays[i] - mean);
    }
    if (dev == 0){ FibDelFlag = 1; }
    sum = std::accumulate( trigger_delays.begin(), trigger_delays.end(), 0.0 );
    mean = sum / trigger_delays.size();
    dev = 0;
    for (size_t i=0; i<trigger_delays.size(); i++) {
      dev += abs(trigger_delays[i] - mean);
    }
    if (dev == 0){ TrigDelFlag = 1; }
    logFile_Flags << "Trigger delay dev flag: " << TrigDelFlag << endl;
    logFile_Flags << "Fibre delay dev flag: " << FibDelFlag << endl;

    // 8: direct hits - time
    double HitsPeak = fDirectHits->GetMean();
    if (HitsPeak > (val_hitPeak-val_hitPeakDev) & HitsPeak < (val_hitPeak+val_hitPeakDev)){ hitPeakFlag = 1; }
    logFile_Flags << "Hit peak flag: " << hitPeakFlag << endl;
    fprintf(logFile, "Hit peak: %f\n", HitsPeak);

    // 9: direct hits - stability over time
    int xBins = fTimevsEvent->GetXaxis()->GetNbins();
    int yBins = fTimevsEvent->GetYaxis()->GetNbins();
    int timeCheck = 0;
    std::vector<double> colTimes;

    for (size_t i=1; i<=xBins; i++) {
      int countColEls = 0;
      double sumCol = 0;
      for (size_t j=1; j<=yBins; j++) {
        sumCol += fTimevsEvent->GetBinContent(i,j)*fTimevsEvent->GetYaxis()->GetBinCenter(j);
        countColEls += fTimevsEvent->GetBinContent(i,j);
      }
      double meanCol = sumCol / countColEls;
      colTimes.push_back( meanCol );
      fTimevsEventFit->SetPoint(i, fTimevsEvent->GetXaxis()->GetBinCenter(i), meanCol);
    }
    double timeMean = std::accumulate( colTimes.begin(), colTimes.end(), 0.0 ) / colTimes.size();
    fprintf(logFile, "Mean of the time over event distribution: %f\n", timeMean);
    double tempMean = 0;
    double x,y;
    for (size_t i=1; i<=fTimevsEventFit->GetN(); i++) {
      fTimevsEventFit->GetPoint(i, x, y);
      tempMean = y;
      if (tempMean > (timeMean-val_timeDev) & tempMean < (timeMean+val_timeDev)){ timeRollFlag += 1; }
      else { cout << "BAD TIME: " << tempMean << endl;
             fprintf(logFile, "BAD TIME: %i\t%f\n", i, timeMean);}
    }
    cout << "fit bins: " << fTimevsEventFit->GetN() << endl;
    cout << "roll checks: " << timeRollFlag << endl;
    if ( timeRollFlag > (fTimevsEventFit->GetN() - 10) ) { timeRollFlag = 1; }
    else { timeRollFlag = 0; }
    logFile_Flags << "Time over event flag: " << timeRollFlag << endl;

    // 10: direct hits - # peaks
    TSpectrum *s = new TSpectrum();
    Int_t nfound = s->Search(fDirectHits);
    double peakTime = *s->GetPositionX();
    if (peakTime > (val_hitPeak-val_hitPeakDev) & peakTime < (val_hitPeak+val_hitPeakDev) & nfound==1){ peakFlag = 1; }
    logFile_Flags << "Peak check flag: " << peakFlag << endl;
    fprintf(logFile, "Direct hits, peaks: %i\n", nfound);
    fprintf(logFile, "Direct hits, peak time: %f\n", peakTime);

    // 11: # PMTs in beamspot
    pmtsBSP = double(pmtsBS)/double(NPMTS)*100;
    pmtsBSgoodOccupP = double(pmtsBSgoodOccup)/double(NPMTS)*100;
    double goodPMTSP = double(pmtsBSgoodOccup)/double(pmtsBS)*100;
    if ( (pmtsBS >= val_pmtsBS) & (pmtsBSP >= val_pmtsBSP) & (pmtsBSgoodOccup >= val_pmtsBSGO) & (pmtsBSgoodOccupP >= val_pmtsBSGOP) & (goodPMTSP >= val_pmtsBSGO_ratio) ){ PMTsFlag = 1; }
    logFile_Flags << "Check on PMTs (beamspot / occupancy) flag: " << PMTsFlag << endl;
    fprintf(logFile, "PMTs in BS: %i\t%f\n", pmtsBS, pmtsBSP);
    fprintf(logFile, "PMTs in BS, good occup: %i\t%f\n", pmtsBSgoodOccup, pmtsBSgoodOccupP);
    fprintf(logFile, "Good PMTs: %f\n", double(goodPMTSP));

    // 12: run length
    EXTAlength = (lastEXTA - firstEXTA);
    if ( EXTAlength <= val_runtime ){ runtimeFlag = 1; }
    logFile_Flags << "Runtime flag: " << runtimeFlag << endl;
    fprintf(logFile, "EXTAlength: %f\n", double(EXTAlength));

    // 13: events in subruns
    for (size_t i=0; i<n_subruns; i++) {
      if ( EXTAtimes10[i].size() >= (val_evsSub-val_evsSubDev*val_evsSub) & EXTAtimes10[i].size() <= (val_evsSub+val_evsSubDev*val_evsSub) ){ evsSubFlag += 1; }
    }
    if ( evsSubFlag == val_nsubs ){ evsSubFlag = 1; }
    logFile_Flags << "Check on evens in subruns: " << evsSubFlag << endl;

    // 14: frequency - from length
    int totalEXTA = 0;
    for (size_t i=0; i<EXTAtimes10.size(); i++) {
      totalEXTA += EXTAtimes10[i].size();
    }
    long double freqLenght = double(totalEXTA) / double(EXTAlength);
    if ( freqLenght >= (val_freq-val_freqDev*val_freq) & freqLenght <= (val_freq+val_freqDev*val_freq) ){ freqRuntimeFlag = 1; }
    logFile_Flags << "Frequency flag: " << freqRuntimeFlag << endl;
    fprintf(logFile, "Frequency from length: %f\n", double(freqLenght));

    // 15: frequency - from delay between events
    cout.precision(10);
    for (size_t i=0; i<EXTAtimes10.size(); i++) {
      for (size_t j=0; j<EXTAtimes10[i].size(); j++) {
        fEXTADelta10->Fill( EXTAtimes10[i][j]/double(1e7)*1000. );
        fEXTADelta50->Fill( EXTAtimes50[i][j]/double(5e7)*1000. );
        // stolen trigs here
        if ( EXTAtimes10[i][j]/double(1e7)*1000. >= 2*((val_freq-val_freqDev*val_freq)/1000) & EXTAtimes10[i][j]/double(1e7)*1000. <= 2*((val_freq+val_freqDev*val_freq)/1000) ){
           fEXTAStolen->Fill(EXTAtimes10[i][j]/double(1e7)*1000.);
        }
      }
    }

    // 16: Occupancy check: All - integral 1-5
    double integAll = fOccup->Integral( fOccup->FindBin(1), fOccup->FindBin(5) );
    double integAll_all = fOccup->Integral( fOccup->FindBin(0), fOccup->FindBin(fOccup->GetNbinsX()) );
    double integALLR = integAll/integAll_all*100;
    if ( integALLR >= val_integ_all ){ integAllFlag = 1; }
    logFile_Flags << "Occupancy all flag: " << integAllFlag << endl;
    fprintf(logFile, "Occup 1-5 to all: %f\n", integALLR);

    // 17: Occupancy check: Beamspot - integral 1-5
    double integBS = fOccupBeamSpot->Integral( fOccupBeamSpot->FindBin(1), fOccupBeamSpot->FindBin(5) );
    double integBS_all = fOccupBeamSpot->Integral( fOccupBeamSpot->FindBin(0), fOccupBeamSpot->FindBin(fOccupBeamSpot->GetNbinsX()) );
    double integBSR = integBS/integBS_all*100;
    if ( integBSR >= val_integ_bs ){ integBSFlag = 1; }
    logFile_Flags << "Occupancy beamspot flag: " << integBSFlag << endl;
    fprintf(logFile, "Occup BM 1-5 to all: %f\n", integBSR);

    // 18: Occupancy check: Beamspot cut - number of PMTs
    int entriesBSCut = fOccupCut->GetEntries();
    if ( entriesBSCut >= val_bs_pmts ){ BSpmtsFlag = 1; }
    logFile_Flags << "Beamspot PMTs flag: " << BSpmtsFlag << endl;
    fprintf(logFile, "Good PMTs BS 1-5 occ: %i\n", entriesBSCut);

    // 19: PIN RMS
    double pinRMS = fPINHist->GetRMS();
    if ( pinRMS <= val_pin_rms ){ PINFlag = 1; }
    logFile_Flags << "PIN RMS flag: " << PINFlag << endl;
    fprintf(logFile, "PIN RMS: %f\n", pinRMS);

    // 20: PIN vs NHit
    double covariance = fPINvsNHit->GetCovariance();
    double correlationF = fPINvsNHit->GetCorrelationFactor();
    if ( covariance >= val_cov & correlationF >= val_corf ){ PINNHitFlag = 1; }
    logFile_Flags << "PIN-NHit cov and corr flag: " << PINNHitFlag << endl;
    fprintf(logFile, "Covariance: %f\n", covariance);
    fprintf(logFile, "Correlation factor: %f\n", correlationF);
    fprintf(logFile, "\n");

    // Add flags to log file
    fprintf(logFile, logFile_Flags.str().c_str());


    // Store plots
    // ROOT file
    // FileName.str("");
    // FileName << run_id;
    // FileName << ".root";
    // outname = FileName.str();
    // const char * FileNameC = outname.c_str();
    // RootFile = new TFile(FileNameC,"RECREATE");
    // RootFile->cd();

    fOccupC->cd();
    fOccup->SetTitle("fOccup");
    fOccup->GetXaxis()->SetTitle("Occupancy [%]");
    fOccup->GetYaxis()->SetTitle("count");
    fOccup->GetYaxis()->SetTitleOffset(1.3);
    fOccupC->SetLogy();
    fOccup->Draw();
    fOccupC->SaveAs(Form("%i_val1_fOccup.png",run_id));
    fOccupC->Close();

    fOccupBeamSpotC->cd();
    fOccupBeamSpot->SetTitle("fOccupBeamSpot");
    fOccupBeamSpot->GetXaxis()->SetTitle("Occupancy [%]");
    fOccupBeamSpot->GetYaxis()->SetTitle("count");
    fOccupBeamSpot->GetYaxis()->SetTitleOffset(1.3);
    fOccupBeamSpot->Draw();
    fOccupBeamSpotC->SaveAs(Form("%i_val1_fOccupBeamSpot.png",run_id));
    fOccupBeamSpotC->Close();

    fOccupCutC->cd();
    fOccupCut->SetTitle("fOccupCut");
    fOccupCut->GetXaxis()->SetTitle("Occupancy [%]");
    fOccupCut->GetYaxis()->SetTitle("count");
    fOccupCut->GetYaxis()->SetTitleOffset(1.3);
    fOccupCut->Draw();
    fOccupCutC->SaveAs(Form("%i_val1_fOccupCut.png",run_id));
    fOccupCutC->Close();

    fDirectHitsC->cd();
    fDirectHits->SetTitle("fDirectHits");
    fDirectHits->GetXaxis()->SetTitle("hit time [ns]");
    fDirectHits->GetYaxis()->SetTitle("count");
    fDirectHits->GetYaxis()->SetTitleOffset(1.3);
    fDirectHits->Draw();
    fDirectHitsC->SetLogy();
    fDirectHitsC->SaveAs(Form("%i_val1_fDirectHits.png",run_id));
    fDirectHitsC->Close();

    fNHitHistC->cd();
    fNHitHist->SetTitle("fNHitHist");
    fNHitHist->GetXaxis()->SetTitle("NHit");
    fNHitHist->GetYaxis()->SetTitle("count");
    fNHitHist->GetYaxis()->SetTitleOffset(1.3);
    fNHitHist->SetLineColorAlpha(kBlue, 0.2);
    fCNHitHist->SetLineColorAlpha(kRed, 0.2);
    fSNHitHists->Add(fNHitHist);
    fSNHitHists->Add(fCNHitHist);
    fSNHitHists->Draw("nostack");
    fSNHitHists->SetTitle("NHit distribution");
    fSNHitHists->GetXaxis()->SetTitle("NHit");
    fSNHitHists->GetYaxis()->SetTitle("count");
    fSNHitHists->GetYaxis()->SetTitleOffset(1.2);
    TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
    legend->AddEntry(fNHitHist, "Raw NHit");
    legend->AddEntry(fCNHitHist, "Cleaned NHit");
    legend->SetTextSize(0.04);
    legend->Draw("same");
    fNHitHistC->SaveAs(Form("%i_val1_fNHitHist.png",run_id));
    fNHitHistC->Close();

    fNHitvsEventC->cd();
    fNHitvsEvent->SetTitle("fNHitvsEvent");
    fNHitvsEvent->GetXaxis()->SetTitle("event");
    fNHitvsEvent->GetYaxis()->SetTitle("NHit");
    fNHitvsEvent->GetYaxis()->SetTitleOffset(1.3);
    fNHitvsEvent->Fit("pol0");
    fNHitvsEvent->GetXaxis()->SetRangeUser(0,200000);
    fNHitvsEvent->Draw("PZCOL");
    smooth->SetLineWidth(2);
    //smooth->SetLineColor(2);
    smooth->Draw("same");
    fNHitvsEventC->SaveAs(Form("%i_val1_fNHitvsEvent.png",run_id));
    fNHitvsEventC->Close();

    fNHitvsEventZOOMC->cd();
    fNHitvsEventZOOM->SetTitle("fNHitvsEventZOOM");
    fNHitvsEventZOOM->GetXaxis()->SetTitle("event");
    fNHitvsEventZOOM->GetYaxis()->SetTitle("NHit");
    fNHitvsEventZOOM->GetYaxis()->SetTitleOffset(1.3);
    fNHitvsEventZOOM->SetMarkerStyle(5);
    fNHitvsEventZOOM->SetMarkerColorAlpha(kBlack, 0.8);
    //fNHitvsEventZOOM->Fit("pol0");
    //fNHitvsEventZOOM->GetXaxis()->SetRangeUser(0,3000);
    fNHitvsEventZOOM->Draw("AP");
    fNHitvsEventZOOMC->SaveAs(Form("%i_val1_fNHitvsEventZOOM.png",run_id));
    fNHitvsEventZOOMC->Close();

    fTimevsAngleC->cd();
    fTimevsAngle->SetTitle("fTimevsAngle");
    fTimevsAngle->GetXaxis()->SetTitle("Angle [deg]");
    fTimevsAngle->GetYaxis()->SetTitle("Direct hit time [ns]");
    fTimevsAngle->GetYaxis()->SetTitleOffset(1.3);
    fTimevsAngle->Draw("PZCOL");
    fTimevsAngleC->SaveAs(Form("%i_val1_fTimevsAngle.png",run_id));
    fTimevsAngleC->Close();

    fTimevsEventC->cd();
    fTimevsEvent->SetTitle("fTimevsEvent");
    fTimevsEvent->GetXaxis()->SetTitle("event");
    fTimevsEvent->GetYaxis()->SetTitle("Direct hit time [ns]");
    fTimevsEvent->GetYaxis()->SetTitleOffset(1.3);
    fTimevsEvent->GetYaxis()->SetRangeUser(HitsPeak-30, HitsPeak+30);
    fTimevsEvent->Draw("PZCOL");
    fTimevsEventFit->SetLineWidth(2);
    fTimevsEventFit->SetLineColor(kRed);
    fTimevsEventFit->Draw("C same");
    fTimevsEventC->SaveAs(Form("%i_val1_fTimevsEvent.png",run_id));
    fTimevsEventC->Close();

    fFibreDelaysC->cd();
    fFibreDelaysC->SetGrid();
    fFibreDelays->SetTitle("fFibreDelays");
    fFibreDelays->GetXaxis()->SetTitle("subrun #");
    fFibreDelays->GetYaxis()->SetTitle("fibre delay [ns]");
    fFibreDelays->GetYaxis()->SetTitleOffset(1.3);
    fFibreDelays->SetMarkerStyle(5);
    fFibreDelays->Draw("AP");
    fFibreDelays->GetXaxis()->SetRangeUser(-1,val_nsubs);
    fFibreDelays->GetYaxis()->SetRangeUser(*min_element(fibre_delays.begin(), fibre_delays.end())-1,*max_element(fibre_delays.begin(), fibre_delays.end())+1);
    fFibreDelaysC->SaveAs(Form("%i_val1_fFibreDelays.png",run_id));
    fFibreDelaysC->Close();

    fTriggerDelaysC->cd();
    fTriggerDelaysC->SetGrid();
    fTriggerDelays->SetTitle("fTriggerDelays");
    fTriggerDelays->GetXaxis()->SetTitle("subrun #");
    fTriggerDelays->GetYaxis()->SetTitle("trigger delay [ns]");
    fTriggerDelays->GetYaxis()->SetTitleOffset(1.3);
    fTriggerDelays->SetMarkerStyle(5);
    fTriggerDelays->Draw("AP");
    fTriggerDelays->GetXaxis()->SetRangeUser(-1,val_nsubs);
    fTriggerDelays->GetYaxis()->SetRangeUser(*min_element(trigger_delays.begin(), trigger_delays.end())-1,*max_element(trigger_delays.begin(), trigger_delays.end())+1);
    fTriggerDelaysC->SaveAs(Form("%i_val1_fTriggerDelays.png",run_id));
    fTriggerDelaysC->Close();

    fPINHistC->cd();
    fPINHist->SetTitle("fPINHist");
    fPINHist->GetXaxis()->SetTitle("PIN");
    fPINHist->GetYaxis()->SetTitle("count");
    fPINHist->GetYaxis()->SetTitleOffset(1.3);
    fPINHist->Draw();
    fPINHistC->SaveAs(Form("%i_val1_fPINHist.png",run_id));
    fPINHistC->Close();

    fPINC->cd();
    fPINC->SetGrid();
    fPIN->SetTitle("fPIN");
    fPIN->GetXaxis()->SetTitle("subrun #");
    fPIN->GetYaxis()->SetTitle("PIN");
    fPIN->GetYaxis()->SetTitleOffset(1.3);
    fPIN->SetMarkerStyle(5);
    fPIN->GetXaxis()->SetRangeUser(-1,val_nsubs);
    fPIN->Draw("AP");
    fPINC->SaveAs(Form("%i_val1_fPIN.png",run_id));
    fPINC->Close();

    favgNHitHistC->cd();
    favgNHitHist->SetTitle("favgNHitHist");
    favgNHitHist->GetXaxis()->SetTitle("avg NHit");
    favgNHitHist->GetYaxis()->SetTitle("count");
    favgNHitHist->GetYaxis()->SetTitleOffset(1.3);
    favgNHitHist->Draw();
    favgNHitHistC->SaveAs(Form("%i_val1_favgNHitHist.png",run_id));
    favgNHitHistC->Close();

    favgNHitC->cd();
    favgNHitC->SetGrid();
    favgNHit->SetTitle("favg_NHit");
    favgNHit->GetXaxis()->SetTitle("subrun #");
    favgNHit->GetYaxis()->SetTitle("NHit");
    favgNHit->GetYaxis()->SetTitleOffset(1.3);
    favgNHit->SetMarkerStyle(5);
    favgNHit->GetXaxis()->SetRangeUser(-1,val_nsubs);
    favgNHit->Draw("AP");
    favgNHitC->SaveAs(Form("%i_val1_favgNHit.png",run_id));
    favgNHitC->Close();

    fPINvsNHitC->cd();
    fPINvsNHitC->SetGrid();
    fPINvsNHit->SetTitle("fPINvsNHit");
    fPINvsNHit->GetXaxis()->SetTitle("NHit");
    fPINvsNHit->GetYaxis()->SetTitle("PIN");
    fPINvsNHit->GetYaxis()->SetTitleOffset(1.3);
    fPINvsNHit->SetMarkerStyle(5);
    fPINvsNHit->Draw("AP");
    fPINvsNHitC->SaveAs(Form("%i_val1_fPINvsNHit.png",run_id));
    fPINvsNHitC->Close();

    fXYThetaC->cd();
    fXYThetaC->SetMargin(0.05,0.15,0.05,0.05);
    fXYTheta->SetMarkerStyle(8);
    fXYTheta->SetMarkerSize(1.0);
    fXYTheta->GetXaxis()->SetLabelSize(0.03);
    fXYTheta->GetYaxis()->SetLabelSize(0.03);
    fXYTheta->SetTitle("fXYTheta");
    view->RotateView(270, 0);
    fXYTheta->Draw("AZCOLPCOL");
    fXYThetaC->SaveAs(Form("%i_val1_fXYTheta.png",run_id));
    fXYThetaC->Close();

    fXYThetaAtPMTC->cd();
    fXYThetaAtPMTC->SetMargin(0.05,0.15,0.05,0.05);
    fXYThetaAtPMT->SetMarkerStyle(8);
    fXYThetaAtPMT->SetMarkerSize(1.0);
    fXYThetaAtPMT->GetXaxis()->SetLabelSize(0.03);
    fXYThetaAtPMT->GetYaxis()->SetLabelSize(0.03);
    fXYThetaAtPMT->SetTitle("fXYThetaAtPMT");
    view->RotateView(270, 0);
    fXYThetaAtPMT->Draw("AZCOLPCOL");
    fXYThetaAtPMTC->SaveAs(Form("%i_val1_fXYThetaAtPMT.png",run_id));
    fXYThetaAtPMTC->Close();

    fXYOccupAllC->cd();
    fXYOccupAllC->SetMargin(0.05,0.15,0.05,0.05);
    fXYOccupAll->SetMarkerStyle(8);
    fXYOccupAll->SetMarkerSize(1.0);
    fXYOccupAll->GetXaxis()->SetLabelSize(0.03);
    fXYOccupAll->GetYaxis()->SetLabelSize(0.03);
    fXYOccupAll->SetTitle("fXYOccupAll");
    fXYOccupAll->Draw("AZCOLPCOL");
    view->RotateView(270, 0);
    fXYOccupAllC->SaveAs(Form("%i_val1_fXYOccupAll.png",run_id));
    fXYOccupAllC->Close();

    fXYOccupBeamSpotC->cd();
    fXYOccupBeamSpotC->SetMargin(0.05,0.15,0.05,0.05);
    fXYOccupBeamSpot->SetMarkerStyle(8);
    fXYOccupBeamSpot->SetMarkerSize(1.0);
    fXYOccupBeamSpot->GetXaxis()->SetLabelSize(0.03);
    fXYOccupBeamSpot->GetYaxis()->SetLabelSize(0.03);
    fXYOccupBeamSpot->SetTitle("fXYOccupBeamSpot");
    fXYOccupBeamSpot->Draw("AZCOLPCOL");
    view->RotateView(270, 0);
    fXYOccupBeamSpotC->SaveAs(Form("%i_val1_fXYOccupBeamSpot.png",run_id));
    fXYOccupBeamSpotC->Close();

    fEvsSubrunC->cd();
    fEvsSubrun->SetTitle("fEvsSubrun");
    fEvsSubrun->GetXaxis()->SetTitle("Evs in subrun");
    fEvsSubrun->GetYaxis()->SetTitle("count");
    fEvsSubrun->GetYaxis()->SetTitleOffset(1.3);
    fEvsSubrun->Draw();
    fEvsSubrunC->SaveAs(Form("%i_val1_fEvsSubrun.png",run_id));
    fEvsSubrunC->Close();

    fEXTADeltaC->cd();
    fEXTADeltaC->SetLogy();
    fEXTADelta10->SetLineColor(3);
    fEXTADelta50->SetLineColor(4);
    fEXTAStolen->SetLineColor(2);
    fEXTADelta10->SetFillStyle(3145);
    fEXTADelta50->SetFillStyle(3154);
    fEXTAStolen->SetFillStyle(3305);
    fEXTADelta10->SetFillColorAlpha(3, 0.8);
    fEXTADelta50->SetFillColorAlpha(4, 0.8);
    fEXTAStolen->SetFillColorAlpha(2, 0.8);
    fEXTADelta->Add(fEXTADelta10);
    fEXTADelta->Add(fEXTADelta50);
    fEXTADelta->Add(fEXTAStolen);
    fEXTADelta->Draw("nostack");
    fEXTADelta->SetTitle("fEXTADelta");
    fEXTADelta->GetXaxis()->SetTitle("Event separation [ms]");
    fEXTADelta->GetYaxis()->SetTitle("count");
    fEXTADelta->GetYaxis()->SetTitleOffset(1.2);
    TLegend* legend2 = new TLegend(0.65,0.7,0.9,0.9);
    legend2->AddEntry(fEXTADelta10, "10 MHz clock");
    legend2->AddEntry(fEXTADelta50, "50 MHz clock");
    legend2->AddEntry(fEXTAStolen, "Stolen evs");
    legend2->SetTextSize(0.04);
    legend2->Draw("same");
    fEXTADeltaC->SaveAs(Form("%i_val1_fEXTADelta.png",run_id));
    fEXTADeltaC->Close();

    //RootFile->Close();

    // Delete
    delete fNHitvsEvent;
    delete fTimevsEvent;
    delete fOccup;
    delete fOccupBeamSpot;
    delete fOccupCut;
    delete fDirectHits;
    delete fNHitHist;
    delete fCNHitHist;
    delete fTimevsAngle;
    delete fPINHist;
    delete fEvsSubrun;
    delete fEXTADelta10;
    delete fEXTADelta50;
    delete fEXTAStolen;
    delete favgNHitHist;
    fibre_delays.clear();
    trigger_delays.clear();
    EXTAtimes10.clear();
    EXTAtimes50.clear();
    NHits.clear();
    NHits_cleaned.clear();
    subruns.clear();
    PINs.clear();
    PIN_errs.clear();
    x_errors.clear();
    avgNHits.clear();
    avgNHitsE.clear();
    PIN_RMSes.clear();
    fPMTs.clear();
    delete smooth; delete nhitGraph;
    delete fPIN;
    delete fPINvsNHit;

  }

  Processor::Result MyUserProc::DSEvent(DS::Run& run, DS::Entry& ds) {

    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {
      const RAT::DS::EV& ev = ds.GetEV(iEv);

      count++;

      // Trigger type
      trig = ev.GetTrigType();
      if (!(trig & 0x8000)){ CNEXTA++; continue; }    // EXT trigger only

      allEvs++;

      int subrun = ds.GetSubRunID();

      // fix for old tellie runs 
      if (allEvs == 1) {
        if (subrun == 1) {
          flag = 1;
          cout << "old tellie run, applying fix for subruns..." << endl;
        }
        else if (subrun == 0) {
          flag = 0;
        }
      }
      if (flag == 1){ subrun = subrun - 1; }
      //
      if (subrun == 40){ continue; } // this should not happen, but sometimes does ?!

      if (subrun != prevSubrun){ firstEXTAinSubrun = 1; }
      prevSubrun = subrun;

      avgNHit[subrun][0][0]++;

      pmt_hits_passed = 0;
      const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();   // calib PMTs
      for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
        pmt_hits++;
        RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
        int pmtID = pmt.GetID();

        const RAT::DU::PMTCalStatus& pmtStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
        const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
        const RAT::DU::PMTInfo& pmtinfo_loop = RAT::DU::Utility::Get()->GetPMTInfo();
        unsigned int status = pmtStatus.GetHitStatus(pmt);
        if(status & (1<<pmtStatus.kCHSBit)){ CCHS++; continue; }
        if(status & (1<<pmtStatus.kECABit)){ CECA++; continue; }
        if(status & (1<<pmtStatus.kPCABit)){ CPCA++; continue; }
        if(status & (1<<pmtStatus.kXTalkBit)){ CXT++; continue; }
        if ( !chs.IsTubeOnline(pmtID) ){ COFF++; continue; }
        if ( !chs.IsEnabled() ){ CE++; continue; }
        if ( !chs.IsChannelOnline(pmtID) ){ CCO++; continue; }
        if ( !chs.IsDAQEnabled(pmtID) ){ CDAQ++; continue; }
        if ( pmtinfo_loop.GetType(pmtID) != 1 ){ CNN++; continue; }

        // Get info for this PMT
        const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
        TVector3 pmtPos = pmtinfo.GetPosition(pmtID);       // position [mm]
        TVector3 pmtDir = pmtinfo.GetDirection(pmtID);      // direction
        if (pmtPos.Mag()==0) { CMAG++; continue; }

        double pmtTime = pmt.GetTime();                     // hit time [ns]
        if (pmtTime > 0){
          if (pmtTime < time_min){time_min = pmtTime;}
          if (pmtTime > time_max){time_max = pmtTime;}
        }

        // Get light travel time
        lpc.CalcByPosition(fibrepos, pmtPos, energy, LOCALITY);

        // LPC checks
        if (lpc.GetTIR() == 1) { CTIR++; continue; }           // total internal reflection
        if (lpc.GetPathValid() == 0) { CPV++; continue; }      // check whether path parameters are valid
        if (lpc.GetResvHit() == 1) { CRH++; continue; }        // whether end point was within locality

        // Get light emission angle
        TVector3 startDir = lpc.GetInitialLightVec();         // start direction at fibre
        double theta = startDir.Angle(fibredir)*180./pi;      // emission angle at fibre

        if (lpc.GetTotalDist() <= 12000){ CDIST++; continue;}        // this rejects near reflections
        if (lpc.GetDistInInnerAV() <= 7000){ CDAV++; continue;}      // this rejects other weird paths

        // ToF
        double distInInnerAV = lpc.GetDistInInnerAV();
        double distInAV = lpc.GetDistInAV();
        double distInWater = lpc.GetDistInWater();
        double lightTravelTime = gv.CalcByDistance(distInInnerAV, distInAV, distInWater, energy);

        TVector3 endDir = lpc.GetIncidentVecOnPMT();        // end direction at PMT
        double thetaAtPMT = endDir.Angle(pmtDir)*180./pi;   // incident angle with bucket face

        // Angular cut here
        fPMTs[pmtID][7][0]++;
        if (fPMTs[pmtID][5][0] == -5){
          fPMTs[pmtID][5][0] = theta;         // at fibre
          fPMTs[pmtID][6][0] = thetaAtPMT;    // at pmt
        }
        if ( (theta > 12) || (theta < 0) ) { CANG++; continue; }

        // Store hits
        fPMTs[pmtID][3].push_back(pmtTime);
        fTimevsEvent->Fill(allEvs, pmtTime);

        if (fPMTs[pmtID][4][0] == 0){
          fPMTs[pmtID][4][0] = 1;             // is in beam spot
        }

        avgNHit[subrun][1][0]++;

        pmt_hits_passed++;
        pmt_hits_passed_total++;
      } // pmt loop

      if (pmt_hits_passed > 0){
        // Set first EXTA event time here
        if (passedEvs == 0) {
          firstEXTA = ev.GetClockCount10()/10000000;
        }
        lastEXTA = ev.GetClockCount10()/10000000;

        // event dt here
        if (firstEXTAinSubrun != 1) {
          EXTAtimes10[subrun].push_back( ev.GetClockCount10() - prevEXTA10 );
          EXTAtimes50[subrun].push_back( ev.GetClockCount50() - prevEXTA50 );
        }
        prevEXTA10 = ev.GetClockCount10();
        prevEXTA50 = ev.GetClockCount50();
        firstEXTAinSubrun = 0;

        passedEvs++;
        // Get NHit data
        NHit = ev.GetNhits();
        NHit_cleaned = ev.GetNhitsCleaned(); //not doing anything with this yet
        NHits.push_back(NHit);
        NHits_cleaned.push_back(NHit_cleaned);
        fNHitvsEvent->Fill(allEvs, NHit_cleaned);
        avgNHit[subrun][2].push_back( pmt_hits_passed );
        if (allEvs < 2500){
          fNHitvsEventZOOM->SetPoint(allEvs, allEvs, NHit_cleaned);
        }
      }

    } // entries
  }

  /// Display vector as a string
  string MyUserProc::printVector(const TVector3& v) {
    string out;
    if (v.Mag() < 10) out = Form("(%.3f, %.3f, %.3f)", v.X(),  v.Y(), v.Z());
    else              out = Form("(%.1f, %.1f, %.1f)", v.X(),  v.Y(), v.Z());
    return out.c_str();
  }

  void MyUserProc::Smooth(TGraph* g, TGraph* h, int width) {
    int n = h->GetN();
    double* y = h->GetY();
    double sum;
    int cnt, min, max, b;
    for (int i=0; i<n; i++) {
      if (width==0) {
        g->SetPoint(i,i,y[i]);
        continue;
      }
      min = (i<width) ? 0 : i-width;
      max = (i>n-width) ? n : i+width;
      sum = 0;
      cnt = 0;
      for (b=min; b<max; b++) {
        sum += y[b];
        cnt++;
      }
      g->SetPoint(i,i,sum/cnt);
    }
  }
