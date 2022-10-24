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

  virtual double AngularSystematic(double b, double angle);

protected:

  // Fit results
  string fits_folder;
  string fits_file_dir;
  string fits_file_ang;
  string fits_file_offset;
  string delim;
  stringstream dir_fit_file;
  stringstream ang_fit_file;
  stringstream offset_fit_file;
  std::vector<string> fibre_ang_name;
  std::vector<double> fibre_ang_B;
  std::vector<double> fibre_ang_B_err;
  int temp_count;
  std::vector<string> fibre_dir_name;
  std::vector<double> fibre_dir_X;
  std::vector<double> fibre_dir_Y;
  std::vector<double> fibre_dir_Z;
  double AngB;
  std::vector<string> fibre_pcaoffset_name;
  std::vector<double> fibre_pcaoffset;
  std::vector<double> fibre_fibre_pcaoffset_er;

  // LPC
  double energy;
  double fLEDWavelength;
  double LOCALITY;

  // RAT stuff
  RAT::DU::GroupVelocity gv;
  RAT::DU::LightPathCalculator lpc;
  int NPMTS;
  int run_id;

  // Event loop
  int goodHit;
  int trig;
  std::vector<double> t_ToF;
  std::vector<double> t_bucket;
  std::vector<double> t_angsys;
  std::vector<double> dir_hits;
  std::vector<double> resid_hits;

  // Database
  RAT::DBLinkPtr tellie_run_data;
  json::Value sub_run_info;
  std::string fibre_db;
  RAT::DBLinkPtr fibre_data;
  std::vector<int> IPWs;
  std::vector<double> freq;
  TVector3 fibrepos, fibredir, fitfibredir;

  // ROOT
  TFile *RootFile;
  stringstream FileName;
  std::string outname;

  // Validation
  double val_tof_mean;
  double val_tof_rms;
  double val_buck_mean;
  double val_buck_rms;
  double val_as_mean_min;
  double val_as_mean_max;
  double val_as_rms;
  double val_cor_dev;
  double val_tof_max;
  double val_tof_min;
  double val_buc_max;
  double val_buc_min;
  double val_as_max;
  double val_as_min;
  double val_direct;
  double val_direct_dev;
  double val_resid;
  double val_resid_dev;
  double val_DirResid_par0;
  double val_DirResid_par1;
  double val_fitSlope;
  double val_fitSlope_dev;

  int tofFlag;
  int buckFlag;
  int asFlag;
  int tofDistFlag;
  int bucDistFlag;
  int asDistFlag;
  int tofDecreaseFlag;
  int asIncreaseFlag;
  int residHitsFlag;
  int DirResidFlag;
  int DirAngFlag;
  int ResAngFlag;

  TCanvas *fToFHC;
  TCanvas *fBucketHC;
  TCanvas *fAngSysHC;
  TCanvas *fToFC;
  TCanvas *fBucketC;
  TCanvas *fAngSysC;
  TCanvas *fHitTimesC;
  TCanvas *fHitTimesRelC;
  TCanvas *fDirHitTimesAngC;
  TCanvas *fResidHitTimesAngC;

  TH1D *fToFH;
  TH1D *fBucketH;
  TH1D *fAngSysH;
  TGraph *fToF;
  TGraph *fBucket;
  TGraph *fAngSys;
  THStack *fTimes;
  TH1D *fHitTimes;
  TH1D *fResidTimes;
  TH2D *fHitTimesRel;
  TH2D *fHitTimesAng;
  TH2D *fResidTimesAng;

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

    // Load env variables
    LOCALITY = atoi(getenv("LOCALITY"));
    fLEDWavelength = 506.0*1e-6;
    energy = util::WavelengthToEnergy(fLEDWavelength);
    fits_folder = getenv("FITS_FOLDER");
    fits_file_ang = getenv("ANG_FIT_FILE");
    fits_file_dir = getenv("DIR_FIT_FILE");
    fits_file_offset = getenv("OFFSET_FIT_FILE");
    ang_fit_file << fits_file_ang;
    dir_fit_file << fits_file_dir;
    offset_fit_file << fits_file_offset;
    delim = "\t";

    // Validation
    val_tof_mean = atof(getenv("TOF_MEAN"));
    val_tof_rms = atof(getenv("TOF_RMS"));
    val_buck_mean = atof(getenv("BUCK_MEAN"));
    val_buck_rms = atof(getenv("BUCK_RMS"));
    val_as_mean_min = atof(getenv("AS_MEAN_MIN"));
    val_as_mean_max = atof(getenv("AS_MEAN_MAX"));
    val_as_rms = atof(getenv("AS_RMS"));
    val_cor_dev = atof(getenv("COR_DEV"))/100;
    val_tof_max = atof(getenv("TOF_MAX"));
    val_tof_min = atof(getenv("TOF_MIN"));
    val_buc_max = atof(getenv("BUC_MAX"));
    val_buc_min = atof(getenv("BUC_MIN"));
    val_as_max = atof(getenv("AS_MAX"));
    val_as_min = atof(getenv("AS_MIN"));
    val_direct = atof(getenv("HIT_PEAK"));
    val_direct_dev = atof(getenv("HIT_PEAK_DEV"));
    val_resid = atof(getenv("RESID_PEAK"));
    val_resid_dev = atof(getenv("RESID_PEAK_DEV"));
    val_DirResid_par0 = atof(getenv("DIR_RESID_0"));
    val_DirResid_par1 = atof(getenv("DIR_RESID_1"));
    val_fitSlope = atof(getenv("FIT_SLOPE"));
    val_fitSlope_dev = atof(getenv("FIT_SLOPE_DEV"));

    // Load angular fit results
    cout << "Opening ang fit file: " << ang_fit_file.str() << endl;
    std::ifstream fit_file( ang_fit_file.str().c_str() );
    for ( string line; getline(fit_file, line); ){
      char *token = strtok( (char *)line.c_str(), (char *)delim.c_str());
      temp_count = 0;
      while (token != NULL){
        switch (temp_count){
          case 0:
            fibre_ang_name.push_back( token );
            break;
          case 3:
            fibre_ang_B.push_back( atof(token) );
            break;
          case 4:
            fibre_ang_B_err.push_back( atof(token) );
            break;
        }
        token = strtok(NULL, (char *)delim.c_str());
        temp_count += 1;
      }
    }

    // Load direction fit results
    cout << "Opening dir fit file: " << dir_fit_file.str() << endl;
    std::ifstream fit_file2( dir_fit_file.str().c_str() );
    for ( string line; getline(fit_file2, line); ){
      char *token = strtok( (char *)line.c_str(), (char *)delim.c_str());
      temp_count = 0;
      while (token != NULL){
        switch (temp_count){
          case 0:
            fibre_dir_name.push_back( token );
            break;
          case 1:
            fibre_dir_X.push_back( atof(token) );
            break;
          case 2:
            fibre_dir_Y.push_back( atof(token) );
            break;
          case 3:
            fibre_dir_Z.push_back( atof(token) );
            break;
        }
        token = strtok(NULL, (char *)delim.c_str());
        temp_count += 1;
      }
    }

    // Load pca offset fit results
    cout << "Opening pcaoffset fit file: " << offset_fit_file.str() << endl;
    std::ifstream fit_file3( offset_fit_file.str().c_str() );
    for ( string line; getline(fit_file3, line); ){
      char *token = strtok( (char *)line.c_str(), (char *)delim.c_str());
      temp_count = 0;
      while (token != NULL){
        switch (temp_count){
          case 0:
            fibre_pcaoffset_name.push_back( token );
            break;
          case 1:
            fibre_pcaoffset.push_back( atof(token) );
            break;
          case 2:
            fibre_fibre_pcaoffset_er.push_back( atof(token) );
            break;
        }
        token = strtok(NULL, (char *)delim.c_str());
        temp_count += 1;
      }
    }

  }

  MyUserProc::~MyUserProc() {
  }

  void MyUserProc::BeginOfRun( DS::Run& run ) {

    // Initialise counters
    AngB = 0.0;
    goodHit = 0;

    // Validation
    tofFlag = 0;
    buckFlag = 0;
    asFlag = 0;
    tofDistFlag = 0;
    bucDistFlag = 0;
    asDistFlag = 0;
    tofDecreaseFlag = 0;
    asIncreaseFlag = 0;
    residHitsFlag = 0;
    DirResidFlag = 0;
    DirAngFlag = 0;
    ResAngFlag = 0;

    // Initialise RAT vars
    run_id = run.GetRunID();
    gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    gv.BeginOfRun();
    lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lpc.BeginOfRun();
    lpc.SetELLIEEvent(true);
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
    NPMTS = pmtinfo.GetCount();
    cout << "Initialised RAT. Number of PMTs is " << NPMTS << "." << endl;

    // Load run fibre details
    try {
      tellie_run_data = DB::Get()->GetLink("TELLIE_RUN");
      sub_run_info = tellie_run_data->GetJSON("sub_run_info");

      // Access first subtun to get stuff
      fibre_db = sub_run_info[0]["fibre"].getString();
    }
    catch (DBNotFoundError& e) {
      std::cout << "Proc(): Couldn't load TELLIE Run details from ratdb!" << std::endl;
    }

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

    // Set fitted direction
    for (size_t i=0; i<fibre_dir_name.size(); i++) {
      if (fibre_dir_name[i] == fibre_db){
        fitfibredir.SetXYZ( fibre_dir_X[i], fibre_dir_Y[i], fibre_dir_Z[i] );
      }
    }
    // Set ang sys parameter
    for (size_t i=0; i<fibre_ang_name.size(); i++) {
      if (fibre_ang_name[i] == fibre_db){
        AngB = fibre_ang_B[i];
      }
    }
    cout << "Set fibre as: fibre " << fibre_db << ", pos " << printVector(fibrepos) << ", dir " << printVector(fitfibredir) << endl;
    cout << "ang b parameter: " << AngB << endl;

    // Init event loop graphs
    fToF = new TGraph();
    fBucket = new TGraph();
    fAngSys = new TGraph();
    fHitTimesAng = new TH2D ("fHitTimesAng", "", 100, 0, 12.5, 100, val_direct-25, val_direct+25 );
    fResidTimesAng = new TH2D ("fResidTimesAng", "", 100, 0, 12.5, 100, val_resid-25, val_resid+25 );

    // Reset vectors
    t_ToF.clear();
    t_bucket.clear();
    t_angsys.clear();
    dir_hits.clear();
    resid_hits.clear();

    // Output file for stats
    logFile_namess.str("");
    logFile_Flags.str("");
    logFile_namess << run_id << "_val2.log";
    logFile_name = logFile_namess.str();

  }

  void MyUserProc::EndOfRun( DS::Run& run ) {

    // Store plots
    // ROOT file
    // FileName.str("");
    // FileName << run_id;
    // FileName << ".root";
    // outname = FileName.str();
    // const char * FileNameC = outname.c_str();
    // RootFile = new TFile(FileNameC,"RECREATE");

    // ROOT stuff
    fToFHC = new TCanvas ("fToFHC", "fToFHC", 1024, 768);
    fBucketHC = new TCanvas ("fBucketHC", "fBucketHC", 1024, 768);
    fAngSysHC = new TCanvas ("fAngSysHC", "fAngSysHC", 1024, 768);
    fToFC = new TCanvas ("fToFC", "fToFC", 1024, 768);
    fBucketC = new TCanvas ("fBucketC", "fBucketC", 1024, 768);
    fAngSysC = new TCanvas ("fAngSysC", "fAngSysC", 1024, 768);
    fHitTimesC = new TCanvas ("fHitTimesC", "fHitTimesC", 1024, 768);
    fHitTimesRelC = new TCanvas ("fHitTimesRelC", "fHitTimesRelC", 1024, 768);
    fDirHitTimesAngC = new TCanvas ("fDirHitTimesAngC", "fDirHitTimesAngC", 1024, 768);
    fResidHitTimesAngC = new TCanvas ("fResidHitTimesAngC", "fResidHitTimesAngC", 1024, 768);

    fToFH = new TH1D ("fToFH", "", 100, *min_element(t_ToF.begin(), t_ToF.end())-1, *max_element(t_ToF.begin(), t_ToF.end())+1 );
    fBucketH = new TH1D ("fBucketH", "", 100, *min_element(t_bucket.begin(), t_bucket.end())-0.01, *max_element(t_bucket.begin(), t_bucket.end())+0.01 );
    fAngSysH = new TH1D ("fAngSysH", "", 100, *min_element(t_angsys.begin(), t_angsys.end())-0.01, *max_element(t_angsys.begin(), t_angsys.end())+0.01 );
    fTimes = new THStack("fTimes", "");
    fHitTimes = new TH1D ("fHitTimes", "", 100, *min_element(dir_hits.begin(), dir_hits.end())-2, *max_element(dir_hits.begin(), dir_hits.end())+2 );
    fResidTimes = new TH1D ("fResidTimes", "", 100, *min_element(resid_hits.begin(), resid_hits.end())-2, *max_element(resid_hits.begin(), resid_hits.end())+2 );
    fHitTimesRel = new TH2D ("fHitTimesRel", "", 100, val_direct-15, val_direct+15, 100, val_resid-15, val_resid+15 );

    // Prepare plots
    for (size_t i=0; i<t_ToF.size(); i++) {
      fToFH->Fill(t_ToF[i]);
      fBucketH->Fill(t_bucket[i]);
      fAngSysH->Fill(t_angsys[i]);
      fHitTimes->Fill(dir_hits[i]);
      fResidTimes->Fill(resid_hits[i]);
      fHitTimesRel->Fill(dir_hits[i], resid_hits[i]);
    }

    // Validation checks
    logFile = fopen(logFile_name.c_str(), "w");
    fprintf(logFile, "### VALIDATION ###\n");
    logFile_Flags << "### FLAGS ###" << endl;
    // 0: ToF mean and rms
    double tof_mean = fToFH->GetMean();
    double tof_rms = fToFH->GetRMS();
    if( tof_mean >= val_tof_mean-(val_cor_dev*val_tof_mean) & tof_mean <= val_tof_mean+(val_cor_dev*val_tof_mean) & tof_rms <= val_tof_rms )
    { tofFlag = 1; }
    logFile_Flags << "TOF flag: " << tofFlag << endl;
    fprintf(logFile, "TOF mean: %f\n", tof_mean);
    fprintf(logFile, "TOF RMS: %f\n", tof_rms);

    // 1: bucket mean and rms
    double buc_mean = fBucketH->GetMean();
    double buc_rms = fBucketH->GetRMS();
    if( buc_mean >= val_buck_mean-(val_cor_dev*val_buck_mean) & buc_mean <= val_buck_mean+(val_cor_dev*val_buck_mean) & buc_rms <= val_buck_rms )
    { buckFlag = 1; }
    logFile_Flags << "Bucket flag: " << buckFlag << endl;
    fprintf(logFile, "Bucket mean: %f\n", buc_mean);
    fprintf(logFile, "Bucket RMS: %f\n", buc_rms);

    // 2: angsys mean and rms. For ang sys the deviations are nominally bigger. Therefore we allow for a wiiiiide range. Should be tuned later.
    double angsys_mean = fAngSysH->GetMean();
    double angsys_rms = fAngSysH->GetRMS();
    if( angsys_mean >= val_as_mean_min & angsys_mean <= val_as_mean_max & angsys_rms <= val_as_rms ){ asFlag = 1; }
    logFile_Flags << "Ang sys flag: " << asFlag << endl;
    fprintf(logFile, "Ang sys mean: %f\n", angsys_mean);
    fprintf(logFile, "Ang sys RMS: %f\n", angsys_rms);

    // 3: mean and max of ToF
    double TOFmax = TMath::MaxElement( fToF->GetN(), fToF->GetY() );
    double TOFmin = TMath::MinElement( fToF->GetN(), fToF->GetY() );
    if ( TOFmax <=  val_tof_max & TOFmin >= val_tof_min ){ tofDistFlag = 1; }
    logFile_Flags << "TOF min-max flag: " << tofDistFlag << endl;
    fprintf(logFile, "TOF min: %f\n", TOFmin);
    fprintf(logFile, "TOF max: %f\n", TOFmax);

    // 4: mean and max of bucket
    double BUCmax = TMath::MaxElement( fBucket->GetN(), fBucket->GetY() );
    double BUCmin = TMath::MinElement( fBucket->GetN(), fBucket->GetY() );
    if ( BUCmax <=  val_buc_max & BUCmin >= val_buc_min ){ bucDistFlag = 1; }
    logFile_Flags << "Bucket min-max flag: " << bucDistFlag << endl;
    fprintf(logFile, "Bucket min: %f\n", BUCmin);
    fprintf(logFile, "Bucket max: %f\n", BUCmax);

    // 5: mean and max of ang sys
    double ASmax = TMath::MaxElement( fAngSys->GetN(), fAngSys->GetY() );
    double ASmin = TMath::MinElement( fAngSys->GetN(), fAngSys->GetY() );
    if ( ASmax <=  val_as_max & ASmin >= val_as_min ){ asDistFlag = 1; }
    logFile_Flags << "Ang sys min-max flag: " << asDistFlag << endl;
    fprintf(logFile, "Ang sys min: %f\n", ASmin);
    fprintf(logFile, "Ang sys max: %f\n", ASmax);

    // 6: TOF evaluate
    for (size_t i=1; i<12; i+=2) {
      double now = fToF->Eval(i);
      double prev = fToF->Eval(i-1);
      if (i == 1){
        if (now <= prev+1.5) { tofDecreaseFlag += 1; }
      } else {
        if (now <= prev+0.5) { tofDecreaseFlag += 1; }
      }
      double now2 = fAngSys->Eval(i);
      double prev2 = fAngSys->Eval(i-1);
      if (now2 >= prev2) { asIncreaseFlag += 1; }
    }
    if ( tofDecreaseFlag == 6 ){ tofDecreaseFlag = 1; }
    else { tofDecreaseFlag = 0; }
    logFile_Flags << "TOF progression flag: " << tofDecreaseFlag << endl;

    // 7: AS evaluate
    if ( asIncreaseFlag == 6 ){ asIncreaseFlag = 1; }
    else { asIncreaseFlag = 0; }
    logFile_Flags << "Ang sys progression flag: " << asIncreaseFlag << endl;

    // 8: hit time residuals - peak
    TSpectrum *s = new TSpectrum();
    Int_t nfound = s->Search(fHitTimes);
    double peakTime = *s->GetPositionX();
    TSpectrum *sR = new TSpectrum();
    Int_t nfoundR = sR->Search(fResidTimes);
    double peakTimeR = *sR->GetPositionX();
    if ( nfound == 1 & nfoundR == 1 & peakTimeR >= val_resid-val_resid_dev & peakTimeR <= val_resid+val_resid_dev ){ residHitsFlag = 1; }
    logFile_Flags << "Resid peak flag: " << residHitsFlag << endl;
    fprintf(logFile, "Direct peaks: %i\n", nfound);
    fprintf(logFile, "Direct peaktime: %f\n", peakTime);
    fprintf(logFile, "Residual peaks: %i\n", nfoundR);
    fprintf(logFile, "Residual peaktime: %f\n", peakTimeR);

    // 9: direct vs residual times
    TF1 *fit = new TF1("fit", "pol1");
    fHitTimesRel->Fit("fit");
    TF1 *fitinfo = fHitTimesRel->GetFunction("fit");
    double fit_el1 = fitinfo->GetParameter(0);
    double fit_el2 = fitinfo->GetParameter(1);
    if ( fit_el1 <= val_DirResid_par0 & fit_el2 <= val_DirResid_par1 ){ DirResidFlag = 1; }
    logFile_Flags << "Dir-Resid fit flag: " << DirResidFlag << endl;
    fprintf(logFile, "Dir-Resid fit: %f\t%f\n", fit_el1, fit_el2);

    // 10: Direct hit times f angle
    TF1 *fit2 = new TF1("fit2", "pol1");
    fHitTimesAng->Fit("fit2");
    TF1 *fitinfo2 = fHitTimesAng->GetFunction("fit2");
    double fit2_el1 = fitinfo2->GetParameter(0);
    double fit2_el2 = fitinfo2->GetParameter(1);
    if ( fit2_el1 >= val_direct-val_direct_dev & fit2_el1 <= val_direct+val_direct_dev & abs(fit2_el2) > val_fitSlope-val_fitSlope_dev & abs(fit2_el2) < val_fitSlope+val_fitSlope_dev ){ DirAngFlag = 1; }
    logFile_Flags << "Direct vs Angle fit flag: " << DirAngFlag << endl;
    fprintf(logFile, "Direct vs Angle fit: %f\t%f\n", fit2_el1, fit2_el2);

    // 11: Residual hit times f angle
    TF1 *fit3 = new TF1("fit3", "pol1");
    fResidTimesAng->Fit("fit3");
    TF1 *fitinfo3 = fResidTimesAng->GetFunction("fit3");
    double fit3_el1 = fitinfo3->GetParameter(0);
    double fit3_el2 = fitinfo3->GetParameter(1);
    if ( fit3_el1 >= val_resid-val_resid_dev & fit3_el1 <= val_resid+val_resid_dev & abs(fit3_el2) > val_fitSlope-val_fitSlope_dev & abs(fit3_el2) < val_fitSlope+val_fitSlope_dev ){ ResAngFlag = 1; }
    logFile_Flags << "Residual vs Angle fit flag: " << ResAngFlag << endl;
    fprintf(logFile, "Residual vs Angle fit: %f\t%f\n", fit3_el1, fit3_el2);

    // Add flags to log file
    fprintf(logFile, "\n");
    fprintf(logFile, logFile_Flags.str().c_str());
    fclose(logFile);

    // Plot
    //RootFile->cd();

    fToFHC->cd();
    fToFH->SetTitle("fToFH");
    fToFH->GetXaxis()->SetTitle("Time of Flight [ns]");
    fToFH->GetYaxis()->SetTitle("count");
    fToFH->GetYaxis()->SetTitleOffset(1.3);
    fToFH->Draw();
    fToFHC->SaveAs(Form("%i_val2_fToFH.png",run_id));
    fToFHC->Close();

    fBucketHC->cd();
    fBucketH->SetTitle("fBucketH");
    fBucketH->GetXaxis()->SetTitle("Bucket time [ns]");
    fBucketH->GetYaxis()->SetTitle("count");
    fBucketH->GetYaxis()->SetTitleOffset(1.3);
    fBucketH->Draw();
    fBucketHC->SaveAs(Form("%i_val2_fBucketH.png",run_id));
    fBucketHC->Close();

    fAngSysHC->cd();
    fAngSysH->SetTitle("fAngSysH");
    fAngSysH->GetXaxis()->SetTitle("Angular systematic correction [ns]");
    fAngSysH->GetYaxis()->SetTitle("count");
    fAngSysH->GetYaxis()->SetTitleOffset(1.3);
    fAngSysH->Draw();
    fAngSysHC->SaveAs(Form("%i_val2_fAngSysH.png",run_id));
    fAngSysHC->Close();

    fToFC->cd();
    fToF->SetTitle("fToF");
    fToF->GetXaxis()->SetTitle("angle");
    fToF->GetYaxis()->SetTitle("Time of Flight [ns]");
    fToF->GetYaxis()->SetTitleOffset(1.3);
    fToF->SetMarkerStyle(5);
    fToF->SetMarkerColorAlpha(kBlack, 0.8);
    fToF->Draw("AP");
    fToFC->SaveAs(Form("%i_val2_fToF.png",run_id));
    fToFC->Close();

    fBucketC->cd();
    fBucket->SetTitle("fBucket");
    fBucket->GetXaxis()->SetTitle("angle");
    fBucket->GetYaxis()->SetTitle("Bucket time [ns]");
    fBucket->GetYaxis()->SetTitleOffset(1.3);
    fBucket->SetMarkerStyle(5);
    fBucket->SetMarkerColorAlpha(kBlack, 0.8);
    fBucket->Draw("AP");
    fBucketC->SaveAs(Form("%i_val2_fBucket.png",run_id));
    fBucketC->Close();

    fAngSysC->cd();
    fAngSys->SetTitle("fAngSys");
    fAngSys->GetXaxis()->SetTitle("angle");
    fAngSys->GetYaxis()->SetTitle("Angular systematic correction [ns]");
    fAngSys->GetYaxis()->SetTitleOffset(1.3);
    fAngSys->SetMarkerStyle(5);
    fAngSys->SetMarkerColorAlpha(kBlack, 0.8);
    fAngSys->Draw("AP");
    fAngSysC->SaveAs(Form("%i_val2_fAngSys.png",run_id));
    fAngSysC->Close();

    fHitTimesC->cd();
    fHitTimes->SetLineColorAlpha(kBlue, 0.2);
    fResidTimes->SetLineColorAlpha(kRed, 0.2);
    fTimes->Add(fHitTimes);
    fTimes->Add(fResidTimes);
    fTimes->Draw("nostack");
    fTimes->SetTitle("NHit distribution");
    fTimes->GetXaxis()->SetTitle("NHit");
    fTimes->GetYaxis()->SetTitle("count");
    fTimes->GetYaxis()->SetTitleOffset(1.2);
    TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
    legend->AddEntry(fHitTimes, "Direct hits");
    legend->AddEntry(fResidTimes, "Residual hits");
    legend->SetTextSize(0.04);
    legend->Draw("same");
    fHitTimesC->SetLogy();
    fHitTimesC->SaveAs(Form("%i_val2_fTimes.png",run_id));
    fHitTimesC->Close();

    fHitTimesRelC->cd();
    fHitTimesRel->SetTitle("fHitTimesRel");
    fHitTimesRel->GetXaxis()->SetTitle("Direct time [ns]");
    fHitTimesRel->GetYaxis()->SetTitle("Residual time [ns]");
    fHitTimesRel->GetYaxis()->SetTitleOffset(1.3);
    fHitTimesRel->Draw("PZCOL");
    fHitTimesRelC->SaveAs(Form("%i_val2_fHitTimesRel.png",run_id));
    fHitTimesRelC->Close();

    fDirHitTimesAngC->cd();
    fHitTimesAng->SetTitle("fHitTimesAng");
    fHitTimesAng->GetXaxis()->SetTitle("angle [deg]");
    fHitTimesAng->GetYaxis()->SetTitle("Direct time [ns]");
    fHitTimesAng->GetYaxis()->SetTitleOffset(1.3);
    fHitTimesAng->Draw("PZCOL");
    fDirHitTimesAngC->SaveAs(Form("%i_val2_fHitTimesAng.png",run_id));
    fDirHitTimesAngC->Close();

    fResidHitTimesAngC->cd();
    fResidTimesAng->SetTitle("fResidTimesAng");
    fResidTimesAng->GetXaxis()->SetTitle("angle [deg]");
    fResidTimesAng->GetYaxis()->SetTitle("Residual time [ns]");
    fResidTimesAng->GetYaxis()->SetTitleOffset(1.3);
    fResidTimesAng->Draw("PZCOL");
    fResidHitTimesAngC->SaveAs(Form("%i_val2_fResidTimesAng.png",run_id));
    fResidHitTimesAngC->Close();

    //RootFile->Close();

    delete fHitTimes; delete fResidTimes; delete fHitTimesRel;
    delete fHitTimesAng; delete fResidTimesAng;
    delete fToFH; delete fBucketH; delete fAngSysH;

  }

  Processor::Result MyUserProc::DSEvent(DS::Run& run, DS::Entry& ds) {

    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {
      const RAT::DS::EV& ev = ds.GetEV(iEv);

      // Trigger type
      trig = ev.GetTrigType();
      if (!(trig & 0x8000)){ continue; }    // EXT trigger only

      const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();   // calib PMTs
      for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
        RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
        int pmtID = pmt.GetID();

        const RAT::DU::PMTCalStatus& pmtStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
        const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
        const RAT::DU::PMTInfo& pmtinfo_loop = RAT::DU::Utility::Get()->GetPMTInfo();
        unsigned int status = pmtStatus.GetHitStatus(pmt);
        if(status & (1<<pmtStatus.kCHSBit)){ continue; }
        if(status & (1<<pmtStatus.kECABit)){ continue; }
        if(status & (1<<pmtStatus.kPCABit)){ continue; }
        if(status & (1<<pmtStatus.kXTalkBit)){ continue; }
        if ( !chs.IsTubeOnline(pmtID) ){ continue; }
        if ( !chs.IsEnabled() ){ continue; }
        if ( !chs.IsChannelOnline(pmtID) ){ continue; }
        if ( !chs.IsDAQEnabled(pmtID) ){ continue; }
        if ( pmtinfo_loop.GetType(pmtID) != 1 ){ continue; }

        // Get info for this PMT
        const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
        TVector3 pmtPos = pmtinfo.GetPosition(pmtID);       // position [mm]
        TVector3 pmtDir = pmtinfo.GetDirection(pmtID);      // direction
        if (pmtPos.Mag()==0) { continue; }

        double pmtTime = pmt.GetTime();                     // hit time [ns]

        // Get light travel time
        //lpc.CalcByPosition(fibrepos, pmtPos, energy, LOCALITY);
        lpc.CalcByPositionPartial( fibrepos, pmtPos, energy, LOCALITY ); // partial

        // LPC checks
        if (lpc.GetTIR() == 1) { continue; }            // total internal reflection
        if (lpc.GetPathValid() == 0) { continue; }      // check whether path parameters are valid
        if (lpc.GetResvHit() == 1) { continue; }        // whether end point was within locality

        // Get light emission angle
        TVector3 startDir = lpc.GetInitialLightVec();         // start direction at fibre
        double theta = startDir.Angle(fitfibredir)*180./pi;      // emission angle at fibre

        // Angular cut here
        if ( (theta > 12) || (theta < 0) ) { continue; }

        if (lpc.GetTotalDist() <= 12000){ continue;}        // this rejects near reflections
        //if (lpc.GetDistInInnerAV() <= 7000){ continue;}      // this rejects other weird paths

        // ToF
        //double distInInnerAV = lpc.GetDistInInnerAV();
        double distInInnerAV = lpc.GetDistInUpperTarget(); // partial
        double distInAV = lpc.GetDistInAV();
        //double distInWater = lpc.GetDistInWater();
        double distInWater = lpc.GetDistInWater() + lpc.GetDistInLowerTarget();  //partial
        double lightTravelTime = gv.CalcByDistance(distInInnerAV, distInAV, distInWater, energy);

        TVector3 endDir = lpc.GetIncidentVecOnPMT();            // end direction at PMT
        double thetaAtPMT = endDir.Angle(pmtDir)*180./pi;       // incident angle with bucket face
        double lightBucketTime = gv.PMTBucketTime(thetaAtPMT);  // DocDB 3138

        // get ang sys correction
        double thetaRads = theta/180*M_PI;
        double additionalDelayFromAngSys = AngularSystematic(AngB, thetaRads);

        // Store hits
        dir_hits.push_back( pmtTime );
        resid_hits.push_back( pmtTime - lightTravelTime - lightBucketTime - additionalDelayFromAngSys );
        t_ToF.push_back( lightTravelTime );
        t_bucket.push_back( lightBucketTime );
        t_angsys.push_back( additionalDelayFromAngSys );
        fToF->SetPoint( goodHit, theta, lightTravelTime );
        fBucket->SetPoint( goodHit, thetaAtPMT, lightBucketTime );
        fAngSys->SetPoint( goodHit, theta, additionalDelayFromAngSys );
        fHitTimesAng->Fill( theta, pmtTime );
        fResidTimesAng->Fill( theta, pmtTime - lightTravelTime - lightBucketTime - additionalDelayFromAngSys );
        goodHit++;

      } // pmt loop
    } // entries

  }

  /// Display vector as a string
  string MyUserProc::printVector(const TVector3& v) {
    string out;
    if (v.Mag() < 10) out = Form("(%.3f, %.3f, %.3f)", v.X(),  v.Y(), v.Z());
    else              out = Form("(%.1f, %.1f, %.1f)", v.X(),  v.Y(), v.Z());
    return out.c_str();
  }

  // Angular systematic fit function (the angle is in radians)
  double MyUserProc::AngularSystematic(double b, double angle){
    double y = 0.0;
    // currently not using the ang a value (set to 0)
    y = 0 + b * ( ( 1./cos(angle) ) - 1 );
    return y;
  }
