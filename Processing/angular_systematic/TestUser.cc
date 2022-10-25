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

#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>

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

  virtual string printVector(const TVector3& v);

  virtual Processor::Result DSEvent(DS::Run& run, DS::Entry& ds);

  virtual void FitPromptPeaks(TH2D htime, int NPMTS, std::vector<double> occupancy, std::vector<double> pmtangs, TGraph2DErrors &result);

  virtual void printProgress(int it, int n);

protected:

  // Fit results
  string fits_folder;
  string fits_file;
  string delim;
  stringstream dir_fit_file;
  std::vector<string> fibre_dir_name;
  std::vector<double> fibre_dir_X;
  std::vector<double> fibre_dir_Y;
  std::vector<double> fibre_dir_Z;
  int temp_count;
  string fpath;
  string fname;
  int NCOL;
  TFile *extfile;
  TH2D *hfineD;
  std::vector<TGraph*> gDir;
  TGraph *pFibDir;
  TGraph *pWgtDir;
  TGraph *pFitDir;
  TVector3 *fitDir;
  TVector3 *maxima;
  float nearmax;
  double val_hitPeak;
  double val_residPeak;

  // RAT stuff
  RAT::DU::GroupVelocity gv;
  RAT::DU::LightPathCalculator lpc;
  int NPMTS;

  // Database
  RAT::DBLinkPtr tellie_run_data;
  std::string fibre_db;
  json::Value sub_run_info;
  RAT::DBLinkPtr fibre_data;
  TVector3 fibrepos, fibredir, fitfibredir;
  DS::Run temp_run;
  DB *ratdb;
  RAT::DBLinkPtr belly_fibres;
  std::vector<string> belly_fibres_list;
  std::vector<string> belly_fibres_list_index;
  bool IS_BELLY_FIBRE;

  // LPC
  double energy;
  double fLEDWavelength;
  double LOCALITY;

  // Event loop
  int count;
  int allEvs;
  int trig;
  TVector3 pmtPos;
  TVector3 pmtDir;
  double pmtTime;
  int MAXANG;
  int NBINS;
  std::vector<double> pmtOccup;
  std::vector<double> pmtAngle;

  // ROOT
  TFile *rootfile;
  TF1 *fitResult;
  TH2D *htime;
  TH1D *hpmtseg;
  TGraph2DErrors *pmtfits;
  TGraphErrors *gpmts;
  TH1D *hpeak;
  TH1D *hpeakerr;
  TH1D *hmean;
  TH1D *hmeanerr;
  TH1D *hwdth;
  TH1D *hwdtherr;
  TH2D *hprofile;
  TH1D *hpmtgood;
  double width_min;
  double width_max;
  TCanvas *fHitsC;
  TH1D *fHits;
  TCanvas *fHitsRMSC;
  TH1D *fHitsRMS;

  // Global
  int run_id;
  bool MORE_OUTPUT;
  string imgfile;

  // File for fit
  FILE *out;

  // File for results
  stringstream logFile_namess;
  string logFile_name;
  FILE *logFile;

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
    fits_file = getenv("DIR_FIT_FILE");
    dir_fit_file << fits_file;
    delim = "\t";
    MAXANG = atoi(getenv("ANG_SYS_ANG"));
    NBINS = 24;
    MORE_OUTPUT = true;
    NCOL = 20;
    nearmax = -1;
    val_hitPeak = atof(getenv("HIT_PEAK"));
    val_residPeak = atof(getenv("RESID_PEAK"));
    out = fopen("angular_fit.txt","a");

    // Load direction fit results
    cout << "Opening fit file: " << dir_fit_file.str() << endl;
    std::ifstream fit_file( dir_fit_file.str().c_str() );
    for ( string line; getline(fit_file, line); ){
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

    // Retrieve rotated detector view from fibre validation document
    /// this is obsolete now
    // gDir.resize(NCOL+2);
    // fpath = "/home/michal/Automation/Processing/position_fit";
    // //fname = Form("%s/PCA_%s.root",fpath.c_str(),fibre.c_str());
    // fname = Form("%s/TEST.root",fpath.c_str());
    // ifstream ext(fname.c_str());
    // if (ext.good()) {
    //   extfile = new TFile(fname.c_str(),"READ");
    //   hfineD = (TH2D*)extfile->Get("hfineD");
    //   //gDir2D = (TGraph2D*)extfile->Get("gDir2D");
    //   for(int s=0; s<NCOL+2; s++) {
    //     string ss = Form("[%d]",s);
    //     gDir[s] = (TGraph*)extfile->Get(("gDir"+ss).c_str());
    //   }
    //   pFibDir = (TGraph*)extfile->Get("pFibDir");
    //   pWgtDir = (TGraph*)extfile->Get("pWgtDir");
    //   pFitDir = (TGraph*)extfile->Get("pFitDir");
    //   fitDir = (TVector3*)extfile->Get("fitDir");
    //   maxima = (TVector3*)extfile->Get("maxima");
    //   nearmax = maxima->X();
    //   cout << "Loaded detector view plots." << endl;
    // } else {
    //   cout << "*** WARNING: File containing 2D graphs not found!" << endl;
    // }

    // Load belly plate fibres from ratdb file
    ratdb = RAT::DB::Get();
    temp_run.SetRunID(275000);      // this allows to load custom table outside BeginOfRun for all runs
    ratdb->BeginOfRun(temp_run);    // the temp number is random
    try {
      belly_fibres = ratdb->GetLink("TELLIE_BELLY_PLATES");
      belly_fibres_list = belly_fibres->GetSArray("list");
      belly_fibres_list_index = belly_fibres->GetSArray("fibre_index");
    }
    catch (DBNotFoundError& e) {
      std::cout << "Couldn't load belly fibres from ratdb!" << std::endl;
      cout << "Reason: " << e.what() << endl;
    }
    cout << "Loaded belly plates data: " << belly_fibres_list.size() << " fibres affected." << endl;

    /*
    for (size_t i=0; i<fibre_dir_name.size(); i++) {
      cout << fibre_dir_name[i] << " " << fibre_dir_X[i] << " " << fibre_dir_Y[i] << " " << fibre_dir_Z[i] << endl;
    }
    */

  }

  MyUserProc::~MyUserProc() {
  }

  void MyUserProc::BeginOfRun( DS::Run& run ) {

    // Initialise counters
    count = 0;
    allEvs = 0;
    run_id = run.GetRunID();
    width_min = 5000;
    width_max = 0;

    // Initialise RAT vars
    gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    gv.BeginOfRun();
    lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lpc.BeginOfRun();
    lpc.SetELLIEEvent(true);
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    NPMTS = pmtinfo.GetCount();
    cout << "Initialised RAT. Number of PMTs is " << NPMTS << "." << endl;

    // Reset vectors
    pmtOccup.clear();
    pmtAngle.clear();
    pmtOccup.resize(NPMTS);
    pmtAngle.resize(NPMTS);

    // Load run fibre details
    try {
      tellie_run_data = DB::Get()->GetLink("TELLIE_RUN");
      sub_run_info = tellie_run_data->GetJSON("sub_run_info");
      fibre_db = sub_run_info[0]["fibre"].getString();
    }
    catch (DBNotFoundError& e) {
      std::cout << "Proc(): Couldn't load TELLIE Run details from ratdb!" << std::endl;
    }
    cout << "FIBRE: " << fibre_db << endl;

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
    cout << "Set fibre as: fibre " << fibre_db << ", pos " << printVector(fibrepos) << ", dir " << printVector(fitfibredir) << endl;

    // Find out if direct light from fibre is affected by belly plates
    IS_BELLY_FIBRE = false;
    for (size_t i=0; i<belly_fibres_list.size(); i++) {
      if ( belly_fibres_list[i] == fibre_db ){
        IS_BELLY_FIBRE = true;
      }
    }
    cout << "Is affected by belly plate? " << IS_BELLY_FIBRE << endl;

    // Initialise ROOT stuff
    htime = new TH2D("htime","",NPMTS,0,NPMTS,500,0,500);
    hpmtseg = new TH1D("hpmtseg",fibre_db.c_str(),NBINS,0,MAXANG);
    gpmts = new TGraphErrors();
    fHitsC = new TCanvas ("fHitsC", "fHitsC", 1024, 768);
    fHitsRMSC = new TCanvas ("fHitsRMSC", "fHitsRMSC", 1024, 768);
    fHits = new TH1D ("fHits", "", 100, val_residPeak-40, val_residPeak+40 );
    fHitsRMS = new TH1D ("fHitsRMS", "", 100, 0, 0.7 );

    // Output file for stats
    logFile_namess.str("");
    logFile_namess << run_id << "_as.log";
    logFile_name = logFile_namess.str();

  }

  void MyUserProc::EndOfRun( DS::Run& run ) {

    // Turn array of PMTs into percentage (occupancy)
    cout << "Number of events was " << count << endl;
    cout << "Number of EXTA events was " << allEvs << endl;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      pmtOccup[iPMT] /= (float)allEvs;
    }

    // Fit PMT hit times (Gaussian prompt peak)
    pmtfits = new TGraph2DErrors();
    FitPromptPeaks(*htime, NPMTS, pmtOccup, pmtAngle, *pmtfits);
    if (!pmtfits) { cout << "*** ERROR *** Could not fit PMT prompt peaks!" << endl; }
    if (pmtfits->GetN()==0) { cout << "*** ERROR *** Graph contains zero points!" << endl; }

    double *tx = pmtfits->GetX();   // Constant (A)
    double *ty = pmtfits->GetY();   // Mean value (mu)
    double *tz = pmtfits->GetZ();   // Uncertainty (sigma)
    double *tex = pmtfits->GetEX();
    double *tey = pmtfits->GetEY();
    double *tez = pmtfits->GetEZ();

    // Fill graph with time vs angle (use only PMTs of interest)
    gpmts = new TGraphErrors();
    int npmts=0;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      if (tx[iPMT]+ty[iPMT]+tz[iPMT]==0) continue;
      if (ty[iPMT] > 1e3) {
        cout << "Fitted mean hit time for PMT #" << iPMT << " is (" << ty[iPMT] << " +- " << tey[iPMT] << ") ns! Setting to 0..." << endl;
        tx[iPMT]=0; ty[iPMT]=0; tz[iPMT]=0;
        continue;
      }
      else {
        cout << "Fitted time for PMT: " << iPMT << " = " << ty[iPMT] << " +- " << tey[iPMT] << endl;
        fHits->Fill(ty[iPMT]);
        fHitsRMS->Fill(tey[iPMT]);
      }
      gpmts->SetPoint(npmts, pmtAngle[iPMT], ty[iPMT]);
      gpmts->SetPointError(npmts, 0, tey[iPMT]);
      npmts++;
    }
    cout << endl;

    // Mean and RMS of entire graph
    double meanhittime = gpmts->GetMean(2);
    double rmshittime = gpmts->GetRMS(2);
    cout << "Mean hit time = " << meanhittime << " ns, RMS = " << rmshittime << " ns." << endl;
    logFile = fopen(logFile_name.c_str(), "w");
    fprintf(logFile, "Mean hit time and RMS: %f\t%f\n", meanhittime, rmshittime);

    // Investigate PMTs with unusual offsets w.r.t. mean hit time
    // TODO: add check here for any meantime that's + 20 of average

    // Fill fit results into histograms
    hpeak = new TH1D("hpeak","",NBINS,0,NBINS);         // mean intensity
    hmean = new TH1D("hmean","",NBINS,0,NBINS);         // mean hit time
    hwdth = new TH1D("hwdth","",NBINS,0,NBINS);         // mean signal width
    hpeakerr = new TH1D("hpeakerr","",NBINS,0,NBINS);   // error on mean intensity
    hmeanerr = new TH1D("hmeanerr","",NBINS,0,NBINS);   // error on mean hit time
    hwdtherr = new TH1D("hwdtherr","",NBINS,0,NBINS);   // error on mean signal width
    int tb;
    std::vector<int> tn;
    std::vector<double> tavgx;
    std::vector<double> tavgy;
    std::vector<double> tavgz;
    std::vector<double> terrx;
    std::vector<double> terry;
    std::vector<double> terrz;
    tn.resize(NBINS);
    tavgx.resize(NBINS);
    tavgy.resize(NBINS);
    tavgz.resize(NBINS);
    terrx.resize(NBINS);
    terry.resize(NBINS);
    terrz.resize(NBINS);
    // Loop over all PMTs in the graph and bin the data
    for (int i=0; i<NPMTS; i++) {
      if (tx[i]+ty[i]+tz[i]==0) continue;     // not good PMT
      tb = (int)pmtAngle[i];                  // bin index (1 degree bins)
      tn[tb]++;                               // number of PMTs in each bin
      tavgx[tb] += tx[i]/(tex[i]*tex[i]);     // sum of weighted intensities
      terrx[tb] += 1./(tex[i]*tex[i]);        // sum of errors squared
      tavgy[tb] += ty[i]/(tey[i]*tey[i]);     // sum of weighted hit times
      terry[tb] += 1./(tey[i]*tey[i]);        // sum of errors squared
      tavgz[tb] += tz[i]/(tez[i]*tez[i]);     // sum of weighted signal widths
      terrz[tb] += 1./(tez[i]*tez[i]);        // sum of errors squared
    }
    // Loop over each bin to get the average and propagated error
    for (int j=0; j<NBINS; j++) {
      if (tn[j]==0) continue;                 // no PMTs in this bin
      // Mean intensity
      tavgx[j] = tavgx[j]/terrx[j];           // weighted mean intensity
      terrx[j] = sqrt(1./terrx[j]);           // error on weighted mean
      hpeak->SetBinContent(j+1,tavgx[j]);
      hpeak->SetBinError(j+1,terrx[j]);
      hpeakerr->SetBinContent(j+1,terrx[j]);
      // Mean hit time
      tavgy[j] = tavgy[j]/terry[j];           // weighted mean hit time
      terry[j] = sqrt(1./terry[j]);           // error on weighted mean
      hmean->SetBinContent(j+1,tavgy[j]);
      hmean->SetBinError(j+1,terry[j]);
      hmeanerr->SetBinContent(j+1,terry[j]);
      // Mean signal width
      tavgz[j] = tavgz[j]/terrz[j];           // weighted mean signal width
      terrz[j] = sqrt(1./terrz[j]);           // error on weighted mean
      hwdth->SetBinContent(j+1,tavgz[j]);
      hwdth->SetBinError(j+1,terrz[j]);
      hwdtherr->SetBinContent(j+1,terrz[j]);
      if (tavgz[j] > width_max){width_max = tavgz[j];}
      if (tavgz[j] < width_min){width_min = tavgz[j];}
    }

    // Get central hit time and initialise 2D histogram
    int centraltime = (int)round(meanhittime + 5/2.);
    centraltime -= centraltime % 5;           // rounded to the nearest multiple of 5
    hprofile = new TH2D("hprofile",fibre_db.c_str(),NBINS,0,MAXANG,30,centraltime-15,centraltime+15);
    hpmtgood = new TH1D("hpmtgood",fibre_db.c_str(),NBINS,0,MAXANG);

    for (int i=0; i<NPMTS; i++) {
      if (tx[i]+ty[i]+tz[i]==0) continue;     // not good PMT
      hpmtgood->Fill(pmtAngle[i]);
      // Only fill histogram if pmttime is within prompt peak (90% area ~ 1.645 sigma)
      for (int j=0; j<=htime->GetNbinsY()+1; j++) {
        float time = htime->GetYaxis()->GetBinCenter(j);
        if (fabs(time - ty[i]) > 1.645*tz[i]) continue; // not in prompt peak
        for (int k=0; k<(int)htime->GetBinContent(i+1,j+1); k++)
          hprofile->Fill(pmtAngle[i], time);  // Angle of PMT w.r.t. fitted fibre direction [deg], time [ns]
      }
    }

    // Normalise angular slices to number of PMTs in slice
    double tmp, goodpmts;
    for (int binx=0; binx<=hprofile->GetNbinsX()+1; binx++) {     // include underflow/overflow
      goodpmts = hpmtgood->GetBinContent(binx);
      for (int biny=0; biny<=hprofile->GetNbinsY()+1; biny++) {   // include underflow/overflow
        tmp = hprofile->GetBinContent(binx,biny);
        if (goodpmts<=0) hprofile->SetBinContent(binx,biny,0);
        else hprofile->SetBinContent(binx,biny,tmp/goodpmts);
      }
    }

    // Write all objects to file (histograms included automatically, graphs not)
    /*
    rootfile = new TFile("TEST.root","RECREATE");
    pmtfits->Write("pmtfits");
    gpmts->Write("gpmts");
    rootfile->Write();
    cout << "Wrote output to file " << "TEST.root" << endl;
    */

    double meanhittime2 = gpmts->GetMean(2);
    int centraltime2 = (int)round(meanhittime2 + 5/2.);
    centraltime2 -= centraltime2 % 5;           // rounded to the nearest multiple of 5
    TH2D *hprofile2 = new TH2D("hprofile2",fibre_db.c_str(),NBINS,0,MAXANG,30,centraltime2-15,centraltime2+15);
    for (int i=0; i<=hprofile->GetNbinsX()+1; i++) {
      for (int j=0; j<=hprofile->GetNbinsY()+1; j++) {
        hprofile2->SetBinContent(i,j,hprofile->GetBinContent(i,j));
      }
    }

    // Fit angular systematic: y = a - b + b/cos(x)
    TCanvas *c1 = new TCanvas("c1","",800,600);
    TF1 *fitResult = new TF1("fitResult", "[0] - [1] + [1]/cos(x/180.*pi)", 0, 24);
    fitResult->SetParameter(0, gpmts->GetMean(2));  // value at zero angle
    fitResult->SetParameter(1, 0);                  // assume flat line a priori
    gpmts->Fit("fitResult","R,q");                  // force range, quiet mode
    c1->Close(); delete c1;

    // Fill histograms with time residual
    TH1D *hresid = new TH1D("hresid","",30,-3,3);
    //TH1D *hpulls = new TH1D("hpulls","",40,-20,20);
    double x,y,ey,y0;
    for (int i=0; i<gpmts->GetN(); i++) {
      gpmts->GetPoint(i,x,y);
      //ey = gpmts->GetErrorY(i);
      y0 = fitResult->Eval(x);
      hresid->Fill(y-y0);
      //hpulls->Fill((y-y0)/ey);
    }
    c1 = new TCanvas("c1","",800,600);
    TF1 *fitresid = new TF1("fitresid","gaus",-3,3);
    fitresid->SetParameters(hresid->GetMaximum(),hresid->GetMean(),hresid->GetRMS());
    hresid->Fit("fitresid","R,q");
    double sigma_resid = fitresid->GetParameter(2);
    c1->Close(); delete c1;

    // ******************
    //  PLOTTING SECTION
    // ******************
    int minvalY = (int)round(fitResult->GetParameter(0));

    // Plotting options
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.6);
    gStyle->SetStatY(0.88);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.18);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadBorderSize(0);

    // *****
    // Define canvas and pads
    TCanvas *pad0 = new TCanvas("pad0","",1024,768);
    TCanvas *pad1 = new TCanvas("pad1","",1024,768);
    TCanvas *pad2 = new TCanvas("pad2","",1024,768);
    TCanvas *pad3 = new TCanvas("pad3","",1024,768);
    TCanvas *pad4 = new TCanvas("pad4","",1024,768);
    TCanvas *pad5 = new TCanvas("pad5","",1024,768);
    TCanvas *pad6 = new TCanvas("pad6","",1024,768);

    // *****
    // Text boxes to display fit results
    TBox *tbox = new TBox(0.6,minvalY+4.45,15.5,minvalY+6.85);
    tbox->SetLineColor(2);
    tbox->SetFillColor(kYellow-9);
    TLatex *tfit[4] = {NULL};
    tfit[0] = new TLatex(1, minvalY+6.7, "Fit function: y = a #plus b#left(#frac{1}{cos(x)} #minus 1#right)");
    tfit[1] = new TLatex(1, minvalY+5.9, Form(" #Rightarrow a = ( %.3lf #pm %.3lf ) ns",
                                              fitResult->GetParameter(0), fitResult->GetParError(0)));
    tfit[2] = new TLatex(1, minvalY+5.4, Form(" #Rightarrow b = ( %.3lf #pm %.3lf ) ns",
                                              fitResult->GetParameter(1), fitResult->GetParError(1)));
    tfit[3] = new TLatex(1, minvalY+4.95, Form(" #Rightarrow #chi^{2}/ndf = %d / %d",
                                               (int)round(fitResult->GetChisquare()), fitResult->GetNDF()));
    for (int l=0; l<4; l++) {
      tfit[l]->SetTextAlign(13);
      tfit[l]->SetTextFont(62);
      tfit[l]->SetTextSize(0.03);
    }

    // *****
    // PMT hit time offsets and parametrised angular systematic
    pad0->cd();
    pad0->SetGrid();
    //pad0->SetLeftMargin(0.12);    // for label
    //pad0->SetBottomMargin(0.12);  // for label
    string tstr = Form("TELLIE angular systematic (run %d);Angle of PMT w.r.t. fitted fibre direction [deg];Offset in PMT hit time [ns]", run_id);
    gpmts->SetTitle(tstr.c_str());
    gpmts->SetMarkerColor(4);
    gpmts->SetMarkerStyle(6);
    gpmts->Draw("AP");
    gpmts->GetXaxis()->SetTitleOffset(1.2);
    gpmts->GetYaxis()->SetTitleOffset(1.3);
    gpmts->GetXaxis()->SetLimits(0,24);
    gpmts->GetYaxis()->SetRangeUser(minvalY-3,minvalY+7); // suppresses outliers!
    //tbox->Draw("L same");
    //for (int l=0; l<4; l++) tfit[l]->Draw("same");

    //Create a histogram to hold the confidence intervals
    TH1D *hint = new TH1D("hint", "Fitted function with error band", 10*NBINS, 0, MAXANG);
    TH1D *hint2 = new TH1D("hint2", "Fitted function with error band", 10*NBINS, 0, MAXANG);
    TH1D *hint3 = new TH1D("hint3", "Fitted function with error band", 10*NBINS, 0, MAXANG);
    //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
    //Now the "hint" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors
    for(int i=0;i<=hint->GetNbinsX()+1;i++) {
      float x = (float)MAXANG/hint->GetNbinsX()*(i-1);
      hint->SetBinContent(i,fitResult->Eval(x));
      hint->SetBinError(i,sigma_resid);
      hint2->SetBinContent(i,fitResult->Eval(x));
      hint2->SetBinError(i,2*sigma_resid);
      hint3->SetBinContent(i,fitResult->Eval(x));
      hint3->SetBinError(i,3*sigma_resid);
    }
    hint->SetStats(kFALSE);
    hint->SetLineColor(0);
    hint->SetFillColorAlpha(2, 0.2);
    hint->Draw("e3 same");
    hint2->SetStats(kFALSE);
    hint2->SetLineColor(0);
    hint2->SetFillColorAlpha(2, 0.2);
    hint2->Draw("e3 same");
    hint3->SetStats(kFALSE);
    hint3->SetLineColor(0);
    hint3->SetFillColorAlpha(2, 0.2);
    hint3->Draw("e3 same");
    gpmts->Draw("P same");

    stringstream temp;
    temp << run_id << "_pad0.png";
    pad0->SaveAs(temp.str().c_str());
    pad0->Close();

    // *****
    // Normalised intensity profile (time vs angle)
    pad1->cd();
    pad1->SetGrid();
    //pad1->SetRightMargin(0.15);   // for TH2D color scale
    hprofile2->SetTitle(Form("Normalised intensity profile (%s);Angle of PMT w.r.t. fitted fibre direction [deg];Offset in PMT hit time [ns]",fibre_db.c_str()));
    hprofile2->SetStats(0);
    hprofile2->Draw("colz");
    hprofile2->GetXaxis()->SetTitleOffset(1.2);
    hprofile2->GetYaxis()->SetTitleOffset(1.3);
    stringstream temp1;
    temp1 << run_id << "_pad1.png";
    pad1->SaveAs(temp1.str().c_str());
    pad1->Close();

    // *****
    // Mean hit times (binned PMTs)
    pad2->cd();
    pad2->SetGrid();
    hmean->SetTitle("Signal mean (#mu);Angle of PMT w.r.t. fitted fibre direction [deg];Mean PMT signal offset [ns]");
    hmean->SetLineWidth(2);
    hmean->Draw();
    hmean->GetXaxis()->SetTitleOffset(1.2);
    hmean->GetYaxis()->SetTitleOffset(1.3);
    hmean->GetYaxis()->SetRangeUser(206, *max_element(tavgy.begin(), tavgy.end())+2);
    //hmeanerr->GetYaxis()->SetRangeUser(0,0.16);
    stringstream temp2;
    temp2 << run_id << "_pad2.png";
    pad2->SaveAs(temp2.str().c_str());
    pad2->Close();

    // *****
    // Mean signal widths
    pad3->cd();
    pad3->SetGrid();
    hwdth->SetTitle("Signal width (#sigma);Angle of PMT w.r.t. fitted fibre direction [deg];Mean PMT signal width [ns]");
    hwdth->SetLineWidth(2);
    hwdth->Draw();
    hwdth->GetXaxis()->SetTitleOffset(1.2);
    hwdth->GetYaxis()->SetTitleOffset(1.3);
    hwdth->GetYaxis()->SetRangeUser(width_min-0.25, width_max+0.25);
    stringstream temp3;
    temp3 << run_id << "_pad3.png";
    pad3->SaveAs(temp3.str().c_str());
    pad3->Close();

    // *****
    // Mean intensities and number of PMTs in angular segment
    pad4->cd();
    pad4->SetGrid();
    //pad4->SetLeftMargin(0.12);
    hpeak->SetTitle("Signal amplitude (A);Angle of PMT w.r.t. fitted fibre direction [deg];Mean PMT signal amplitude [a.u.]");
    hpeak->SetLineWidth(2);
    hpeak->Draw();
    hpeak->GetXaxis()->SetTitleOffset(1.2);
    hpeak->GetYaxis()->SetTitleOffset(1.3);
    hpeak->GetYaxis()->SetRangeUser(0, *max_element(tavgx.begin(), tavgx.end())+200);
    // PMTs in that bin
    hpmtseg->SetLineWidth(2);
    hpmtseg->SetLineColor(2);
    hpmtseg->SetTitle("# All PMTs");
    hpmtseg->Draw("same");
    // PMTs with occupancy >1%
    hpmtgood->SetLineWidth(2);
    hpmtgood->SetLineColor(3);
    hpmtgood->SetTitle("# Good PMTs");
    hpmtgood->Draw("same");
    // Legend
    TLegend *leg = pad4->BuildLegend();
    stringstream temp4;
    temp4 << run_id << "_pad4.png";
    pad4->SaveAs(temp4.str().c_str());
    pad4->Close();

    // *****
    // Mean errors (calculated with weighted arithmetic mean method)
    pad5->cd();
    pad5->SetGrid();
    //pad5->SetLeftMargin(0.12);
    hmeanerr->SetTitle("Error on mean (s_{#mu});Angle of PMT w.r.t. fitted fibre direction [deg];Error on mean PMT signal offset [ns]");
    hmeanerr->SetLineWidth(2);
    hmeanerr->Draw();
    hmeanerr->GetXaxis()->SetTitleOffset(1.2);
    hmeanerr->GetYaxis()->SetTitleOffset(1.3);
    hmeanerr->GetYaxis()->SetRangeUser(*min_element(terry.begin(), terry.end())-0.005, *max_element(terry.begin(), terry.end())+0.005);
    stringstream temp5;
    temp5 << run_id << "_pad5.png";
    pad5->SaveAs(temp5.str().c_str());
    pad5->Close();

    // *****
    // Time residuals
    pad6->cd();
    pad6->SetGrid();
    hresid->SetTitle("PMT time residuals;Time residual [ns];#PMTs/bin");
    hresid->SetLineWidth(2);
    hresid->Draw();
    hresid->GetXaxis()->SetTitleOffset(1.2);
    hresid->GetYaxis()->SetTitleOffset(1.3);
    hresid->GetYaxis()->SetRangeUser(0, hresid->GetMaximum()+75);
    stringstream temp6;
    temp6 << run_id << "_pad6.png";
    pad6->SaveAs(temp6.str().c_str());
    pad6->Close();

    // View from direct light spot (fitted)
    // this is also obsolete now
    // pad7->cd()->SetGrid();
    // ifstream ext(fname.c_str());
    // if (ext.good()) {
    //   /*
    //   TH3F *hempty = new TH3F("hempty","",10,-10,10,10,-10,10,10,0,1000+1);
    //   hempty->Draw("");             // empty histogram for plot range
    //   pad7->SetTheta(90-0.001);      // view from above
    //   pad7->SetPhi(0+0.001);         // no x-y rotation
    //   gDir2D->Draw("pcol,same");
    //   */
    //   hfineD->SetTitle("Direct light spot;X' [m];Y' [m]");
    //   hfineD->GetXaxis()->SetTitleOffset(1.3);
    //   hfineD->GetYaxis()->SetTitleOffset(1.4);
    //   hfineD->Draw("scat");
    //   hfineD->SetStats(0);
    //   for(int s=0;s<NCOL+2;s++) {
    //     if(!gDir[s]) continue;
    //     if(gDir[s]->GetN()==0) continue;
    //     gDir[s]->SetMarkerStyle(6);
    //     gDir[s]->Draw("P same");
    //   }
    //   pFibDir->Draw("P same");
    //   pWgtDir->Draw("P same");
    //   TGraph *pFitDir2 = (TGraph*)pFitDir->Clone();
    //   double fitD_rotX = fitDir->X()/1e3;
    //   double fitD_rotY = fitDir->Y()/1e3;
    //   if (IS_BELLY_FIBRE) pFitDir2->SetMarkerStyle(28); // open cross to indicate Gaussian fit is not used
    //   pFitDir2->DrawGraph(1,&fitD_rotX,&fitD_rotY,"P same");
    //   // Indicate possible belly plate effect
    //   TLatex *txtB = new TLatex(-9.25,9.5,"BELLY PLATE");
    //   txtB->SetTextAlign(13);
    //   txtB->SetTextFont(102);
    //   txtB->SetTextSize(0.04);
    //   if (IS_BELLY_FIBRE) txtB->Draw();
    //   // Indicate colour scale
    //   TLatex *txtD = new TLatex(9.5,9.5,Form("Occup. #leq %.2f%%",100.*nearmax));
    //   txtD->SetTextAlign(33);
    //   txtD->SetTextFont(82);
    //   txtD->SetTextSize(0.04);
    //   txtD->Draw();
    // }

    // *****
    // Save canvas and close
    //imgfile = Form("angular_%s",fibre_db.c_str());
    //c0->Print(Form("%s.png",imgfile.c_str()));
    //c0->Print(Form("%s.pdf",imgfile.c_str()));
    //c0->Close();

    fHitsC->cd();
    fHits->SetTitle("fHits");
    fHits->GetXaxis()->SetTitle("Residual time [ns]");
    fHits->GetYaxis()->SetTitle("count");
    fHits->Draw();
    stringstream hitspng;
    hitspng << run_id << "_hits.png";
    fHitsC->SaveAs(hitspng.str().c_str());
    fHitsC->Close();

    fHitsRMSC->cd();
    fHitsRMS->SetTitle("fHitsRMS");
    fHitsRMS->GetXaxis()->SetTitle("Residual time - RMS [ns]");
    fHitsRMS->GetYaxis()->SetTitle("count");
    fHitsRMS->Draw();
    stringstream hitsrmspng;
    hitsrmspng << run_id << "_hitsrms.png";
    fHitsRMSC->SaveAs(hitsrmspng.str().c_str());
    fHitsRMSC->Close();

    // Close files and free memory
    //rootfile->Close();

    // *****
    // Clear ROOT graphs&memory
    delete hprofile; delete hprofile2; delete hpmtgood;
    delete htime; delete hpmtseg; delete gpmts;
    delete hpeak; delete hmean; delete hwdth; delete hpeakerr; delete hmeanerr; delete hwdtherr;
    delete hresid; delete hint; delete hint2; delete hint3;
    delete pmtfits;
    delete fHits; delete fHitsRMS;

    // Store to file for next steps
    stringstream temp_output;
    temp_output << fibre_db << "\t" << fitResult->GetParameter(0) << "\t" << fitResult->GetParError(0) << "\t" << fitResult->GetParameter(1) << "\t" << fitResult->GetParError(1) << "\n";
    out = fopen("angular_fit.txt","a");
    fprintf(out, temp_output.str().c_str());
    fclose(out);

    fprintf(logFile, "Ang a: %f\t%f\n", fitResult->GetParameter(0), fitResult->GetParError(0));
    fprintf(logFile, "Ang b: %f\t%f\n", fitResult->GetParameter(1), fitResult->GetParError(1));
    fclose(logFile);

    return;

  }

  Processor::Result MyUserProc::DSEvent(DS::Run&, DS::Entry& ds) {

    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {
      const RAT::DS::EV& ev = ds.GetEV(iEv);

      count++;

      // Trigger type
      trig = ev.GetTrigType();
      if (!(trig & 0x8000)) continue;    // EXT trigger only

      allEvs++;

      const RAT::DS::CalPMTs& pmts = ev.GetCalPMTs();   // calib PMTs
      for(int iPMT=0; iPMT<pmts.GetNormalCount(); iPMT++) {
        RAT::DS::PMTCal pmt = pmts.GetNormalPMT(iPMT);
        int pmtID = pmt.GetID();

        const RAT::DU::PMTCalStatus& pmtStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
        const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
        const RAT::DU::PMTInfo& pmtinfo_loop = RAT::DU::Utility::Get()->GetPMTInfo();
        unsigned int status = pmtStatus.GetHitStatus(pmt);
        if(status & (1<<pmtStatus.kCHSBit)) continue;
        if(status & (1<<pmtStatus.kECABit)) continue;
        if(status & (1<<pmtStatus.kPCABit)) continue;
        if(status & (1<<pmtStatus.kXTalkBit)) continue;
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

        if (lpc.GetTIR() == 1) { continue; }           // total internal reflection
        if (lpc.GetPathValid() == 0) { continue; }      // check whether path parameters are valid
        if (lpc.GetResvHit() == 1) { continue; }        // whether end point was within locality

        // additional path checks
        //if (lpc.GetTotalDist() <= 12000){ continue;}        // this rejects near reflections
        //if (lpc.GetDistInInnerAV() <= 7000){ continue;}      // this rejects other weird paths
        if (lpc.GetTotalDistPartial() <= 6000){ CDIST++; continue;}  //partial
        if (lpc.GetTotalDistPartial() == lpc.GetDistInWater()){ CDIST++; continue;} //partial

        //double distInInnerAV = lpc.GetDistInInnerAV();
        double distInInnerAV = lpc.GetDistInUpperTarget(); // partial
        double distInAV = lpc.GetDistInAV();
        //double distInWater = lpc.GetDistInWater();
        double distInWater = lpc.GetDistInWater() + lpc.GetDistInLowerTarget();  //partial
        double lightTravelTime = gv.CalcByDistance(distInInnerAV, distInAV, distInWater, energy);

        // Get light bucket time
        TVector3 endDir = lpc.GetIncidentVecOnPMT();        // end direction at PMT
        double thetaAtPMT = endDir.Angle(pmtDir)*180./pi;   // incident angle with bucket face
        double lightBucketTime = gv.PMTBucketTime(thetaAtPMT);  // DocDB 3138

        // Get residual time after correcting for all offsets
        double emissionTime = pmtTime - lightTravelTime - lightBucketTime;

        // Get light emission angle
        TVector3 startDir = lpc.GetInitialLightVec();       // start direction at fibre
        double theta = startDir.Angle(fitfibredir)*180./pi;      // emission angle at fibre

        /*
        cout << "hit time: " << pmtTime << endl;
        cout << "ToF: " << lightTravelTime << endl;
        cout << "bucket: " << lightBucketTime << endl;
        cout << "injection: " << emissionTime << endl;
        cout << "angle: " << theta << endl << endl;
        */

        // Fill histograms/arrays
        if (theta > MAXANG) continue;
        htime->Fill(pmtID, emissionTime);                   // time residual vs PMT ID
        if (pmtOccup[pmtID] == 0) {                   // fill angles only once per PMT
          hpmtseg->Fill(theta);
          pmtAngle[pmtID] = theta;
          //cout << "setting angle: pmtID: " << pmtID << ", angle: " << theta << endl;
        }
        pmtOccup[pmtID]++;

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

  // -----------------------------------------------------------------------------
  /// Fit the prompt peaks (direct light) of the individual PMTs hit time
  /// Input:  run number, time-vs-pmtID histogram, number of PMTs, PMT occupancy,
  ///         PMT angle w.r.t. fibre, array for overall fit parameters
  /// Output: 3D graph containing fit results for each PMT
  void MyUserProc::FitPromptPeaks(TH2D htime, int NPMTS, std::vector<double> occupancy, std::vector<double> pmtangs, TGraph2DErrors &result) {

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> ex;
    std::vector<double> ey;
    std::vector<double> ez;
    x.resize(NPMTS);
    y.resize(NPMTS);
    z.resize(NPMTS);
    ex.resize(NPMTS);
    ey.resize(NPMTS);
    ez.resize(NPMTS);

    std::vector<int> saved_pmt;
    saved_pmt.resize(24);
    string histname[24] = {""};
    string histtitle[24] = {""};
    string fitname[24] = {""};
    std::vector<TF1*> plotted_fit;
    std::vector<TH1D*> plotted_hist;
    plotted_fit.resize(24);
    plotted_hist.resize(24);
    TCanvas *c1 = NULL;
    TH1D *temp = NULL;
    TF1 *fitPMT = NULL;

    cout << "Fitting PMT prompt peaks..." << endl;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      //printProgress(iPMT, NPMTS);

      // Reject PMTs outside ROI
      if (occupancy[iPMT]<0.01) continue; // only consider PMTs with >=1% occupancy
      if (pmtangs[iPMT]>24.0) continue;   // only consider PMTs within 2x nominal aperture (12 deg)
      temp = htime.ProjectionY("temp",iPMT+1,iPMT+1,""); // histogram bins in [1,N]
      if (!temp) { cout << "*** WARNING: Could not obtain projection for PMT #" << iPMT << "!" << endl; continue; }
      if (temp->GetEntries()==0) { cout << "*** WARNING: Projection for PMT #" << iPMT << " is empty!" << endl; continue; }
      if (temp->GetMean()>1e3) { cout << "*** WARNING: Mean value for PMT #" << iPMT << " is " << temp->GetMean() << "ns! Skipping..." << endl; continue; }

      // Define prompt peak range (>20% of max. intensity)
      int lobin = temp->GetMaximumBin();
      int hibin = temp->GetMaximumBin();
      while(temp->GetBinContent(lobin) > 0.2*temp->GetMaximum()) lobin--;
      while(temp->GetBinContent(hibin) > 0.2*temp->GetMaximum()) hibin++;

      // Fit 1D Gaussian to each slice to get prompt hit time
      c1 = new TCanvas("c1","NULL",800,600);  // suppress default Canvas creation
      fitPMT = new TF1("fitPMT", "gaus", temp->GetBinLowEdge(lobin), temp->GetBinLowEdge(hibin+1));
      fitPMT->SetParameters(temp->GetMaximum(), temp->GetBinCenter(htime.GetMaximumBin()), temp->GetRMS());
      temp->Fit("fitPMT","R,q");  // force range, quiet mode
      c1->Close();
      c1 = NULL;

      // Fill fitted hit times vs. angle into graph
      x[iPMT]  = fitPMT->GetParameter(0);
      y[iPMT]  = fitPMT->GetParameter(1);
      z[iPMT]  = fitPMT->GetParameter(2);
      ex[iPMT] = fitPMT->GetParError(0);
      ey[iPMT] = fitPMT->GetParError(1);
      ez[iPMT] = fitPMT->GetParError(2);

      // Save fit results for individual PMT (for proof of concept)
      int thisbin = (int)floor(pmtangs[iPMT]);
      if (MORE_OUTPUT && saved_pmt[thisbin]==0) {
        // Histogram for this bin
        histname[thisbin] = Form("hangle_%d",thisbin);
        histtitle[thisbin] = Form("PMT #%d (%d deg)",iPMT,thisbin);
        plotted_hist[thisbin] = (TH1D*)temp->Clone();
        plotted_hist[thisbin]->SetDirectory(0); // required to avoid segfault (?)
        plotted_hist[thisbin]->SetName(histname[thisbin].c_str());
        // Fit result for this bin
        fitname[thisbin] = Form("hfit_%d",thisbin);
        plotted_fit[thisbin] = (TF1*)fitPMT->Clone();
        plotted_fit[thisbin]->SetName(fitname[thisbin].c_str());
        saved_pmt[thisbin] = 1;
      }

      if (fitPMT) delete fitPMT;
      if (temp) delete temp;
      fitPMT = NULL;
      temp = NULL;

    } // PMT loop

    if (c1) delete c1;
    if (fitPMT) delete fitPMT;
    if (temp) delete temp;

    // Plot fit results for individual PMT (for proof of concept)
    TCanvas *c0 = NULL;
    if (MORE_OUTPUT) {
      TCanvas *c0 = new TCanvas("c0","Single PMT fits",3000,2000);
      c0->Divide(6,4);
      for (int thisbin=0; thisbin<24; thisbin++) {
        if (saved_pmt[thisbin]==0 || plotted_hist[thisbin]==NULL) continue;
        //cout << "PMT at angle " << thisbin << " deg has " << plotted_hist[thisbin]->GetEntries() << " entries around " << plotted_hist[thisbin]->GetMean() << endl;
        c0->cd(thisbin+1)->SetGrid();
        //c0->cd(thisbin+1)->SetLogy();
        plotted_hist[thisbin]->SetStats(0);
        plotted_hist[thisbin]->Draw();
        plotted_hist[thisbin]->Fit(fitname[thisbin].c_str(),"R,q");
        plotted_hist[thisbin]->SetTitle(histtitle[thisbin].c_str());
        plotted_hist[thisbin]->GetXaxis()->SetTitle("Hit time offset [ns]");
        plotted_hist[thisbin]->GetYaxis()->SetTitle("Number of events");
        float tmpmean = plotted_hist[thisbin]->GetMean();
        float tmprms  = plotted_hist[thisbin]->GetRMS();
        plotted_hist[thisbin]->GetXaxis()->SetRangeUser(val_residPeak-50,val_residPeak+50);
        plotted_hist[thisbin]->GetYaxis()->SetRangeUser(0,1.5*plotted_hist[thisbin]->GetMaximum());
        plotted_hist[thisbin]->GetXaxis()->SetTitleOffset(1.3);
        plotted_hist[thisbin]->GetYaxis()->SetTitleOffset(1.6);
      }
      stringstream tempx;
      tempx << run_id << "_angular_singlePMTs.png";
      c0->SaveAs(tempx.str().c_str());
      //c0->Print((sglstr+".pdf").c_str());
      c0->Close();
    }
    if (c0) delete c0;
    saved_pmt.clear();
    for (size_t i=0; i<24; i++) {
      delete plotted_hist[i];
      delete plotted_fit[i];
    }
    plotted_hist.clear();
    plotted_fit.clear();

    // Return 2D graph with fit results (by reconstructing object at given address)
    result = TGraph2DErrors(NPMTS,&x[0],&y[0],&z[0],&ex[0],&ey[0],&ez[0]);

    return;
  }

  // -----------------------------------------------------------------------------
  /// Display progress bar within a loop (can't have any other output in loop!)
  void MyUserProc::printProgress(int it, int n) {
    int div = (n - (n % 100))/100;              // 1% division
    if (it % div != 0 && it != n-1) return;
    float prog = (float)it/n;
    int barWidth = 70;
    cout << "[";
    int pos = barWidth * prog;
    for (int i=0; i<barWidth; ++i) {
      if (i < pos)                cout << "=";  // processed
      else if (pos+1 == barWidth) cout << "=";  // reached 100%
      else if (i == pos)          cout << ">";  // processing
      else                        cout << " ";  // not yet processed
    }
    int perc = (int)round(100.*prog);
    if (it < n-1) cout << "] " << perc << "%\r" << flush;
    else          cout << "] " << perc << "%\r" << endl;
    return;
  }
