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

//RAT includes
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
#include <RAT/json.hh>
#include <RAT/DBLink.hh>
#include <RAT/DBTable.hh>
#include <RAT/BitManip.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DataQualityProc.hh>

//ROOT includes
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
#include <TGaxis.h>
#include <TGraph2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TVector2.h>
#include <TTree.h>
#include <TSpectrum.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TKey.h>
#include <TImage.h>
#include <TView3D.h>
#include <TPaveText.h>
#include <TPaveStats.h>

//C++ includes
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <Riostream.h>
#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>

using namespace std;
using namespace RAT;

// Class declaration
// (if this gets big, you should move it to its own header file)
class MyUserProc : public Processor {
public:

  MyUserProc();

  virtual ~MyUserProc();

  virtual void BeginOfRun( DS::Run& run );

  /// @param[in] run DS::Run location in data structure
  virtual void EndOfRun( DS::Run& run );

  virtual Processor::Result DSEvent(DS::Run& run, DS::Entry& ds);

  virtual void SetS( const std::string& param, const std::string& value );

protected:

  std::string GraphType;

  // RATDB, Utils
  //RAT::DB *db = RAT::DB::Get();
  DBLinkPtr dblink;
  RAT::DataQualityProc *dq;

  // RAT stuff for time of flight
  RAT::DU::GroupVelocity gv;
  RAT::DU::LightPathCalculator lpc;

  // Bucket times
  TVector3 endDir;
  double thetaAtPMT;
  double lightBucketTime;

  // Laserball
  double x0, y0, z0;
  TVector3 lbpos;

  int NEvents;
  int nhits;
  int nhits_pass;

  // Events
  int trig;
  double energy;
  double fLEDWavelength;
  double LOCALITY;

  // Define ECA and PCA masks and corresponding labels
  unsigned int fECAMask;
  unsigned int fPCAMask;

  // data
  std::vector< std::vector < std::vector <double> > >  fPMTs;
  int NPMTs;
  TVector3 pmtPosAll;
  TVector2 flat;
  int pmt_id;

  // Gaussians + Occup
  double occup;
  double maxValue;
  double PMTHitTimesMean;
  double PMTHitTimesRMS;
  string temp_name;
  string temp_name2;
  TF1 * PMTHitTimesPeakFit;
  TF1 * PMTResidualTimesPeakFit;
  double maxValueR;
  double PMTHitTimesMeanR;
  double PMTHitTimesRMSR;
  double PMTHitTimesMedianR;
  string temp_nameR;
  string temp_name2R;
  TF1 * PMTHitTimesPeakFitR;
  TF1 * PMTResidualTimesPeakFitR;

  double minTime;
  double maxTime;
  double minTimeR;
  double maxTimeR;

  // ROOT file and plots
  std::string FileName;
  TFile *File;

  TCanvas *fPMTcoverageC;
  TCanvas *fPMTcoverage2C;
  TCanvas *lPMTHitTimesC;
  TCanvas *lPMTHitTimesZOOMC;
  TCanvas *AllPMTHitTimesC;
  TCanvas *AllPMTGaussiansC;
  TCanvas *lPMTHitTimesRC;
  TCanvas *lPMTHitTimesZOOMRC;
  TCanvas *AllPMTHitTimesRC;
  TCanvas *AllPMTHitTimesRLogC;
  TCanvas *AllPMTHitTimesRLogC2;
  TCanvas *AllPMTHitTimesRMedianC;
  TCanvas *AllPMTHitTimesRLogMedianC;
  TCanvas *AllPMTHitTimesRLogMedianC2;
  TCanvas *AllPMTGaussiansRC;
  TCanvas *HistPMTGaussiansC;
  TCanvas *HistPMTGaussiansRC;
  TCanvas *HistPMTGaussiansLogC;
  TCanvas *HistPMTGaussiansRLogC;
  TCanvas *AllPMTHitTimesRMedian2C;
  TCanvas *AllPMTGaussiansC2;
  TCanvas *AllPMTGaussiansRC2;
  TCanvas *AllPMTGaussiansRC3;
  TCanvas *AllPMTMedianRC;
  TCanvas *AllPMTMedianRC2;
  TCanvas *PMTHitsHistC;
  TCanvas *PMTHits2DC;

  TGraph2D *fPMTcoverage;
  TH1F *lPMTHitTimes;
  TH1F *lPMTHitTimesZOOM;
  TH1F *AllPMTHitTimes;
  TGraphErrors *AllPMTGaussians;
  TH1F *lPMTHitTimesR;
  TH1F *lPMTHitTimesZOOMR;
  TH1F *AllPMTHitTimesR;
  TH1F *AllPMTHitTimesRLog;
  TH1F *AllPMTHitTimesRLog2;
  TH1F *AllPMTHitTimesRMedian;
  TH1F *AllPMTHitTimesRLogMedian;
  TH1F *AllPMTHitTimesRLogMedian2;
  TGraphErrors *AllPMTGaussiansR;
  TH1F *HistPMTGaussians;
  TH1F *HistPMTGaussiansR;
  TH1F *AllPMTHitTimesRMedian2;
  TGraphErrors *AllPMTGaussians2;
  TGraphErrors *AllPMTGaussiansR2;
  TGraphErrors *AllPMTMedianR;
  TGraphErrors *AllPMTMedianR2;
  TH1F *PMTHitsHist;
  TGraph *PMTHits2D;

  // Final global plots
  std::vector<double> HitTimesID;
  std::vector<double> HitTimesMean;
  std::vector<double> HitTimesRMS;
  std::vector<double> HitTimesMeanR;
  std::vector<double> HitTimesRMSR;
  std::vector<double> HitTimesMedianR;

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

  }

  MyUserProc::~MyUserProc() {
  }

  void MyUserProc::BeginOfRun( DS::Run& run ) {

    // Initialise RAT vars
    gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    gv.BeginOfRun();
    lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lpc.BeginOfRun();
    lpc.SetELLIEEvent(true);
    dq = new RAT::DataQualityProc("");
    RAT::DS::BitMask PMTCalStatus = RAT::DS::BitMask();
    const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();

    // other vars
    x0 = 0.0, y0 = 0.0, z0 = 0.0;
    NEvents = 0;
    nhits = 0;
    nhits_pass = 0;
    fECAMask = 0xB00008E;
    fPCAMask = 0x3C;
    minTime = 1000;
    maxTime = 0;
    minTimeR = 1000;
    maxTimeR = 0;

    // Laserball position
    dblink = RAT::DB::Get()->GetLink("CALIB_COMMON_RUN_LEVEL","MANIP");

    try {
     std::vector <double> lb_db = dblink->GetDArray("position");
     x0 = lb_db[0];
     y0 = lb_db[1];
     z0 = lb_db[2];
     lbpos.SetX(x0);
     lbpos.SetY(y0);
     lbpos.SetZ(z0);
     cout << "LB position: " << x0 << ", " << y0 << ", " << z0 << endl;
    } catch (...) {
     std::cerr << "Failed to extract LB position !!! " << std::endl;

     // set to default for now
     lbpos.SetX(0.0);
     lbpos.SetY(-254.0);
     lbpos.SetZ(644.0);
    }
    cout << "LB position set as: x=" << lbpos.X() << ", y=" << lbpos.Y() << ", z=" << lbpos.Z() << endl;

    // PMT data
    NPMTs = pmtinfo.GetCount();
    std::cout << "There is " << NPMTs << " PMTs!" << endl;
    fPMTs.resize(NPMTs);

    for (int h = 0; h < NPMTs; h++) {
      fPMTs[h].resize(9);
      pmtPosAll = pmtinfo.GetPosition(h);
      int face;
      flat = dq->IcosProject(pmtPosAll, face);
      for (int g = 0; g < 9; g++){
        fPMTs[h][g].resize(1);
      }
      fPMTs[h][0][0] = 0;            // ID
      fPMTs[h][1][0] = flat.X();     // x pos
      fPMTs[h][2][0] = flat.Y();     // y pos
      fPMTs[h][3][0] = 0;            // direct hits array
      fPMTs[h][4][0] = 0;            // residual time array
      fPMTs[h][5][0] = 0;            // direct gaussian mean
      fPMTs[h][6][0] = 0;            // direct gaussian rms
      fPMTs[h][7][0] = 0;            // residual gaussian mean
      fPMTs[h][8][0] = 0;            // residual gaussian rms
    }

    //ROOT stuff
    const char * FileNameC = FileName.c_str();
    File = new TFile(FileNameC, "RECREATE");

    fPMTcoverageC = new TCanvas ("fPMTcoverage", "fPMTcoverage", 1024, 768);
    fPMTcoverage2C = new TCanvas ("fPMTcoverage2", "fPMTcoverage2", 1024, 768);
    AllPMTHitTimesC = new TCanvas ("AllPMTHitTimes", "AllPMTHitTimes", 1024, 768);
    AllPMTGaussiansC = new TCanvas ("AllPMTGaussians", "AllPMTGaussians", 1024, 768);
    AllPMTHitTimesRC = new TCanvas ("AllPMTHitTimesR", "AllPMTHitTimesR", 1024, 768);
    AllPMTHitTimesRLogC = new TCanvas ("AllPMTHitTimesRLog", "AllPMTHitTimesRLog", 1024, 768);
    AllPMTHitTimesRLogC2 = new TCanvas ("AllPMTHitTimesRLog2", "AllPMTHitTimesRLog2", 1024, 768);
    AllPMTHitTimesRMedianC = new TCanvas ("AllPMTHitTimesRMedian", "AllPMTHitTimesRMedian", 1024, 768);
    AllPMTHitTimesRMedian2C = new TCanvas ("AllPMTHitTimesRMedian2", "AllPMTHitTimesRMedian2", 1024, 768);
    AllPMTHitTimesRLogMedianC = new TCanvas ("AllPMTHitTimesRLogMedian", "AllPMTHitTimesRLogMedian", 1024, 768);
    AllPMTHitTimesRLogMedianC2 = new TCanvas ("AllPMTHitTimesRLogMedian2", "AllPMTHitTimesRLogMedian2", 1024, 768);
    AllPMTGaussiansRC = new TCanvas ("AllPMTGaussiansR", "AllPMTGaussiansR", 1024, 768);
    HistPMTGaussiansC = new TCanvas ("HistPMTGaussians", "HistPMTGaussians", 1024, 768);
    HistPMTGaussiansRC = new TCanvas ("HistPMTGaussiansR", "HistPMTGaussiansR", 1024, 768);
    HistPMTGaussiansLogC = new TCanvas ("HistPMTGaussiansLogC", "HistPMTGaussiansLogC", 1024, 768);
    HistPMTGaussiansRLogC = new TCanvas ("HistPMTGaussiansRLogC", "HistPMTGaussiansRLogC", 1024, 768);
    AllPMTGaussiansC2 = new TCanvas ("AllPMTGaussians2", "AllPMTGaussians2", 1024, 768);
    AllPMTGaussiansRC2 = new TCanvas ("AllPMTGaussiansR2", "AllPMTGaussiansR2", 1024, 768);
    AllPMTMedianRC = new TCanvas ("AllPMTMedianR", "AllPMTMedianR", 1024, 768);
    AllPMTMedianRC2 = new TCanvas ("AllPMTMedianRC", "AllPMTMedianRC", 1024, 768);
    PMTHitsHistC = new TCanvas ("PMTHitsHist", "PMTHitsHist", 1024, 768);
    PMTHits2DC = new TCanvas ("PMTHits2D", "PMTHits2D", 1024, 768);

    fPMTcoverage = new TGraph2D();
    AllPMTGaussians = new TGraphErrors();
    AllPMTGaussiansR = new TGraphErrors();
    AllPMTGaussians2 = new TGraphErrors();
    AllPMTGaussiansR2 = new TGraphErrors();
    AllPMTMedianR = new TGraphErrors();
    AllPMTMedianR2 = new TGraphErrors();
    PMTHits2D = new TGraph();

    gStyle->SetOptFit();

    cout << "Looping through events ..." << endl;
  }

  Processor::Result MyUserProc::DSEvent(DS::Run&, DS::Entry& ds) {

    for(int iEv=0; iEv<ds.GetEVCount(); iEv++) {
        const RAT::DS::EV& ev = ds.GetEV(iEv);

        // Trigger type
        trig = ev.GetTrigType();
        if (!(trig & 0x8000)){ continue; }    // EXT trigger only

        const RAT::DS::CalPMTs& calibratedPMTs = ev.GetCalPMTs();
        for( int ipmt = 0; ipmt < calibratedPMTs.GetNormalCount(); ipmt++ ) {

          const RAT::DS::PMTCal& calPMT = calibratedPMTs.GetNormalPMT(ipmt);

          // Get ECA and PCA statuses
          unsigned int eca_status = (unsigned int)calPMT.GetStatus().GetBits(0,32);
          unsigned int pca_status = (unsigned int)calPMT.GetStatus().GetBits(32,32);

          nhits ++;

          int pmtID = calPMT.GetID();
          const RAT::DU::PMTCalStatus& pmtStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
          const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();
          const RAT::DU::PMTInfo& pmtinfo_loop = RAT::DU::Utility::Get()->GetPMTInfo();
          unsigned int status = pmtStatus.GetHitStatus(calPMT);
          if(status & (1<<pmtStatus.kCHSBit)){ continue; }
          if(status & (1<<pmtStatus.kECABit)){ continue; }
          if(status & (1<<pmtStatus.kPCABit)){ continue; }
          if(status & (1<<pmtStatus.kXTalkBit)){ continue; }
          if ( !chs.IsTubeOnline(pmtID) ){ continue; }
          if ( !chs.IsEnabled() ){ continue; }
          if ( !chs.IsChannelOnline(pmtID) ){ continue; }
          if ( !chs.IsDAQEnabled(pmtID) ){ continue; }
          if ( pmtinfo_loop.GetType(pmtID) != 1 ){ continue; }

          nhits_pass ++;

          // Get time of flight info
          const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
          TVector3 pmtpos = pmtinfo.GetPosition(calPMT.GetID());
          TVector3 pmtdir = pmtinfo.GetDirection(calPMT.GetID());
          if (pmtpos.Mag()==0) { continue; }
          int face;
          flat = dq->IcosProject(pmtpos, face);
          lpc.CalcByPosition( lbpos, pmtpos, energy, LOCALITY);	// water

          if (lpc.GetTIR() == 1) { continue; }            // total internal reflection
          if (lpc.GetPathValid() == 0) { continue; }      // check whether path parameters are valid
          if (lpc.GetResvHit() == 1) { continue; }        // whether end point was within locality

          //lpc.CalcByPositionPartial( lbpos, pmtpos);	// partial
          double distInScint = lpc.GetDistInInnerAV();	// water
          //double distInScint = lpc.GetDistInUpperTarget() + lpc.GetDistInLowerTarget();	// partial
          double distInAV = lpc.GetDistInAV();
          double distInWater = lpc.GetDistInWater();
          double tof = gv.CalcByDistance( distInScint, distInAV, distInWater );

          double hitTime = calPMT.GetTime();

          if (hitTime <= 0){continue;}

          // Bucket time correction
          endDir = lpc.GetIncidentVecOnPMT();                   // end direction at PMT
          thetaAtPMT = endDir.Angle(pmtdir)*180./M_PIl;         // incident angle with bucket face
          lightBucketTime = gv.PMTBucketTime(thetaAtPMT);

          double pmtTimeRes = hitTime - tof - lightBucketTime;

          //cout << "hit: " << hitTime << ", rhit: " << pmtTimeRes << endl;
          if (hitTime < minTime){minTime = hitTime;}
          if (hitTime > maxTime){maxTime = hitTime;}
          if (pmtTimeRes < minTimeR){minTimeR = pmtTimeRes;}
          if (pmtTimeRes > maxTimeR){maxTimeR = pmtTimeRes;}
          //cout << "min: " << minTime << ", max: " << maxTime << ", rmin: " << minTimeR << ", rmax: " << maxTimeR << endl;

          //store data
          for (int lPMT = 0; lPMT < fPMTs.size(); lPMT++){
            if (flat.X() == fPMTs[lPMT][1][0]){
              fPMTs[lPMT][0][0] = pmtID;
              fPMTs[lPMT][3].push_back(hitTime);
              fPMTs[lPMT][4].push_back(pmtTimeRes);
            }   // position check
          }     // lPMT loop
        }       // calibrated PMTs loop
        NEvents++;
    }
    return Processor::OK;
  }

  void MyUserProc::EndOfRun( DS::Run& run ){

    cout << "Events: " << NEvents << endl;
    cout << "Hits: " << nhits << endl;
    cout << "Hits that pass checks: " << nhits_pass << endl;
    cout << "Hit times MIN: " << minTime << ", MAX: " << maxTime << endl;
    cout << "Residual hit times MIN: " << minTimeR << ", MAX: " << maxTimeR << endl;

    double range = (maxTimeR-minTimeR)/2;
    cout << "range: " << range << endl;

    // for now change to 250
    range = 300;

    AllPMTHitTimesRLog = new TH1F("AllPMTHitTimesRLog", "", 100, minTimeR, maxTimeR);
    AllPMTHitTimes = new TH1F("AllPMTHitTimes", "", 100, minTime, maxTime);
    AllPMTHitTimesR = new TH1F("AllPMTHitTimesR", "", 100, minTime, maxTime);
    PMTHitsHist = new TH1F("PMTHitsHistC", "", 100, 0, 33500);
    HistPMTGaussians = new TH1F("HistPMTGaussians", "", 100, range+75-50, range+75+50);
    HistPMTGaussiansR = new TH1F("HistPMTGaussiansR", "", 100, range-50, range+50);

    // ?
    AllPMTHitTimesRLogMedian = new TH1F("AllPMTHitTimesRLogMedian", "", 100, range-15, range+15);
    AllPMTHitTimesRLogMedian2 = new TH1F("AllPMTHitTimesRLogMedian2", "", 100, range-15, range+15);
    AllPMTHitTimesRMedian = new TH1F("AllPMTHitTimesRMedian", "", 100, range-15, range+15);
    AllPMTHitTimesRMedian2 = new TH1F("AllPMTHitTimesRMedian2", "", 100, range-15, range+15);

    // Analysing + plots
    cout << "Analysing stored data ..." << endl;

    File->cd();
    for (int lPMT = 0; lPMT < fPMTs.size(); lPMT++){
      if (fPMTs[lPMT][0][0] != 0){            // check that ID is not 0
        cout << lPMT << " " << fPMTs[lPMT][3].size() << endl;
        if (fPMTs[lPMT][3].size() > 25){       // real pmt with hits (25 is arbitrary here, had issue with PMT with 1 entry [421])

          //erase first 0 element
          fPMTs[lPMT][3].erase(fPMTs[lPMT][3].begin());
          fPMTs[lPMT][4].erase(fPMTs[lPMT][4].begin());

          maxValue = 0.0;
          PMTHitTimesMean = 0.0;
          PMTHitTimesRMS = 0.0;
          maxValueR = 0.0;
          PMTHitTimesMeanR = 0.0;
          PMTHitTimesRMSR = 0.0;
          //lPMTHitTimes->Reset();
          //lPMTHitTimesZOOM->Reset();

          // hit times plot setup
          ostringstream oss;
          oss.str("");
          oss << "PMT-" << fPMTs[lPMT][0][0];
          temp_name = oss.str();
          oss << "-Log";
          temp_name2 = oss.str();
          lPMTHitTimesC = new TCanvas (temp_name.c_str(), temp_name.c_str(), 1024, 768);
          lPMTHitTimesZOOMC = new TCanvas (temp_name2.c_str(), temp_name2.c_str(), 1024, 768);
          lPMTHitTimes = new TH1F ( temp_name.c_str(), temp_name.c_str(), 50, minTime-10, maxTime+10 );

          ostringstream ossR;
          ossR.str("");
          ossR << "PMT-" << fPMTs[lPMT][0][0];
          ossR << "-Residual";
          temp_nameR = ossR.str();
          ossR << "-Log";
          temp_name2R = ossR.str();
          lPMTHitTimesR = new TH1F ( temp_nameR.c_str(), temp_nameR.c_str(), 50, minTimeR-10, maxTimeR+10 );
          lPMTHitTimesRC = new TCanvas (temp_nameR.c_str(), temp_nameR.c_str(), 1024, 768);
          lPMTHitTimesZOOMRC = new TCanvas (temp_name2R.c_str(), temp_name2R.c_str(), 1024, 768);

          for (int ltimes = 0; ltimes < fPMTs[lPMT][3].size(); ltimes++){
            lPMTHitTimes->Fill(fPMTs[lPMT][3][ltimes]);
            AllPMTHitTimes->Fill(fPMTs[lPMT][3][ltimes]);
            lPMTHitTimesR->Fill(fPMTs[lPMT][4][ltimes]);
            AllPMTHitTimesR->Fill(fPMTs[lPMT][4][ltimes]);
            AllPMTHitTimesRLog->Fill(fPMTs[lPMT][4][ltimes]);
          }

          // get peaks here for better peak fit
          int peakBin = lPMTHitTimes->GetMaximumBin();
          double peakBinVal = lPMTHitTimes->GetXaxis()->GetBinCenter(peakBin);
          int peakBinR = lPMTHitTimesR->GetMaximumBin();
          double peakBinValR = lPMTHitTimesR->GetXaxis()->GetBinCenter(peakBinR);
          //cout << "peak bin val: " << peakBinVal << endl;
          //cout << "peak bin val R: " << peakBinValR << endl;
          lPMTHitTimesZOOMR = new TH1F ( temp_name2R.c_str(), temp_name2R.c_str(), 100, peakBinValR-15, peakBinValR+15 );
          lPMTHitTimesZOOM = new TH1F ( temp_name2.c_str(), temp_name2.c_str(), 100, peakBinVal-15, peakBinVal+15 );
          for (int ltimes = 0; ltimes < fPMTs[lPMT][3].size(); ltimes++){
            lPMTHitTimesZOOM->Fill(fPMTs[lPMT][3][ltimes]);
            lPMTHitTimesZOOMR->Fill(fPMTs[lPMT][4][ltimes]);
          }

          // plot times for this PMT
          lPMTHitTimesC->cd();
          lPMTHitTimesC->Clear();
          lPMTHitTimesC->Update();
          lPMTHitTimes->GetXaxis()->SetTitle("hit time [ns]");
          lPMTHitTimes->GetYaxis()->SetTitle("count");
          lPMTHitTimes->Draw();
          lPMTHitTimesC->Write();
          lPMTHitTimesC->Close();

          lPMTHitTimesRC->cd();
          lPMTHitTimesRC->Clear();
          lPMTHitTimesRC->Update();
          lPMTHitTimesR->GetXaxis()->SetTitle("residual time [ns]");
          lPMTHitTimesR->GetYaxis()->SetTitle("count");
          lPMTHitTimesR->Draw();
          lPMTHitTimesRC->Write();
          lPMTHitTimesRC->Close();

          // gaussian fit here
          cout << "PMT: " << lPMT << ", hits: " << fPMTs[lPMT][3].size() << endl;
          maxValue = lPMTHitTimesZOOM->GetMaximum();      // get max
          lPMTHitTimesZOOM->SetMinimum(0.1 * maxValue);   // get rid of low 10%
          lPMTHitTimesZOOM->Fit("gaus");                  // get mean and rms
          PMTHitTimesPeakFit = lPMTHitTimesZOOM->GetFunction("gaus");
          PMTHitTimesMean = PMTHitTimesPeakFit->GetParameter(1);
          PMTHitTimesRMS = PMTHitTimesPeakFit->GetParameter(2);
          fPMTs[lPMT][5][0] = PMTHitTimesMean;
          fPMTs[lPMT][6][0] = PMTHitTimesRMS;
          HistPMTGaussians->Fill(PMTHitTimesMean);

          lPMTHitTimesZOOMC->cd();
          lPMTHitTimesZOOMC->Clear();
          lPMTHitTimesZOOMC->Update();
          lPMTHitTimesZOOMC->SetLogy();
          lPMTHitTimesZOOM->GetXaxis()->SetTitle("hit time [ns]");
          lPMTHitTimesZOOM->GetYaxis()->SetTitle("count");
          lPMTHitTimesZOOM->Draw();
          lPMTHitTimesZOOMC->Write();
          lPMTHitTimesZOOMC->Close();

          HitTimesID.push_back(fPMTs[lPMT][0][0]);
          HitTimesMean.push_back(fPMTs[lPMT][5][0]);
          HitTimesRMS.push_back(fPMTs[lPMT][6][0]);

          // residuals
          maxValueR = lPMTHitTimesZOOMR->GetMaximum();      // get max
          lPMTHitTimesZOOMR->SetMinimum(0.1 * maxValueR);   // get rid of low 10%
          lPMTHitTimesZOOMR->Fit("gaus");                  // get mean and rms
          PMTHitTimesPeakFitR = lPMTHitTimesZOOMR->GetFunction("gaus");
          PMTHitTimesMeanR = PMTHitTimesPeakFitR->GetParameter(1);
          PMTHitTimesRMSR = PMTHitTimesPeakFitR->GetParameter(2);
          fPMTs[lPMT][7][0] = PMTHitTimesMeanR;
          fPMTs[lPMT][8][0] = PMTHitTimesRMSR;
          HistPMTGaussiansR->Fill(PMTHitTimesMeanR);
          cout << "gauss: " << PMTHitTimesMeanR << endl;

          lPMTHitTimesZOOMRC->cd();
          lPMTHitTimesZOOMRC->Clear();
          lPMTHitTimesZOOMRC->Update();
          lPMTHitTimesZOOMRC->SetLogy();
          lPMTHitTimesZOOMR->GetXaxis()->SetTitle("residual time [ns]");
          lPMTHitTimesZOOMR->GetYaxis()->SetTitle("count");
          lPMTHitTimesZOOMR->Draw();
          lPMTHitTimesZOOMRC->Write();
          lPMTHitTimesZOOMRC->Close();

          HitTimesMeanR.push_back(fPMTs[lPMT][7][0]);
          HitTimesRMSR.push_back(fPMTs[lPMT][8][0]);

          Double_t x, q;
          q = 0.5; // 0.5 for "median"
          lPMTHitTimesZOOMR->ComputeIntegral(); // just a precaution
          lPMTHitTimesZOOMR->GetQuantiles(1, &x, &q);
          std::cout << "median: " << x << std::endl;
          HitTimesMedianR.push_back(x);

          AllPMTHitTimesRMedian->Fill(x);
          AllPMTHitTimesRLogMedian->Fill(x);
          AllPMTHitTimesRLogMedian2->Fill(x);
          AllPMTHitTimesRMedian2->Fill(x);

          // hits per PMT plots
          int PMTHits = fPMTs[lPMT][3].size();
          PMTHitsHist->Fill(PMTHits);
          PMTHits2D->SetPoint(fPMTs[lPMT][0][0], fPMTs[lPMT][0][0], PMTHits);

        } // real hits (occup)

        // coverage
        occup = ((double)fPMTs[lPMT][3].size()/(double)NEvents) * 100;
        fPMTcoverage->SetPoint(fPMTs[lPMT][0][0], fPMTs[lPMT][1][0], fPMTs[lPMT][2][0], occup );
      } // ID check
    } // PMT loop

    // Storing plots
    cout << "Storing plots ..." << endl;
    //File->cd();

    fPMTcoverageC->cd();
    fPMTcoverage->SetTitle("fPMTcoverage");
    fPMTcoverage->SetMarkerStyle(8);
    fPMTcoverage->SetMarkerSize(1.0);
    fPMTcoverage->Draw("AP");
    fPMTcoverageC->Write();
    fPMTcoverageC->Close();

    fPMTcoverage2C->cd();
    //TGaxis::SetMaxDigits(4);
    fPMTcoverage2C->SetMargin(0.05,0.15,0.05,0.05);
    fPMTcoverage->SetMarkerStyle(8);
    fPMTcoverage->SetMarkerSize(1.0);
    //gStyle->SetLabelSize(0.02);
    fPMTcoverage->GetXaxis()->SetLabelSize(0.03);
    fPMTcoverage->GetYaxis()->SetLabelSize(0.03);
    fPMTcoverage->SetTitle("fPMTcoverage_final");
    fPMTcoverage->Draw("AZCOLPCOL");
    fPMTcoverage2C->Write();
    fPMTcoverage2C->Close();

    AllPMTHitTimesC->cd();
    AllPMTHitTimes->SetTitle("AllPMTHitTimes");
    AllPMTHitTimes->GetXaxis()->SetTitle("hit time [ns]");
    AllPMTHitTimes->GetYaxis()->SetTitle("count");
    AllPMTHitTimes->Draw();
    AllPMTHitTimesC->Write();
    AllPMTHitTimesC->Close();

    AllPMTHitTimesRC->cd();
    AllPMTHitTimesR->SetTitle("AllPMTHitTimesR");
    AllPMTHitTimesR->GetXaxis()->SetTitle("residual time [ns]");
    AllPMTHitTimesR->GetYaxis()->SetTitle("count");
    AllPMTHitTimesR->Draw();
    AllPMTHitTimesRC->Write();
    AllPMTHitTimesRC->Close();

    AllPMTHitTimesRLogC->cd();
    AllPMTHitTimesRLogC->SetLogy();
    AllPMTHitTimesRLog->SetTitle("AllPMTHitTimesRLog");
    AllPMTHitTimesRLog->GetXaxis()->SetTitle("residual time [ns]");
    AllPMTHitTimesRLog->GetYaxis()->SetTitle("count");
    AllPMTHitTimesRLog->Draw();
    AllPMTHitTimesRLogC->Write();
    AllPMTHitTimesRLogC->Close();

    // get prompt peak and edges here
    int binmax = AllPMTHitTimesRLog->GetMaximumBin();
    // scan to left, get first increasing
    double previousB = 1;
    double thisB = 0;
    int i = 1;
    while (previousB > thisB){
      thisB = AllPMTHitTimesRLog->GetBinContent(binmax-i);
      previousB = AllPMTHitTimesRLog->GetBinContent(binmax-i+1);
      i++;
    }
    double lowerBound = AllPMTHitTimesRLog->GetXaxis()->GetBinCenter(binmax-i+1);
    // scan to right, get first increasing
    double previousB2 = 1;
    double thisB2 = 0;
    int ii = 1;
    while (previousB2 > thisB2){
      thisB2 = AllPMTHitTimesRLog->GetBinContent(binmax+ii);
      previousB2 = AllPMTHitTimesRLog->GetBinContent(binmax+ii-1);
      ii++;
    }
    double upperBound = AllPMTHitTimesRLog->GetXaxis()->GetBinCenter(binmax+ii-1);
    cout << "prompt time peak range: " << lowerBound << " - " << upperBound << endl;

    //AllPMTHitTimesRLog2 = new TH1F("AllPMTHitTimesRLog2", "", 100, lowerBound+14, upperBound-14);
    AllPMTHitTimesRLog2 = new TH1F("AllPMTHitTimesRLog2", "", 100, 287, 307);
    for (int lPMT = 0; lPMT < fPMTs.size(); lPMT++){
      if (fPMTs[lPMT][0][0] != 0){
        if (fPMTs[lPMT][3].size() > 25){
          for (int ltimes = 0; ltimes < fPMTs[lPMT][3].size(); ltimes++){
            AllPMTHitTimesRLog2->Fill(fPMTs[lPMT][4][ltimes]);
          }
        }
      }
    }

    AllPMTHitTimesRLogC2->cd();
    AllPMTHitTimesRLogC2->SetLogy();
    AllPMTHitTimesRLog2->SetTitle("AllPMTHitTimesRLog2");
    AllPMTHitTimesRLog2->GetXaxis()->SetTitle("residual time [ns]");
    AllPMTHitTimesRLog2->GetYaxis()->SetTitle("count");
    AllPMTHitTimesRLog2->Fit("gaus");
    AllPMTHitTimesRLog2->Draw();
    AllPMTHitTimesRLogC2->Write();
    AllPMTHitTimesRLogC2->Close();

    AllPMTHitTimesRMedianC->cd();
    AllPMTHitTimesRMedian->SetTitle("AllPMTHitTimesRMedian");
    AllPMTHitTimesRMedian->GetXaxis()->SetTitle("median time [ns]");
    AllPMTHitTimesRMedian->GetYaxis()->SetTitle("count");
    AllPMTHitTimesRMedian->Draw();
    AllPMTHitTimesRMedianC->Write();
    AllPMTHitTimesRMedianC->Close();

    AllPMTHitTimesRMedian2C->cd();
    AllPMTHitTimesRMedian2->SetTitle("AllPMTHitTimesRMedian2");
    AllPMTHitTimesRMedian2->GetXaxis()->SetTitle("median time [ns]");
    AllPMTHitTimesRMedian2->GetYaxis()->SetTitle("count");
    AllPMTHitTimesRMedian2->Fit("gaus");
    AllPMTHitTimesRMedian2->Draw();
    AllPMTHitTimesRMedian2C->Write();
    AllPMTHitTimesRMedian2C->Close();

    AllPMTHitTimesRLogMedianC->cd();
    AllPMTHitTimesRLogMedianC->SetLogy();
    AllPMTHitTimesRLogMedian->SetTitle("AllPMTHitTimesRLogMedian");
    AllPMTHitTimesRLogMedian->GetXaxis()->SetTitle("median time [ns]");
    AllPMTHitTimesRLogMedian->GetYaxis()->SetTitle("count");
    AllPMTHitTimesRLogMedian->Draw();
    AllPMTHitTimesRLogMedianC->Write();
    AllPMTHitTimesRLogMedianC->Close();

    AllPMTHitTimesRLogMedianC2->cd();
    AllPMTHitTimesRLogMedianC2->SetLogy();
    AllPMTHitTimesRLogMedian2->SetTitle("AllPMTHitTimesRLogMedian2");
    AllPMTHitTimesRLogMedian2->GetXaxis()->SetTitle("median time [ns]");
    AllPMTHitTimesRLogMedian2->GetYaxis()->SetTitle("count");
    AllPMTHitTimesRLogMedian2->Fit("gaus");
    AllPMTHitTimesRLogMedian2->Draw();
    AllPMTHitTimesRLogMedianC2->Write();
    AllPMTHitTimesRLogMedianC2->Close();

    HistPMTGaussiansC->cd();
    HistPMTGaussians->SetTitle("HistPMTGaussians");
    HistPMTGaussians->GetXaxis()->SetTitle("hit time [ns]");
    HistPMTGaussians->GetYaxis()->SetTitle("count");
    HistPMTGaussians->Fit("gaus");
    HistPMTGaussians->Draw();
    HistPMTGaussiansC->Write();
    HistPMTGaussiansC->Close();

    HistPMTGaussiansRC->cd();
    HistPMTGaussiansR->SetTitle("HistPMTGaussiansR");
    HistPMTGaussiansR->GetXaxis()->SetTitle("residual time [ns]");
    HistPMTGaussiansR->GetYaxis()->SetTitle("count");
    HistPMTGaussiansR->Fit("gaus");
    HistPMTGaussiansR->Draw();
    HistPMTGaussiansRC->Write();
    HistPMTGaussiansRC->Close();

    HistPMTGaussiansLogC->cd();
    HistPMTGaussiansLogC->SetLogy();
    HistPMTGaussians->SetTitle("HistPMTGaussiansLog");
    HistPMTGaussians->GetXaxis()->SetTitle("hit time [ns]");
    HistPMTGaussians->GetYaxis()->SetTitle("count");
    HistPMTGaussians->Draw();
    HistPMTGaussiansLogC->Write();
    HistPMTGaussiansLogC->Close();

    HistPMTGaussiansRLogC->cd();
    HistPMTGaussiansRLogC->SetLogy();
    HistPMTGaussiansR->SetTitle("HistPMTGaussiansRLog");
    HistPMTGaussiansR->GetXaxis()->SetTitle("residual time [ns]");
    HistPMTGaussiansR->GetYaxis()->SetTitle("count");
    HistPMTGaussiansR->Fit("gaus");
    HistPMTGaussiansR->Draw();
    HistPMTGaussiansRLogC->Write();
    HistPMTGaussiansRLogC->Close();

    int No = HitTimesID.size();
    std::vector<double> x_errors (No, 0);
    double *HitTimesID_array = &HitTimesID[0];
    double *HitTimesMean_array = &HitTimesMean[0];
    double *x_errors_array = &x_errors[0];
    double *HitTimesRMS_array = &HitTimesRMS[0];
    AllPMTGaussians = new TGraphErrors( No, HitTimesID_array, HitTimesMean_array, x_errors_array, HitTimesRMS_array);
    AllPMTGaussiansC->cd();
    AllPMTGaussians->SetTitle("AllPMTGaussians");
    AllPMTGaussians->GetXaxis()->SetTitle("PMT ID");
    AllPMTGaussians->GetYaxis()->SetTitle("fitted hit time [ns]");
    AllPMTGaussians->Draw("A*");
    AllPMTGaussiansC->Write();
    AllPMTGaussiansC->Close();

    AllPMTGaussians2 = new TGraphErrors( No, HitTimesID_array, HitTimesMean_array, x_errors_array, x_errors_array);
    AllPMTGaussiansC2->cd();
    AllPMTGaussians2->SetTitle("AllPMTGaussians2");
    AllPMTGaussians2->GetXaxis()->SetTitle("PMT ID");
    AllPMTGaussians2->GetYaxis()->SetTitle("fitted hit time [ns]");
    AllPMTGaussians2->SetMarkerStyle(5);
    AllPMTGaussians2->Draw("AP");
    AllPMTGaussiansC2->Write();
    AllPMTGaussiansC2->Close();

    double *HitTimesMeanR_array = &HitTimesMeanR[0];
    double *HitTimesRMSR_array = &HitTimesRMSR[0];
    AllPMTGaussiansR = new TGraphErrors( No, HitTimesID_array, HitTimesMeanR_array, x_errors_array, HitTimesRMSR_array);
    AllPMTGaussiansRC->cd();
    AllPMTGaussiansR->SetTitle("AllPMTGaussiansR");
    AllPMTGaussiansR->GetXaxis()->SetTitle("PMT ID");
    AllPMTGaussiansR->GetYaxis()->SetTitle("fitted residual time [ns]");
    AllPMTGaussiansR->Draw("A*");
    AllPMTGaussiansRC->Write();
    AllPMTGaussiansRC->Close();

    AllPMTGaussiansR2 = new TGraphErrors( No, HitTimesID_array, HitTimesMeanR_array, x_errors_array, x_errors_array);
    AllPMTGaussiansRC2->cd();
    AllPMTGaussiansR2->SetTitle("AllPMTGaussiansR2");
    AllPMTGaussiansR2->GetXaxis()->SetTitle("PMT ID");
    AllPMTGaussiansR2->GetYaxis()->SetTitle("fitted residual time [ns]");
    AllPMTGaussiansR2->SetMarkerStyle(5);
    AllPMTGaussiansR2->Draw("AP");
    AllPMTGaussiansRC2->Write();
    AllPMTGaussiansRC2->Close();

    double *HitTimesMedianR_array = &HitTimesMedianR[0];
    AllPMTMedianR = new TGraphErrors( No, HitTimesID_array, HitTimesMedianR_array, x_errors_array, HitTimesRMSR_array);
    AllPMTMedianRC->cd();
    AllPMTMedianR->SetTitle("AllPMTMedianR");
    AllPMTMedianR->GetXaxis()->SetTitle("PMT ID");
    AllPMTMedianR->GetYaxis()->SetTitle("fitted residual time [ns]");
    AllPMTMedianR->Draw("A*");
    AllPMTMedianRC->Write();
    AllPMTMedianRC->Close();

    AllPMTMedianR2 = new TGraphErrors( No, HitTimesID_array, HitTimesMedianR_array, x_errors_array, x_errors_array);
    AllPMTMedianRC2->cd();
    AllPMTMedianR2->SetTitle("AllPMTMedianR2");
    AllPMTMedianR2->GetXaxis()->SetTitle("PMT ID");
    AllPMTMedianR2->GetYaxis()->SetTitle("fitted residual time [ns]");
    AllPMTMedianR2->SetMarkerStyle(5);
    AllPMTMedianR2->Draw("AP");
    AllPMTMedianRC2->Write();
    AllPMTMedianRC2->Close();

    PMTHitsHistC->cd();
    PMTHitsHist->SetTitle("PMTHitsHist");
    PMTHitsHist->GetXaxis()->SetTitle("PMT hits");
    PMTHitsHist->GetYaxis()->SetTitle("count");
    PMTHitsHist->Draw();
    PMTHitsHistC->Write();
    PMTHitsHistC->Close();

    PMTHits2DC->cd();
    PMTHits2D->SetTitle("PMTHits2D");
    PMTHits2D->GetXaxis()->SetTitle("PMT");
    PMTHits2D->GetYaxis()->SetTitle("Hits");
    PMTHits2D->SetMarkerStyle(5);
    PMTHits2D->Draw("AP");
    PMTHits2DC->Write();
    PMTHits2DC->Close();

    File->Close();

  }

  void MyUserProc::SetS( const std::string& param, const std::string& value ){
    if( param == "file" ){
      FileName = value;
      cout << "Output set to: " << FileName << endl;
    }
    else{
      throw ParamUnknown( param );
    }
  }
