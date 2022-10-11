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

  virtual Processor::Result DSEvent(DS::Run& run, DS::Entry& ds);

  virtual string printVector(const TVector3& v);

  virtual double AngularSystematic(double b, double angle);

protected:

  // Fit results
  string fits_folder;
  string fits_file_dir;
  string fits_file_ang;
  string delim;
  stringstream dir_fit_file;
  stringstream ang_fit_file;
  std::vector<string> fibre_ang_name;
  std::vector<double> fibre_ang_B;
  std::vector<double> fibre_ang_B_err;
  int temp_count;
  std::vector<string> fibre_dir_name;
  std::vector<double> fibre_dir_X;
  std::vector<double> fibre_dir_Y;
  std::vector<double> fibre_dir_Z;
  double AngB;

  // Database
  RAT::DBLinkPtr tellie_run_data;
  std::string fibre_db;
  json::Value sub_run_info;
  RAT::DBLinkPtr fibre_data;
  TVector3 fibrepos, fibredir, fitfibredir;

  // LPC
  double energy;
  double fLEDWavelength;
  double LOCALITY;

  // RAT stuff
  RAT::DU::GroupVelocity gv;
  RAT::DU::LightPathCalculator lpc;
  int NPMTS;

  // Event loop
  int count;
  int allEvs;
  int trig;
  std::vector<double> pmtOccup;
  std::vector<double> pmtAngle;
  std::vector< std::vector<double> > pmtResidTimes;
  double residMin;
  double residMax;
  int nearL;
  int farL;

  // Analysis
  int PMTBeamSpot;
  double maxPeak;
  int maxBin;
  double maxPeak2;
  int maxBin2;
  double lowHalfMax;
  double highHalfMax;
  double lowHalfMax2;
  double highHalfMax2;
  TF1 *peakFit;
  TF1 *peakFit2;
  TF1 *peakFitTwo;
  TF1 *peakFitTwo2;
  double pca_peak;
  double pca_peak_err;
  double pca_width;
  double pca_width_err;
  double pca_peak2;
  double pca_peak_err2;
  double pca_width2;
  double pca_width_err2;

  // ROOT
  TFile *rootfile;
  TH1D *ResidualHits;
  TH1D *ResidualHits2;

  // Global
  int run_id;

  // File for fits
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
    fits_folder = getenv("FITS_FOLDER");
    fits_file_ang = getenv("ANG_FIT_FILE");
    fits_file_dir = getenv("DIR_FIT_FILE");
    ang_fit_file << fits_file_ang;
    dir_fit_file << fits_file_dir;
    delim = "\t";
    LOCALITY = atoi(getenv("LOCALITY"));
    fLEDWavelength = 506.0*1e-6;
    energy = util::WavelengthToEnergy(fLEDWavelength);
    out = fopen("pcaoffset_fit.txt","a");

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

    /*
    cout << fibre_ang_name.size() << endl;
    cout << fibre_ang_B.size() << " " << fibre_ang_B_err.size() << endl;
    cout << fibre_dir_name.size() << endl;
    cout << fibre_dir_X.size() << " " << fibre_dir_Y.size() << " " << fibre_dir_Z.size() << endl;
    for (size_t i=0; i<fibre_ang_name.size(); i++) {
      cout << i << ", fibre: " << fibre_ang_name[i] << ", a: " << fibre_ang_B[i] << ", a_err: " << fibre_ang_B_err[i] << endl;
    }
    for (size_t i=0; i<fibre_dir_name.size(); i++) {
      cout << i << ", fibre: " << fibre_dir_name[i] << ", x: " << fibre_dir_X[i] << ", y: " << fibre_dir_Y[i] << ", z: " << fibre_dir_Z[i] << endl;
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
    residMin = 500;
    residMax = 0;
    PMTBeamSpot = 0;
    nearL = 0;
    farL = 0;
    lowHalfMax = 0;
    highHalfMax = 0;
    lowHalfMax2 = 0;
    highHalfMax2 = 0;
    AngB = 0.0;
    maxPeak = 0;
    maxBin = 0;
    maxPeak2 = 0;
    maxBin2 = 0;
    lowHalfMax = 0;
    highHalfMax = 0;
    lowHalfMax2 = 0;
    highHalfMax2 = 0;

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
    pmtResidTimes.clear();
    pmtOccup.resize(NPMTS);
    pmtAngle.resize(NPMTS);
    pmtResidTimes.resize(NPMTS);

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
    // Set ang sys parameter
    for (size_t i=0; i<fibre_ang_name.size(); i++) {
      if (fibre_ang_name[i] == fibre_db){
        AngB = fibre_ang_B[i];
      }
    }
    cout << "Set fibre as: fibre " << fibre_db << ", pos " << printVector(fibrepos) << ", dir " << printVector(fitfibredir) << endl;
    cout << "ang b parameter: " << AngB << endl;

    // Output file for stats
    logFile_namess.str("");
    logFile_namess << run_id << "_pca.log";
    logFile_name = logFile_namess.str();
    logFile = fopen(logFile_name.c_str(), "w");

  }

  void MyUserProc::EndOfRun( DS::Run& run ) {

    // Turn array of PMTs into percentage (occupancy)
    cout << "Number of events was " << count << endl;
    cout << "Number of EXTA events was " << allEvs << endl;
    cout << "Near to far hits perc: " << double(nearL)/double(farL)*100 << endl;
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      pmtOccup[iPMT] /= (float)allEvs;
      //if (pmtOccup[iPMT] > 0.01 && pmtOccup[iPMT] < 0.05) {cout << iPMT << ", occ: " << pmtOccup[iPMT] << endl;}
    }

    // Fill residual hits histogram
    ResidualHits = new TH1D ("ResidualHits", "", 10000, residMin-1, residMax+1);
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      if (pmtAngle[iPMT] <= 12.0){                                      // angular cut (beamspot only)
        if (pmtOccup[iPMT] >= 0.01 && pmtOccup[iPMT] <= 0.05){          // occupancy cut
          PMTBeamSpot++;
          for (int hit=0; hit<(pmtResidTimes[iPMT].size()); hit++){     // goes through valid hits here
            ResidualHits->Fill( pmtResidTimes[iPMT][hit] );
          }
        }
      }
    }
    cout << "PMT in beamspot with acceptable occupancy: " << PMTBeamSpot << endl;

    // Fits
    maxPeak = ResidualHits->GetMaximum();
    maxBin = ResidualHits->GetMaximumBin();
    //Fit a gaussian to the peak defined as the function between both half maximum points
    for(int i=1; i<maxBin; i++){
        if(ResidualHits->GetBinContent(i) > 0.2*maxPeak){     // cut data < 20% of peak
            lowHalfMax = ResidualHits->GetBinCenter(i-2);
            break;
        }
    }
    for(int i=maxBin; i<=ResidualHits->GetXaxis()->GetNbins(); i++){
        if(ResidualHits->GetBinContent(i) < 0.2*maxPeak){
            highHalfMax = ResidualHits->GetBinCenter(i);
            break;
        }
    }
    ResidualHits->Fit("gaus","","",lowHalfMax,highHalfMax);
    peakFit = ResidualHits->GetFunction("gaus");
    pca_peak = peakFit->GetParameter(1);
    pca_peak_err = peakFit->GetParError(1);
    pca_width = peakFit->GetParameter(2);
    pca_width_err = peakFit->GetParError(2);

    // One more loop for more precise peak fit
    ResidualHits2 = new TH1D ("ResidualHits2", "", 1000, pca_peak-3, pca_peak+2);
    for (int iPMT=0; iPMT<NPMTS; iPMT++) {
      if (pmtAngle[iPMT] <= 12.0){                                      // angular cut (beamspot only)
        if (pmtOccup[iPMT] >= 0.01 && pmtOccup[iPMT] <= 0.05){          // occupancy cut
          PMTBeamSpot++;
          for (int hit=0; hit<(pmtResidTimes[iPMT].size()); hit++){     // goes through valid hits here
            if ( (pmtResidTimes[iPMT][hit] > (pca_peak-3)) && (pmtResidTimes[iPMT][hit] < (pca_peak+2)) ){
              if (pmtResidTimes[iPMT][hit] > 0.2*pca_peak){
                ResidualHits2->Fill( pmtResidTimes[iPMT][hit] );
              }
            }
          }
        }
      }
    }

    // Fits2
    maxPeak2 = ResidualHits2->GetMaximum();
    maxBin2 = ResidualHits2->GetMaximumBin();
    //Fit a gaussian to the peak defined as the function between both half maximum points
    for(int i=1; i<maxBin2; i++){
        if(ResidualHits2->GetBinContent(i) > 0.01*maxPeak2){
            lowHalfMax2 = ResidualHits2->GetBinCenter(i-2);
            break;
        }
    }
    for(int i=maxBin2; i<=ResidualHits2->GetXaxis()->GetNbins(); i++){
        if(ResidualHits2->GetBinContent(i) < 0.01*maxPeak2){
            highHalfMax2 = ResidualHits2->GetBinCenter(i);
            break;
        }
    }

    ResidualHits2->Fit("gaus","","",lowHalfMax2,highHalfMax2);
    peakFit2 = ResidualHits2->GetFunction("gaus");
    pca_peak2 = peakFit2->GetParameter(1);
    pca_peak_err2 = peakFit2->GetParError(1);
    pca_width2 = peakFit2->GetParameter(2);
    pca_width_err2 = peakFit2->GetParError(2);

    cout << ResidualHits->GetMean() << " " << ResidualHits->GetRMS() << endl;
    cout << ResidualHits2->GetMean() << " " << ResidualHits2->GetRMS() << endl;
    cout << pca_peak << " " << pca_width << endl;
    cout << pca_peak2 << " " << pca_width2 << endl;
    cout << maxPeak << " " << maxPeak2 << endl;

    // Write all objects to ROOT file
    // TODO: save this as png for each fibre

    // rootfile = new TFile("TEST.root","RECREATE");
    // rootfile->cd();
    // ResidualHits->Write();
    // ResidualHits2->Write();
    // rootfile->Write();
    // rootfile->Close();
    // cout << "Wrote output to file " << "TEST.root" << endl;

    TCanvas *pad0 = new TCanvas("pad0","",1024,768);
    TCanvas *pad1 = new TCanvas("pad1","",1024,768);

    pad0->cd();
    ResidualHits->GetXaxis()->SetTitle("Residual hit time [ns]");
    ResidualHits->GetYaxis()->SetTitle("Count");
    ResidualHits->GetYaxis()->SetTitleOffset(1.3);
    ResidualHits->Draw();
    stringstream temp;
    temp << run_id << "_res1.png";
    pad0->SaveAs(temp.str().c_str());
    pad0->Close();

    pad1->cd();
    ResidualHits2->GetXaxis()->SetTitle("Residual hit time [ns]");
    ResidualHits2->GetYaxis()->SetTitle("Count");
    ResidualHits2->GetYaxis()->SetTitleOffset(1.3);
    ResidualHits2->Draw();
    stringstream temp2;
    temp2 << run_id << "_res2.png";
    pad1->SaveAs(temp2.str().c_str());

    cout << "Saved plots succesfully." << endl;

    // Store to file for next steps
    stringstream temp_output;
    temp_output << fibre_db << "\t" << run_id << "\t" << ResidualHits2->GetMean() << "\t" << ResidualHits2->GetRMS() << "\n";
    fprintf(out, temp_output.str().c_str());

    fprintf(logFile, "Injection time: %f\t%f\n", ResidualHits2->GetMean(), ResidualHits2->GetRMS());

    // Delete plots for next iter
    delete ResidualHits; delete ResidualHits2;

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
        lpc.CalcByPosition(fibrepos, pmtPos, energy, LOCALITY);

        // LPC checks
        if (lpc.GetTIR() == 1) { continue; }            // total internal reflection
        if (lpc.GetPathValid() == 0) { continue; }      // check whether path parameters are valid
        if (lpc.GetResvHit() == 1) { continue; }        // whether end point was within locality

        // Get light emission angle
        TVector3 startDir = lpc.GetInitialLightVec();            // start direction at fibre
        double theta = startDir.Angle(fitfibredir)*180./pi;      // emission angle at fibre

        // Angular cut here
        if ( (theta > 12) || (theta < 0) ) { continue; }

        if (lpc.GetTotalDist() <= 12000){nearL++;} else {farL++;}
        if (lpc.GetTotalDist() <= 12000){continue;}         // this rejects near reflections
        if (lpc.GetDistInInnerAV() <= 7000){continue;}      // this rejects other weird paths

        // ToF
        double distInInnerAV = lpc.GetDistInInnerAV();
        double distInAV = lpc.GetDistInAV();
        double distInWater = lpc.GetDistInWater();
        double lightTravelTime = gv.CalcByDistance(distInInnerAV, distInAV, distInWater, energy);

        // Get light bucket time
        TVector3 endDir = lpc.GetIncidentVecOnPMT();            // end direction at PMT
        double thetaAtPMT = endDir.Angle(pmtDir)*180./pi;       // incident angle with bucket face
        double lightBucketTime = gv.PMTBucketTime(thetaAtPMT);  // DocDB 3138

        // get ang sys correction
        double thetaRads = theta/180*M_PI;
        double additionalDelayFromAngSys = AngularSystematic(AngB, thetaRads);

        // Get residual time after correcting for all offsets
        double emissionTime = pmtTime - lightTravelTime - lightBucketTime - additionalDelayFromAngSys;
        //cout << theta << " " << pmtTime << " " << lightTravelTime << " " << lightBucketTime << " " << additionalDelayFromAngSys << " " << emissionTime << endl;
        pmtResidTimes[pmtID].push_back( emissionTime );
        if (emissionTime > residMax) residMax = emissionTime;
        if (emissionTime < residMin) residMin = emissionTime;

        if (pmtOccup[pmtID] == 0) {   // fill angles only once per PMT
          pmtAngle[pmtID] = theta;
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

  // Angular systematic fit function (the angle is in radians)
  double MyUserProc::AngularSystematic(double b, double angle){
    double y = 0.0;
    // currently not using the ang a value (set to 0)
    y = 0 + b * ( ( 1./cos(angle) ) - 1 );
    return y;
  }
