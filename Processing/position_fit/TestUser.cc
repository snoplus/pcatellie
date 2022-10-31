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

  virtual void GetLimits(std::vector<double> occupancy, int NPMTS, double &HOTLIMIT, double &COLDLIMIT);

  virtual void BinLog(TAxis *axis, const Double_t non0start=-999.);

  virtual void GetMaxColVal(const TVector3& center, std::vector<double> occupancy, int NPMTS, float &nearval, float &farval, const RAT::DU::PMTInfo& pmtinfo);

  virtual void FitLightSpot(TGraph2D* graph, double radius, double cone, double* params);

  virtual void GetRotationAngles(const TVector3& vec, double &rot_Z, double &rot_X);

  virtual void FillHemisphere(const TVector3& center, std::vector<double> occupancy, int NPMTS, std::vector<TGraph> &dots, TGraph2D* graph, int NCOL, float MAXVAL, const RAT::DU::PMTInfo& pmtinfo);

  virtual void DrawCircle(const TVector3& center, double angle, std::vector<TVector3> &dots, int NDOTS);

  virtual void CorrectFibre(const TVector3& fit_origin, double min_dist);

  virtual std::vector<string> ListFibres();

  virtual Processor::Result DSEvent(DS::Run& run, DS::Entry& ds);

protected:

  // Database
  DS::Run temp_run;
  DB *ratdb;
  RAT::DBLinkPtr belly_fibres;
  std::vector<string> belly_fibres_list;
  std::vector<string> belly_fibres_list_index;
  RAT::DBLinkPtr tellie_run_data;
  std::string fibre_db;
  json::Value sub_run_info;
  RAT::DBLinkPtr fibre_data;
  std::vector<string> fibresInUse;

  // ENV variables
  int dir_light_ang;
  int ref_light_ang;
  int NCOL;
  double min_dist;

  // Fits for light spots
  TVector3 *dirfit;
  TVector3 *reffit;
  bool IS_BELLY_FIBRE;
  TVector3 fibrepos, fibredir, lightpos;
  TVector3 lightOrigin;
  int correctFibre;

  // RAT variables
  int NPMTs;
  int run_id;

  // Event loop
  int trig;
  int count;        // actual all events
  int allEvs;       // just EXTA events
  int allNHit;
  std::vector<int> PMTHits;
  std::vector<double> PMTOccup;

  // For PMT limits
  Double_t EPSILON;
  double HOTLIMIT;
  double COLDLIMIT;
  float nearmax;
  float farmax;

  // Fiting faces
  RAT::DataQualityProc *dq;

  // Final spot fits
  double rot_Z1;
  double rot_X1;
  double rot_Z2;
  double rot_X2;
  double DIRANG;
  double REFANG;

  // Plots
  int NDOTS;             // number of points in circle
  TFile *rootfile;
  TH2D *hicos;
  TH2D *hcoarse;
  TH2D *hfineD;
  TH2D *hfineR;
  TH1F *hoccup;
  TH1F *hoccupoff;
  TH1F *hoccuplo;
  TH1F *hoccuphi;
  TGraph2D *gDir2D;
  TGraph2D *gRef2D;
  TGraph *pcontD;
  TGraph *pcontR;
  TGraph *pFibPos;
  TGraph *pFibDir;
  TGraph *pWgtDir;
  TGraph *pWgtRef;
  TGraph *pFitDir;
  TGraph *pFitRef;
  TGraph *pcircD1;
  TGraph *pcircD2ULL;
  TGraph *pcircR1;
  TGraph *pcircR2;
  TVector2 *angles;
  TVector2 *maxima;
  TVector2 *limits;
  TVector3 *fitDir;
  TVector3 *fitRef;
  TVector3 *fitResD;
  TVector3 *fitResR;
  std::vector<TGraph*> icos; std::vector<TGraph> gDir; std::vector<TGraph> gRef;

  // LPC
  double energy;
  double fLEDWavelength;
  RAT::DU::LightPathCalculator lpc;
  double LOCALITY;
  TVector3 fittedDir;
  double dirDiff;

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
    dir_light_ang = atoi(getenv("DIR_LIGHT_ANG"));
    ref_light_ang = atoi(getenv("REF_LIGHT_ANG"));
    LOCALITY = atoi(getenv("LOCALITY"));
    min_dist = atoi(getenv("MIN_DIST"));
    NDOTS = 360;
    fLEDWavelength = 506.0*1e-6;
    energy = util::WavelengthToEnergy(fLEDWavelength);

    // TODO: make sure this does not load bad PCA,ECA,CSS tables...
    ratdb = RAT::DB::Get();
    temp_run.SetRunID(275000);      // this allows to load custom table outside BeginOfRun for all runs
    ratdb->BeginOfRun(temp_run);    // the temp number is random

    // Load belly plate fibres from ratdb file
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

  }

  MyUserProc::~MyUserProc() {

  }

  void MyUserProc::BeginOfRun( DS::Run& run ) {

    // Initialise vars
    count = 0;
    allEvs = 0;
    HOTLIMIT = -1;
    COLDLIMIT = -1;
    NCOL = 20;
    EPSILON = 1e-10;
    nearmax = -1;
    farmax = -1;
    allNHit = 0;
    PMTHits.clear();
    PMTOccup.clear();
    dq = new RAT::DataQualityProc("");
    gDir2D = new TGraph2D();
    gRef2D = new TGraph2D();
    run_id = run.GetRunID();
    rot_Z1 = -1;
    rot_X1 = -1;
    rot_Z2 = -1;
    rot_X2 = -1;
    DIRANG = -1;
    REFANG = -1;
    correctFibre = 0;

    icos.clear(); gDir.clear(); gRef.clear();
    icos.resize(NCOL+2); gDir.resize(NCOL+2); gRef.resize(NCOL+2);

    fibrepos.SetXYZ(0,0,0);
    fibredir.SetXYZ(0,0,0);
    lightpos.SetXYZ(0,0,0);

    // lpc
    lpc = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lpc.BeginOfRun();
    lpc.SetELLIEEvent(true);
    dirDiff = 0;

    // Load basic fibre details
    try {
      tellie_run_data = DB::Get()->GetLink("TELLIE_RUN");
      sub_run_info = tellie_run_data->GetJSON("sub_run_info");
      fibre_db = sub_run_info[0]["fibre"].getString();
    }
    catch (DBNotFoundError& e) {
      std::cout << "Proc(): Couldn't load TELLIE Run details from ratdb!" << std::endl;
    }
    cout << "FIBRE: " << fibre_db << endl;

    // Find out if direct light from fibre is affected by belly plates
    IS_BELLY_FIBRE = false;
    for (size_t i=0; i<belly_fibres_list.size(); i++) {
      if ( belly_fibres_list[i] == fibre_db ){
        IS_BELLY_FIBRE = true;
      }
    }
    cout << "Is affected by belly plate? " << IS_BELLY_FIBRE << endl;

    // Load fibre position and direction
    try {
      fibre_data = DB::Get()->GetLink("FIBRE", fibre_db);
      fibrepos.SetXYZ(fibre_data->GetD("x"), fibre_data->GetD("y"), fibre_data->GetD("z"));  // position
      fibredir.SetXYZ(fibre_data->GetD("u"), fibre_data->GetD("v"), fibre_data->GetD("w"));  // direction
      lightpos = fibrepos + 2*fibrepos.Mag()*fibredir;                        // projected light spot centre
    }
    catch (DBNotFoundError& e) {
      std::cout << "Proc(): Couldn't load FIBRE details from ratdb!" << std::endl;
    }
    cout << "RATDB: fibre " << fibre_db << ", pos " << printVector(fibrepos) << ", dir " << printVector(fibredir) << endl;

    // Initialise RAT stuff
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
    NPMTs = pmtinfo.GetCount();
    cout << "Initialised RAT. Number of PMTs is " << NPMTs << "." << endl;
    PMTOccup.resize(NPMTs);
    PMTHits.resize(NPMTs);

    // Load current set of fibres (patched fibres in use)
    fibresInUse = ListFibres();
    cout << "We have " << fibresInUse.size() << " fibres in use!" << endl;

    // Output file for stats
    logFile_namess.str("");
    logFile_namess << run_id << "_pos.log";
    logFile_name = logFile_namess.str();

  }

  void MyUserProc::EndOfRun( DS::Run& run ) {

    // Init RAT stuff (again :/)
    const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();

    // Initialise histograms
    hicos = new TH2D("hicos","PMT positions",200,0,1,200,0,1); // icosahedral
    hcoarse = new TH2D("hcoarse","PMT positions",20,-1,1.001,10,-1.001,0); // units of pi
    hfineD = new TH2D("hfineD","PMT positions",1000,-10,10,1000,-10,10); // fine grained
    hfineR = new TH2D("hfineR","PMT positions",1000,-10,10,1000,-10,10); // fine grained
    hoccup = new TH1F("hoccup","PMT hit count",60,0,1); // occupancy
    BinLog(hoccup->GetXaxis(),1e-6); // minimum for log axis

    // Convert hits on PMTs to Occupancy
    for (size_t i=0; i<NPMTs; i++) {
      PMTOccup[i] = double(PMTHits[i])/double(allEvs);
      if(PMTOccup[i]==0) hoccup->Fill(2e-6); // visible on log scale
      else hoccup->Fill(PMTOccup[i]);
    }

    // Basic stats
    cout << "Run: " << run_id << endl;                   // run number
    cout << "Events: " << allEvs << endl;                // events with EXTA trigger
    cout << "TotalNHit: " << allNHit << endl;            // calibrated PMT hits with removed crosstalk
    cout << "NPMTS: " << NPMTs << endl;                  // number of PMTs

    // Get PMT limits
    GetLimits(PMTOccup, NPMTs, HOTLIMIT, COLDLIMIT);

    cout << "HOTLIMIT: " << HOTLIMIT << endl;
    cout << "COLDLIMIT: " << COLDLIMIT << endl;

    // ********************************************************************
    // Get icosahedral projection of detector
    // ********************************************************************
    TVector2 icospos;
    TVector3 pmtpos;
    std::vector<int> icosN;
    icosN.resize(NCOL+2);
    std::vector< std::vector<double> > icosX;
    std::vector< std::vector<double> > icosY;
    icosX.resize(NCOL+2);
    icosY.resize(NCOL+2);
    for (size_t i=0; i<NCOL+2; i++) {
      icosX[i].resize(NPMTs);
      icosY[i].resize(NPMTs);
    }
    std::vector<int> pmtface;
    pmtface.resize(NPMTs);
    int nfacepmts[20] = {0}, nfacegoodpmts[20] = {0};
    TVector3 facecenter[20], faceweight[20];
    cout << "Creating icosahedral projection..." << endl;
    for(int id=0; id<NPMTs; id++) {
      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue; // not a valid PMT
      if (pmtinfo.GetType(id) != 1) continue; // not a normal PMT

      // Obtain PSUP face
      int face; // side of icosahedron containing PMT
      icospos = dq->IcosProject(pmtpos,face);
      pmtface[id] = face; // face of PMT (1-20)
      facecenter[face-1] += pmtpos;
      nfacepmts[face-1]++;

      // Get correct bin for colour scale
      int step;
      if (PMTOccup[id] == 0) step=0;                  // off PMT
      else if (PMTOccup[id] < COLDLIMIT)  step=1;     // cold PMT
      else if (PMTOccup[id] > HOTLIMIT) step=NCOL+1;  // hot PMT
      // logarithmic colour scale
      else step = (int)TMath::Ceil((log10(PMTOccup[id])-log10(COLDLIMIT)) / ((log10(HOTLIMIT)-log10(COLDLIMIT))/(NCOL-1))) + 1;

      // PMT checks
      icosX[step][icosN[step]] = icospos.X();
      icosY[step][icosN[step]] = icospos.Y();
      icosN[step]++;
      if (step<2) continue; // cold or off PMT
      if (step>NCOL) continue; // hot PMT

      // Calculate PSUP face heat
      faceweight[face-1] += pmtpos*PMTOccup[id];
      nfacegoodpmts[face-1]++;
    }
    int maxface=-1;
    double faceheat[20]={0}, maxfaceheat=-1;
    for(int fc=0; fc<20; fc++) {
      if (nfacepmts[fc]==0 || nfacegoodpmts[fc]==0) continue;
      facecenter[fc] *= 1./nfacepmts[fc];     // all PMTs weighted equally
      faceweight[fc] *= 1./nfacegoodpmts[fc]; // good PMTS weighted by intensity
      faceheat[fc] = faceweight[fc].Mag();
      cout << "Face #" << fc+1 << " has intensity " << faceheat[fc] << endl;
      if (faceheat[fc]>maxfaceheat) { maxfaceheat=faceheat[fc]; maxface=fc+1; }
    }
    printf("Hot faces: ");
    TVector3 bestguess(0,0,0);
    int hotfaces=0;
    for(int fc=0; fc<20; fc++) {
      if(faceheat[fc]<maxfaceheat/5.) continue; // face is "hot" if intensity is >20% of hottest face
      printf("%d ",fc+1);
      bestguess += faceweight[fc];
      hotfaces++;
    }
    if (hotfaces>0) { bestguess *= 1./hotfaces; }
    else { cout << "Something went wrong! No hot faces found." << endl; }
    printf(" -> best guess direction: %s\n",printVector(bestguess.Unit()).c_str());

    // Make graphs for icosahedral projection
    for (int s=0; s<NCOL+2; s++) {
      int col;
      if      (s==0) col = 16; // off PMTs (grey)
      else if (s==1) col = 51; // cold PMTs (purple)
      else if (s==21) col = 1; // hot PMTs (black)
      else col = (int)(50.+s*(50./NCOL)); // scale from 55-100 (s=2-20)
      if (icosN[s]==0) {continue;}
      icos[s] = new TGraph( icosN[s], &icosX[s][0], &icosY[s][0]);
      icos[s]->SetMarkerStyle(7);
      icos[s]->SetMarkerColor(col);
    }

    // ********************************************************************
    // First iteration: Fill coarse histogram with PMT hit positions
    // ********************************************************************
    float hitsum = 0;
    double pmtrad = 0;
    TVector3 pmtsum(0,0,0);
    for(int id=0; id<NPMTs; id++) {
      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue;            // not a valid PMT
      if (pmtinfo.GetType(id) != 1) continue;   // not a normal PMT
      if (PMTOccup[id] == 0) continue;         // off PMTs
      if (PMTOccup[id] > HOTLIMIT) continue;   // hot PMTs
      if (PMTOccup[id] < COLDLIMIT) continue;  // cold PMTs
      pmtrad += pmtpos.Mag()*PMTOccup[id];
      pmtsum += pmtpos*PMTOccup[id];
      hitsum += PMTOccup[id];
    }
    pmtrad /= hitsum;
    if(fabs(lightpos.Mag()/pmtrad-1.) > 0.01) {  // tolerate 1%
      cout << "*** WARNING *** Projected light position deviates from PSUP sphere!" << endl;
      return;
    }
    cout << "PMT radius = " << pmtrad << " mm\n" << endl;

    // Get the maximum bin as guesstimate for light spot
    bestguess.SetMag(pmtrad);
    TVector3 guess_dir =  bestguess;
    TVector3 guess_ref = -bestguess;  // TODO - improve guesstimate for reflected light?

    cout << "End of first itter: " << endl;
    printf(" -> best guess direction: %s\n",printVector(bestguess.Unit()).c_str());

    // ********************************************************************
    // Second iteration: take weighted average around estimated light spots
    // ********************************************************************
    TVector3 direct(0.,0.,0.);
    TVector3 reflected(0.,0.,0.);
    float weight, empty1, empty2;
    GetMaxColVal(guess_dir, PMTOccup, NPMTs, nearmax, empty1, pmtinfo);
    GetMaxColVal(guess_ref, PMTOccup, NPMTs, farmax, empty2, pmtinfo);

    if(nearmax==0) cout << "*** WARNING *** No good PMTs in direct light cone!" << endl;
    if(farmax==0) cout << "*** WARNING *** No good PMTs in reflected light cone!" << endl;

    cout << "Calculating weighted average..." << endl;
    for(int id=0; id<NPMTs; id++) {
      //if (id % int(NPMTS/100.) == 0) printProgress(id, NPMTS);
      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue;                  // not a valid PMT position
      if (pmtinfo.GetType(id) != 1) continue;         // not a normal PMT
      weight = PMTOccup[id];
      if (weight == 0) continue;                      // off PMTs
      if (weight > HOTLIMIT) continue;                // hot PMTs
      if (weight < COLDLIMIT) continue;               // cold PMTs
      if (pmtpos.Angle(guess_dir) < pi*dir_light_ang/180.) // within cone of direct light spot
        direct += pmtpos*(1.*weight/hitsum);
      if (pmtpos.Angle(guess_ref) < pi*ref_light_ang/180.) // within cone of reflected light spot
        reflected += pmtpos*(1.*weight/hitsum);
    }
    if (direct.Mag() == 0) {
      cout << "*** WARNING *** No good PMTs in direct light cone! Setting to z-axis." << endl;
      direct.SetXYZ(0,0,1);
    }
    if (reflected.Mag() == 0) {
      cout << "*** WARNING *** No good PMTs in reflected light cone! Setting to z-axis." << endl;
      reflected.SetXYZ(0,0,1);
    }
    direct.SetMag(pmtrad);
    reflected.SetMag(pmtrad);

    cout << "End of second itter: " << endl;
    printf(" -> best guess direction: %s\n",printVector(direct.Unit()).c_str());

    // Fill graphs with new view, centered on weighted light spot
    FillHemisphere(direct, PMTOccup, NPMTs, gDir, gDir2D, NCOL, nearmax, pmtinfo);
    FillHemisphere(reflected, PMTOccup, NPMTs, gRef, gRef2D, NCOL, farmax, pmtinfo);

    // Get rotation angles for rotating view over weighted light spots
    GetRotationAngles(direct, rot_Z1, rot_X1);
    GetRotationAngles(reflected, rot_Z2, rot_X2);

    // ********************************************************************
    // Third iteration: perform 2D Gaussian fit around weighted light spots
    // ********************************************************************
    // Fit 2D graphs with Gaussian surface
    const int NPARS = 4;
    double parDir[NPARS], parRef[NPARS];
    cout << "Performing 2D Gaussian fits..." << endl;
    FitLightSpot(gDir2D,pmtrad/1e3,dir_light_ang,parDir); // units: dist [m], ang [deg]
    FitLightSpot(gRef2D,pmtrad/1e3,ref_light_ang,parRef);
    cout << "FIT RESULTS:" << endl;
    cout << "- Direct light";
    for (int p=0; p<NPARS; p++) printf(" %.3lf", parDir[p]);
    cout << endl;
    cout << "- Reflected light";
    for (int p=0; p<NPARS; p++) printf(" %.3lf", parRef[p]);
    cout << endl;

    // Translate to point on PSUP sphere
    fitDir = new TVector3(1e3*parDir[1],1e3*parDir[2],sqrt(pmtrad*pmtrad-1e6*parDir[1]*parDir[1]-1e6*parDir[2]*parDir[2]));
    fitRef = new TVector3(1e3*parRef[1],1e3*parRef[2],sqrt(pmtrad*pmtrad-1e6*parRef[1]*parRef[1]-1e6*parRef[2]*parRef[2]));
    // Standard deviation (Gaussian sigma) translated to angle [deg]
    double sigmaAngDir = asin(1e3*parDir[3]/pmtrad)*180./pi;
    double sigmaAngRef = asin(1e3*parRef[3]/pmtrad)*180./pi;

    // Output values (Gaussian fit)
    dirfit = new TVector3(0,0,0);
    dirfit->SetXYZ(fitDir->X(), fitDir->Y(), fitDir->Z());
    dirfit->RotateY(rot_Z1);
    dirfit->RotateZ(rot_X1);
    reffit = new TVector3(0,0,0);
    reffit->SetXYZ(fitDir->X(), fitDir->Y(), fitDir->Z());
    reffit->RotateY(rot_Z2);
    reffit->RotateZ(rot_X2);

    DIRANG = lightpos.Angle(*dirfit);   // angle between expected and fitted light spot
    REFANG = fibrepos.Angle(*reffit);   // angle between fibre position and fitted light spot

    // Overwrite value for direct light if affected by belly plates
    if (IS_BELLY_FIBRE) {
      cout << "Direct light is affected by belly plates. Using weighted average instead of Gaussian fit." << endl;
      DIRANG = lightpos.Angle(direct);         // angle between expected and weighted light spot
      //OVERWRITE DIRECT LIGHT FIT RESULT FOR FIBRES AFFECTED BY BELLY PLATES
      dirfit->SetXYZ(direct.X(), direct.Y(), direct.Z());
    }

    // Check what fibre was firing
    lightOrigin = pmtrad*(reffit->Unit());
    CorrectFibre(lightOrigin, min_dist);
    cout << "Correct fibre? " << correctFibre << endl;

    cout << "Final DIR fit: " << dirfit->X() << " " << dirfit->Y() << " " << dirfit->Z() << endl;
    cout << "DIR fit angular deviation from expected: " << DIRANG << endl;
    cout << "Final REF fit: " << reffit->X() << " " << reffit->Y() << " " << reffit->Z() << endl;
    cout << "REF fit angular deviation from expected: " << REFANG << endl;

    // ********************************************************************
    // Getting the new fibre directions here
    // ********************************************************************

    lpc.CalcByPosition(fibrepos, *dirfit, energy, LOCALITY);
    //lpc.CalcByPositionPartial( fibrepos, *dirfit, energy, LOCALITY ); // partial
    fittedDir = lpc.GetInitialLightVec();
    dirDiff = fibredir.Angle(fittedDir)*180./M_PIl;

    cout << "New fitted fibre direction: " << fittedDir.X() << " " << fittedDir.Y() << " " << fittedDir.Z() << endl;
    cout << "Angular difference: " << dirDiff << endl;

    // Store to file for next steps
    // Currently storing the direct light position and the direction
    stringstream temp_output;
    temp_output << fibre_db << "\t" << fittedDir.X() << "\t" << fittedDir.Y() << "\t" << fittedDir.Z();
    temp_output << "\t" << dirfit->X() << "\t" << dirfit->Y() << "\t" << dirfit->Z() << "\n";
    out = fopen("direction_fit.txt","a");
    fprintf(out, temp_output.str().c_str());
    fclose(out);

    // ********************************************************************
    // Create histograms and graphs for output
    // ********************************************************************

    // Get contours around fitted light spots in (phi,theta)
    std::vector<TVector3> dotsD;
    std::vector<TVector3> dotsR;
    dotsD.resize(NDOTS);
    dotsR.resize(NDOTS);
    if (IS_BELLY_FIBRE) DrawCircle(direct, sigmaAngDir, dotsD, NDOTS); // weighted method (TODO - use fitted sigma?)
    else DrawCircle(*dirfit, sigmaAngDir, dotsD, NDOTS);               // Gaussian fit
    DrawCircle(*reffit, sigmaAngRef, dotsR, NDOTS);

    std::vector<double> pcontDx;
    std::vector<double> pcontDy;
    std::vector<double> pcontRx;
    std::vector<double> pcontRy;
    pcontDx.resize(NDOTS);
    pcontDy.resize(NDOTS);
    pcontRx.resize(NDOTS);
    pcontRy.resize(NDOTS);
    for (int d=0; d<NDOTS; d++) {
      int dummy;
      icospos = dq->IcosProject(dotsD[d],dummy);
      pcontDx[d]=icospos.X();
      pcontDy[d]=icospos.Y();
      icospos = dq->IcosProject(dotsR[d],dummy);
      pcontRx[d]=icospos.X();
      pcontRy[d]=icospos.Y();
    }
    pcontD = new TGraph(NDOTS,&pcontDx[0],&pcontDy[0]);
    pcontR = new TGraph(NDOTS,&pcontRx[0],&pcontRy[0]);
    pcontD->SetMarkerStyle(6);
    pcontR->SetMarkerStyle(6);

    int weightMarker = 2;
    int fitMarker = 34;
    int trueMarker = 24;

    // Create markers for relevant points (icosahedral view)
    int nope;
    double ptX, ptY;
    if (IS_BELLY_FIBRE) icospos = dq->IcosProject(direct,nope); // use weighted average
    else icospos = dq->IcosProject(*dirfit,nope);               // use Gaussian fit
    ptX = icospos.X();
    ptY = icospos.Y();
    TGraph *pFitDir = new TGraph(1,&ptX,&ptY);          // direct light (fitted)
    pFitDir->SetMarkerStyle(fitMarker);
    pFitDir->SetMarkerColor(1);
    pFitDir->SetMarkerSize(2);

    icospos = dq->IcosProject(*reffit,nope);
    ptX = icospos.X();
    ptY = icospos.Y();
    TGraph *pFitRef = new TGraph(1,&ptX,&ptY);          // reflected light (fitted)
    pFitRef->SetMarkerStyle(fitMarker);
    pFitRef->SetMarkerColor(1);
    pFitRef->SetMarkerSize(1.5);

    // Create markers for relevant points (rotated view)
    lightpos.RotateZ(-rot_X1);
    lightpos.RotateY(-rot_Z1);
    ptX = lightpos.X()/1e3;
    ptY = lightpos.Y()/1e3;
    TGraph *pFibDir = new TGraph(1,&ptX,&ptY);      // expected light position
    pFibDir->SetMarkerStyle(trueMarker);
    pFibDir->SetMarkerColor(1);
    pFibDir->SetMarkerSize(2);

    fibrepos.RotateZ(-rot_X2);
    fibrepos.RotateY(-rot_Z2);
    ptX = fibrepos.X()/1e3;
    ptY = fibrepos.Y()/1e3;
    TGraph *pFibPos = new TGraph(1,&ptX,&ptY);      // expected fibre position
    pFibPos->SetMarkerStyle(trueMarker);
    pFibPos->SetMarkerColor(1);
    pFibPos->SetMarkerSize(1.5);

    // View is centered over weighted light spot
    const int nil = 0;
    TGraph *pWgtDir = new TGraph(1,&nil,&nil);         // direct light (weighted)
    if (IS_BELLY_FIBRE) pWgtDir->SetMarkerStyle(fitMarker);
    else pWgtDir->SetMarkerStyle(weightMarker);
    pWgtDir->SetMarkerColor(1);
    pWgtDir->SetMarkerSize(2);

    TGraph *pWgtRef = new TGraph(1,&nil,&nil);         // reflected light (weighted)
    pWgtRef->SetMarkerStyle(weightMarker);
    pWgtRef->SetMarkerColor(1);
    pWgtRef->SetMarkerSize(1.5);

    std::vector< std::vector<double> > pcircDx;
    std::vector< std::vector<double> > pcircDy;
    std::vector< std::vector<double> > pcircRx;
    std::vector< std::vector<double> > pcircRy;
    pcircDx.resize(2);
    pcircDy.resize(2);
    pcircRx.resize(2);
    pcircRy.resize(2);
    for (size_t i=0; i<2; i++) {
      pcircDx[i].resize(NDOTS);
      pcircDy[i].resize(NDOTS);
      pcircRx[i].resize(NDOTS);
      pcircRy[i].resize(NDOTS);
    }
    int nfrontD=0, nbackD=0, nfrontR=0, nbackR=0;
    for (int d=0; d<NDOTS; d++) {
      // Rotate circle around direct spot
      dotsD[d].RotateZ(-rot_X1);
      dotsD[d].RotateY(-rot_Z1);
      if(dotsD[d].Z() > 0) {
        pcircDx[0][nfrontD]=dotsD[d].X()/1e3;
        pcircDy[0][nfrontD]=dotsD[d].Y()/1e3;
        nfrontD++;
      } else {
        pcircDx[1][nbackD]=dotsD[d].X()/1e3;
        pcircDy[1][nbackD]=dotsD[d].Y()/1e3;
        nbackD++;
      }
      // Rotate circle around reflected spot
      dotsR[d].RotateZ(-rot_X2);
      dotsR[d].RotateY(-rot_Z2);
      if(dotsR[d].Z() > 0) {
        pcircRx[0][nfrontR]=dotsR[d].X()/1e3;
        pcircRy[0][nfrontR]=dotsR[d].Y()/1e3;
        nfrontR++;
      } else {
        pcircRx[1][nbackR]=dotsR[d].X()/1e3;
        pcircRy[1][nbackR]=dotsR[d].Y()/1e3;
        nbackR++;
      }
    }
    TGraph *pcircD1=NULL, *pcircD2=NULL, *pcircR1=NULL, *pcircR2=NULL;
    if(nfrontD>0) { pcircD1 = new TGraph(nfrontD,&pcircDx[0][0],&pcircDy[0][0]);
                    pcircD1->SetMarkerStyle(6); }
    if(nbackD>0)  { pcircD2 = new TGraph(nbackD,&pcircDx[1][0],&pcircDy[1][0]);
                    pcircD2->SetMarkerStyle(1); }
    if(nfrontR>0) { pcircR1 = new TGraph(nfrontR,&pcircRx[0][0],&pcircRy[0][0]);
                    pcircR1->SetMarkerStyle(6); }
    if(nbackR>0)  { pcircR2 = new TGraph(nbackR,&pcircRx[1][0],&pcircRy[1][0]);
                    pcircR2->SetMarkerStyle(1); }

    // ********************************************************************
    // Storing output
    // ********************************************************************

    // We do not need this, but leaving if we decide we want more output
    // rootfile = new TFile("TEST.root","RECREATE");
    // rootfile->cd();
    //
    // // Write all objects to file (histograms included automatically, graphs not)
    // angles = new TVector2(DIRANG,REFANG);
    // angles->Write("angles");
    // TVector2 *maxima = new TVector2(nearmax,farmax);
    // maxima->Write("maxima");
    // TVector2 *limits = new TVector2(HOTLIMIT,COLDLIMIT);
    // limits->Write("limits");
    // fitDir->Write("fitDir");
    // fitRef->Write("fitRef");
    // dirfit->Write("dirfit");
    // reffit->Write("reffit");
    // gDir2D->Write("gDir2D");
    // for (int s=0; s<NCOL+2; s++) {
    //   string ss = Form("[%d]",s);
    //   if(icos[s]){ icos[s]->Write(("icos"+ss).c_str()); }
    //   if(gDir[s].GetN() > 0){ gDir[s].Write(("gDir"+ss).c_str()); }
    //   if(gRef[s].GetN() > 0) gRef[s].Write(("gRef"+ss).c_str());
    // }
    // pcontD->Write("pcontD");
    // pcontR->Write("pcontR");
    // pFibPos->Write("pFibPos");
    // pFibDir->Write("pFibDir");
    // pWgtDir->Write("pWgtDir");
    // pWgtRef->Write("pWgtRef");
    // pFitDir->Write("pFitDir");
    // pFitRef->Write("pFitRef");
    // if(nfrontD>0) pcircD1->Write("pcircD1");
    // if(nbackD>0)  pcircD2->Write("pcircD2");
    // if(nfrontR>0) pcircR1->Write("pcircR1");
    // if(nbackR>0)  pcircR2->Write("pcircR2");

    // ********************************************************************
    // Plotting section
    // ********************************************************************

    // No rotation required, fit was performed in rotated view
    double fitD_rotX = fitDir->X()/1e3;
    double fitD_rotY = fitDir->Y()/1e3;
    double fitR_rotX = fitRef->X()/1e3;
    double fitR_rotY = fitRef->Y()/1e3;

    TCanvas *c0 = new TCanvas("c0","",1024,768);
    TCanvas *c1 = new TCanvas("c1","",1024,768);
    TCanvas *c2 = new TCanvas("c2","",1024,768);
    TCanvas *c3 = new TCanvas("c3","",1024,768);
    TCanvas *c4 = new TCanvas("c4","",1024,768);

    gStyle->SetTitleOffset(1.2,"xyz");

    // ----------
    // Run summary
    c0->cd();
    TLatex *title, *t[6], *v[6];
    float avgnhit = (float)allNHit/allEvs;
    char* text[6] = { Form(" - Run number:"),
                      Form(" - Fibre name:"),
                      Form(" - Number of events:"),
                      Form(" - Average NHit:"),
                      Form(" - Fit deviation (dir.):"),
                      Form(" - Fit deviation (refl.):")
                    };
    char* vals[6] = { Form("%d",run_id),               // avgnhit calc, add exta evs
                      Form("%s",fibre_db.c_str()),
                      Form("%d",allEvs),
                      Form("%.2f",avgnhit),
                      Form("%.2f#circ",DIRANG/pi*180),
                      Form("%.2f#circ",REFANG/pi*180)
                    };

    // Highlight unsatisfactory values in bold
    int nexpected = 2e5;
    if (count < 0.99*nexpected) vals[2] = Form("#bf{%s}",vals[2]);     // lost more than 1% of triggers
    if (avgnhit<25 || avgnhit>45) vals[3] = Form("#bf{%s}",vals[3]);   // nhit far outside optimal range (30-40)
    if (DIRANG/pi*180 >= 10.) vals[4] = Form("#bf{%s}",vals[4]); // bad direct light fit
    if (REFANG/pi*180 >= 10.) vals[5] = Form("#bf{%s}",vals[5]); // bad reflected light fit

    // Place text in pad
    title = new TLatex(0.05,0.9,"SNO+ TELLIE PCA data");
    title->SetTextAlign(12);
    title->SetTextFont(62);
    title->SetTextSize(0.08);
    title->Draw();
    for (int l=0; l<6; l++) {
      t[l] = new TLatex(0.05,0.7-0.12*l,text[l]);
      t[l]->SetTextAlign(11);
      t[l]->SetTextFont(82);
      t[l]->SetTextSize(0.04);
      t[l]->Draw();
      v[l] = new TLatex(0.75,0.7-0.12*l,vals[l]);
      v[l]->SetTextAlign(11);
      v[l]->SetTextFont(82);
      v[l]->SetTextSize(0.04);
      v[l]->Draw();
    }

    // Indicate logarithmic scale
    TLatex *txtL = new TLatex(0.01,0.9,"LOG SCALE");
    txtL->SetTextAlign(13);
    txtL->SetTextFont(102);
    txtL->SetTextSize(0.025);

    // Indicate possible belly plate effect
    TLatex *txtB = new TLatex(-9.25,9.5,"BELLY PLATE");
    txtB->SetTextAlign(13);
    txtB->SetTextFont(102);
    txtB->SetTextSize(0.02);

    // Indicate colour scale
    TLatex *txtD = new TLatex(9.5,9.5,Form("Occup. #leq %.2f%%",100.*nearmax));
    txtD->SetTextAlign(33);
    txtD->SetTextFont(82);
    txtD->SetTextSize(0.02);
    TLatex *txtR = new TLatex(9.5,9.5,Form("Occup. #leq %.2f%%",100.*farmax));
    txtR->SetTextAlign(33);
    txtR->SetTextFont(82);
    txtR->SetTextSize(0.02);

    stringstream temp;
    temp << run_id << "_c0.png";
    c0->SaveAs(temp.str().c_str());
    c0->Close();

    // ----------
    // PMT hit count histogram
    c1->cd();
    c1->SetGrid();
    c1->SetLogx();
    c1->SetLogy();
    hoccup->SetMinimum(0.5);
    hoccup->SetMaximum(5e3);
    hoccup->SetTitle("PMT occupancy;Occupancy;NPMTs");
    hoccup->SetLineWidth(2);
    hoccup->SetLineColor(1);
    hoccup->GetXaxis()->SetTitleOffset(1.3);
    //hoccup->GetYaxis()->SetTitleOffset(1.3);
    hoccup->Draw();

    // Draw off, cold and hot PMT bins in their respective colour
    hoccupoff = (TH1F*)hoccup->Clone("hoccupoff");
    hoccupoff->SetFillColor(16);
    hoccupoff->GetXaxis()->SetRangeUser(1e-6,3e-6);
    hoccupoff->Draw("same");
    hoccuplo = (TH1F*)hoccup->Clone("hoccuplo");
    hoccuplo->SetFillColor(51);
    double coldlimedge;
    for (int b=hoccuplo->GetNbinsX(); b>0; b--) {
      if(hoccuplo->GetBinCenter(b)>COLDLIMIT) continue;
      coldlimedge = hoccuplo->GetBinCenter(b);
      break;
    }
    hoccuplo->GetXaxis()->SetRangeUser(3e-6,coldlimedge);
    hoccuplo->Draw("same");
    hoccuphi = (TH1F*)hoccup->Clone("hoccuphi");
    hoccuphi->SetFillColor(1);
    double hotlimedge;
    for (int b=0; b<hoccuphi->GetNbinsX(); b++) {
      if(hoccuphi->GetBinCenter(b)<HOTLIMIT) continue;
      hotlimedge = hoccuphi->GetBinCenter(b);
      break;
    }
    hoccuphi->GetXaxis()->SetRangeUser(hotlimedge,1.);
    hoccuphi->Draw("same");

    // Draw line indicating "hot PMT" limit
    TLine *lhot = new TLine(HOTLIMIT,0,HOTLIMIT,5e3);
    lhot->SetLineWidth(2);
    lhot->SetLineColor(2);
    lhot->Draw("same");

    // Draw line indicating "cold PMT" limit
    TLine *lcold = new TLine(COLDLIMIT,0,COLDLIMIT,5e3);
    lcold->SetLineWidth(2);
    lcold->SetLineColor(4);
    lcold->Draw("same");

    stringstream temp1;
    temp1 << run_id << "_c1.png";
    c1->SaveAs(temp1.str().c_str());
    c1->Close();

    // ----------
    // Icosahedral projection of detector display
    c2->cd();
    c2->SetGrid();
    //hicos->SetTitle("Detector display (PMT hit sum)");
    //hicos->Draw();
    for(int s=0;s<NCOL+2;s++) {
      if(!icos[s]) continue;
      if(icos[s]->GetN()==0) continue;
      icos[s]->Draw("P same");
    }
    //phot->Draw("P same");
    pFitDir->Draw("P same");
    pFitRef->Draw("P same");
    pcontD->Draw("P same");
    pcontR->Draw("P same");
    txtL->Draw();

    stringstream temp2;
    temp2 << run_id << "_c2.png";
    c2->SaveAs(temp2.str().c_str());
    c2->Close();

    // ----------
    // View from direct light spot (fitted)
    c3->cd();
    c3->SetGrid();
    hfineD->SetTitle("Direct light (PMT hit sum);X' [m];Y' [m]");
    hfineD->GetXaxis()->SetTitleOffset(1.3);
    hfineD->GetYaxis()->SetTitleOffset(1.4);
    hfineD->Draw("scat");
    hfineD->SetStats(0);
    for(int s=0;s<NCOL+2;s++) {;
      if(gDir[s].GetN()==0) continue;
      gDir[s].Draw("P same");
    }
    pFibDir->Draw("P same");
    pWgtDir->Draw("P same");
    TGraph *pFitDir2 = (TGraph*)pFitDir->Clone();
    if (IS_BELLY_FIBRE) pFitDir2->SetMarkerStyle(28); // open cross = Gaussian fit not used
    pFitDir2->DrawGraph(1,&fitD_rotX,&fitD_rotY,"P same");
    if(pcircD1) pcircD1->Draw("P same");
    if(pcircD2) pcircD2->Draw("P same");
    txtD->Draw();
    if (IS_BELLY_FIBRE) txtB->Draw();
    //hfineD->Write("hfineD");

    stringstream temp3;
    temp3 << run_id << "_c3.png";
    c3->SaveAs(temp3.str().c_str());
    c3->Close();

    // ----------
    // View from reflected light spot (fitted)
    c4->cd();
    c4->SetGrid();
    hfineR->SetTitle("Reflected light (PMT hit sum);X'' [m];Y'' [m]");
    hfineR->GetXaxis()->SetTitleOffset(1.3);
    hfineR->GetYaxis()->SetTitleOffset(1.4);
    hfineR->Draw("scat");
    hfineR->SetStats(0);
    for(int s=0;s<NCOL+2;s++) {
      if(gRef[s].GetN()==0) continue;
      gRef[s].Draw("P same");
    }
    pFibPos->Draw("P same");
    pWgtRef->Draw("P same");
    pFitRef->DrawGraph(1,&fitR_rotX,&fitR_rotY,"P same");
    if(pcircR1) pcircR1->Draw("P same");
    if(pcircR2) pcircR2->Draw("P same");
    txtR->Draw();

    stringstream temp4;
    temp4 << run_id << "_c4.png";
    c4->SaveAs(temp4.str().c_str());
    c4->Close();

    // ----------
    // Save canvas and close
    // stringstream file_name;
    // file_name << fibre_db;
    // string outfile = file_name.str();
    // c0->SaveAs(Form("%s.png",outfile.c_str()));
    // //c0->Print(Form("%s.pdf",outfile.c_str()));
    // c0->Close();

    // Delete pointers created with 'new'
    if(hicos) delete hicos;
    if(hcoarse) delete hcoarse;
    if(hfineD) delete hfineD;
    if(hfineR) delete hfineR;
    if(hoccup) delete hoccup;
    if(hoccupoff) delete hoccupoff;
    if(hoccuplo) delete hoccuplo;
    if(hoccuphi) delete hoccuphi;
    //if(gDir2D) delete gDir2D;
    //if(gRef2D) delete gRef2D;
    for(int s=0; s<NCOL+2; s++) delete icos[s];
    //if (icos) { for(int s=0; s<NCOL+2; s++) delete icos[s]; }
    //if (gDir) { for(int s=0; s<NCOL+2; s++) delete gDir[s]; }
    //if (gRef) { for(int s=0; s<NCOL+2; s++) delete gRef[s]; }
    if(pcontD) delete pcontD;
    if(pcontR) delete pcontR;
    if(pFibPos) delete pFibPos;
    if(pFibDir) delete pFibDir;
    if(pWgtDir) delete pWgtDir;
    if(pWgtRef) delete pWgtRef;
    if(pcircD1) delete pcircD1;
    if(pcircD2) delete pcircD2;
    if(pcircR1) delete pcircR1;
    if(pcircR2) delete pcircR2;
    if(c0) delete c0;
    if(title) delete title;
    if(t) for(int l=0;l<6;l++) delete t[l];
    if(v) for(int l=0;l<6;l++) delete v[l];
    if(txtD) delete txtD;
    if(txtR) delete txtR;
    if(lhot) delete lhot;

    /*
    cout << "Wrote output to file " << "TEST.root" << endl;
    rootfile->Close();
    */

    // Print to log File
    logFile = fopen(logFile_name.c_str(), "w");
    fprintf(logFile, "Run: %i\n", run_id);
    fprintf(logFile, "Fibre: %s\n", fibre_db.c_str());
    fprintf(logFile, "Is affected by belly plate?: %i\n", IS_BELLY_FIBRE);
    fprintf(logFile, "RATDB pos: %f\t%f\t%f\n", fibrepos.X(), fibrepos.Y(), fibrepos.Z());
    fprintf(logFile, "RATDB dir: %f\t%f\t%f\n", fibredir.X(), fibredir.Y(), fibredir.Z());
    fprintf(logFile, "Correct fibre?: %i\n", correctFibre);
    fprintf(logFile, "Final DIR fit: %f\t%f\t%f\n", dirfit->X(), dirfit->Y(), dirfit->Z());
    fprintf(logFile, "DIR fit angular deviation from expected: %f\n", DIRANG/pi*180);
    fprintf(logFile, "Final REF fit: %f\t%f\t%f\n", reffit->X(), reffit->Y(), reffit->Z());
    fprintf(logFile, "REF fit angular deviation from expected: %f\n", REFANG/pi*180);
    fprintf(logFile, "New fitted fibre direction: %f\t%f\t%f\n", fittedDir.X(), fittedDir.Y(), fittedDir.Z());
    fprintf(logFile, "Angular difference: %f\n", dirDiff);
    fclose(logFile);

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

        // PMT position check
        const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
        TVector3 pmtPos = pmtinfo.GetPosition(pmtID);       // position [mm]
        if (pmtPos.Mag()==0) { continue; }

        PMTHits[pmtID]++;
        allNHit++;

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

  /// Get intensity limit (hits/PMT) above which PMT is considered a screamer
  void MyUserProc::GetLimits(std::vector<double> occupancy, int NPMTS, double &HOTLIMIT, double &COLDLIMIT) {
    const int NBINS = 60;
    float MAX_NHIT = *std::max_element( occupancy.begin(), occupancy.end() );   // hottest PMT
    TH1F *hOccup = new TH1F("hOccup","",NBINS,0.,1.); // 0-100% occupancy
    BinLog(hOccup->GetXaxis(),1e-6); // minimum for log scale
    for(int id=0; id<NPMTS; id++) {
      if (occupancy[id]==0) continue;
      hOccup->Fill(occupancy[id]);
    }
    int i = hOccup->GetMaximumBin(); // start at max. bin
    int j = hOccup->GetMaximumBin();
    //cout << "Starting at bin=" << i << " with occupancy " << hOccup->GetBinCenter(i) << " and nevents " << hOccup->GetBinContent(i) << endl;
    while (hOccup->GetBinContent(i) >= 10) i++;
    while (hOccup->GetBinContent(j) >= 10) j--;
    HOTLIMIT = hOccup->GetBinLowEdge(i);
    COLDLIMIT = hOccup->GetBinLowEdge(j+1);
    //cout << "Stopping at bin=" << i << " with occupancy " << hOccup->GetBinCenter(i) << " and nevents " << hOccup->GetBinContent(i) << endl;
    if (hOccup) delete hOccup;
  }

  void MyUserProc::BinLog(TAxis *axis, const Double_t non0start) {

    const Int_t bins = axis->GetNbins();

    const Double_t xmin = axis->GetXmin();
    const Double_t xmax = axis->GetXmax();

    Bool_t k0start = kFALSE;
    if (xmin<EPSILON){
      k0start = kTRUE;
      if(non0start<EPSILON){
        printf("*** ERROR - BinLog bad non0start = %f\n", non0start); exit(1);
      }
    }

    Double_t *new_bins = new Double_t[bins + 1];

    const Double_t factor = k0start ? (TMath::Power(xmax/non0start, 1./(bins-1)))
                                    : (TMath::Power(xmax/xmin,      1./bins));

    new_bins[0] = xmin;
    new_bins[1] = k0start ? non0start : (new_bins[0]*factor);

    for (int i = 2; i <= bins; i++) {
      new_bins[i] = factor * new_bins[i-1];
    }
    axis->Set(bins, new_bins);
    delete [] new_bins;
  }

  void MyUserProc::GetMaxColVal(const TVector3& center, std::vector<double> occupancy, int NPMTS, float &nearval, float &farval, const RAT::DU::PMTInfo& pmtinfo) {
    const int NBINS = 60;
    TH1F *hNear = new TH1F("hNear","",NBINS,0,1);
    TH1F *hFar  = new TH1F("hFar","",NBINS,0,1);
    BinLog(hNear->GetXaxis(),1e-6);
    BinLog(hFar->GetXaxis(),1e-6);
    TVector3 pmtpos, newpos;
    for(int id=0; id<NPMTS; id++) {
      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue;                // not a valid PMT position
      if (pmtinfo.GetType(id) != 1) continue;       // not a normal PMT
      if (occupancy[id] > HOTLIMIT) continue;       // hot PMT
      if (occupancy[id] < COLDLIMIT) continue;      // cold PMT
      if (center.Angle(pmtpos)*180./pi <= 24)       // narrow cone around central point
        hNear->Fill(occupancy[id]);
      else if (center.Angle(pmtpos)*180./pi >= 156) // same around opposite side
        hFar->Fill(occupancy[id]);
    }
    int i=hNear->GetMaximumBin();
    int j=hFar->GetMaximumBin();
    //while (hNear->Integral(1,j)/hNear->Integral(1,NBINS) < 0.99) j++;
    //while (hFar->Integral(1,k)/hFar->Integral(1,NBINS) < 0.99) k++;
    while (hNear->GetBinContent(i) >= 1) i++;
    while (hFar->GetBinContent(j) >= 1) j++;
    nearval = hNear->GetXaxis()->GetBinLowEdge(i);
    farval = hFar->GetXaxis()->GetBinLowEdge(j);
    cout << "Nearval is " << nearval << endl; // TODO - replace this function with hard limits? e.g. nearval=NHit/8, farval=NHit/40
    if(hNear) delete hNear;
    if(hFar) delete hFar;
  }

  // -----------------------------------------------------------------------------
  /// Fit 2D Gaussian surface over intensities for PMTs in a plane
  void MyUserProc::FitLightSpot(TGraph2D* graph, double radius, double cone, double* params) {
    // Get input graph entries
    int npts = graph->GetN();
    double *xpts = graph->GetX();   // X projection
    double *ypts = graph->GetY();   // Y projection
    double *zpts = graph->GetZ();   // Intensity

    // First loop: Find maximum value in desired range
    double maxval=-1;
    for (int n=0; n<npts; n++) {
      double rpt = sqrt(xpts[n]*xpts[n] + ypts[n]*ypts[n]);
      if (rpt>radius) continue;
      double ang = asin(rpt/radius);
      if (ang/pi*180. > cone) continue;
      if (maxval<zpts[n]) maxval=zpts[n];
    }

    // Second loop: Fill values above 20% of maximum values into output graph
    TGraph2D *graf = new TGraph2D();
    int pts=0;
    for (int n=0; n<npts; n++) {
      if(zpts[n] < maxval*0.2) continue;
      double rpt = sqrt(xpts[n]*xpts[n] + ypts[n]*ypts[n]);
      if (rpt>radius) continue;
      double ang = asin(rpt/radius);
      if (ang/pi*180. > cone) continue;
      double scale = 1./cos(ang);
      graf->SetPoint(pts,scale*xpts[n],scale*ypts[n],zpts[n]);
      pts++;
    }

    // Fit 2D Gaussian surface to selected points
    TF2 *fit = new TF2("gaus2d","[0]*TMath::Gaus(x,[1],[3])*TMath::Gaus(y,[2],[3])",-10,10,-10,10);
    double aperture = radius*tan(cone/2. / 180.*pi);  // half input angle <=> 12 deg aperture
    fit->SetParameters(0.8*maxval,0.,0.,aperture);    // a priori: 80% max. intensity, centered, nominal aperture
    fit->SetParLimits(0,0.2*maxval,maxval);           // amplitude range [20%, 100%] of max. intensity
    fit->SetParLimits(1,-aperture,aperture);          // x-pos range [-12, +12] degrees from weighted centre
    fit->SetParLimits(2,-aperture,aperture);          // y-pos range [-12, +12] degrees from weighted centre
    fit->SetParLimits(3,0.5*aperture,1.5*aperture);   // sigma range [-50%, +50%] of nominal aperture
    graf->Fit("gaus2d","R,q");                        // fit gaussian (force range, quiet mode)

    // Raw fit parameters
    double amp = fit->GetParameter(0);
    double mux = fit->GetParameter(1);
    double muy = fit->GetParameter(2);
    double sig = fit->GetParameter(3);

    // Scale mean value back to sphere
    double rad = sqrt(mux*mux+muy*muy);
    double sx = mux*cos(atan(rad/radius));
    double sy = muy*cos(atan(rad/radius));

    // Scale (mean +- sigma) values back to sphere
    double xp = (mux+sig);
    double rxp = sqrt(xp*xp+muy*muy);
    double sxp = xp*cos(atan(rxp/radius));
    double xm = (mux-sig);
    double rxm = sqrt(xm*xm+muy*muy);
    double sxm = xm*cos(atan(rxm/radius));
    double yp = (muy+sig);
    double ryp = sqrt(yp*yp+mux*mux);
    double syp = yp*cos(atan(ryp/radius));
    double ym = (muy-sig);
    double rym = sqrt(ym*ym+mux*mux);
    double sym = ym*cos(atan(rym/radius));

    // Average over all scaled sigmas
    double sigma = (fabs(sxp-sx)+fabs(sxm-sx)+fabs(syp-sy)+fabs(sym-sy))/4.;
    //printf("Fit result scaled to sphere:\tp1 = %.3lf\tp2 = %.3lf\tp3 = %.3lf\n",sx,sy,sigma);

    // Set fit results
    params[0] = amp;    // amplitude
    params[1] = sx;     // mu_x
    params[2] = sy;     // mu_y
    params[3] = sigma;  // sigma

    if (graf) delete graf;
    if (fit) delete fit;
  }

  // -----------------------------------------------------------------------------
  /// Get rotation angles to rotate from given direction to xyz-frame
  void MyUserProc::GetRotationAngles(const TVector3& vec, double &rot_Z, double &rot_X) {

    // Get rotation angle from z-axis to central direction
    const TVector3 e1(1,0,0);
    const TVector3 e3(0,0,1);
    rot_Z = acos(e3*(vec.Unit()));         // rotate counter-clockwise, towards x-axis by [0,pi]

    // Get rotation angle from x-axis to central direction (projected to transverse plane)
    int sign = 0;
    if (vec[1] > 0) sign = +1;
    else            sign = -1;
    TVector3 vect(vec[0],vec[1],0);
    rot_X = sign*acos(e1*(vect.Unit()));   // rotate around z-axis by [-pi,pi]
  }

  // -----------------------------------------------------------------------------
  /// Create detector view from above central point
  void MyUserProc::FillHemisphere(const TVector3& center, std::vector<double> occupancy, int NPMTS, std::vector<TGraph> &dots, TGraph2D* graph, int NCOL, float MAXVAL, const RAT::DU::PMTInfo& pmtinfo) {

    std::vector<int> ndot;
    ndot.resize(NCOL+2);
    std::vector< std::vector<double> > dotx;
    std::vector< std::vector<double> > doty;
    dotx.resize(NCOL+2);
    doty.resize(NCOL+2);
    for (size_t i=0; i<NCOL+2; i++) {
      dotx[i].resize(NPMTS);
      doty[i].resize(NPMTS);
    }

    // Find rotation angles for this frame
    double rot_X, rot_Z;
    GetRotationAngles(center,rot_Z,rot_X);

    // Loop over all PMTs
    int counter=0;
    TVector3 pmtpos, newpos;
    for(int id=0; id<NPMTS; id++) {

      pmtpos = pmtinfo.GetPosition(id);
      if (pmtpos.Mag()==0) continue;                 // not a valid PMT position
      if (center.Angle(pmtpos) > pi/2.) continue;    // not in same hemisphere as central point

      // Rotate particles into correct frame
      newpos = pmtpos;
      newpos.RotateZ(-rot_X);
      newpos.RotateY(-rot_Z);

      // Get correct bin for colour scale
      int step;
      if (occupancy[id] == 0) step = 0;                 // off PMT
      else if (occupancy[id] < COLDLIMIT) step = 1;     // cold PMT
      else if (occupancy[id] > HOTLIMIT) step = NCOL+1; // hot PMT
      else if (occupancy[id] > MAXVAL) step = NCOL;     // cap color range
      // linear colour scale
      else step = (int)TMath::Ceil(occupancy[id]/(1.*MAXVAL/(NCOL-1)))+1;

      // Fill array of 1D graphs (bad practice) - TODO: replace this (currently plotted)
      dotx[step][ndot[step]]=newpos.X()/1e3;
      doty[step][ndot[step]]=newpos.Y()/1e3;
      ndot[step]++;

      // Fill 2D graph - TODO: more effective? (currently used for fit)
      if (newpos.Z() <= 0) continue;                 // not in hemisphere (safety check)
      if (pmtinfo.GetType(id) != 1) continue;        // not a normal PMT (remove OWLEs)
      else if (occupancy[id] == 0) continue;         // off PMT
      else if (occupancy[id] < COLDLIMIT) continue;  // cold PMT
      else if (occupancy[id] > HOTLIMIT) continue;   // hot PMT
      else if (occupancy[id] > MAXVAL) occupancy[id]=MAXVAL; // cap range
      double xpt = newpos.X()/1e3;//*cos(pi/4)/cos(newpos.Theta());
      double ypt = newpos.Y()/1e3;//*cos(pi/4)/cos(newpos.Theta());
      graph->SetPoint(counter, xpt, ypt, occupancy[id]);
      counter++;
    }

    for (int s=0; s<NCOL+2; s++) {
      int col;
      if      (s==0)  col=16; // off PMTs (grey)
      else if (s==1)  col=51; // cold PMTs (violet)
      else if (s==21) col=1;  // hot PMTs (black)
      else            col=(int)(50.+s*(50./NCOL)); // scale from 55-100 (s=2-20)
      if (ndot[s]==0) {continue;}
      dots[s] = TGraph(ndot[s],&dotx[s][0],&doty[s][0]);
      dots[s].SetMarkerStyle(7);
      dots[s].SetMarkerColor(col);
    }

  }

  // -----------------------------------------------------------------------------
  /// Draw a circle around a point in a plane orthogonal to an angle with N dots
  void MyUserProc::DrawCircle(const TVector3& center, double angle, std::vector<TVector3> &dots, int NDOTS) {
    TVector3 dot(0,0,0);
    for (int i=0; i<NDOTS; i++) {
      dot.SetMagThetaPhi(center.Mag(), pi*angle/180., i*2.*pi/NDOTS);
      dot.RotateUz(center.Unit());   // rotate into frame of input vector
      dots[i].SetXYZ(dot.X(),dot.Y(),dot.Z());
    }
  }

  // -----------------------------------------------------------------------------
  /// Calculate what fibre was firing
  void MyUserProc::CorrectFibre(const TVector3& fit_origin, double min_dist) {
    int maxNumFibre = fibresInUse.size();
    std::vector<string> likelyFibre;
    for ( int i = 0; i < maxNumFibre; i++) {
      string fibNum = fibresInUse[i];
      TVector3 thisFibrePos(0,0,0);
      try {
        fibre_data = DB::Get()->GetLink("FIBRE", fibNum);
        thisFibrePos.SetXYZ(fibre_data->GetD("x"), fibre_data->GetD("y"), fibre_data->GetD("z"));
      }
      catch (DBNotFoundError& e) {
        std::cout << "Proc(): Couldn't load FIBRE details from ratdb!" << std::endl;
      }
      double dist = (fit_origin-thisFibrePos).Mag();
      if (dist < min_dist) {
        likelyFibre.push_back( fibNum );
      }
    }
    for (size_t i=0; i<likelyFibre.size(); i++) {
      cout << likelyFibre[i] << endl;
      if (likelyFibre[i] == fibre_db){ correctFibre = 1; }
    }
  }

  // -----------------------------------------------------------------------------
  /// Gets a list of fibres currently in use (patched)
  std::vector<string> MyUserProc::ListFibres() {
    std::vector<string> fibres;
    try {
      fibre_data = DB::Get()->GetLink("TELLIE_PATCH_MAPPING");
      fibres = fibre_data->GetSArray("fibres");
    }
    catch (DBNotFoundError& e) {
      std::cout << "Proc(): Couldn't load FIBRE details from ratdb!" << std::endl;
    }
    return fibres;
  }
