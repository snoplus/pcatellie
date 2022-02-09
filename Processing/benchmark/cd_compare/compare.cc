//C++ includes
#include "Riostream.h"
#include "sys/stat.h"
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>

//ROOT includes
#include "TGraph.h"
#include "TGraph2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TKey.h"
#include "TImage.h"
#include "TView3D.h"
#include "TVectorD.h"
#include "TPaveText.h"

//RAT includes
#include <RAT/DB.hh>
#include <RAT/DBLink.hh>
#include <RAT/DBTable.hh>
#include <RAT/BitManip.hh>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/PCABits.hh>

using namespace std;
using namespace RAT;
using namespace RAT::DS;

TVector2 TransformCoord( const TVector3& V1, const TVector3& V2, const TVector3& V3, const TVector2& A1, const TVector2& A2, const TVector2& A3,const TVector3& P ){
  TVector3 xV = V2 - V1;
  TVector3 yV = ( ( V3 - V1 ) + ( V3 - V2 ) ) * 0.5;
  TVector3 zV = xV.Cross( yV ).Unit();

  double planeD = V1.Dot( zV );

  double t = planeD / P.Dot( zV );

  TVector3 localP = t*P - V1;

  TVector2 xA = A2 - A1;
  TVector2 yA = ( ( A3 - A1 ) +( A3 - A2 ) ) * 0.5;

  double convUnits = xA.Mod() / xV.Mag();

  TVector2 result;
  result = localP.Dot( xV.Unit() ) * xA.Unit() * convUnits;
  result += localP.Dot( yV.Unit() ) * yA.Unit() * convUnits + A1;
  return result;
}

TVector2 IcosProject( TVector3 pmtPos , int &segment ){

    double fa; ///<constant for psup projection
    double fb; ///<constant for psup projection

    TVector2 *fA12a; ///<vector for psup projection
    TVector2 *fA12b; ///<vector for psup projection
    TVector2 *fA12c; ///<vector for psup projection
    TVector2 *fA12d; ///<vector for psup projection
    TVector2 *fA12e; ///<vector for psup projection
    TVector2 *fA2a; ///<vector for psup projection
    TVector2 *fA2b; ///<vector for psup projection
    TVector2 *fA17a; ///<vector for psup projection
    TVector2 *fA17b; ///<vector for psup projection
    TVector2 *fA51a; ///<vector for psup projection
    TVector2 *fA51b; ///<vector for psup projection
    TVector2 *fA51c; ///<vector for psup projection
    TVector2 *fA51d; ///<vector for psup projection
    TVector2 *fA51e; ///<vector for psup projection
    TVector2 *fA27; ///<vector for psup projection
    TVector2 *fA46; ///<vector for psup projection
    TVector2 *fA31; ///<vector for psup projection
    TVector2 *fA6; ///<vector for psup projection
    TVector2 *fA37; ///<vector for psup projection
    TVector2 *fA33; ///<vector for psup projection
    TVector2 *fA58; ///<vector for psup projection
    TVector2 *fA54; ///<vector for psup projection

    fa = 1.0 / 5.5;
    fb = fa * sqrt( 3.0 ) / 2.0;

    fA12a = new TVector2( fa / 2.0, 0.0 );
    fA12b = new TVector2( 3.0 * fa / 2.0, 0.0 );
    fA12c = new TVector2( 5.0 * fa / 2.0, 0.0 );
    fA12d = new TVector2( 7.0 *fa / 2.0, 0.0 );
    fA12e = new TVector2( 9.0 * fa / 2.0, 0.0 );
    fA2a = new TVector2( 0.0, fb );
    fA2b = new TVector2( 5.0 * fa, fb );
    fA17a = new TVector2( fa / 2.0 , 2.0 * fb );
    fA17b = new TVector2( 11.0 * fa / 2.0 , 2.0 * fb );
    fA51a = new TVector2( fa, 3.0 * fb );
    fA51b = new TVector2( 2.0 * fa, 3.0 * fb );
    fA51c = new TVector2( 3.0 * fa, 3.0 * fb );
    fA51d = new TVector2( 4.0 * fa, 3.0 * fb );
    fA51e = new TVector2( 5.0 * fa, 3.0 * fb );
    fA27 = new TVector2( 4.0 * fa, fb );
    fA46 = new TVector2( 3.0 * fa, fb );
    fA31 = new TVector2( 2.0 * fa, fb );
    fA6 = new TVector2( fa, fb );
    fA37 = new TVector2( 9.0 * fa / 2.0 , 2.0 * fb );
    fA33 = new TVector2( 3.0 * fa / 2.0 , 2.0 * fb );
    fA58 = new TVector2( 5.0 * fa / 2.0 , 2.0 * fb );
    fA54 = new TVector2( 7.0 * fa / 2.0 , 2.0 * fb );

    TVector3 pointOnSphere( pmtPos.X(), pmtPos.Y(), pmtPos.Z() );
    pointOnSphere = pointOnSphere.Unit();
    pointOnSphere.RotateX( -45.0 );
    // From http://www.rwgrayprojects.com/rbfnotes/polyhed/PolyhedraData/Icosahedron/Icosahedron.pdf
    const double t = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
    const TVector3 V2 = TVector3( t * t, 0.0, t * t * t ).Unit();
    const TVector3 V6 = TVector3( -t * t, 0.0, t * t * t ).Unit();
    const TVector3 V12 = TVector3( 0.0, t * t * t, t * t ).Unit();
    const TVector3 V17 = TVector3( 0.0, -t * t * t, t * t ).Unit();
    const TVector3 V27 = TVector3( t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V31 = TVector3( -t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V33 = TVector3( -t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V37 = TVector3( t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V46 = TVector3( 0.0, t * t * t, -t * t ).Unit();
    const TVector3 V51 = TVector3( 0.0, -t * t * t, -t * t ).Unit();
    const TVector3 V54 = TVector3( t * t, 0.0, -t * t * t ).Unit();
    const TVector3 V58 = TVector3( -t * t, 0.0, -t * t * t ).Unit();

    // Faces {{ 2, 6, 17}, { 2, 12, 6}, { 2, 17, 37}, { 2, 37, 27}, { 2, 27, 12}, {37, 54, 27},
    // {27, 54, 46}, {27, 46, 12}, {12, 46, 31}, {12, 31, 6}, { 6, 31, 33}, { 6, 33, 17},
    // {17, 33, 51}, {17, 51, 37}, {37, 51, 54}, {58, 54, 51}, {58, 46, 54}, {58, 31, 46},
    // {58, 33, 31}, {58, 51, 33}}

    vector<TVector3> IcosahedralCentres;
    IcosahedralCentres.push_back( ( V2 + V6 + V17 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V12 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V17 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V37 + V27 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V27 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V54 + V27 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V27 + V54 + V46 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V27 + V46 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V46 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V31 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V31 + V33 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V33 + V17 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V17 + V33 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V17 + V51 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V51 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V54 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V46 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V31 + V46 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V33 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V51 + V33 ) * ( 1.0 / 3.0 ) );

    vector<double> distFromCentre;
    unsigned int uLoop;
    for( uLoop = 0; uLoop < IcosahedralCentres.size(); uLoop++ ){
      distFromCentre.push_back( ( IcosahedralCentres[uLoop] - pointOnSphere ).Mag() );
    }
    const int face = min_element( distFromCentre.begin(), distFromCentre.end() ) - distFromCentre.begin() + 1;
    TVector2 resultPosition;
    switch(face){
    case 1://{ 2, 6, 17}
      resultPosition = TransformCoord( V2, V6, V17, *fA2a, *fA6, *fA17a, pointOnSphere );
      break;
    case 2://{ 2, 12, 6}
      resultPosition = TransformCoord( V2, V12, V6, *fA2a, *fA12a, *fA6, pointOnSphere );
      break;
    case 3://{ 2, 17, 37}
      resultPosition = TransformCoord( V2, V17, V37, *fA2b, *fA17b, *fA37, pointOnSphere );
      break;
    case 4://{ 2, 37, 27}
      resultPosition = TransformCoord( V2, V37, V27, *fA2b, *fA37, *fA27, pointOnSphere );
      break;
    case 5://{ 2, 27, 12}
      resultPosition = TransformCoord( V2, V27, V12, *fA2b, *fA27, *fA12e, pointOnSphere );
      break;
    case 6://{37, 54, 27}
      resultPosition = TransformCoord( V37, V54, V27, *fA37, *fA54, *fA27, pointOnSphere );
      break;
    case 7://{27, 54, 46}
      resultPosition = TransformCoord( V27, V54, V46, *fA27, *fA54, *fA46, pointOnSphere );
      break;
    case 8://{27, 46, 12}
      resultPosition = TransformCoord( V27, V46, V12, *fA27, *fA46, *fA12d, pointOnSphere );
      break;
    case 9://{12, 46, 31}
      resultPosition = TransformCoord( V12, V46, V31, *fA12c, *fA46, *fA31, pointOnSphere );
      break;
    case 10://{12, 31, 6}
      resultPosition = TransformCoord( V12, V31, V6, *fA12b, *fA31, *fA6, pointOnSphere );
      break;
    case 11://{ 6, 31, 33}
      resultPosition = TransformCoord( V6, V31, V33, *fA6, *fA31, *fA33, pointOnSphere );
      break;
    case 12://{ 6, 33, 17}
      resultPosition = TransformCoord( V6, V33, V17, *fA6, *fA33, *fA17a, pointOnSphere );
      break;
    case 13://{17, 33, 51}
      resultPosition = TransformCoord( V17, V33, V51, *fA17a, *fA33, *fA51a, pointOnSphere );
      break;
    case 14://{17, 51, 37}
      resultPosition = TransformCoord( V17, V51, V37, *fA17b, *fA51e, *fA37, pointOnSphere );
      break;
    case 15://{37, 51, 54}
      resultPosition = TransformCoord( V37, V51, V54, *fA37, *fA51d, *fA54, pointOnSphere );
      break;
    case 16://{58, 54, 51}
      resultPosition = TransformCoord( V58, V54, V51, *fA58, *fA54, *fA51c, pointOnSphere );
      break;
    case 17://{58, 46, 54}
      resultPosition = TransformCoord( V58, V46, V54, *fA58, *fA46, *fA54, pointOnSphere );
      break;
    case 18://{58, 31, 46}
      resultPosition = TransformCoord( V58, V31, V46, *fA58, *fA31, *fA46, pointOnSphere );
      break;
    case 19://{58, 33, 31}
      resultPosition = TransformCoord( V58, V33, V31, *fA58, *fA33, *fA31, pointOnSphere );
      break;
    case 20://{58, 51, 33}
      resultPosition = TransformCoord( V58, V51, V33, *fA58, *fA51b, *fA33, pointOnSphere );
      break;
    }
    // output face, if needed
    segment = face;
    // 1 - x,y pos to project the same as node map
    return TVector2(1.0 - resultPosition.X(), 1.0 - 2.0 * resultPosition.Y() );
}


int main(int argc, char* argv[]) {

  if (argc != 5){
    std::cout<<" Usage: Compare TW1.ratdb TW2.ratdb run1 run2"<<std::endl;
    exit(0);
  }

  string TWFilename1, TWFilename2;
  int runid1, runid2;
  TWFilename1 = argv[1];
  TWFilename2 = argv[2];
  runid1 = atoi(argv[3]);
  runid2 = atoi(argv[4]);

  DB *ratdb = RAT::DB::Get();
  ratdb->LoadDefaults();

  std::cout<<" Loading first constants "<<std::endl;

  DS::Run run1;
  run1.SetRunID(runid1);

  try{
    ratdb->Load(TWFilename1,true);
    ratdb->BeginOfRun(run1);
  }
  catch(...){
    std::cout<< " Problem getting constants. Exit!" << std::endl;
    exit(0);
  }

  RAT::DU::Utility::LoadDBAndBeginRun();
  const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();


  PCABits fPCABits;
  fPCABits.BeginOfRun(run1);

  TView3D *view;
  std::vector<double> data;

  DBLinkPtr TWbank1 = DB::Get()->GetLink("PCA_TW");
  std::cout << " Analyzing PCA TW calibration: run" << TWbank1->GetIArray("run_range")[0] << std::endl;
  std::vector<int> TWstatus1 = TWbank1->GetIArray("PCATW_status");
  std::vector<int> TWpoints1 = TWbank1->GetIArray("twinter");
  std::vector<float> TWpointErrs1 = TWbank1->GetFArrayFromD("twinterrms");
  //std::vector<float> TWpointErrs1mean = TWbank1->GetFArrayFromD("meanerr");
  ratdb->EndOfRun(run1);
  std::cout<<" Done "<<std::endl;
  ratdb->Clear();
  std::cout << " Clear RATDB " << std::endl;

  std::cout<<" Loading second constants "<<std::endl;

  DS::Run run2;
  run2.SetRunID(runid2);

  try{
    ratdb->Load(TWFilename2);
    ratdb->BeginOfRun(run2);
  }
  catch(...){
    std::cout<< " Problem getting constants. Exit!" << std::endl;
    exit(0);
  }
  DBLinkPtr TWbank2 = DB::Get()->GetLink("PCA_TW");
  std::cout << " Analyzing PCA TW calibration: run" << TWbank2->GetIArray("run_range")[0] << std::endl;
  std::vector<int> TWstatus2 = TWbank2->GetIArray("PCATW_status");
  std::vector<int> TWpoints2 = TWbank2->GetIArray("twinter");
  std::vector<float> TWpointErrs2 = TWbank2->GetFArrayFromD("twinterrms");
  //std::vector<float> TWpointErrs2mean = TWbank2->GetFArrayFromD("meanerr");
  std::cout<<" Done "<<std::endl;

  cout << TWstatus2.size() << " " << TWpoints2.size() << " " << TWpointErrs2.size() << endl;

  // First unpack the ocnstants
  // unpack the data and store in new arrays,this is a (q,t) interpolation point
  double thisr, thisi, thisi_old, thisr_old;
  int loc_word_fit;
  int packresult;

  double min_cable = 9000;
  double max_cable = 0;
  double min_cable2 = 9000;
  double max_cable2 = 0;
  double min_diff = 9000;
  double max_diff = 0;

  // Plots
  TCanvas *cableC = new TCanvas("cableC","cableC",800,600);
  TCanvas *cable2C = new TCanvas("cable2C","cable2C",800,600);
  TCanvas *cableCardC = new TCanvas("cableCardC","cableCardC",800,600);
  TCanvas *cableCard2C = new TCanvas("cableCard2C","cableCard2C",800,600);
  TCanvas *cableDiffC = new TCanvas("cableDiffC","cableDiffC",800,600);
  TCanvas *cableDiffCardC = new TCanvas("cableDiffCardC","cableDiffCardC",800,600);
  TCanvas *cableDiffHistC = new TCanvas("cableDiffHistC","cableDiffHistC",800,600);
  TCanvas *cableFlatC = new TCanvas("cableFlatC","cableFlatC",800,600);
  TCanvas *cableFlatZoomC = new TCanvas("cableFlatZoomC","cableFlatZoomC",800,600);
  TCanvas *cableFlat2C = new TCanvas("cableFlat2C","cableFlat2C",800,600);
  TCanvas *cableFlatDiffC = new TCanvas("cableFlatDiffC","cableFlatDiffC",800,600);
  TCanvas *cableFlatDiffZoomC = new TCanvas("cableFlatDiffZoomC","cableFlatDiffZoomC",800,600);
  TH2D* cable = new TH2D("cable", "cable", 32*16*19, 0, 32*16*19, 500, -250, 250);
  TH2D* cable2 = new TH2D("cable2", "cable2", 32*16*19, 0, 32*16*19, 500, -250, 250);
  TH2D* cableCard = new TH2D("cableCard", "cableCard", 16, 0, 16, 500, -250, 250);
  TH2D* cableCard2 = new TH2D("cableCard2", "cableCard2", 16, 0, 16, 500, -250, 250);
  TH2D* cableDiff = new TH2D("cableDiff", "cableDiff", 32*16*19, 0, 32*16*19, 500, -250, 250);
  TH2D* cableDiffCard = new TH2D("cableDiffCard", "cableDiffCard", 16, 0, 16, 500, -250, 250);
  TH1D* cableDiffHist = new TH1D("cableDiffHist", "cableDiffHist", 100, -10, 10);
  TGraph2D* cableFlat = new TGraph2D();
  TGraph2D* cableFlatZoom = new TGraph2D();
  TGraph2D* cableFlat2 = new TGraph2D();
  TGraph2D* cableFlatDiff = new TGraph2D();
  TGraph2D* cableFlatDiffZoom = new TGraph2D();
  TCanvas* fSummaryC = new TCanvas("fSummaryC","fSummaryC",800,600);
  TPaveText *fSummary = new TPaveText(.01,.01,.99,.99);

  int wOff = 0;
  int wZO = 0;
  int wLO = 0;
  int goodCal = 0;
  int otherPMTs = 0;

  std::vector<int> out_ID;
  std::vector<double> temp_diff;
  std::vector< std::vector<double> > out_diff;
  std::vector<int> nowGood;
  std::vector<int> nowBad;

  for (int i = 0; i < 9728; i++) {
    int card = ( i%(16*32) )/32 ;

    int face;
    TVector3 pmtPos = pmtinfo.GetPosition(i);
    TVector2 flat = IcosProject(pmtPos, face);
    double Xpos = flat.X();
    double Ypos = flat.Y();
    if ((Xpos == 0) | (Ypos == 0)){continue;}
    //std::cout << i << ": " << Xpos << " " << Ypos << endl;

    if (pmtinfo.GetType(i) != 1){
      otherPMTs++;
    }

    loc_word_fit = i * 11 + 10; // Gets last valuestored for each PMT (assuming 11 are stored for each)
    packresult = fPCABits.PCAPacker(fPCABits.UNPACK, fPCABits.GRADIENT_INTERCEPT, TWpoints1[loc_word_fit], thisi, thisr);
    //cout << i << " " << thisi << endl;
    std::string binary = std::bitset<32>(TWstatus1[i]).to_string();
    //cout << binary << endl;
    int isOff = atoi(&binary[30]);
    int isZO = atoi(&binary[29]);
    int isLO = atoi(&binary[28]);
    int statNow = 0;
    if (isOff != 1){ wOff++;
      if (pmtinfo.GetType(i) != 1){
        otherPMTs++;
        cableFlat->SetPoint(i, Xpos, Ypos, -17);
        cableFlatZoom->SetPoint(i, Xpos, Ypos, -17);
      }
      cableFlat->SetPoint(i, Xpos, Ypos, -20);
      cableFlatZoom->SetPoint(i, Xpos, Ypos, -20);
    }
    if (isZO != 1){ wZO++; }
    if (isLO != 1){ wLO++; }
    if ( (isOff == 1) && (isZO == 1) && (isLO == 1) ){ goodCal++; statNow = 1;}

    if(statNow == 1){
      cable->Fill(i, thisi);
      cableCard->Fill(card, thisi);
      if (thisi > max_cable) max_cable = thisi;
      if (thisi < min_cable) min_cable = thisi;
      cableFlat->SetPoint(i, Xpos, Ypos, thisi);
      if (thisi < 10){
        cableFlatZoom->SetPoint(i, Xpos, Ypos, thisi);
      } else { cout << "PMT outlier!!! ID: " << i << ", delay: " << thisi << endl; }
    }

    packresult = fPCABits.PCAPacker(fPCABits.UNPACK, fPCABits.GRADIENT_INTERCEPT, TWpoints2[loc_word_fit], thisi_old, thisr_old);
    std::string binary_old = std::bitset<32>(TWstatus2[i]).to_string();
    int isOff_old = atoi(&binary_old[30]);
    int isZO_old = atoi(&binary_old[29]);
    int isLO_old = atoi(&binary_old[28]);
    int statOld = 0;
    if ( (isOff_old == 1) && (isZO_old == 1) && (isLO_old == 1) ){ statOld = 1;}
    if (isOff_old != 1){
      if (pmtinfo.GetType(i) != 1){
        cableFlat2->SetPoint(i, Xpos, Ypos, -17);
      }
      cableFlat2->SetPoint(i, Xpos, Ypos, -20);
    }
    if(statOld == 1){
      cable2->Fill(i, thisi_old);
      cableCard2->Fill(card, thisi_old);
      if (thisi_old > max_cable2) max_cable2 = thisi_old;
      if (thisi_old < min_cable2) min_cable2 = thisi_old;
      cableFlat2->SetPoint(i, Xpos, Ypos, thisi_old);
    };

    if ( (statNow == 1) && (statOld == 1) ) {
      if (thisi-thisi_old > max_diff) max_diff = thisi-thisi_old;
      if (thisi-thisi_old < min_diff) min_diff = thisi-thisi_old;
      cableDiff->Fill(i, thisi-thisi_old);
      cableDiffCard->Fill(card, thisi-thisi_old);
      cableDiffHist->Fill(thisi-thisi_old);
      cableFlatDiff->SetPoint(i, Xpos, Ypos, thisi-thisi_old);
      if (thisi < 10){
        cableFlatDiffZoom->SetPoint(i, Xpos, Ypos, thisi-thisi_old);
      }

      // New
      if ( (thisi-thisi_old > 10.) | (thisi-thisi_old < -10.) ){
        cout << "WARNING: difference is greater than 10, PMT ID = " << i << endl;
        out_ID.push_back( i );
        temp_diff.push_back( thisi );
        temp_diff.push_back( thisi_old );
        out_diff.push_back( temp_diff );
        temp_diff.clear();
      }
    } else {
      if (pmtinfo.GetType(i) != 1){
        cableFlatDiff->SetPoint(i, Xpos, Ypos, -17);
        cableFlatDiffZoom->SetPoint(i, Xpos, Ypos, -17);
      } else {
      cableFlatDiff->SetPoint(i, Xpos, Ypos, -20);
      cableFlatDiffZoom->SetPoint(i, Xpos, Ypos, -20);
      }
    }

    if ( (statNow == 1) & (statOld == 0) ){
      nowGood.push_back( i );
    } else if ( (statNow == 0) & (statNow == 1) ) {
      nowBad.push_back( i );
    }

  }

  cout << "Summary" << endl;
  cout << "word OFF: " << wOff << endl;
  cout << "word ZERO OCC: " << wZO << endl;
  cout << "word LOW OCC: " << wLO << endl;
  cout << "word goodCal: " << goodCal << endl;
  cout << "otherPMTs: " << otherPMTs << endl;
  cout << "nowGood: " << nowGood.size() << endl;
  cout << "nowBad: " << nowBad.size() << endl;

  cout << min_diff << " " << max_diff << endl;

  // Store to log file
  stringstream logFile_namess;
  string logFile_name;
  FILE *logFile;
  logFile_namess.str("");
  logFile_namess << runid1 << "-" << runid2 << "_cd_compare.log";
  logFile_name = logFile_namess.str();
  logFile = fopen(logFile_name.c_str(), "w");

  fprintf(logFile, "Run1: %i\n", runid1);
  fprintf(logFile, "Run2: %i\n", runid2);
  fprintf(logFile, "PMTOff: %i\n", wOff);
  fprintf(logFile, "PMTZeroOc: %i\n", wZO);
  fprintf(logFile, "PMTLowOc: %i\n", wLO);
  fprintf(logFile, "PMTGoodCal: %i\n", goodCal);
  fprintf(logFile, "minDiff: %f\n", min_diff);
  fprintf(logFile, "maxDiff: %f\n", max_diff);

  if (out_ID.size() > 1){
    fprintf(logFile, "out_ID: ");
    for (size_t i=0; i<out_ID.size()-1; i++) {
      fprintf(logFile, "%i, ", out_ID[i]);
    }
    fprintf(logFile, "%i\n", out_ID[out_ID.size()-1]);
  } else if (out_ID.size() == 1){
    fprintf(logFile, "out_ID: %i\n", out_ID[0]);
  } else {
    fprintf(logFile, "out_ID: 0\n");
  }

  if (out_diff.size() > 1){
    fprintf(logFile, "out_diff: ");
    for (size_t i=0; i<out_diff.size()-1; i++) {
      fprintf(logFile, "[%f, %f], ", out_diff[i][0], out_diff[i][1]);
    }
    fprintf(logFile, "[%f, %f]\n", out_diff[out_diff.size()-1][0], out_diff[out_diff.size()-1][1]);
  } else if (out_diff.size() == 1){
    fprintf(logFile, "out_diff: %i\n", out_diff[0]);
  } else {
    fprintf(logFile, "out_diff: 0\n");
  }

  if (nowGood.size() > 1){
    fprintf(logFile, "nowGood: ");
    //cout << nowGood.size() << endl;
    for (size_t i=0; i<nowGood.size()-1; i++) {
      fprintf(logFile, "%i, ", nowGood[i]);
    }
    fprintf(logFile, "%i\n", nowGood[nowGood.size()-1]);
  } else if (nowGood.size() == 1){
    fprintf(logFile, "nowGood: %i\n", nowGood[0]);
  } else {
    fprintf(logFile, "nowGood: 0\n");
  }

  if (nowBad.size() > 1){
    fprintf(logFile, "nowBad: ");
    for (size_t i=0; i<nowBad.size()-1; i++) {
      fprintf(logFile, "%i, ", nowBad[i]);
    }
    fprintf(logFile, "%i\n", nowBad[nowBad.size()-1]);
  } else if (nowBad.size() == 1){
    fprintf(logFile, "nowBad: %i\n", nowBad[0]);
  } else {
    fprintf(logFile, "nowBad: 0\n");
  }

  //TFile *ofile = new TFile(Form("PCA_cable_comparison_%i_%i.root",runid1,runid2),"RECREATE");
  //ofile->cd();

  cableC->cd();
  cable->SetStats(0);
  cable->GetXaxis()->SetTitle("Cable delays");
  cable->GetXaxis()->SetTitleOffset(1.0);
  cable->GetYaxis()->SetTitleOffset(1.0);
  cable->GetYaxis()->SetTitle("Delay [ns]");
  cable->GetYaxis()->SetRangeUser(min_cable-0.2, max_cable+0.2);
  cable->Draw();
  cableC->SaveAs("cd_comp_cable.png");
  cableC->Close();

  cable2C->cd();
  cable2->SetStats(0);
  cable2->GetXaxis()->SetTitle("Cable delays 2");
  cable2->GetXaxis()->SetTitleOffset(1.0);
  cable2->GetYaxis()->SetTitleOffset(1.0);
  cable2->GetYaxis()->SetTitle("Delay [ns]");
  cable2->GetYaxis()->SetRangeUser(min_cable2-0.2, max_cable2+0.2);
  cable2->Draw();
  cable2C->SaveAs("cd_comp_cable2.png");
  cable2C->Close();

  cableCardC->cd();
  cableCard->SetStats(0);
  cableCard->GetXaxis()->SetTitle("Cable delays - Card");
  cableCard->GetXaxis()->SetTitleOffset(1.0);
  cableCard->GetYaxis()->SetTitleOffset(1.0);
  cableCard->GetYaxis()->SetTitle("Delay [ns]");
  cableCard->GetYaxis()->SetRangeUser(min_cable-0.2, max_cable+0.2);
  cableCard->Draw();
  cableCardC->SaveAs("cd_comp_cableCard.png");
  cableCardC->Close();

  cableCard2C->cd();
  cableCard2->SetStats(0);
  cableCard2->GetXaxis()->SetTitle("Cable delays - Card");
  cableCard2->GetXaxis()->SetTitleOffset(1.0);
  cableCard2->GetYaxis()->SetTitleOffset(1.0);
  cableCard2->GetYaxis()->SetTitle("Delay [ns]");
  cableCard2->GetYaxis()->SetRangeUser(min_cable2-0.2, max_cable2+0.2);
  cableCard2->Draw();
  cableCard2C->SaveAs("cd_comp_cableCard2.png");
  cableCard2C->Close();

  cableDiffC->cd();
  //cableDiff->SetStats(0);
  cableDiff->GetXaxis()->SetTitle("Cable delays diff");
  cableDiff->GetXaxis()->SetTitleOffset(1.0);
  cableDiff->GetYaxis()->SetTitleOffset(1.0);
  cableDiff->GetYaxis()->SetTitle("Delay diff [ns]");
  cableDiff->GetYaxis()->SetRangeUser(min_diff-0.1, max_diff+0.1);
  cableDiff->Draw();
  cableDiffC->SaveAs("cd_comp_cableDiff.png");
  cableDiffC->Close();

  cableDiffCardC->cd();
  //cableDiffCard->SetStats(0);
  cableDiffCard->GetXaxis()->SetTitle("Cable delays diff - Card");
  cableDiffCard->GetXaxis()->SetTitleOffset(1.0);
  cableDiffCard->GetYaxis()->SetTitleOffset(1.0);
  cableDiffCard->GetYaxis()->SetTitle("Delay diff [ns]");
  cableDiffCard->GetYaxis()->SetRangeUser(min_diff-0.1, max_diff+0.1);
  cableDiffCard->Draw();
  cableDiffCardC->SaveAs("cd_comp_cableDiffCard.png");
  cableDiffCardC->Close();

  cableDiffHistC->cd();
  //cableDiffHist->SetStats(0);
  cableDiffHist->GetXaxis()->SetTitle("Cable delays diff - Card");
  cableDiffHist->GetXaxis()->SetTitleOffset(1.0);
  cableDiffHist->GetYaxis()->SetTitleOffset(1.0);
  cableDiffHist->GetYaxis()->SetTitle("Delay diff [ns]");
  //cableDiffHist->GetYaxis()->SetRangeUser(min_diff-0.1, max_diff+0.1);
  cableDiffHist->Draw();
  cableDiffHistC->SaveAs("cd_comp_cableDiffHist.png");
  cableDiffHistC->Close();

  view = new TView3D();
  cableFlatC->cd();
  cableFlatC->SetMargin(0.05,0.15,0.05,0.05);
  cableFlat->GetZaxis()->SetTitle("cable delay");
  cableFlat->SetMarkerStyle(8);
  cableFlat->SetMarkerSize(0.7);
  cableFlat->GetXaxis()->SetLabelSize(0.03);
  cableFlat->GetYaxis()->SetLabelSize(0.03);
  cableFlat->SetTitle("Flat Map of delays");
  view->RotateView(270, 0);
  gStyle->SetNumberContours(200);
  cableFlat->Draw("AZCOLPCOL");
  cableFlatC->SaveAs("cd_comp_cableFlat.png");
  cableFlatC->Close();

  view = new TView3D();
  cableFlatZoomC->cd();
  cableFlatZoomC->SetMargin(0.05,0.15,0.05,0.05);
  cableFlatZoom->GetZaxis()->SetTitle("cable delay");
  cableFlatZoom->SetMarkerStyle(8);
  cableFlatZoom->SetMarkerSize(0.7);
  cableFlatZoom->GetXaxis()->SetLabelSize(0.03);
  cableFlatZoom->GetYaxis()->SetLabelSize(0.03);
  cableFlatZoom->SetTitle("Flat Map of delays");
  view->RotateView(270, 0);
  gStyle->SetNumberContours(200);
  cableFlatZoom->Draw("AZCOLPCOL");
  cableFlatZoomC->SaveAs("cd_comp_cableFlatZoom.png");
  cableFlatZoomC->Close();

  cableFlat2C->cd();
  cableFlat2C->SetMargin(0.05,0.15,0.05,0.05);
  cableFlat2->GetZaxis()->SetTitle("cable delay");
  cableFlat2->SetMarkerStyle(8);
  cableFlat2->SetMarkerSize(0.7);
  cableFlat2->GetXaxis()->SetLabelSize(0.03);
  cableFlat2->GetYaxis()->SetLabelSize(0.03);
  cableFlat2->SetTitle("Flat Map of delays");
  view->RotateView(270, 0);
  gStyle->SetNumberContours(200);
  cableFlat2->Draw("AZCOLPCOL");
  cableFlat2C->SaveAs("cd_comp_cableFlat2.png");
  cableFlat2C->Close();

  cableFlatDiffC->cd();
  cableFlatDiffC->SetMargin(0.05,0.15,0.05,0.05);
  cableFlatDiff->GetZaxis()->SetTitle("cable delay");
  cableFlatDiff->SetMarkerStyle(8);
  cableFlatDiff->SetMarkerSize(0.7);
  cableFlatDiff->GetXaxis()->SetLabelSize(0.03);
  cableFlatDiff->GetYaxis()->SetLabelSize(0.03);
  cableFlatDiff->SetTitle("Flat Map of delays");
  view->RotateView(270, 0);
  gStyle->SetNumberContours(200);
  cableFlatDiff->Draw("AZCOLPCOL");
  cableFlatDiffC->SaveAs("cd_comp_cableFlatDiff.png");
  cableFlatDiffC->Close();

  cableFlatDiffZoomC->cd();
  cableFlatDiffZoomC->SetMargin(0.05,0.15,0.05,0.05);
  cableFlatDiffZoom->GetZaxis()->SetTitle("cable delay");
  cableFlatDiffZoom->SetMarkerStyle(8);
  cableFlatDiffZoom->SetMarkerSize(0.7);
  cableFlatDiffZoom->GetXaxis()->SetLabelSize(0.03);
  cableFlatDiffZoom->GetYaxis()->SetLabelSize(0.03);
  cableFlatDiffZoom->SetTitle("Flat Map of delays");
  view->RotateView(270, 0);
  gStyle->SetNumberContours(200);
  cableFlatDiffZoom->Draw("AZCOLPCOL");
  cableFlatDiffZoomC->SaveAs("cd_comp_cableFlatDiffZoom.png");
  cableFlatDiffZoomC->Close();

  // data.resize(4);
  // data[0] = offPMTs;
  // data[1] = Calibrated;
  // data[2] = UNCalibrated;
  // data[3] = CANTbe;
  // //ofile->WriteObject(&data, "data");
  //
  // fSummaryC->cd();
  // stringstream data0_ss; data0_ss << "Offline PMTs: " << data[0]; string data0_str = data0_ss.str(); const char *data0_s = data0_str.c_str();
  // stringstream data1_ss; data1_ss << "Calibrated: " << data[1]; string data1_str = data1_ss.str(); const char *data1_s = data1_str.c_str();
  // stringstream data2_ss; data2_ss << "UNCalibrated: " << data[2]; string data2_str = data2_ss.str(); const char *data2_s = data2_str.c_str();
  // stringstream data3_ss; data3_ss << "CANT be: " << data[3]; string data3_str = data3_ss.str(); const char *data3_s = data3_str.c_str();
  // fSummary->AddText(data0_s);
  // fSummary->AddText(data1_s);
  // fSummary->AddText(data2_s);
  // fSummary->AddText(data3_s);
  // fSummary->Draw();
  // //fSummaryC->SaveAs("cd_comp_fSummary.png");
  // fSummaryC->Close();

  //ofile->Close();

}
