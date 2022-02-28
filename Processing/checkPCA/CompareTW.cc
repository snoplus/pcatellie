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


int main(int argc, char* argv[]) {

  if (argc != 7){
    std::cout<<" Usage: CompareTW TW1.ratdb TW2.ratdb t_min t_max run1 run2"<<std::endl;
    exit(0);
  }

  string TWFilename1, TWFilename2;
  double mintime, maxtime;
  int runid1, runid2;
  TWFilename1 = argv[1];
  TWFilename2 = argv[2];
  mintime = atoi(argv[3]);
  maxtime = atoi(argv[4]);
  runid1 = atoi(argv[5]);
  runid2 = atoi(argv[6]);

  DB *ratdb = RAT::DB::Get();
  ratdb->LoadDefaults();

  std::vector<vector <float> > Time1;
  std::vector<vector <float> > Charge1;
  std::vector<vector <float> > TimeErr1;
  std::vector<vector <float> > TimeErr1mean;
  std::vector<int> Flag1;
  std::vector<float> Intercept1;
  std::vector<float> Gradient1;

  std::vector<vector <float> > Time2;
  std::vector<vector <float> > Charge2;
  std::vector<vector <float> > TimeErr2;
  std::vector<vector <float> > TimeErr2mean;
  std::vector<int> Flag2;
  std::vector<float> Intercept2;
  std::vector<float> Gradient2;

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

  PCABits PCABitPack;
  PCABitPack.BeginOfRun(run1);

  DBLinkPtr TWbank1 = DB::Get()->GetLink("PCA_TW");
  std::cout << " Analyzing PCA TW calibration: run" << TWbank1->GetIArray("run_range")[0] << std::endl;
  std::vector<int> TWstatus1 = TWbank1->GetIArray("PCATW_status");
  std::vector<int> TWpoints1 = TWbank1->GetIArray("twinter");
  std::vector<float> TWpointErrs1 = TWbank1->GetFArrayFromD("twinterrms");
  std::vector<float> TWpointErrs1mean = TWbank1->GetFArrayFromD("meanerr");
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
  std::vector<float> TWpointErrs2mean = TWbank2->GetFArrayFromD("meanerr");
  std::cout<<" Done "<<std::endl;

  // Make a few comparison plots
  TImage *img = TImage::Create();
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  char imagename[30];

  // First unpack the ocnstants
  // unpack the data and store in new arrays,this is a (q,t) interpolation point
  double thisq, thist, thisr, thisi;
  Time1.resize(9728);
  TimeErr1.resize(9728);
  TimeErr1mean.resize(9728);
  Charge1.resize(9728);
  Intercept1.resize(9728);
  Gradient1.resize(9728);
  Time2.resize(9728);
  TimeErr2.resize(9728);
  TimeErr2mean.resize(9728);
  Charge2.resize(9728);
  Intercept2.resize(9728);
  Gradient2.resize(9728);
  int loc_word_points;
  int loc_word_fit;
  int packresult;

  for (int i = 0; i < 9728; i++) {
    for (int j = 0; j < 10; j++) {
      thisq = 0.0;
      thist = 0.0;
      loc_word_points = i * 11 + j;  // the location of the first inter point
      packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.QTINTERPOLATION, TWpoints1[loc_word_points], thisq, thist);
      if (thisq != 0) {
        Time1[i].push_back(thist);
        Charge1[i].push_back(thisq);
      } else {
        Time1[i].push_back(-9999);
        Charge1[i].push_back(0);
      }
      packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.QTINTERPOLATION, TWpoints2[loc_word_points], thisq, thist);
      if (thisq != 0) {
        Time2[i].push_back(thist);
        Charge2[i].push_back(thisq);
      } else {
        Time2[i].push_back(-9999);
        Charge2[i].push_back(0);
      }
    }

  loc_word_fit = i * 11 + 10;
  packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.GRADIENT_INTERCEPT, TWpoints1[loc_word_fit], thisi, thisr);
  Intercept1[i] = thisi;
  Gradient1[i] = thisr;
  packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.GRADIENT_INTERCEPT, TWpoints2[loc_word_fit], thisi, thisr);
  Intercept2[i] = thisi;
  Gradient2[i] = thisr;
}
for (int i = 0; i < 9728; i++) {
  for (int j = 0; j < 10; j++) {
    loc_word_points = i * 10 + j;
    TimeErr1[i].push_back(TWpointErrs1[loc_word_points]);
    TimeErr1mean[i].push_back(TWpointErrs1mean[loc_word_points]);
    TimeErr2[i].push_back(TWpointErrs2[loc_word_points]);
    TimeErr2mean[i].push_back(TWpointErrs2mean[loc_word_points]);
   }
}
//TFile *ofile = new TFile(Form("PCA_TW_comparison_%i_%i.root",runid1,runid2),"RECREATE");
//ofile->cd();
TH2D* check = new TH2D("", "", 100, -50, 50, 100, -5, 5);
TH1D* tw = new TH1D("", "", 100, -5, 5);
int nchannels = 9728;
for (int i = 0; i < nchannels; i++) {
  //cout << i << " " << TWstatus1[i] << " " << TWstatus2[i] << endl;
  // std::string binary = std::bitset<32>(TWstatus1[i]).to_string();
  // std::string binary2 = std::bitset<32>(TWstatus2[i]).to_string();
  // int isOff = atoi(&binary[30]);
  // int isOff2 = atoi(&binary2[30]);
  // if ( isOff == 1 & isOff2 == 1 ){ continue; }

  double chargearray1[10];
  double timearray1[10];
  double chargearrayerr1[10];
  double timearrayerr1[10];
  double timearrayerr1mean[10];
  double chargearray2[10];
  double timearray2[10];
  double chargearrayerr2[10];
  double timearrayerr2mean[10];
  double timearrayerr2[10];

  for (int h = 0; h < 10; h++) {
    chargearray1[h] = Charge1[i][h];
    timearray1[h] = Time1[i][h];
    timearrayerr1[h] = TimeErr1[i][h];
    timearrayerr1mean[h] = TimeErr1mean[i][h];
    chargearrayerr1[h] = 0.0;
    chargearray2[h] = Charge2[i][h];
    timearray2[h] = Time2[i][h];
    timearrayerr2[h] = TimeErr2[i][h];
    timearrayerr2mean[h] = TimeErr2mean[i][h];
    chargearrayerr2[h] = 0.0;
    if (h < 5) {
    check->Fill(Charge1[i][h]-Charge2[i][h], Time1[i][h]-Time2[i][h]);
    tw->Fill(Time1[i][h]-Time2[i][h]);
    }
  }
  if (i%100==0) {
    std::cout << i <<"/"<<nchannels<<" channels"<< std::endl;
  }
  // Time for some plots
  TF1 *FitLine1 = new TF1("f", "[0]+[1]*x", 0, 1000);
  TGraphErrors *interpolation1 = new TGraphErrors(10, chargearray1, timearray1, chargearrayerr1, timearrayerr1);
  TGraphErrors *interpolation1mean = new TGraphErrors(10, chargearray1, timearray1, chargearrayerr1, timearrayerr1mean);
  FitLine1->SetParameter(0, Intercept1[i]);
  FitLine1->SetParameter(1, Gradient1[i]);
  TF1 *FitLine2 = new TF1("f", "[0]+[1]*x", 0, 1000);
  TGraphErrors *interpolation2 = new TGraphErrors(10, chargearray2, timearray2, chargearrayerr2, timearrayerr2);
  TGraphErrors *interpolation2mean = new TGraphErrors(10, chargearray2, timearray2, chargearrayerr2, timearrayerr2mean);
  FitLine2->SetParameter(0, Intercept2[i]);
  FitLine2->SetParameter(1, Gradient2[i]);

  if ( interpolation1mean->GetMean(2) == -9999 || interpolation2mean->GetMean(2) == -9999 ){ continue; }

  interpolation1mean->SetMaximum(15.0);
  interpolation1mean->SetMinimum(-15.0);
  interpolation1mean->SetMarkerColor(4);
  interpolation1mean->SetLineColor(4);
  interpolation1mean->SetLineWidth(2);
  interpolation1mean->SetMarkerStyle(5);
  interpolation1mean->SetMarkerSize(2.0);
  interpolation1mean->GetXaxis()->SetTitle("QHS [cap]");
  interpolation1mean->GetYaxis()->SetTitle("Time [ns]");
  interpolation1mean->SetTitle( Form("PMT-%i", i) );
  interpolation1mean->Draw("A*");
  //interpolation1->SetMarkerStyle(5);
  //interpolation1->Draw("SAME*");
  FitLine1->SetLineWidth(1);
  FitLine1->SetLineColor(4);
  FitLine1->SetLineStyle(2);
  FitLine1->Draw("SAME");

  interpolation2mean->SetMarkerColor(2);
  interpolation2mean->SetLineColor(2);
  interpolation2mean->SetLineWidth(2);
  interpolation2mean->SetMarkerStyle(5);
  interpolation2mean->SetMarkerSize(2.0);
  interpolation2mean->Draw("SAME*");
  //interpolation2->SetMarkerStyle(5);
  //interpolation2->Draw("SAME*");
  FitLine2->SetLineWidth(1);
  FitLine2->SetLineColor(2);
  FitLine2->SetLineStyle(2);
  FitLine2->Draw("SAME");

  TLegend* leg = new TLegend(.7,.7,.9,.9);
  leg->AddEntry(interpolation1mean,Form("Run %i",runid1),"L");
  leg->AddEntry(interpolation2mean,Form("Run %i",runid2),"L");
  leg->SetFillColor(0);
  leg->Draw("same");

  img->Clear();
  img->FromPad(c1);
  //c1->Write( Form("PMT-%i", i) );
  c1->SaveAs( Form("TW_PMT-%i.png", i) );
  delete leg;
  delete interpolation1;
  delete interpolation1mean;
  delete FitLine1;
  delete interpolation2;
  delete interpolation2mean;
  delete FitLine2;
}

//ofile->Close();

}
