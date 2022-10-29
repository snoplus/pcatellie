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

  std::vector<vector <float> > Time1;
  std::vector<vector <float> > Charge1;
  std::vector<vector <float> > TimeErr1;
  std::vector<vector <float> > TimeErr1mean;
  std::vector<int> Flag1;

  std::vector<vector <float> > Time2;
  std::vector<vector <float> > Charge2;
  std::vector<vector <float> > TimeErr2;
  std::vector<vector <float> > TimeErr2mean;
  std::vector<int> Flag2;

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
  //const RAT::DU::PMTInfo& pmtinfo = RAT::DU::Utility::Get()->GetPMTInfo();
  //const RAT::DU::ChanHWStatus& chs = RAT::DU::Utility::Get()->GetChanHWStatus();

  PCABits PCABitPack;
  PCABitPack.BeginOfRun(run1);

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

  // First unpack the ocnstants
  // unpack the data and store in new arrays,this is a (q,t) interpolation point
  double thisq, thist, thisr, thisi;
  double thisq_old, thist_old, thisr_old, thisi_old;
  Time1.resize(9728);
  TimeErr1.resize(9728);
  TimeErr1mean.resize(9728);
  Charge1.resize(9728);
  Time2.resize(9728);
  TimeErr2.resize(9728);
  TimeErr2mean.resize(9728);
  Charge2.resize(9728);
  int loc_word_points;
  int loc_word_fit;
  int packresult;
  std::vector<int> IDs1;
  std::vector<int> IDs2;
  std::vector<double> Intercept1;
  std::vector<double> Intercept2;
  std::vector<double> Gradient1;
  std::vector<double> Gradient2;
  int off1 = 0;
  int off2 = 0;

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
      packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.QTINTERPOLATION, TWpoints2[loc_word_points], thisq_old, thist_old);
      if (thisq != 0) {
        Time2[i].push_back(thist_old);
        Charge2[i].push_back(thisq_old);
      } else {
        Time2[i].push_back(-9999);
        Charge2[i].push_back(0);
      }
    }

    loc_word_fit = i * 11 + 10;
    packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.GRADIENT_INTERCEPT, TWpoints1[loc_word_fit], thisi, thisr);
    std::string binary = std::bitset<32>(TWstatus1[i]).to_string();
    int isOff = atoi(&binary[30]);
    if (isOff != 1){
      IDs1.push_back( i );
      Intercept1.push_back( thisi );
      Gradient1.push_back( thisr );
    } else { off1++; }
    packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.GRADIENT_INTERCEPT, TWpoints2[loc_word_fit], thisi_old, thisr_old);
    std::string binary2 = std::bitset<32>(TWstatus2[i]).to_string();
    int isOff_old = atoi(&binary2[30]);
    if (isOff_old != 1){
      IDs2.push_back( i );
      Intercept2.push_back( thisi_old );
      Gradient2.push_back( thisr_old );
    } else { off2++; }
  }

  for (int i = 0; i < 9728; i++) {
    for (int j = 0; j < 10; j++) {
      loc_word_points = i * 10 + j;
      TimeErr1[i].push_back(TWpointErrs1[loc_word_points]);
      //TimeErr1mean[i].push_back(TWpointErrs1mean[loc_word_points]);
      TimeErr2[i].push_back(TWpointErrs2[loc_word_points]);
      //TimeErr2mean[i].push_back(TWpointErrs2mean[loc_word_points]);
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  double min_grad;
  double max_grad;
  double min_inter;
  double max_inter;
  int goodCount = 0;
  int badCount = 0;

  std::vector<double> grad_diff;
  std::vector<double> inter_diff;
  for (size_t i=0; i<IDs1.size(); i++) {
    int flag = 0;
    for (size_t j=0; j<IDs2.size(); j++) {
      if (IDs1[i] == IDs2[j]){
        flag = 1;
        goodCount++;
        grad_diff.push_back( Gradient1[i] - Gradient2[j] );
        inter_diff.push_back( Intercept1[i] - Intercept2[j] );
      }
    }
    if (flag == 0){
      badCount++;
    }
  }

  min_grad = *min_element(grad_diff.begin(), grad_diff.end());
  max_grad = *max_element(grad_diff.begin(), grad_diff.end());
  min_inter = *min_element(inter_diff.begin(), inter_diff.end());
  max_inter = *max_element(inter_diff.begin(), inter_diff.end());

  TH1D* grad = new TH1D("", "", 50, min_grad, max_grad);
  TH1D* inter = new TH1D("", "", 50, min_inter, max_inter);

  for (size_t i=0; i<goodCount; i++) {
    grad->Fill(grad_diff[i]);
    inter->Fill(inter_diff[i]);
  }

  c1->cd();
  c1->SetLogy();
  grad->GetXaxis()->SetTitle("TW fit: gradient");
  grad->GetYaxis()->SetTitle("Count");
  grad->SetTitle( "TimeWalk Fit - Comparison" );
  grad->Draw();
  c1->SaveAs("tw_comp_grad.png");
  c1->Close();

  c2->cd();
  c2->SetLogy();
  inter->GetXaxis()->SetTitle("TW fit: intercept");
  inter->GetYaxis()->SetTitle("Count");
  inter->SetTitle( "TimeWalk Fit - Comparison" );
  inter->Draw();
  c2->SaveAs("tw_comp_inter.png");
  c2->Close();

  //ofile->Close();

  // Store to log file
  stringstream logFile_namess;
  string logFile_name;
  FILE *logFile;
  logFile_namess.str("");
  logFile_namess << runid1 << "-" << runid2 << "_tw_compare.log";
  logFile_name = logFile_namess.str();
  logFile = fopen(logFile_name.c_str(), "w");

  fprintf(logFile, "Run1: %i\n", runid1);
  fprintf(logFile, "Run2: %i\n", runid2);
  fprintf(logFile, "IDs1: %i\n", IDs1.size());
  fprintf(logFile, "IDs2: %i\n", IDs2.size());
  fprintf(logFile, "off1: %i\n", off1);
  fprintf(logFile, "off2: %i\n", off2);
  fprintf(logFile, "goodCount: %i\n", goodCount);
  fprintf(logFile, "badCount: %i\n", badCount);
  fprintf(logFile, "min_grad: %f\n", min_grad);
  fprintf(logFile, "max_grad: %f\n", max_grad);
  fprintf(logFile, "min_inter: %f\n", min_inter);
  fprintf(logFile, "max_inter: %f\n", max_inter);

}
