//C++ includes
#include "Riostream.h"
#include "sys/stat.h"
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>

//ROOT includes
#include "TGaxis.h"
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

//#include "PCABits.hh"

using namespace std;
using namespace RAT;
using namespace RAT::DS;

int main(int argc, char* argv[]) {

  string GFFilename, TWFilename, GFFilename_old, TWFilename_old;
  double mintime, maxtime;
  int runidthis = 0;
  int runidold = 0;
  double offset = 0.0;
  // Get the arguments
  if (!(argc == 11)) {   // We print argv[0] assuming it is the program name
    cout << "usage: " << argv[0]
	 << " <.root filename>" << " <.ratdb TW name>" << " <.ratdb GF name> "
	 << "Min time [ns] " << "Max time [ns] " << "run number "
	 << "time global offset" << " <.old ratdb TW name>" << " <.old ratdb GF name>" 
   << " old run number" <<"\n";
    return -1;
  } else {
    TWFilename = argv[2];
    GFFilename = argv[3];
    mintime = atoi(argv[4]);
    maxtime = atoi(argv[5]);
    runidthis = atoi(argv[6]);
    offset = atof(argv[7]);
    TWFilename_old = argv[8];
    GFFilename_old = argv[9];
    runidold = atoi(argv[10]);
  }

  DB *ratdb = RAT::DB::Get();
  ratdb->LoadDefaults();
  
  try{
    ratdb->Load(TWFilename_old,true);
    ratdb->Load(GFFilename_old,true);
    cout << "Loaded local tables!" << endl;
  }
  catch(...){
    std::cout<< " PROBLEM LOADING RATDB! COMPARISON WITH PREVIOUS RUN WON'T BE CORRECT!" << std::endl;
  }

  DS::Run run;
  run.SetRunID(runidold + 1);

  try{
    ratdb->BeginOfRun(run);
  }
  catch(...){
    std::cout<< " PROBLEM LOADING RATDB! RUN WON'T BE CORRECT!" << std::endl;
  }

  PCABits PCABitPack;
  PCABitPack.BeginOfRun(run);

  DBLinkPtr fPmtCalib = DB::Get()->GetLink("PMTCALIB");
  DBLinkPtr TWbank = DB::Get()->GetLink("PCA_TW");

  std::cout << " Comparing to PCA GF calibration: run" << fPmtCalib->GetIArray("run_range")[0] << std::endl;
  std::cout << " Comparing to PCA TW calibration: run" << TWbank->GetIArray("run_range")[0] << std::endl;

  std::vector<int> TWpoints_old = TWbank->GetIArray("twinter");
  std::vector<float> fQhshhp = fPmtCalib->GetFArrayFromD("QHS_hhp");
  std::vector<float> fQhspeak = fPmtCalib->GetFArrayFromD("QHS_peak");
  std::vector<float> fQhsthresh = fPmtCalib->GetFArrayFromD("QHS_threshold");
  std::vector<float> fQhlhhp = fPmtCalib->GetFArrayFromD("QHL_hhp");
  std::vector<float> fQhlpeak = fPmtCalib->GetFArrayFromD("QHL_peak");
  std::vector<float> fQhlthresh = fPmtCalib->GetFArrayFromD("QHL_threshold");
  float fMeanhhpnew;
  float fMinhhpnew;
  float fMaxhhpnew;
  std::vector<float> fQhshhpnew;
  std::vector<float> fQhspeaknew;
  std::vector<float> fQhsthreshnew;
  std::vector<float> fQhlhhpnew;
  std::vector<float> fQhlpeaknew;
  std::vector<float> fQhlthreshnew;
  std::vector<vector <float> > Time;
  std::vector<vector <float> > Charge;
  std::vector<vector <float> > TimeErr;
  std::vector<int> Flag;
  std::vector<float> Intercept;
  std::vector<float> Gradient;
  ratdb->EndOfRun(run);
  std::cout << " EndOfRun " << std::endl;

  // Look at the extracted calibration values,
  // that are in the file that is given as an argument
  ratdb->SetAirplaneModeStatus(true);
  ratdb->Clear();
  std::cout << " Clear " << std::endl;
  run.SetRunID(runidthis + 1); // new constants
  std::cout << " SetRunID " << std::endl;
  ratdb->Load(GFFilename,true);
  ratdb->Load(TWFilename,true);
  std::cout << " Loads " << std::endl;
  try{
    ratdb->BeginOfRun(run);
  }
  catch(...){
    std::cout<< " PROBLEM LOADING RATDB! COMPARISON WITH PREVIOUS RUN WON'T BE CORRECT!" << std::endl;
  }

  std::cout << " Analyzing PCA GF calibration: run" << fPmtCalib->GetIArray("run_range")[0] << std::endl;

  DBLinkPtr fPmtCalibnew = DB::Get()->GetLink("PMTCALIB");
  fMeanhhpnew = fPmtCalibnew->GetD("mean_hhp");
  fMinhhpnew = fPmtCalibnew->GetD("min_hhp");
  fMaxhhpnew = fPmtCalibnew->GetD("max_hhp");
  fQhshhpnew = fPmtCalibnew->GetFArrayFromD("QHS_hhp");
  fQhspeaknew = fPmtCalibnew->GetFArrayFromD("QHS_peak");
  fQhsthreshnew = fPmtCalibnew->GetFArrayFromD("QHS_threshold");
  fQhlhhpnew = fPmtCalibnew->GetFArrayFromD("QHL_hhp");
  fQhlpeaknew = fPmtCalibnew->GetFArrayFromD("QHL_peak");
  fQhlthreshnew = fPmtCalibnew->GetFArrayFromD("QHL_threshold");

  // get the list of keys from the rootfile
  TFile *f1;
  f1 = new TFile(argv[1]);
  TIter next(f1->GetListOfKeys());
  TKey *key;

  ratdb->Clear();
  ratdb->Load(TWFilename);

  std::cout << " Analyzing PCA TW calibration: run" << TWbank->GetIArray("run_range")[0] << std::endl;

  std::vector<int> TWpoints = TWbank->GetIArray("twinter");
  std::vector<float> TWpointErrs = TWbank->GetFArrayFromD("twinterrms");
  //std::vector<float> TWpointMErrs = TWbank->GetFArrayFromD("MeanErr");
  // std::vector<int> TWFlags=TWbank->GetIArray("PCATW_status");
  TH2D *Cable = new TH2D("cable", "", 32*16*19, 0, 32*16*19, 500, -250, 250);
  TH2D *CableDiff = new TH2D("cableDiff", "", 32*16*19, 0, 32*16*19, 200, -5, 5);
  TH2D *CableDiffCard = new TH2D("cableDiffCard", "", 16, 0, 16, 200, -5, 5);
  TH1D *CableDiffHisto = new TH1D("CableDiffHisto", "", 40, -4, 4);

  // unpack the data and store in new arrays,
  // this is a (q,t) interpolation point
  double thisq, thist, thisr, thisi, thisi_old, thisr_old;
  Time.resize(9728);
  TimeErr.resize(9728);
  Charge.resize(9728);
  Intercept.resize(9728);
  Gradient.resize(9728);
  int loc_word_points;
  int loc_word_fit;
  int packresult;
  TImage *img = TImage::Create();
  TH1F *intersec = new TH1F("intersec", "", 1000, -100, 600);
  char imagename[30];
  TView3D *view = new TView3D();
  TCanvas *c1 = new TCanvas("c1", "c1", 1024, 768);

  std::cout << " Getting PMT Coverage plot " << std::endl;

  // Fiber delays from Bootstrap
  TGraphErrors *ledOffset =  (TGraphErrors*)f1->Get("ledOffsets");
  if(ledOffset != NULL) {
    c1->cd();
    gStyle->SetLabelSize(0.02);
    ledOffset->GetXaxis()->SetLabelSize(0.03);
    ledOffset->GetYaxis()->SetLabelSize(0.03);
    ledOffset->SetTitle("");
    ledOffset->Draw("APE1");
    img->Clear();
    img->FromPad(c1);
    img->WriteImage("ledoffsets.png");
  }

  // Fiber delay diff Bootstrap-RATDB
  TGraphErrors *fiberDelayDiff =  (TGraphErrors*)f1->Get("fiberDelayDiff");
  if(fiberDelayDiff != NULL) {
    c1->cd();
    gStyle->SetLabelSize(0.02);
    fiberDelayDiff->GetXaxis()->SetLabelSize(0.03);
    fiberDelayDiff->GetYaxis()->SetLabelSize(0.03);
    fiberDelayDiff->SetTitle("");
    fiberDelayDiff->Draw("APE1");
    img->Clear();
    img->FromPad(c1);
    img->WriteImage("fiberDelayDiff.png");
  }

  // General PMT coverage plot
  TGraph2D *coverage =  (TGraph2D*)f1->Get("PMT Coverage");
  c1->cd();
  TGaxis::SetMaxDigits(4);
  c1->SetMargin(0.05,0.15,0.05,0.05);
  coverage->SetMarkerStyle(8);
  coverage->SetMarkerSize(0.5);
  gStyle->SetLabelSize(0.02);
  coverage->GetXaxis()->SetLabelSize(0.03);
  coverage->GetYaxis()->SetLabelSize(0.03);
  coverage->SetTitle("");
  //coverage->SetMaximum(80000);
  coverage->Draw("AZCOLPCOL");
  view->RotateView(270, 0);
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("HitCoverage.png");

  // General PMT coverage plot -- CUT
  TGraph2D *coverageCut =  (TGraph2D*)f1->Get("PMT Coverage Cut");
  c1->cd();
  TGaxis::SetMaxDigits(4);
  c1->SetMargin(0.05,0.15,0.05,0.05);
  coverageCut->SetMarkerStyle(8);
  coverageCut->SetMarkerSize(0.5);
  gStyle->SetLabelSize(0.02);
  coverageCut->GetXaxis()->SetLabelSize(0.03);
  coverageCut->GetYaxis()->SetLabelSize(0.03);
  coverageCut->SetTitle("");
  coverageCut->Draw("AZCOLPCOL");
  view->RotateView(270, 0);
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("HitCoverageCut.png");

  std::cout << " Getting LED Coverage plot " << std::endl;

  // General LED coverage plot
  TGraph2D *ledCoverage =  (TGraph2D*)f1->Get("ledCoverage");
  if(ledCoverage != NULL){
    c1->cd();
    TGaxis::SetMaxDigits(4);
    c1->SetMargin(0.05,0.15,0.05,0.05);
    ledCoverage->SetMarkerStyle(8);
    ledCoverage->SetMarkerSize(0.5);
    gStyle->SetLabelSize(0.02);
    ledCoverage->GetXaxis()->SetLabelSize(0.03);
    ledCoverage->GetYaxis()->SetLabelSize(0.03);
    ledCoverage->SetTitle("");
    ledCoverage->Draw("AZCOLPCOL");
    view->RotateView(270, 0);
    img->Clear();
    img->FromPad(c1);
    img->WriteImage("ledcoverage.png");
  }

  std::cout << " Getting Best Fibre plot " << std::endl;

  // General Best Fibre plot
  TGraph2D *BestFibre =  (TGraph2D*)f1->Get("Best Fibre");
  if(BestFibre != NULL){
    c1->cd();
    TGaxis::SetMaxDigits(4);
    c1->SetMargin(0.05,0.15,0.05,0.05);
    BestFibre->SetMarkerStyle(8);
    BestFibre->SetMarkerSize(0.5);
    gStyle->SetLabelSize(0.02);
    gStyle->SetNumberContours(95);
    BestFibre->GetXaxis()->SetLabelSize(0.03);
    BestFibre->GetYaxis()->SetLabelSize(0.03);
    BestFibre->SetTitle("");
    BestFibre->Draw("AZCOLPCOL");
    view->RotateView(270, 0);
    img->Clear();
    img->FromPad(c1);
    img->WriteImage("BestFibre.png");
  } else {std::cout << "did not load best fibre from root" << std::endl;}

  std::cout << " Getting ledOffsets plot " << std::endl;

  std::cout << " Getting Offline PMTs plot " << std::endl;

  // Plot offline channels
  TGraph2D *offline = (TGraph2D*)f1->Get("Offline PMTs");
  c1->cd();
  TGaxis::SetMaxDigits(4);
  c1->SetMargin(0.05,0.15,0.05,0.05);
  offline->SetMarkerStyle(8);
  offline->SetMarkerSize(0.5);
  gStyle->SetLabelSize(0.02);
  offline->GetXaxis()->SetLabelSize(0.03);
  offline->GetYaxis()->SetLabelSize(0.03);
  offline->SetTitle("");
  offline->Draw("AZCOLPCOL");
  view->RotateView(270, 0);
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("OfflinePMTs.png");

  std::cout << " Getting HHP flatmaps " << std::endl;

  // The HHP QHL flatmap
  TGraph2D *hhpqhl = (TGraph2D*)f1->Get("HHP QHL flatmap");
  c1->cd();
  TGaxis::SetMaxDigits(4);
  c1->SetMargin(0.05,0.15,0.05,0.05);
  hhpqhl->SetMarkerStyle(8);
  hhpqhl->SetMarkerSize(0.5);
  gStyle->SetLabelSize(0.02);
  hhpqhl->GetXaxis()->SetLabelSize(0.03);
  hhpqhl->GetYaxis()->SetLabelSize(0.03);
  hhpqhl->SetTitle("");
  hhpqhl->Draw("AZCOLPCOL");
  view->RotateView(270, 0);
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("HHPQHLflat.png");

  // The HHP QHS flatmap
  TGraph2D *hhpqhs = (TGraph2D*)f1->Get("HHP QHS flatmap");
  c1->cd();
  TGaxis::SetMaxDigits(4);
  c1->SetMargin(0.05,0.15,0.05,0.05);
  hhpqhs->SetMarkerStyle(8);
  hhpqhs->SetMarkerSize(0.5);
  gStyle->SetLabelSize(0.02);
  hhpqhs->GetXaxis()->SetLabelSize(0.03);
  hhpqhs->GetYaxis()->SetLabelSize(0.03);
  hhpqhs->SetTitle("");
  hhpqhs->Draw("AZCOLPCOL");
  view->RotateView(270,0);
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("HHPQHSflat.png");

  std::cout << " Getting PMT plots. " << std::endl;

  // Time for the PMT plots!
  for (int i = 0; i < 9728; i++) {
    int card = ( i%(16*32) )/32 ;

    for (int j = 0; j < 10; j++) {
      thisq = 0.0;
      thist = 0.0;
      loc_word_points = i * 11 + j;
      packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.QTINTERPOLATION, TWpoints[loc_word_points], thisq, thist);
      if (thisq !=0) {
        Time[i].push_back(thist + offset);
        Charge[i].push_back(thisq);
      } else {
        Time[i].push_back(-9999);
        Charge[i].push_back(0);
      }
    }
    loc_word_fit = i * 11 + 10;
    packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.GRADIENT_INTERCEPT, TWpoints[loc_word_fit], thisi, thisr);
    if(thisi == -9999.) continue;
    Intercept[i] = thisi + offset;
    Gradient[i] = thisr;
    Cable->Fill(i, thisi);
    double iterat = thisr * 100 + thisi;
    packresult = PCABitPack.PCAPacker(PCABitPack.UNPACK, PCABitPack.GRADIENT_INTERCEPT, TWpoints_old[loc_word_fit], thisi_old, thisr_old);
    if(thisi_old == -9999.) continue;
    //std::cout << i << ", ";
    //std::cout << thisi << ", ";
    //std::cout << thisi_old << ", ";
    CableDiff->Fill(i, thisi - thisi_old);
    CableDiffCard->Fill(card, thisi - thisi_old);
    CableDiffHisto->Fill(thisi - thisi_old);
    intersec->Fill(iterat);
  }

  // Get the errors as well
  for (int i = 0; i < 9728; i++) {
    for (int j = 0; j < 10; j++) {
      loc_word_points = i * 10 + j;
      TimeErr[i].push_back(TWpointErrs[loc_word_points]);
    }
  }

  std::cout << " DONE " << std::endl;

  c1->cd();
  c1->SetMargin(0.1,0.1,0.1,0.1);
  c1->Clear();
  intersec->GetXaxis()->SetTitle("Delay [ns]");
  intersec->GetXaxis()->SetTitleOffset(1.0);
  intersec->GetYaxis()->SetTitleOffset(1.0);
  intersec->GetYaxis()->SetTitle("Counts");
  //intersec->SaveAs("intersec.root");
  intersec->Draw();
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("intersec.png");
  c1->cd();
  c1->Clear();

  Cable->SetStats(0);
  Cable->GetXaxis()->SetTitle("LCN");
  Cable->GetXaxis()->SetTitleOffset(1.0);
  Cable->GetYaxis()->SetTitleOffset(1.0);
  Cable->GetYaxis()->SetTitle("Delay [ns]");
  //Cable->SaveAs("cable.root");
  Cable->GetYaxis()->SetRangeUser(-20, 20);
  Cable->Draw();
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("cable.png");
  c1->cd();
  c1->Clear();

  CableDiff->SetStats(0);
  CableDiff->GetXaxis()->SetTitle("LCN");
  CableDiff->GetXaxis()->SetTitleOffset(1.0);
  CableDiff->GetYaxis()->SetTitleOffset(1.0);
  CableDiff->GetYaxis()->SetTitle("New - old delay [ns]");
  //CableDiff->SaveAs("cableDiff.root");
  CableDiff->Draw("");
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("cableDiff.png");
  c1->cd();
  c1->Clear();

  CableDiffCard->SetStats(0);
  CableDiffCard->GetXaxis()->SetTitle("Card");
  CableDiffCard->GetXaxis()->SetTitleOffset(1.0);
  CableDiffCard->GetYaxis()->SetTitleOffset(1.0);
  CableDiffCard->GetYaxis()->SetTitle("New - old delay [ns]");
  //CableDiffCard->SaveAs("cableDiffCard.root");
  CableDiffCard->Draw("");
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("cableDiffCard.png");
  c1->cd();
  c1->Clear();

  //CableDiffHisto->SetStats(0);
  CableDiffHisto->GetXaxis()->SetTitle("Channel delay - Difference [ns]");
  CableDiffHisto->GetXaxis()->SetTitleOffset(1.0);
  CableDiffHisto->GetYaxis()->SetTitleOffset(1.0);
  CableDiffHisto->GetYaxis()->SetTitle("Count");
  //CableDiffHisto->SaveAs("cableDiffHisto.root");
  CableDiffHisto->Draw("");
  img->Clear();
  img->FromPad(c1);
  img->WriteImage("cableDiffHisto.png");
  c1->cd();
  c1->Clear();


  std::cout << " Getting Diff plots. " << std::endl;
  TCanvas *c2 = new TCanvas("c2", "c2", 1024, 768);

  //Now get the general plots
  img->Clear();
  gPad->SetLogy(0);
  TH2D *PCAdiffPEAKQHS = (TH2D*)f1->Get("fPCAdiffPEAKQHS");
  PCAdiffPEAKQHS->SetTitle("");
  PCAdiffPEAKQHS->GetXaxis()->SetTitle("PMT ID");
  PCAdiffPEAKQHS->GetXaxis()->SetTitleColor(1);
  PCAdiffPEAKQHS->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffPEAKQHS.png");
  //c2->SaveAs("PCAdiffPEAKQHS.root");

  img->Clear();
  TH2D *PCAdiffTHRESHQHS = (TH2D*)f1->Get("fPCAdiffTHRESHQHS");
  PCAdiffTHRESHQHS->SetTitle("");
  PCAdiffTHRESHQHS->GetXaxis()->SetTitle("PMT ID");
  PCAdiffTHRESHQHS->GetXaxis()->SetTitleColor(1);
  PCAdiffTHRESHQHS->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffTHRESHQHS.png");
  //c2->SaveAs("PCAdiffTHRESHQHS.root");

  img->Clear();
  TH2D *PCAdiffHHPQHS = (TH2D*)f1->Get("fPCAdiffHHPQHS");
  PCAdiffHHPQHS->SetTitle("");
  PCAdiffHHPQHS->GetXaxis()->SetTitle("PMT ID");
  PCAdiffHHPQHS->GetXaxis()->SetTitleColor(1);

  PCAdiffHHPQHS->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffHHPQHS.png");
  //c2->SaveAs("PCAdiffHHPQHS.root");

  img->Clear();
  TH2D *PCAdiffPEAKQHL = (TH2D*)f1->Get("fPCAdiffPEAKQHL");
  PCAdiffPEAKQHL->SetTitle("");
  PCAdiffPEAKQHL->GetXaxis()->SetTitle("PMT ID");
  PCAdiffPEAKQHL->GetXaxis()->SetTitleColor(1);
  PCAdiffPEAKQHL->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffPEAKQHL.png");
  //c2->SaveAs("PCAdiffPEAKQHL.root");

  img->Clear();
  TH2D *PCAdiffTHRESHQHL = (TH2D*)f1->Get("fPCAdiffTHRESHQHL");
  PCAdiffTHRESHQHL->SetTitle("");
  PCAdiffTHRESHQHL->GetXaxis()->SetTitle("PMT ID");
  PCAdiffTHRESHQHL->GetXaxis()->SetTitleColor(1);
  PCAdiffTHRESHQHL->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffTHRESHQHL.png");
  //c2->SaveAs("PCAdiffTHRESHQHL.root");

  img->Clear();
  TH2D *PCAdiffHHPQHL = (TH2D*)f1->Get("fPCAdiffHHPQHL");
  PCAdiffHHPQHL->SetTitle("");
  PCAdiffHHPQHL->GetXaxis()->SetTitle("PMT ID");
  PCAdiffHHPQHL->GetXaxis()->SetTitleColor(1);
  PCAdiffHHPQHL->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffHHPQHL.png");
  //c2->SaveAs("PCAdiffHHPQHL.root");

  gPad->SetLogy();
  img->Clear();
  TH1F *PCAdiffPEAKQHS1d = (TH1F*)f1->Get("fPCAdiffPEAKQHS1d");
  PCAdiffPEAKQHS1d->SetTitle("");
  PCAdiffPEAKQHS1d->GetXaxis()->SetTitle("Difference in QHS peak [cap]");
  PCAdiffPEAKQHS1d->GetXaxis()->SetTitle("QHS (peak_{prev}-peak_{new}) [cap]");
  PCAdiffPEAKQHS1d->GetXaxis()->SetTitleColor(1);
  PCAdiffPEAKQHS1d->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffPEAKQHS1d.png");
  //c2->SaveAs("PCAdiffPEAKQHS1d.root");

  img->Clear();
  TH1F *PCAdiffTHRESHQHS1d = (TH1F*)f1->Get("fPCAdiffTHRESHQHS1d");
  PCAdiffTHRESHQHS1d->SetTitle("");
  PCAdiffTHRESHQHS1d->GetXaxis()->SetTitle("QHS (thesh_{prev}-thresh_{new}) [cap]");
  PCAdiffTHRESHQHS1d->GetXaxis()->SetTitleColor(1);
  PCAdiffTHRESHQHS1d->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffTHRESHQHS1d.png");
  //c2->SaveAs("PCAdiffTHRESHQHS1d.root");


  img->Clear();
  TH1F *PCAdiffHHPQHS1d = (TH1F*)f1->Get("fPCAdiffHHPQHS1d");
  PCAdiffHHPQHS1d->SetTitle("");
  PCAdiffHHPQHS1d->GetXaxis()->SetTitle("QHS (hhp_{prev}-hhp_{new}) [cap]");
  PCAdiffHHPQHS1d->GetXaxis()->SetTitleColor(1);
  PCAdiffHHPQHS1d->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffHHPQHS1d.png");
  //c2->SaveAs("PCAdiffHHPQHS1d.root");

  img->Clear();
  TH1F *PCAdiffPEAKQHL1d = (TH1F*)f1->Get("fPCAdiffPEAKQHL1d");
  PCAdiffPEAKQHL1d->SetTitle("");
  PCAdiffPEAKQHL1d->GetXaxis()->SetTitle("QHL (peak_{prev}-peak_{new}) [cap]");
  PCAdiffPEAKQHL1d->GetXaxis()->SetTitleColor(1);
  PCAdiffPEAKQHL1d->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffPEAKQHL1d.png");
  //c2->SaveAs("PCAdiffPEAKQHL1d.root");


  img->Clear();
  TH1F *PCAdiffTHRESHQHL1d = (TH1F*)f1->Get("fPCAdiffTHRESHQHL1d");
  PCAdiffTHRESHQHL1d->SetTitle("");
  PCAdiffTHRESHQHL1d->GetXaxis()->SetTitle("QHL (thresh_{prev}-thresh_{new}) [cap]");
  PCAdiffTHRESHQHL1d->GetXaxis()->SetTitleColor(1);
  PCAdiffTHRESHQHL1d->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffTHRESHQHL1d.png");
  //c2->SaveAs("PCAdiffTHRESHQHL1d.root");

  img->Clear();
  TH1F *PCAdiffHHPQHL1d = (TH1F*)f1->Get("fPCAdiffHHPQHL1d");
  PCAdiffHHPQHL1d->SetTitle("");
  PCAdiffHHPQHL1d->GetXaxis()->SetTitle("QHL (hhp_{prev}-hhp_{new}) [cap]");
  PCAdiffHHPQHL1d->GetXaxis()->SetTitleColor(1);
  PCAdiffHHPQHL1d->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAdiffHHPQHL1d.png");
  //c2->SaveAs("PCAdiffHHPQHL1d.root");

  img->Clear();
  TH1I *PCAhitPerPMT = (TH1I*)f1->Get("fPCAhitPerPMT");
  PCAhitPerPMT->SetTitle("");
  PCAhitPerPMT->GetXaxis()->SetTitle("#hits/PMT");
  PCAhitPerPMT->GetXaxis()->SetTitleColor(1);
  PCAhitPerPMT->GetXaxis()->SetTitleOffset(1.2);
  PCAhitPerPMT->GetYaxis()->SetTitle("Counts");
  PCAhitPerPMT->Rebin(100);
  PCAhitPerPMT->Draw();
  img->FromPad(c2);
  img->WriteImage("PCAhitPerPMT.png");
  PCAhitPerPMT->SetTitle("");
  //c2->SaveAs("PCAhitPerPMT.root");

  img->Clear();
  TH1I *PCATWRMS = (TH1I*)f1->Get("fRMSHist");
  PCATWRMS->SetTitle("");
  PCATWRMS->GetXaxis()->SetTitle("TW point RMS");
  PCATWRMS->GetXaxis()->SetTitleColor(1);
  PCATWRMS->GetXaxis()->SetTitleOffset(1.2);
  PCATWRMS->GetYaxis()->SetTitle("Counts");
  PCATWRMS->Draw();
  img->FromPad(c2);
  img->WriteImage("PCATWRMS.png");
  PCATWRMS->SetTitle("");
  //c2->SaveAs("PCATWRMS.root");

  std::cout << " DONE " << std::endl;

  std::cout << " Get general plots. " << std::endl;

  TH1F * QHSpeakhistsnew= new TH1F("qhspeakhistsnew","qhspeakhistsnew",220,-20,200);
  TH1F * QHSthreshhistsnew= new TH1F("qhsthreshhistsnew","qhstheshhistsnew",220,-20,200);
  TH1F * QHShhphistsnew= new TH1F("qhshhphistsnew","qhshhphistsnew",220,-20,200);
  TH1F * QHSpeakhistsold= new TH1F("qhspeakhistsold","qhspeakhistsold",220,-20,200);
  TH1F * QHSthreshhistsold= new TH1F("qhsthreshhistsold","qhstheshhistsold",220,-20,200);
  TH1F * QHShhphistsold= new TH1F("qhshhphistsold","qhshhphistsold",220,-20,200);

  TH1F * QHLpeakhistsnew= new TH1F("qhlpeakhistsnew","qhlpeakhistsnew",220,-20,200);
  TH1F * QHLthreshhistsnew= new TH1F("qhlthreshhistsnew","qhltheshhistsnew",220,-20,200);
  TH1F * QHLhhphistsnew= new TH1F("qhlhhphistsnew","qhlhhphistsnew",220,-20,200);
  TH1F * QHLpeakhistsold= new TH1F("qhlpeakhistsold","qhlpeakhistsold",220,-20,200);
  TH1F * QHLthreshhistsold= new TH1F("qhlthreshhistsold","qhltheshhistsold",220,-20,200);
  TH1F * QHLhhphistsold= new TH1F("qhlhhphistsold","qhlhhphistsold",220,-20,200);

  TH2F * peakQHS=new TH2F("peakQHS","peakQHS",60,-10,50,60,-10,50);
  TH2F * threshQHS=new TH2F("threshQHS","threshQHS",60,-10,50,60,-10,50);
  TH2F * hhpQHS=new TH2F("hhpQHS","hhpQHS",110,-10,100,110,-10,100);
  TH2F * peakQHL=new TH2F("peakQHL","peakQHL",60,-10,50,60,-10,50);
  TH2F * threshQHL=new TH2F("threshQHL","threshQHL",60,-10,50,60,-10,50);
  TH2F * hhpQHL=new TH2F("hhpQHL","hhpQHL",110,-10,100,110,-10,100);

  //Fill some of the comparison summary plots
  for(int i=0;i<fQhshhp.size();i++){
    QHSpeakhistsold->Fill(fQhspeak[i]);
    QHSthreshhistsold->Fill(fQhsthresh[i]);
    QHShhphistsold->Fill(fQhshhp[i]);
    QHSpeakhistsnew->Fill(fQhspeaknew[i]);
    QHSthreshhistsnew->Fill(fQhsthreshnew[i]);
    QHShhphistsnew->Fill(fQhshhpnew[i]);
    peakQHS->Fill(fQhspeaknew[i], fQhspeak[i]);
    threshQHS->Fill(fQhsthreshnew[i], fQhsthresh[i]);
    hhpQHS->Fill(fQhshhpnew[i], fQhshhp[i]);
    QHLpeakhistsold->Fill(fQhlpeak[i]);
    QHLthreshhistsold->Fill(fQhlthresh[i]);
    QHLhhphistsold->Fill(fQhlhhp[i]);
    QHLpeakhistsnew->Fill(fQhlpeaknew[i]);
    QHLthreshhistsnew->Fill(fQhlthreshnew[i]);
    QHLhhphistsnew->Fill(fQhlhhpnew[i]);
    peakQHL->Fill(fQhlpeaknew[i], fQhlpeak[i]);
    threshQHL->Fill(fQhlthreshnew[i], fQhlthresh[i]);
    hhpQHL->Fill(fQhlhhpnew[i], fQhlhhp[i]);
  }

  TLine *guideeyes=new TLine(-10,-10,50,50);
  TLine *guideeyes2=new TLine(-10,-10,100,100);
  guideeyes->SetLineColor(2);
  guideeyes2->SetLineColor(2);

  c1->Divide(1,3);
  c1->cd(1);
  threshQHS->GetXaxis()->SetTitle("thesh new [cap]");
  threshQHS->GetYaxis()->SetTitle("thesh prev [cap]");
  threshQHS->Draw();
  guideeyes->Draw("SAME");
  c1->cd(2);
  peakQHS->GetXaxis()->SetTitle("peak new [cap]");
  peakQHS->GetYaxis()->SetTitle("peak prev [cap]");
  peakQHS->Draw();
  guideeyes->Draw("SAME");
  c1->cd(3);
  hhpQHS->GetXaxis()->SetTitle("hhp new [cap]");
  hhpQHS->GetYaxis()->SetTitle("hhp prev [cap]");
  hhpQHS->Draw();
  guideeyes2->Draw("SAME");
  img->FromPad(c1);
  img->WriteImage("peakSummaryQHS.png");

  c1->cd(1);

  threshQHL->GetXaxis()->SetTitle("thesh new [cap]");
  threshQHL->GetYaxis()->SetTitle("thesh prev [cap]");
  threshQHL->Draw();
  guideeyes->Draw("SAME");
  c1->cd(2);
  peakQHL->GetXaxis()->SetTitle("peak new [cap]");
  peakQHL->GetYaxis()->SetTitle("peak prev [cap]");
  peakQHL->Draw();
  guideeyes->Draw("SAME");
  c1->cd(3);
  hhpQHL->GetXaxis()->SetTitle("hhp new [cap]");
  hhpQHL->GetYaxis()->SetTitle("hhp prev [cap]");
  hhpQHL->Draw();
  guideeyes2->Draw("SAME");
  img->FromPad(c1);
  img->WriteImage("peakSummaryQHL.png");

  TLegend *leg1 = new TLegend(0.35, 0.75, 0.80, 0.90);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);

  c1->cd(1);
  gPad->SetLogy();
  QHSthreshhistsold->GetXaxis()->SetTitle("QHS Threshold [cap]");
  QHSthreshhistsold->GetYaxis()->SetTitle("Counts");
  QHSthreshhistsold->SetTitle("");
  QHSthreshhistsold->Draw();
  QHSthreshhistsnew->SetLineColor(2);
  QHSthreshhistsnew->Draw("SAME");
  leg1->AddEntry(QHSthreshhistsnew->DrawCopy("same"), "new");
  leg1->AddEntry(QHSthreshhistsold->DrawCopy("same"), "prev");
  leg1->Draw();
  c1->cd(2);
  gPad->SetLogy();
  QHSpeakhistsold->GetXaxis()->SetTitle("QHS Peak [cap]");
  QHSpeakhistsold->GetYaxis()->SetTitle("Counts");
  QHSpeakhistsold->SetTitle("");
  QHSpeakhistsold->Draw();
  QHSpeakhistsnew->SetLineColor(2);
  QHSpeakhistsnew->Draw("SAME");
  leg1->Clear();
  leg1->AddEntry(QHSpeakhistsnew->DrawCopy("same"), "new");
  leg1->AddEntry(QHSpeakhistsold->DrawCopy("same"), "prev");
  leg1->Draw();
  c1->cd(3);
  gPad->SetLogy();
  QHShhphistsold->GetXaxis()->SetTitle("QHS hhp [cap]");
  QHShhphistsold->GetYaxis()->SetTitle("Counts");
  QHShhphistsold->SetTitle("");
  QHShhphistsold->Draw();
  QHShhphistsnew->SetLineColor(2);
  QHShhphistsnew->Draw("SAME");
  leg1->Clear();
  leg1->AddEntry(QHShhphistsnew->DrawCopy("same"), "new");
  leg1->AddEntry(QHShhphistsold->DrawCopy("same"), "prev");
  leg1->Draw();
  img->FromPad(c1);
  img->WriteImage("HistsQHS.png");
  //c1->SaveAs("HistsQHS.root");

  c1->cd(1);
  gPad->SetLogy();
  QHLthreshhistsold->GetXaxis()->SetTitle("QHL Threshold [cap]");
  QHLthreshhistsold->GetYaxis()->SetTitle("Counts");
  QHLthreshhistsold->SetTitle("");
  QHLthreshhistsold->Draw();
  QHLthreshhistsnew->SetLineColor(2);
  QHLthreshhistsnew->Draw("SAME");
  leg1->AddEntry(QHLthreshhistsnew->DrawCopy("same"), "new");
  leg1->AddEntry(QHLthreshhistsold->DrawCopy("same"), "prev");
  leg1->Draw();
  c1->cd(2);
  gPad->SetLogy();
  QHLpeakhistsold->GetXaxis()->SetTitle("QHL Peak [cap]");
  QHLpeakhistsold->GetYaxis()->SetTitle("Counts");
  QHLpeakhistsold->SetTitle("");
  QHLpeakhistsold->Draw();
  QHLpeakhistsnew->SetLineColor(2);
  QHLpeakhistsnew->Draw("SAME");
  leg1->Clear();
  leg1->AddEntry(QHLpeakhistsnew->DrawCopy("same"), "new");
  leg1->AddEntry(QHLpeakhistsold->DrawCopy("same"), "prev");
  leg1->Draw();
  c1->cd(3);
  gPad->SetLogy();
  QHLhhphistsold->GetXaxis()->SetTitle("QHL hhp [cap]");
  QHLhhphistsold->GetYaxis()->SetTitle("Counts");
  QHLhhphistsold->SetTitle("");
  QHLhhphistsold->Draw();
  QHLhhphistsnew->SetLineColor(2);
  QHLhhphistsnew->Draw("SAME");
  leg1->Clear();
  leg1->AddEntry(QHLhhphistsnew->DrawCopy("same"), "new");
  leg1->AddEntry(QHLhhphistsold->DrawCopy("same"), "prev");
  leg1->Draw();
  img->FromPad(c1);
  img->WriteImage("HistsQHL.png");
  //c1->SaveAs("HistsQHL.root");

  std::cout << " DONE. " << std::endl;



  std::cout << " Getting histograms in ROOT file. This might take a while... " << std::endl;

  size_t pos1, pos2;
  TLine *lhhpold;
  TLine *lthreshold;
  TLine *lpeakold;
  TLine *lhhpnew;
  TLine *lthreshnew;
  TLine *lpeaknew;
  TLegend *leg4 = new TLegend(0.7, 0.70, 0.89, 0.90);
  std::string chrgqhs("QHS");
  std::string chrgqhl("QHL");
  double max = 0;
  double chargearray[10];
  double timearray[10];
  double chargearrayerr[10];
  double timearrayerr[10];
  std::string str;
  int keycounter = 0;

  // Loop over all histograms in the rootfile.. this will take a while
  while ((key = (TKey*)next())) {
    keycounter++;
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1") && !cl->InheritsFrom("TH2")) continue;

    // First up are the 2D histograms: time vs. charge
    if (cl->InheritsFrom("TH2")) {
      TH2 *h2;
      h2 = (TH2*)key->ReadObj();
      str = h2->GetTitle();
      // Find PMTID for this histogram
      pos1 = str.find("-")+1;
      pos2 = str.find(" ")+1;
      if (pos2 > 9728) continue;
      std::string pmtidstr = str.substr(pos1, pos1);
      std::string chargestr = str.substr(pos2, 3);
      int PMTid = atoi(pmtidstr.c_str());
      if (!(chargestr.compare("TW"))) {
        for (int h = 0; h < 10; h++) {
          chargearray[h] = Charge[PMTid][h];
          timearray[h] = Time[PMTid][h];
          timearrayerr[h] = TimeErr[PMTid][h];
          chargearrayerr[h] = 0.0;
          if(PMTid == 175) std::cout<<h<<": "<<chargearray[h]<<" "<<timearray[h]<<std::endl;
        }
        c2->cd();
        gPad->SetLogy(0);
        h2->GetYaxis()->SetTitle("Tpca [ns]");
        h2->GetXaxis()->SetTitle("QHS [cap]");
        h2->GetYaxis()->SetRangeUser(mintime, maxtime);
        h2->GetXaxis()->SetRangeUser(-20, 400);
        h2->GetYaxis()->SetTitleOffset(1.3);
        h2->SetMarkerColor(15);  // Set points in scatter plot to grey
        h2->Draw();
        TF1 *FitLine = new TF1("f", "[0]+[1]*x", 0, 1000);
        TGraphErrors *interpolation= new TGraphErrors(10, chargearray, timearray, chargearrayerr, timearrayerr);
        FitLine->SetParameter(0, Intercept[PMTid]);
        FitLine->SetParameter(1, Gradient[PMTid]);
        interpolation->SetMaximum(50.0);
        interpolation->SetMarkerColor(4);
        interpolation->SetLineColor(4);
        interpolation->SetLineWidth(2);
        interpolation->SetMarkerStyle(20);
        interpolation->SetMarkerSize(1.0);
        interpolation->Draw("SAME*");
        FitLine->SetLineWidth(1);
        FitLine->SetLineColor(2);
        FitLine->SetLineStyle(2);
        FitLine->Draw("SAME");
        //stringstream kek;
        //kek << PMTid << "-" << chargestr.c_str() << ".root";
        //c2->SaveAs(kek.str().c_str());
        img->Clear();
        img->FromPad(c2);
        sprintf(imagename, "PMT-%i-%s.png", PMTid, chargestr.c_str());
        img->WriteImage(imagename);
      }
      delete h2;
    }

    // Now we have the 1D histograms
    if (cl->InheritsFrom("TH1")) {
      TH1 *h1;
      h1 = (TH1*)key->ReadObj();
      str = h1->GetTitle();
      // Find PMTID for this histogram
      pos1 = str.find("-")+1;
      pos2 = str.find(" ")+1;
      if (pos2 > 9728) continue;
      std::string pmtidstr = str.substr(pos1, pos1);
      std::string chargestr = str.substr(pos2, 3);
      int PMTid = atoi(pmtidstr.c_str());
      if (!(chargestr.compare("Tec"))) {
        c2->cd();
        gPad->SetLogy();
        h1->GetXaxis()->SetTitle("Teca [ns]");
        h1->GetYaxis()->SetTitle("Counts");
        h1->Draw();
        img->Clear();
        img->FromPad(c2);
        sprintf(imagename, "PMT-%i-%s.png", PMTid,chargestr.c_str());
        img->WriteImage(imagename);
      }
      if(!(chargestr.compare("QHL")) || !(chargestr.compare("QHS"))){
        gPad->SetLogy(0);
        if(!(chargestr.compare("QHS"))){
          h1->GetXaxis()->SetTitle("QHS [cap]");
          lhhpold = new TLine(fQhshhp[PMTid], 0, fQhshhp[PMTid], h1->GetMaximum()+5);
          lthreshold = new TLine(fQhsthresh[PMTid], 0, fQhsthresh[PMTid], h1->GetMaximum()+5);
          lpeakold = new TLine(fQhspeak[PMTid], 0,fQhspeak[PMTid], h1->GetMaximum()+5);
          lhhpnew = new TLine(fQhshhpnew[PMTid], 0,fQhshhpnew[PMTid], h1->GetMaximum()+5);
          lthreshnew = new TLine(fQhsthreshnew[PMTid], 0,fQhsthreshnew[PMTid],h1->GetMaximum()+5);
          lpeaknew = new TLine(fQhspeaknew[PMTid], 0,fQhspeaknew[PMTid],h1->GetMaximum()+5);
          max = fQhspeaknew[PMTid];
        }
        if(!(chargestr.compare("QHL"))){
          h1->GetXaxis()->SetTitle("QHL [cap]");
          lhhpold=new TLine(fQhlhhp[PMTid],0,fQhlhhp[PMTid],h1->GetMaximum()+5);
          lthreshold=new TLine(fQhlthresh[PMTid],0,fQhlthresh[PMTid],h1->GetMaximum()+5);
          lpeakold=new TLine(fQhlpeak[PMTid],0,fQhlpeak[PMTid],h1->GetMaximum()+5);
          lhhpnew=new TLine(fQhlhhpnew[PMTid],0,fQhlhhpnew[PMTid],h1->GetMaximum()+5);
          lthreshnew=new TLine(fQhlthreshnew[PMTid],0,fQhlthreshnew[PMTid],h1->GetMaximum()+5);
          lpeaknew=new TLine(fQhlpeaknew[PMTid],0,fQhlpeaknew[PMTid],h1->GetMaximum()+5);
          max=fQhlpeaknew[PMTid];
        }
        lhhpold->SetLineColor(2);
        lpeakold->SetLineColor(1);
        lthreshold->SetLineColor(4);
        lhhpold->SetLineStyle(2);
        lpeakold->SetLineStyle(2);
        lthreshold->SetLineStyle(2);
        h1->GetXaxis()->SetTitleColor(1);
        h1->GetYaxis()->SetTitleColor(1);
        h1->GetYaxis()->SetTitle("Counts");
        h1->SetMaximum(h1->GetMaximum()+5);
        if((h1->GetXaxis()->GetXmax()-h1->GetXaxis()->GetXmin())>400) {
          h1->GetXaxis()->SetRangeUser(max-50,max+200);
        }
        h1->Draw();
        lhhpold->Draw("SAME");
        lthreshold->Draw("SAME");
        lpeakold->Draw("SAME");
        leg4->Clear();
        leg4->SetFillStyle(0);
        leg4->SetBorderSize(0);
        // Plot current GF results, only if ratdb file was provided
        lhhpnew->SetLineColor(2);
        lpeaknew->SetLineColor(1);
        lthreshnew->SetLineColor(4);
        lhhpnew->SetLineStyle(1);
        lpeaknew->SetLineStyle(1);
        lthreshnew->SetLineStyle(1);
        lhhpnew->Draw("SAME");
        lthreshnew->Draw("SAME");
        lpeaknew->Draw("SAME");
        leg4->AddEntry(lthreshnew, "TH", "L");
        leg4->AddEntry(lpeaknew, "PK", "L");
        leg4->AddEntry(lhhpnew, "HHP", "L");
        leg4->AddEntry(lthreshold, "prev TH", "L");
        leg4->AddEntry(lpeakold, "prev PK", "L");
        leg4->AddEntry(lhhpold, "prev HHP", "L");
        leg4->Draw();
        //stringstream kek;
        //kek << PMTid << "-" << chargestr.c_str() << ".root";
        //c2->SaveAs(kek.str().c_str());
        img->Clear();
        img->FromPad(c2);
        sprintf(imagename, "PMT-%i-%s.png", PMTid, chargestr.c_str());
        img->WriteImage(imagename);
      }
      delete h1;
    }
  } //end loop over hists in ROOT file

  std::cout << " DONE. " << std::endl;

}
