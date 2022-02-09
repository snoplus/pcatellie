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
#include "TLegend.h"

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {

  string Filename1, Filename2;
  int runid1, runid2;

  if (argc != 5){
    std::cout<<" Usage: Compare bench1.root bench2.root run1 run2"<<std::endl;
    exit(0);
  } else {
    Filename1 = argv[1];
    Filename2 = argv[2];
    runid1 = atoi(argv[3]);
    runid2 = atoi(argv[4]);
  }

  gROOT->SetBatch(kTRUE);

  // load files
  TFile *lb;
  TFile *t2;

  // benchmark root files
  lb = new TFile(Filename1.c_str());
  t2 = new TFile(Filename2.c_str());

  // new root
  //TFile *MyFile = new TFile("Output.root","RECREATE");

  // cable
  TCanvas *c1= new TCanvas ("LogResidualsC", "LogResidualsC", 800, 600);
  TCanvas *lbRC = (TCanvas*)lb->Get("AllPMTHitTimesRLog");
  TH1F *lbR = (TH1F*)lbRC->GetPrimitive("AllPMTHitTimesRLog");
  TCanvas *t2C = (TCanvas*)t2->Get("AllPMTHitTimesRLog");
  TH1F *t2R = (TH1F*)t2C->GetPrimitive("AllPMTHitTimesRLog");

  c1->cd();
  c1->SetLogy();
  lbR->Draw();
  t2R->SetLineColor(kRed);
  t2R->Draw("same");

  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->AddEntry(lbR, "This","l");
  legend->AddEntry(t2R, "Previous","l");
  legend->Draw();

  c1->SaveAs("peak_comp_1.png");
  c1->Close();

  TCanvas *c2= new TCanvas ("LogResidualsC2", "LogResidualsC2", 800, 600);
  TCanvas *lbRC2 = (TCanvas*)lb->Get("AllPMTHitTimesRLog2");
  TH1F *lbR2 = (TH1F*)lbRC2->GetPrimitive("AllPMTHitTimesRLog2");
  TCanvas *t2C2 = (TCanvas*)t2->Get("AllPMTHitTimesRLog2");
  TH1F *t2R2 = (TH1F*)t2C2->GetPrimitive("AllPMTHitTimesRLog2");

  c2->cd();
  c2->SetLogy();
  lbR2->SetLineColorAlpha(kBlack, 0.9);
  lbR2->Draw();
  t2R2->SetLineColorAlpha(kRed, 0.9);
  t2R2->Draw("same");

  TLegend *legend2 = new TLegend(0.1,0.7,0.3,0.9);
  legend2->AddEntry(lbR2, "This","l");
  legend2->AddEntry(t2R2, "Previous","l");
  legend2->Draw();

  c2->SaveAs("peak_comp_2.png");
  c2->Close();

  TCanvas *c3= new TCanvas ("Compare", "Compare", 800, 600);
  TCanvas *c6= new TCanvas ("Compare2", "Compare2", 800, 600);
  TCanvas *lbRC_new = (TCanvas*)lb->Get("AllPMTHitTimesRLog");
  TH1F *lbOrig = (TH1F*)lbRC_new->GetPrimitive("AllPMTHitTimesRLog");
  TCanvas *t2C_new = (TCanvas*)t2->Get("AllPMTHitTimesRLog");
  TH1F *teOrig = (TH1F*)t2C_new->GetPrimitive("AllPMTHitTimesRLog");
  TGraph *compare;
  compare = new TGraph();

  for (int i=0; i<lbOrig->GetNbinsX(); i++){
    compare->SetPoint(i, lbOrig->GetBinCenter(i), teOrig->GetBinContent(i)-lbOrig->GetBinContent(i) );
  }

  TLine *l = new TLine(293,-140,293,155);
  TLine *l2 = new TLine(310,-140,310,155);
  l->SetLineColor(kRed);
  l2->SetLineColor(kRed);

  c3->cd();
  //c3->SetLogy();
  compare->SetMarkerSize(2);
  compare->SetMarkerStyle(22);
  compare->SetMarkerColorAlpha(kBlack, 0.25);
  compare->GetXaxis()->SetTitle("time [ns]");
  compare->SetTitle("This - Previous");
  //compare->GetYaxis()->SetRangeUser(-140,155);
  compare->Draw("AP");
  l->Draw("same");
  l2->Draw("same");
  c3->SaveAs("peak_comp_3.png");
  c3->Close();

  // c6->cd();
  // c6->SetLogy();
  // compare->SetMarkerSize(3);
  // compare->SetMarkerStyle(5);
  // compare->GetXaxis()->SetTitle("time [ns]");
  // compare->SetTitle("TELLIE - LB");
  // //compare->GetYaxis()->SetRangeUser(-140,155);
  // compare->Draw("AP");
  // l->Draw("same");
  // l2->Draw("same");
  // c6->Write();
  // c6->Close();


  TCanvas *c4= new TCanvas ("Compare3", "Compare3", 800, 600);
  TCanvas *c5= new TCanvas ("Compare4", "Compare4", 800, 600);
  TCanvas *lbRC2_new = (TCanvas*)lb->Get("AllPMTHitTimesRLog2");
  TH1F *lbOrig2 = (TH1F*)lbRC2_new->GetPrimitive("AllPMTHitTimesRLog2");
  TCanvas *t2C2_new = (TCanvas*)t2->Get("AllPMTHitTimesRLog2");
  TH1F *teOrig2 = (TH1F*)t2C2_new->GetPrimitive("AllPMTHitTimesRLog2");
  TGraph *compare2;
  compare2 = new TGraph();

  for (int i=0; i<lbOrig2->GetNbinsX(); i++){
    compare2->SetPoint(i, lbOrig2->GetBinCenter(i), abs(teOrig2->GetBinContent(i)-lbOrig2->GetBinContent(i)) );
  }

  c4->cd();
  compare2->SetMarkerSize(2);
  compare2->SetMarkerStyle(22);
  compare2->SetMarkerColorAlpha(kBlack, 0.25);
  compare2->GetXaxis()->SetTitle("time [ns]");
  compare2->SetTitle("This - Previous, peak");
  //compare2->GetYaxis()->SetRangeUser(-500,100);
  compare2->Draw("AP");
  c4->SaveAs("peak_comp_4.png");
  c4->Close();

  // c5->cd();
  // c5->SetLogy();
  // compare2->SetMarkerSize(3);
  // compare2->SetMarkerStyle(5);
  // compare2->GetXaxis()->SetTitle("time [ns]");
  // compare2->SetTitle("TELLIE - LB, peak");
  // compare2->SetMinimum(10e-2);
  // //compare2->SetMaximum(100);
  // compare2->Draw("AP");
  // c5->Write();
  // c5->Close();

  // TCanvas *c7= new TCanvas ("Compare_peak", "Compare_peak", 800, 600);
  //
  // c7->cd();
  // compare->SetMarkerSize(3);
  // compare->SetMarkerStyle(5);
  // compare->GetXaxis()->SetTitle("time [ns]");
  // compare->SetTitle("TELLIE - LB");
  // compare->GetXaxis()->SetRangeUser(280,320);
  // compare->Draw("AP");
  // c7->Write();
  // c7->Close();

  //MyFile->Close();

  // Store to log file
  std::stringstream logFile_namess;
  string logFile_name;
  FILE *logFile;
  logFile_namess.str("");
  logFile_namess << runid1 << "-" << runid2 << "_peak_compare.log";
  logFile_name = logFile_namess.str();
  logFile = fopen(logFile_name.c_str(), "w");

  fprintf(logFile, "Run1: %i\n", runid1);
  fprintf(logFile, "Run2: %i\n", runid2);
  fprintf(logFile, "root1: %s\n", Filename1.c_str());
  fprintf(logFile, "root2: %s\n", Filename2.c_str());

  return 0;
}
