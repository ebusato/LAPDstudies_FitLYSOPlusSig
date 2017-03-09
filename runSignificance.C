///////////////////////////////////////////////////////////
//  > root -l load.C 'runSignificance.C("input.dat")'
///////////////////////////////////////////////////////////

#if defined EXECUTABLE || defined __CLING__

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

#include <TStopwatch.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TString.h>

#include "OpTHyLiC.h"

using namespace std;
using namespace OTH;

#endif

void runSignificance(const std::string& file1) {
  TStopwatch w;
  w.Start();

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // create OpTHyLiC instance
#if defined CPP11
  OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_mt19937); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_mt19937_64); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_minstd_rand); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_minstd_rand0); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux24_base); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux48_base); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux24); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_ranlux48); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::STD_knuth_b); // Provided by the C++11 standard library
//   OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatGammaHyper,OTH::TR3); // using default TRandom3-based pseudo-random number generator, and automatic seed
//   OpTHyLiC oth(OTH::SystMclimit,OTH::StatNormal); // Provided by the C++11 standard library
#else
  OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatLogN); // using default TRandom3-based pseudo-random number generator, and automatic seed
#endif

  oth.addChannel("ch1",file1);

  //oth.printSamples();
  
  const int Nexp=5e5;
  
  // By default, the signal strength is set to 1. If you want another value, just add a third argument to
  // oth.significance(...) with the desired value. 
  std::pair<double, double> s = oth.significance(OTH::SignifExpectedMed,Nexp);
  const double p=s.first;
  const double z=s.second;

  cout << endl << "Results (cpu time=" << w.CpuTime()<< " sec, real time=" << w.RealTime() << " sec): " << endl;
  cout << " -> p=" << p << endl;
  cout << " -> z=" << z << endl;

  TH1F* hLLRb=(TH1F*)oth.getHistoLLRb();
  hLLRb->GetXaxis()->SetTitleSize(0.05);
  hLLRb->GetXaxis()->SetTitleOffset(0.88);
  hLLRb->GetXaxis()->SetTitle("q_{ #mu}");
  hLLRb->GetYaxis()->SetTitle("");
  hLLRb->Rebin(10);
  hLLRb->SetLineColor(kGreen+3);
  TH1F* hLLRsb=(TH1F*)oth.getHistoLLRsb();
  hLLRsb->Rebin(10);

  TCanvas *c1 = new TCanvas("c1", "c1",668,105,700,500);
  c1->SetLogy();
  hLLRb->Draw();
  hLLRsb->Draw("same");
  double qmuobs=oth.computeLLRdata();
  TArrow* arr = new TArrow(qmuobs,hLLRb->GetMaximum()/5.,qmuobs,0,0.02,"|>");
  arr->SetLineWidth(3);
  arr->SetLineColor(kRed);
  arr->SetFillColor(kRed);
  arr->Draw();
  TLatex latex1;
  latex1.SetTextSize(0.05);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(qmuobs,hLLRb->GetMaximum()/2.,"q_{ #mu}^{obs}");

  w.Stop();
}
