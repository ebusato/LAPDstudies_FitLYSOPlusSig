#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
using namespace RooFit ;


RooDataHist* GetDataHistFromTH1(TTree* t, RooRealVar* var, string TH1Name, string histName)
{
  //RooDataSet* data = new RooDataSet("data", "data", t, var);
  //RooDataHist* hist = data->binnedClone();
  //return hist

  TH1* h = new TH1F(TH1Name.c_str(), TH1Name.c_str(), var->getBinning().numBins(), var->getMin(), var->getMax());
  TString s("E[0]>>");
  s += TH1Name.c_str();
  t->Draw(s.Data(), "", "goff");
  return new RooDataHist(histName.c_str(), histName.c_str(), *var, Import(*h)) ;
}

RooAddPdf* MakeModel(RooRealVar* E, RooDataHist* hist_LYSO)
{
  RooHistPdf* histpdf_LYSO = new RooHistPdf("histpdf_LYSO","histpdf_LYSO", *E,*hist_LYSO,0);
  RooRealVar* lyso_yield = new RooRealVar("lyso_yield", "yield of lyso", 100000, 0, 1000000);
 
  RooRealVar* sig_peak_mean = new RooRealVar("sig_peak_mean", "mean of gaussian for signal peak", 511, 420, 610, "keV");
  RooRealVar* sig_peak_sigma = new RooRealVar("sig_peak_sigma", "width of gaussian for signal peak", 40, 0, 90, "keV");
  RooGaussian* sig_gaussian = new RooGaussian("sig_gaussian", "gaussian for signal peak", *E, *sig_peak_mean, *sig_peak_sigma);
  RooRealVar* sig_yield = new RooRealVar("sig_yield", "yield signal peak", 1000, 0, 1000000);
  
  RooArgList* shapes = new RooArgList();
  RooArgList* yields = new RooArgList();
  shapes->add(*histpdf_LYSO); yields->add(*lyso_yield);
  shapes->add(*sig_gaussian);  yields->add(*sig_yield);
  RooAddPdf* model = new RooAddPdf("totalPdf", "sum of signal and background PDF's", *shapes, *yields);
  return model;
}

void MakeCalculationsSensitivity(RooDataHist* hist_LYSO, RooAddPdf* model, RooRealVar* E, int noEntries)
{
  RooAbsPdf* sig_gaussian = (RooAbsPdf*) model->pdfList().find("sig_gaussian");
  RooAbsReal* sig_yield = (RooAbsReal*) model->coefList().find("sig_yield");
  RooAbsReal* lyso_yield = (RooAbsReal*) model->coefList().find("lyso_yield");
  
  // OTH stuff
  TCanvas* c2 = new TCanvas("c2","c2",800,400);
  c2->Divide(2,1);
  c2->cd(1);
  TH1F* hLYSO_gen = (TH1F*) hist_LYSO->createHistogram("E");
  hLYSO_gen->Sumw2();
  hLYSO_gen->Scale(lyso_yield->getVal()/hLYSO_gen->Integral());
  hLYSO_gen->Draw("e");
  cout << "ERROR=" << hLYSO_gen->GetEntries() << "  " << hLYSO_gen->Integral() << "  " << hLYSO_gen->GetBinContent(hLYSO_gen->GetXaxis()->FindBin(400)) << " " << 
hLYSO_gen->GetBinError(hLYSO_gen->GetXaxis()->FindBin(400)) << endl;
  c2->cd(2);
  RooDataHist* dhSig_gen = sig_gaussian->generateBinned(*E, sig_yield->getVal(), ExpectedData());
  TH1F* hSig_gen = (TH1F*) dhSig_gen->createHistogram("E");
  hSig_gen->Draw("e");
  cout << "hSig_gen: " << hSig_gen->Integral() << endl;
  
  RooAbsPdf* histpdf_LYSO = (RooAbsPdf*) model->pdfList().find("histpdf_LYSO");
  // Analytical stuff
  E->setRange("signalWindow", 420, 600) ;
  RooAbsReal* igx_sig = sig_gaussian->createIntegral(*E,NormSet(*E),Range("signalWindow")) ;
  RooAbsReal* igx_lyso = histpdf_LYSO->createIntegral(*E,NormSet(*E),Range("signalWindow")) ;
  double int_sig_window = igx_sig->getVal();
  double int_lyso_window = igx_lyso->getVal();
  cout << "int[E|signal]_Norm[E] = " << int_sig_window << endl;
  cout << "int[E|LYSO]_Norm[E] = " << int_lyso_window << endl;
  
  double yield_sig_window = int_sig_window * sig_yield->getVal();
  double yield_lyso_window = int_lyso_window * lyso_yield->getVal();
  cout << "sig in window = "<< yield_sig_window << endl;
  cout << "lyso in window = "<< yield_lyso_window << endl;
  cout << "s/sqrt(b) = " << yield_sig_window / sqrt(yield_lyso_window) << endl;
  
  double Z = 3.;
  double N1 = yield_lyso_window;
  double N2 = yield_sig_window;
  double time = float(noEntries)/24; // around 24 Hz, to be adjusted
  double tau = 40e-3; // dead time around 40 ms, to be adjusted
  double m2 = N2/time;
  double poly_a = N2*N2;
  double poly_b = -1*Z*Z*N1*m2*tau;
  double poly_c = Z*Z*N1*(m2*tau - 1);
  double delta = poly_b*poly_b - 4*poly_a*poly_c;
  cout << "delta = " << delta << endl;
  double solplus = (-poly_b + sqrt(delta))/(2*poly_a);
  double solminus = (-poly_b - sqrt(delta))/(2*poly_a);
  cout << "sol+ = " << solplus << endl;
  cout << "sol- = " << solminus << endl;
  
  double m1 = N1/time;
  double m1prime = m1/(1+m2*tau*(solplus-1));
  double m2prime = solplus*m2/(1+m2*tau*(solplus-1));
  double N1prime = m1prime * time;
  double N2prime = m2prime * time;
  cout << "N1prime = " << N1prime << endl;
  cout << "N2prime = " << N2prime << endl;
  cout << "s/sqrt(b) = "<< N2prime/sqrt(N1prime) << endl;
 
}

RooFitResult* FitLYSOPlusSig(string dataFile, string lysoFile)
{
  RooRealVar* E = new RooRealVar("E", "Energy", 0, 1200, "keV");
  E->setBins(150);
  E->setRange("whole", E->getMin(), E->getMax());
  
  TFile* f_LYSO = new TFile(lysoFile.c_str());
  TTree* t_LYSO = (TTree*) f_LYSO->Get("tree");
  RooDataHist* hist_LYSO = GetDataHistFromTH1(t_LYSO, E, "hE_lyso", "dhE_lyso");
  cout << "no entries in RooDataHist LYSO = " << hist_LYSO->sum(false) << endl;
  RooAddPdf* model = MakeModel(E, hist_LYSO);
  
  TFile* f = new TFile(dataFile.c_str());
  TTree* t = (TTree*) f->Get("tree");
  int noEntries = t->GetEntries();

  RooDataHist* hist = GetDataHistFromTH1(t, E, "hE_data", "dhE_data");
  cout << "no entries in RooDataHist data = " << hist->sum(false) << endl;

  //E->setRange("betaContinuum", 700, E->getMax());
  RooFitResult* fitRes = model->fitTo(*hist, Extended()); //,Range("betaContinuum"));
  TCanvas* c1 = new TCanvas();
  RooPlot* frame = E->frame(Bins(100)) ;
  hist->plotOn(frame); //, DrawOption("PX"));
  model->plotOn(frame, Range("whole"));
  model->plotOn(frame, Range("whole"), Components("sig_gaussian"),LineColor(kRed));
  model->plotOn(frame, Range("whole"), Components("histpdf_LYSO"),LineColor(kGreen+2));
  hist->plotOn(frame); //, DrawOption("PX"));
  frame->Draw();
  double yShift = 0.07;
  PutText(0.7, 0.85, kBlack, "LAPD");
  PutText(0.7, 0.85-yShift, kBlack, "LPC");
  PutText(0.7, 0.85-2*yShift, kBlack, "Na22 (16 kBq)");
  
  MakeCalculationsSensitivity(hist_LYSO, model, E, noEntries);
  return 0;
}

void FitLYSOPlusSig()
{
  //FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "~/godaq_rootfiles/analysis_v2.10.0/run79.root");
  FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run63.root", "~/godaq_rootfiles/analysis_v2.10.0/run78.root");
}
