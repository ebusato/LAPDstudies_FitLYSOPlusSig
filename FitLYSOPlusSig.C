#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
using namespace RooFit ;


void FitLYSOPlusSig()
{
  RooRealVar E("E", "Energy", 0, 1200, "keV");
  E.setBins(100);
  
  TFile* f = new TFile("data/run63.root");
  TTree* t = (TTree*) f->Get("tree");
  RooDataSet* data = new RooDataSet("data", "data_LYSOPlusSig", t, E);
  RooDataHist* hist = data->binnedClone();

  int noEntries = data->sumEntries();
  cout << "no entries in dataset = " << noEntries << endl;
  
  TFile* f_LYSO = new TFile("data/run78.root");
  TTree* t_LYSO = (TTree*) f_LYSO->Get("tree");
  RooDataSet* data_LYSO = new RooDataSet("data_LYSO", "data_LYSO", t_LYSO, E);
  RooDataHist* hist_LYSO = data_LYSO->binnedClone();
  RooHistPdf* histpdf_LYSO = new RooHistPdf("histpdf_LYSO","histpdf_LYSO",E,*hist_LYSO,0);
  RooRealVar lyso_yield("lyso_yield", "yield of lyso", 100, 0, 1000000);
  
  RooRealVar sig_peak_mean("sig_peak_mean", "mean of gaussian for signal peak", 511, 420, 610, "keV");
  RooRealVar sig_peak_sigma("sig_peak_sigma", "width of gaussian for signal peak", 40, 0, 90, "keV");
  RooGaussian sig_gaussian("sig_gaussian", "gaussian for signal peak", E, sig_peak_mean, sig_peak_sigma);
  RooRealVar sig_yield("sig_yield", "yield signal peak", 1000, 0, 1000000);
  
  RooArgList shapes;
  RooArgList yields;
  shapes.add(*histpdf_LYSO); yields.add(lyso_yield);
  shapes.add(sig_gaussian);  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(*data, Extended());
  
  // Plot unbinned data and histogram pdf overlaid
  RooPlot* frame = E.frame(Bins(100)) ;
  data->plotOn(frame) ;
  totalPdf.plotOn(frame);
  totalPdf.plotOn(frame, Components("sig_gaussian"),LineColor(kRed));
  totalPdf.plotOn(frame, Components("histpdf_LYSO"),LineColor(kGreen+2));
  //histpdf_LYSO->plotOn(frame);
  frame->Draw();
  PutText(0.65, 0.85, kBlack, "Na22 (16 kBq)");
  
  E.setRange("signal_window", 420, 600) ;
  
  RooAbsReal* igx_sig = sig_gaussian.createIntegral(E,NormSet(E),Range("signal_window")) ;
  RooAbsReal* igx_lyso = histpdf_LYSO->createIntegral(E,NormSet(E),Range("signal_window")) ;
  double int_sig_window = igx_sig->getVal();
  double int_lyso_window = igx_lyso->getVal();
  double yield_sig_window = int_sig_window * sig_yield.getVal();
  double yield_lyso_window = int_lyso_window * lyso_yield.getVal();
  cout << "int[E|signal]_Norm[E] = " << int_sig_window << endl ;
  cout << "int[E|LYSO]_Norm[E] = " << int_lyso_window << endl ;
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
