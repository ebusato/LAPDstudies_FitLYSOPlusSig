#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;


void FitLYSOPlusSig()
{
  RooRealVar E("E", "Energy", 0, 1000);
  E.setBins(100);
  
  TFile* f = new TFile("data/run63.root");
  TTree* t = (TTree*) f->Get("tree");
  RooDataSet* data = new RooDataSet("data", "data_LYSOPlusSig", t, E);
  RooDataHist* hist = data->binnedClone();

  TFile* f_LYSO = new TFile("data/run78.root");
  TTree* t_LYSO = (TTree*) f_LYSO->Get("tree");
  RooDataSet* data_LYSO = new RooDataSet("data_LYSO", "data_LYSO", t_LYSO, E);
  RooDataHist* hist_LYSO = data_LYSO->binnedClone();
  RooHistPdf* histpdf_LYSO = new RooHistPdf("histpdf_LYSO","histpdf_LYSO",E,*hist_LYSO,0) ;
  RooRealVar lyso_yield("lyso_yield", "yield of lyso", 100, 0, 1000000);
  
  RooRealVar sig_peak_mean("sig_peak_mean", "mean of gaussian for signal peak", 511, 420, 610, "keV");
  RooRealVar sig_peak_sigma("sig_peak_sigma", "width of gaussian for signal peak", 40, 0, 90, "keV");
  RooGaussian sig_gaussian("gaussian", "gaussian for signal peak", E, sig_peak_mean, sig_peak_sigma);
  RooRealVar sig_yield("sig_yield", "yield signal peak", 1000, 0, 1000000);
  
  RooArgList shapes;
  RooArgList yields;
  shapes.add(*histpdf_LYSO); yields.add(lyso_yield);
  shapes.add(sig_gaussian);  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(*data);
  
  // Plot unbinned data and histogram pdf overlaid
  RooPlot* frame = E.frame(Bins(100)) ;
  data->plotOn(frame) ;
  totalPdf.plotOn(frame);
  frame->Draw();

}
