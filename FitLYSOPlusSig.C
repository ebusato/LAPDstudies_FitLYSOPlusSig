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

TH1F* MakeHistoInRange(TH1F* h, double min, double max, bool setErrorsToZero=false)
{
  TString newName(h->GetName());
  newName += "_inRange";
  int NbinsRange = h->GetXaxis()->FindBin(max) - h->GetXaxis()->FindBin(min) + 1;
  double startVal = h->GetXaxis()->GetBinCenter(h->GetXaxis()->FindBin(min)) - h->GetXaxis()->GetBinWidth(h->GetXaxis()->FindBin(min))/2.;
  double stopVal = h->GetXaxis()->GetBinCenter(h->GetXaxis()->FindBin(max)) + h->GetXaxis()->GetBinWidth(h->GetXaxis()->FindBin(max))/2.;
  TH1F* hNew = new TH1F(newName.Data(), newName.Data(), NbinsRange, startVal, stopVal);
  hNew->SetBinContent(0, 0);
  hNew->SetBinContent(NbinsRange+1, 0);
  for(int i = 1; i < NbinsRange + 1; ++i) {
	  double val = hNew->GetXaxis()->GetBinCenter(i);
	  int origHistoBin = h->GetXaxis()->FindBin(val);
	  double binContent = h->GetBinContent(origHistoBin);
	  double binError = h->GetBinError(origHistoBin);
	  //cout << i << "  " << origHistoBin<< "  " << val << "  " << binContent << "  " << binError << endl;
	  hNew->SetBinContent(i, binContent);
	  if(!setErrorsToZero) {
		hNew->SetBinError(i, binError);
	  }
  }
  hNew->SetMarkerColor(kRed);
  hNew->SetLineColor(kRed);
  return hNew;
}

void MakeCalculationsSensitivity(RooDataHist* hist_LYSO, RooAddPdf* model, RooRealVar* E, int noEntries)
{
  double signalWindow_min = 450;
  double signalWindow_max = 570;
/*
  double signalWindow_min = 505;
  double signalWindow_max = 510;*/
  E->setRange("signalWindow", signalWindow_min, signalWindow_max) ;
	
  RooAbsPdf* sig_gaussian = (RooAbsPdf*) model->pdfList().find("sig_gaussian");
  RooAbsPdf* histpdf_LYSO = (RooAbsPdf*) model->pdfList().find("histpdf_LYSO");
  RooAbsReal* sig_yield   = (RooAbsReal*) model->coefList().find("sig_yield");
  RooAbsReal* lyso_yield  = (RooAbsReal*) model->coefList().find("lyso_yield");
 
  // Analytical stuff
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
  
  double N1 = lyso_yield->getVal();
  double N2 = sig_yield->getVal();
  double time = float(noEntries)/24; // around 24 Hz, to be adjusted
  double tau = 41.15e-3; // dead time around 40 ms, to be adjusted
  double m2 = N2/time;
  double poly_a = N2*N2;
  double Z = 3.*sqrt(int_lyso_window)/int_sig_window;
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
  cout << "s/sqrt(b) = "<< N2prime/sqrt(N1prime)*int_sig_window/sqrt(int_lyso_window) << endl;

  // OTH stuff
  TH1F* hLYSO_gen = (TH1F*) hist_LYSO->createHistogram("E");
  hLYSO_gen->Sumw2();
  
  //hLYSO_gen->Scale(lyso_yield->getVal()/hLYSO_gen->Integral());
  hLYSO_gen->Scale(N1prime/hLYSO_gen->Integral());
  
  cout << "ERROR=" << hLYSO_gen->GetEntries() << "  " << hLYSO_gen->Integral() << "  " << hLYSO_gen->GetBinContent(hLYSO_gen->GetXaxis()->FindBin(400)) << " " << 
hLYSO_gen->GetBinError(hLYSO_gen->GetXaxis()->FindBin(400)) << endl;
  
  //RooDataHist* dhSig_gen = sig_gaussian->generateBinned(*E, sig_yield->getVal(), ExpectedData());
  RooDataHist* dhSig_gen = sig_gaussian->generateBinned(*E, N2prime, ExpectedData());
  
  TH1F* hSig_gen = (TH1F*) dhSig_gen->createHistogram("E");
  cout << "hSig_gen: " << hSig_gen->Integral() << endl;
    
  TH1F* hLYSO_gen_range = MakeHistoInRange(hLYSO_gen, signalWindow_min, signalWindow_max, true);
  TH1F* hSig_gen_range = MakeHistoInRange(hSig_gen, signalWindow_min, signalWindow_max, true);
  TCanvas* c2 = new TCanvas("c2","c2",1200,400);
  c2->Divide(3,1);
  c2->cd(1);
  hLYSO_gen->Draw("e");
  hLYSO_gen_range->Draw("histsame");
  hLYSO_gen_range->Draw("same");
  c2->cd(2);
  hSig_gen->Draw("e");
  hSig_gen_range->Draw("histsame");
  hSig_gen_range->Draw("same");
  c2->cd(3);
  TH1F* hData_gen_range = (TH1F*) hLYSO_gen_range->Clone("hData_gen_range");
  hData_gen_range->Add(hSig_gen_range);
  hData_gen_range->GetYaxis()->SetRangeUser(0, hData_gen_range->GetMaximum()*1.3);
  hData_gen_range->Draw("e");
  hLYSO_gen_range->Draw("same");
  hSig_gen_range->Draw("same");
  
//   TFile* f = new TFile("OTHinput/histos.root", "recreate");
//   hLYSO_gen_range->Write();
//   hSig_gen_range->Write();
//   hData_gen_range->Write();
//   f->Write();
//   f->Close();
  
  ofstream of("OTHinput/inputYield.dat");
  of << "+bg LYSO " << N1prime*int_lyso_window << endl << endl;
  of << "+sig Sig " << N2prime*int_sig_window << endl << endl;
  of << "+data " << N1prime*int_lyso_window + N2prime*int_sig_window << endl;
  system("root -l -b -q load.C 'runSignificance.C(\"OTHinput/inputYield.dat\")'");
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
  PutText(0.65, 0.85, kBlack, "LAPD");
  PutText(0.65, 0.85-yShift, kBlack, "LPC");
  PutText(0.65, 0.85-2*yShift, kBlack, "Na22 (16 kBq)");
  
  MakeCalculationsSensitivity(hist_LYSO, model, E, noEntries);
  return 0;
}

void FitLYSOPlusSig()
{
  //FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "~/godaq_rootfiles/analysis_v2.10.0/run79.root");
  FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run63.root", "~/godaq_rootfiles/analysis_v2.10.0/run78.root");
  
  //FitLYSOPlusSig("~/Travail/Imaging/serverAvirm/DPGA/DataBackup/godaq_rootfiles/analysis_v2.10.0/run63.root", 
//"~/Travail/Imaging/serverAvirm/DPGA/DataBackup/godaq_rootfiles/analysis_v2.10.0/run78.root");
}
