#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
using namespace RooFit ;

class Data {
public:
	Data(TTree* t, RooRealVar* E);
	RooDataHist* GetDataHistFromTH1(string TH1Name, string histName);
	double RunDuration();
	
	TTree* m_tree;
	TCut m_cut;
	
	RooDataHist* m_dh;
	RooRealVar* m_E;
};

Data::Data(TTree* t, RooRealVar* E) : m_tree(t),
		       m_cut(""),
		       m_dh(0),
		       m_E(E)
{}

RooDataHist* Data::GetDataHistFromTH1(string TH1Name, string histName)
{
  TH1* h0 = new TH1F(TH1Name.c_str(), TH1Name.c_str(), m_E->getBinning().numBins(), m_E->getMin(), m_E->getMax());
  TString s0("E[0]>>");
  s0 += TH1Name.c_str();
  m_tree->Draw(s0.Data(), m_cut, "goff");
  
  m_dh = new RooDataHist(histName.c_str(), histName.c_str(), *m_E, Import(*h0));
  cout << "no entries in RooDataHist " << histName << " = " << m_dh->sum(false) << endl;
  
  return m_dh;
}

// returns run duration in minutes
double Data::RunDuration()
{
  m_tree->Draw(">>evtlist", m_cut);
  TEventList *evtlist = (TEventList*)gDirectory->Get("evtlist");
  int Nevents = evtlist->GetN();
  m_tree->GetEntry(evtlist->GetEntry(Nevents - 1));
  //m_tree->GetEntry(m_tree->GetEntries()-1);
  double timeEnd = m_tree->GetLeaf("TimeStamp")->GetValue()*1/64e6;
  m_tree->GetEntry(0);
  double timeBeg = m_tree->GetLeaf("TimeStamp")->GetValue()*1/64e6;
  return (timeEnd - timeBeg)/60;
}

class Model {
public:
  Model(RooRealVar* E, Data* na22, Data* lyso);
  void MakeCalculationsSensitivity(int noEntries);
  RooFitResult* Fit();
  RooRealVar* m_lyso_yield;
  RooRealVar* m_sig_yield;
  RooRealVar* m_sig_peak_mean;
  RooRealVar* m_sig_peak_sigma;
  
  RooAddPdf* m_model;
  
  RooRealVar* m_E;
  
  Data* m_Na22;
  Data* m_LYSO;
};

Model::Model(RooRealVar* E, Data* na22, Data* lyso) : m_E(E),
						      m_Na22(na22),
	     					      m_LYSO(lyso)
{
  m_lyso_yield = new RooRealVar("lyso_yield", "yield of lyso", 100000, 0, 1000000);
  m_sig_yield = new RooRealVar("sig_yield", "yield signal peak", 1000, 0, 1000000);
  m_sig_peak_mean = new RooRealVar("sig_peak_mean", "mean of gaussian for signal peak", 511, 420, 610, "keV");
  m_sig_peak_sigma = new RooRealVar("sig_peak_sigma", "width of gaussian for signal peak", 40, 0, 90, "keV");
	
  cout << "here1 " << m_LYSO << " " << m_LYSO->m_dh << " " << m_E << endl;
  RooHistPdf* histpdf_LYSO = new RooHistPdf("histpdf_LYSO","histpdf_LYSO", *m_E,*(m_LYSO->m_dh),0);
  cout << "here2" << endl;
  RooGaussian* sig_gaussian = new RooGaussian("sig_gaussian", "gaussian for signal peak", *m_E, *m_sig_peak_mean, *m_sig_peak_sigma);
   cout << "here3" << endl;
   
  RooArgList* shapes = new RooArgList();
  RooArgList* yields = new RooArgList();
  shapes->add(*histpdf_LYSO); yields->add(*m_lyso_yield);
  shapes->add(*sig_gaussian);  yields->add(*m_sig_yield);
  m_sig_peak_mean->setVal(508);
  m_sig_peak_mean->setConstant(1);
  m_sig_peak_sigma->setVal(30);
  m_sig_peak_sigma->setConstant(1);
  m_sig_yield->setVal(1500);
  m_sig_yield->setConstant(1);
  m_model = new RooAddPdf("totalPdf", "sum of signal and background PDF's", *shapes, *yields);
  
}

RooFitResult* Model::Fit() 
{
  RooFitResult* fitRes = m_model->fitTo(*(m_Na22->m_dh), Extended(),Range("range_650_Max"));
  TCanvas* c1 = new TCanvas();
  RooPlot* frame = m_E->frame(Bins(100));
  m_Na22->m_dh->plotOn(frame); //, DrawOption("PX"));
  m_model->plotOn(frame, Range("range_whole"));
  m_model->plotOn(frame, Range("range_whole"), Components("sig_gaussian"),LineColor(kRed));
  m_model->plotOn(frame, Range("range_whole"), Components("histpdf_LYSO"),LineColor(kGreen+2));
  m_Na22->m_dh->plotOn(frame); //, DrawOption("PX"));
  frame->Draw();
  double xText = 0.55;
  double yShift = 0.07;
  PutText(xText, 0.81, kBlack, "LAPD");
 // PutText(xText, 0.85-yShift, kBlack, "LPC");
  PutText(xText, 0.81-1*yShift, kBlack, "^{22}Na (16 kBq)");
  stringstream ss;
  ss.precision(3);
  ss << "Run duration: " << m_Na22->RunDuration() << " min";
  PutText(xText, 0.81-2*yShift, kBlack, ss.str().c_str());
  return fitRes;
}

void MakeOTHinput(TString fileName, double LYSO_yield, double LYSO_yieldErr, double Na22_yield)
{
  ofstream of(fileName.Data());
  of << "+bg LYSO " << LYSO_yield << " " << LYSO_yieldErr << endl << endl;
  of << "+sig Sig " << Na22_yield << endl << endl;
  of << "+data " << LYSO_yield + Na22_yield << endl;
}

void Model::MakeCalculationsSensitivity(int noEntries)
{
  // Analytical stuff
 
  double N1_original = m_LYSO->m_dh->sum(kFALSE);
  double N1_original_err = sqrt(N1_original);
  cout << "N1_original = " << N1_original << " +- " << N1_original_err << endl;
    	
  RooAbsReal* sig_yield   = (RooAbsReal*) m_model->coefList().find("sig_yield");
  RooAbsReal* lyso_yield  = (RooAbsReal*) m_model->coefList().find("lyso_yield");
 
  double N1 = lyso_yield->getVal();
  double N2 = sig_yield->getVal();
  double N1_err = N1_original_err * N1 / N1_original;
  cout << "N1 = " << N1 << " +- " << N1_err << endl;
  cout << "N2 = " << N2 << endl;
  cout << "N2/N1 = " << N2 / N1 << endl;

  ///////////////////////////////////////////////////////////
  // Calculation of z (significance)
  
  // Used to compute selection efficiency of E[0] \in signalWindow
  RooAbsPdf* sig_gaussian = (RooAbsPdf*) m_model->pdfList().find("sig_gaussian");
  RooAbsPdf* histpdf_LYSO = (RooAbsPdf*) m_model->pdfList().find("histpdf_LYSO");
  RooAbsReal* igx_sig = sig_gaussian->createIntegral(*m_E,NormSet(*m_E),Range("signalWindow")) ;
  RooAbsReal* igx_lyso = histpdf_LYSO->createIntegral(*m_E,NormSet(*m_E),Range("signalWindow")) ;
  double int_sig_window = igx_sig->getVal();
  double int_lyso_window = igx_lyso->getVal();
  cout << "int[E|signal]_Norm[E] = " << int_sig_window << endl;
  cout << "int[E|LYSO]_Norm[E] = " << int_lyso_window << endl;
  
  // Alternative selection cut: E[0] > 425 && E[0] < 595 && E[1] > 425 && E[1] < 595 && fabs(T30[0] - T30[1]) < 3.6
  // We assume that this cut selects 50% of the signal (rough estimation, this needs to be estimated on simulation)
  // On run79.root, we estimate that the selection efficiency of this cut on LYSO background is 0.0033...% (without this cut we have 800000 events in E[0]>>h histogram,
  // and with the cut we have 2666 events)
  // Try to improve by also cutting on CTR
   int_sig_window = 0.5;
   int_lyso_window = 0.0034;
  
   double Z = 3.*sqrt(int_lyso_window)/int_sig_window;
  ///////////////////////////////////////////////////////////
  
  double time = float(noEntries)/24; // around 24 Hz, to be adjusted
  double tau = 41.e-3; // dead time around 40 ms, to be adjusted
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
  
  ///////////////////////////////////////////////////////////////
  // not used, just to know the values
  double N2_window = int_sig_window * N2;
  double N1_window = int_lyso_window * N1;
  double N1_window_err = N1_err * int_lyso_window;
  cout << "N2_window = "<< N2_window << endl;
  cout << "N1_window = "<< N1_window << " +- " << N1_window_err << endl;
  cout << "s_window/sqrt(b_window) = " << N2_window / sqrt(N1_window) << endl;
  ///////////////////////////////////////////////////////////////
  
  double m1 = N1/time;
  double m1prime = m1/(1+m2*tau*(solplus-1));
  double m2prime = solplus*m2/(1+m2*tau*(solplus-1));
  double N1prime = m1prime * time;
  double N2prime = m2prime * time;
  cout << "N1prime = " << N1prime << endl;
  cout << "N2prime = " << N2prime << endl;
  double N1prime_err = N1_original_err * N1prime / N1_original;
  
  
  double N1primeWindow = N1prime*int_lyso_window;
  double N2primeWindow = N2prime*int_sig_window;
  double N1primeErrorWindow = N1prime_err * int_lyso_window;
  cout << "N1prime in window = " << N1primeWindow << " +- " << N1primeErrorWindow << endl;
  cout << "N2prime in window = " << N2primeWindow << endl;
  cout << "s_window/sqrt(b_window) = "<< N2primeWindow/sqrt(N1primeWindow) << endl;

  // OTH stuff
  /*
  TH1F* hLYSO_gen = (TH1F*) m_LYSO.m_dh->createHistogram("E");
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
  */
//   TFile* f = new TFile("OTHinput/histos.root", "recreate");
//   hLYSO_gen_range->Write();
//   hSig_gen_range->Write();
//   hData_gen_range->Write();
//   f->Write();
//   f->Close();

  MakeOTHinput("OTHinput/inputYield.dat", N1primeWindow, N1primeErrorWindow, N2primeWindow);
  system("root -l -b -q load.C 'runSignificance.C(\"OTHinput/inputYield.dat\")'");
}



/*
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
*/


RooFitResult* FitLYSOPlusSig(string na22File, string lysoFile)
{
  RooRealVar* E = new RooRealVar("E", "Energy", 0, 1200, "keV");
  E->setBins(150);
  E->setRange("range_whole", E->getMin(), E->getMax());
  E->setRange("range_400_Max", 400, E->getMax());
  E->setRange("range_450_Max", 450, E->getMax());
  E->setRange("range_650_Max", 650, E->getMax());
  double signalWindow_min = 450;
  double signalWindow_max = 570;
//   double signalWindow_min = 900;
//   double signalWindow_max = 1100;
  E->setRange("range_sigWindow", signalWindow_min, signalWindow_max) ;
  
  TFile* f_LYSO = new TFile(lysoFile.c_str());
  TTree* t_LYSO = (TTree*) f_LYSO->Get("tree");
  
  TFile* f_Na22 = new TFile(na22File.c_str());
  TTree* t_Na22 = (TTree*) f_Na22->Get("tree");
 
  Data* dataLYSO = new Data(t_LYSO, E);
  Data* dataNa22 = new Data(t_Na22, E);

  dataLYSO->m_cut = "";
  dataNa22->m_cut = "Evt < 50000 && Sat[0] == 0";
  
  RooDataHist* hist_LYSO = dataLYSO->GetDataHistFromTH1("hE_lyso", "dhE_lyso");
  RooDataHist* hist_Na22 = dataNa22->GetDataHistFromTH1("hE_data", "dhE_data");
  
  Model* model = new Model(E, dataNa22, dataLYSO);
 
  model->Fit();
  
  //int noEntries = t_Na22->GetEntries();
  //MakeCalculationsSensitivity(hist_LYSO, model, E, noEntries);
  return 0;
}

void FitLYSOPlusSig()
{
//   FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "~/godaq_rootfiles/analysis_v2.10.0/run78.root");
  FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "~/godaq_rootfiles/analysis_v2.10.0/run78.root");
  
  //FitLYSOPlusSig("~/Travail/Imaging/serverAvirm/DPGA/DataBackup/godaq_rootfiles/analysis_v2.10.0/run63.root", 
//"~/Travail/Imaging/serverAvirm/DPGA/DataBackup/godaq_rootfiles/analysis_v2.10.0/run78.root");
}
