#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooPlot.h"

#include "Headers.h"

using namespace RooFit ;

// Used to compute selection efficiency of E[0] \in signalWindow
//   RooAbsPdf* sig_gaussian = (RooAbsPdf*) m_model->pdfList().find("sig_gaussian");
//   RooAbsPdf* histpdf_LYSO = (RooAbsPdf*) m_model->pdfList().find("histpdf_LYSO");
//   RooAbsReal* igx_sig = sig_gaussian->createIntegral(*m_E,NormSet(*m_E),Range("signalWindow")) ;
//   RooAbsReal* igx_lyso = histpdf_LYSO->createIntegral(*m_E,NormSet(*m_E),Range("signalWindow")) ;
//   double int_sig_window = igx_sig->getVal();
//   double int_lyso_window = igx_lyso->getVal();
//   cout << "int[E|signal]_Norm[E] = " << int_sig_window << endl;
//   cout << "int[E|LYSO]_Norm[E] = " << int_lyso_window << endl;

// Alternative selection cut: E[0] > 425 && E[0] < 595 && E[1] > 425 && E[1] < 595 && fabs(T30[0] - T30[1]) < 3.6
// We assume that this cut selects 50% of the signal (rough estimation, this needs to be estimated on simulation)
// On run79.root, we estimate that the selection efficiency of this cut on LYSO background is 0.0033... (without this cut we have 800000 events in E[0]>>h histogram,
// and with the cut we have 2666 events)
// Try to improve by also cutting on CTR
double eff_signal = 0.5;
double eff_lyso = 0.0034;

double peakEff = 0.9;
  
Data::Data(TTree* t, RooRealVar* E) : m_tree(t),
		       m_cut(""),
		       m_dh(0),
		       m_E(E)
{}

RooDataHist* Data::GetDataHistFromTH1(string TH1Name, string histName)
{
  TH1F* h0 = new TH1F(TH1Name.c_str(), TH1Name.c_str(), m_E->getBinning().numBins(), m_E->getMin(), m_E->getMax());
  TString s0;
  s0 += "E[0]>>";
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
  return (timeEnd - timeBeg)/60.;
}


Model::Model(RooRealVar* E, Data* lyso, RooAbsPdf* sig_gaussian) : m_E(E), 
								   m_LYSO(lyso)
{
  m_lyso_yield = new RooRealVar("lyso_yield", "yield of lyso", 100000, 0, 1000000);
  m_sig_yield = new RooRealVar("sig_yield", "yield signal peak", 1000, 0, 1000000);
  
  RooDataHist* hist_LYSO = lyso->GetDataHistFromTH1("hE_lyso", "dhE_lyso");
  RooHistPdf* histpdf_LYSO = new RooHistPdf("histpdf_LYSO","histpdf_LYSO", *m_E,*hist_LYSO,0);
   
  RooArgList* shapes = new RooArgList();
  RooArgList* yields = new RooArgList();
  shapes->add(*histpdf_LYSO); yields->add(*m_lyso_yield);
  shapes->add(*sig_gaussian);  yields->add(*m_sig_yield);
  m_model = new RooAddPdf("totalPdf", "sum of signal and background PDF's", *shapes, *yields);
}

double Model::GetSigYield(Data* data) 
{
  double Ndata = data->m_dh->sum(kFALSE);
  cout << "Ndata = " << Ndata << " " << m_lyso_yield->getVal() << endl;
  return Ndata - m_lyso_yield->getVal();
}

Result* Model::Fit(Data* data) 
{
  if(data->m_dh == 0) {
    data->GetDataHistFromTH1("hE_data", "dhE_data");
  }
  
  m_sig_yield->setConstant();
//   RooMsgService::instance().setSilentMode(true);
 m_model->fitTo(*(data->m_dh), Extended(),Range("range_650_Max"), PrintEvalErrors(-1));
 
 Result* res = new Result(data->RunDuration(), data->m_activity, GetSigYield(data), m_lyso_yield->getVal(), m_LYSO->m_dh->sum(kFALSE));
 return res;
}

void Model::Plot(Data* data, Result* res, bool plotData)
{
  m_sig_yield->setVal(peakEff*res->m_Nsig);
  m_lyso_yield->setVal(res->m_Nlyso);
  
  RooPlot* frame = m_E->frame(Bins(100));
  if(plotData) {
    data->m_dh->plotOn(frame); //, DrawOption("PX"));
  }
  m_model->plotOn(frame, Range("range_250_Max"));
  m_model->plotOn(frame, Range("range_250_Max"), Components("sig_gaussian"),LineColor(kRed));
  m_model->plotOn(frame, Range("range_200_Max"), Components("histpdf_LYSO"),LineColor(kGreen+2));
  if(plotData) {
    data->m_dh->plotOn(frame); //, DrawOption("PX"));
  }
  frame->Draw();
 
  double xText = 0.55;
  double yShift = 0.07;
  PutText(xText, 0.81, kBlack, "LAPD");
 // PutText(xText, 0.85-yShift, kBlack, "LPC");
  PutText(xText, 0.81-1*yShift, kBlack, "^{22}Na (16 kBq)");
  stringstream ss;
  ss.precision(3);
  ss << "Run duration: " << res->m_time << " min";
  PutText(xText, 0.81-2*yShift, kBlack, ss.str().c_str());
}


// TH1F* MakeHistoInRange(TH1F* h, double min, double max, bool setErrorsToZero=false)
// {
//   TString newName(h->GetName());
//   newName += "_inRange";
//   int NbinsRange = h->GetXaxis()->FindBin(max) - h->GetXaxis()->FindBin(min) + 1;
//   double startVal = h->GetXaxis()->GetBinCenter(h->GetXaxis()->FindBin(min)) - h->GetXaxis()->GetBinWidth(h->GetXaxis()->FindBin(min))/2.;
//   double stopVal = h->GetXaxis()->GetBinCenter(h->GetXaxis()->FindBin(max)) + h->GetXaxis()->GetBinWidth(h->GetXaxis()->FindBin(max))/2.;
//   TH1F* hNew = new TH1F(newName.Data(), newName.Data(), NbinsRange, startVal, stopVal);
//   hNew->SetBinContent(0, 0);
//   hNew->SetBinContent(NbinsRange+1, 0);
//   for(int i = 1; i < NbinsRange + 1; ++i) {
// 	  double val = hNew->GetXaxis()->GetBinCenter(i);
// 	  int origHistoBin = h->GetXaxis()->FindBin(val);
// 	  double binContent = h->GetBinContent(origHistoBin);
// 	  double binError = h->GetBinError(origHistoBin);
// 	  //cout << i << "  " << origHistoBin<< "  " << val << "  " << binContent << "  " << binError << endl;
// 	  hNew->SetBinContent(i, binContent);
// 	  if(!setErrorsToZero) {
// 		hNew->SetBinError(i, binError);
// 	  }
//   }
//   hNew->SetMarkerColor(kRed);
//   hNew->SetLineColor(kRed);
//   return hNew;
// }

Result::Result(double time, double activity, double Nsig, double Nlyso, double NlysoOrigSample): m_time(time),
									   m_activity(activity),
									   m_Nsig(Nsig),
									   m_Nlyso(Nlyso),
									   m_NlysoOrigSample(NlysoOrigSample),
									   m_timeOrig(time),
									   m_activityOrig(activity),
									   m_NsigOrig(Nsig),
									   m_NlysoOrig(Nlyso)
{
  m_deadTime = 41.15e-3;
  m_NlysoErr = 0;
  m_NlysoErrOrig = 0;
}

void Result::Print() 
{
	cout << "Printing result: "<< endl;
	cout << "  -> Activity = " << m_activity << endl;
	cout << "  -> Time = " << m_time*60 << " sec" << endl;
	cout << "  -> Nsig = " << m_Nsig << endl;
	cout << "  -> Nlyso = " << m_Nlyso << " +- " << m_NlysoErr << endl;
	cout << "  -> NlysoOrig = " << m_NlysoOrig << endl;
	double rateSig = m_Nsig / m_time / 60.;
	double rateLyso = m_Nlyso / m_time / 60.;
	cout << "  -> rate sig = " << rateSig << endl;
	cout << "  -> rate lyso = " << rateLyso << endl;
	cout << "  -> rate tot = " << rateSig+rateLyso << endl; 
}

double Result::SolveForAlpha(double Z) 
{  
  double mSig = SigRate();
  double poly_a = m_Nsig*m_Nsig;
  double poly_b = -1*Z*Z*m_Nlyso*mSig*m_deadTime;
  double poly_c = Z*Z*m_Nlyso*(mSig*m_deadTime - 1);
  double delta = poly_b*poly_b - 4*poly_a*poly_c;
  cout << "delta = " << delta << endl;
  double solplus = (-poly_b + sqrt(delta))/(2*poly_a);
  double solminus = (-poly_b - sqrt(delta))/(2*poly_a);
  cout << "sol+ = " << solplus << endl;
  cout << "sol- = " << solminus << endl;
  
  return solplus;
}

void Result::ApplyEff(double effSig, double effLyso)
{
  m_Nsig *= effSig;
  m_Nlyso *= effLyso;
}

void Result::WriteOTHFile(TString fileName) 
{
  ofstream of(fileName.Data());
//   of << "+bg LYSO " << m_Nlyso << " " << m_NlysoErr << endl << endl;
  of << "+bg LYSO " << m_Nlyso << endl << endl;
  of << "+sig Sig " << m_Nsig << endl << endl;
  of << "+data " << m_Nsig + m_Nlyso << endl;
}

std::pair<double, double> Result::CalcZ(TString fileName)
{
  /*
  TString cmd("root -l -b -q load.C 'runSignificance.C(\"");
  cmd += fileName;
  cmd += "\")'";
  cout << "Command is: " << cmd.Data() << endl;
  system(cmd.Data());
  */
  
  WriteOTHFile("OTHinput/inputYield.dat");
  OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatLogN);
  oth.addChannel("ch1",fileName.Data());
  const int Nexp=5e5;
  std::pair<double, double> s = oth.significance(OTH::SignifExpectedMed,Nexp);
//    std::pair<double, double> s = oth.significance(OTH::SignifObserved,Nexp);
  const double p=s.first;
  const double z=s.second;
  
  return std::make_pair(m_Nsig/sqrt(m_Nlyso), z);
}

void Result::RescaleActivity(double factor) {
  // compute rates in sec-1
  double mLYSO = m_Nlyso / m_time / 60.;
  double mSig = m_Nsig / m_time / 60.;
  
  // compute rescaled quantities
  double mLYSOprime = mLYSO / (1+mSig*m_deadTime*(factor-1));
  double mSigprime = factor*mSig / (1+mSig*m_deadTime*(factor-1));
  
  m_activity = m_activity * factor;
  m_Nlyso = mLYSOprime*m_time*60;
  m_Nsig = mSigprime*m_time*60;
}

void Result::RescaleTime(double factor) {
  m_Nlyso *= factor;
  m_Nsig *= factor;
  m_time *= factor;
}

void Result::Restore()
{
  m_Nlyso = m_NlysoOrig;
  m_NlysoErr = m_NlysoErrOrig;
  m_Nsig = m_NsigOrig;
  m_time = m_timeOrig;
  m_activity = m_activityOrig;
}

void SetGraphStyle(TGraph* g, int markerStyle, int markerSize, int color) {
	g->SetMarkerStyle(markerStyle);
	g->SetMarkerSize(markerSize);
	g->SetMarkerColor(color);
}

std::pair<TGraph*, TGraph*> Result::MakeGraphAroundAlpha(double alpha) 
{
  int n = 7;
  TGraph* gAnal = new TGraph(n);
  TGraph* gOTH = new TGraph(n);
  
  for(int i = 0; i < n; i++) {
    Restore();
    double alphaprime = 0;
    if(i <= n/2) {
      alphaprime = alpha/(n/2.-i+1);
    } else {
	alphaprime = (i-n/2.+1)*alpha;    
    }
    RescaleActivity(alphaprime);
//     Print();
    ApplyEff(eff_signal, eff_lyso);
    std::pair<double, double> z = CalcZ("OTHinput/inputYield.dat");
    cout << "zAnal, zOTH = " << z.first << " " << z.second << endl;
    gAnal->SetPoint(i, alphaprime, z.first);
    gOTH->SetPoint(i, alphaprime, z.second);
  
  }

  SetGraphStyle(gAnal, 8, 2, kBlack);
  SetGraphStyle(gOTH, 4, 1, kRed);
return std::make_pair(gAnal, gOTH);
  
}

TF1* Result::MakeFuncZvsAlpha()
{
  TF1* f = new TF1("f", "x*[0]/sqrt([1])*1/sqrt(1+[2]*[3]*(x-1))", 0, 250);
  f->SetParameters(m_Nsig*eff_signal, m_Nlyso*eff_lyso, SigRate(), m_deadTime);
  return f;
}

void FitLYSOPlusSig(string dataFile, string lysoFile, bool na22FromSimu=false)
{
  RooRealVar* E = new RooRealVar("E", "Energy", 0, 1200, "keV");
  E->setBins(150);
  E->setRange("range_whole", E->getMin(), E->getMax());
  E->setRange("range_250_Max", 250, E->getMax());
  E->setRange("range_300_Max", 300, E->getMax());
  E->setRange("range_400_Max", 400, E->getMax());
  E->setRange("range_450_Max", 450, E->getMax());
  E->setRange("range_650_Max", 650, E->getMax());
  E->setRange("range_650_900", 650, 900);
  double signalWindow_min = 450;
  double signalWindow_max = 570;
//   double signalWindow_min = 900;
//   double signalWindow_max = 1100;
  E->setRange("range_sigWindow", signalWindow_min, signalWindow_max) ;
  
  TFile* f_LYSO = new TFile(lysoFile.c_str());
  TTree* t_LYSO = (TTree*) f_LYSO->Get("tree");
  
  TFile* fData = new TFile(dataFile.c_str());
  TTree* tData = (TTree*) fData->Get("tree");
 
  Data* dataLYSO = new Data(t_LYSO, E);
  Data* data = new Data(tData, E);
  data->m_activity = 16e3*0.90; // this 0.9 is the beta+ intensity (0.1 is going to EC)
  
  dataLYSO->m_cut = "";
  data->m_cut = "Evt < 15000 && Sat[0] == 0";
  
  RooAbsPdf* sig_gaussian;
  if(na22FromSimu) {
    ////////////////////////////////////////////////////////////////////////////////
    // Na22 model from G4
    TH1F* h = new TH1F("h", "h", 150, 0, 1200);
    TFile* fNa22Simu = new TFile("../../Data/G4/MakeCrystalsFromHits_Na22_Plot4_source.root", "read");
    TTree* treeNa22Simu = (TTree*) fNa22Simu->Get("tree");
    TTreeReader reader(treeNa22Simu);
    TTreeReaderArray<Double_t> Edep(reader, "EdepSmeared_Constant");
    while (reader.Next()) {
      if(Edep[0]>0) {// && Edep[0] < 0.600){
       h->Fill(Edep[0]*1e3); 
      }
    }
    h->Draw();
  
    RooDataHist* hist_Na22 = new RooDataHist("dhE_sig", "dhE_sig", *E, Import(*h));
    sig_gaussian = new RooHistPdf("sig_gaussian","sig_gaussian", *E,*hist_Na22,0);
  } else {
    RooRealVar* sig_peak_mean = new RooRealVar("sig_peak_mean", "mean of gaussian for signal peak", 507, "keV");
    RooRealVar* sig_peak_sigma = new RooRealVar("sig_peak_sigma", "width of gaussian for signal peak", 32, "keV");

    sig_gaussian = new RooGaussian("sig_gaussian", "gaussian for signal peak", *E, *sig_peak_mean, *sig_peak_sigma);
  }
  
  Model* model = new Model(E, dataLYSO, sig_gaussian);
  Result* res = model->Fit(data);
  res->Print();
  TCanvas* c1 = new TCanvas();
  model->Plot(data, res);
  
  double alpha = res->SolveForAlpha(3*sqrt(eff_lyso)/eff_signal);

  TF1* f = res->MakeFuncZvsAlpha();
  
  std::pair<TGraph*, TGraph*> graphs = res->MakeGraphAroundAlpha(alpha);
  TMultiGraph* multi = new TMultiGraph();
  multi->Add(graphs.first);
  multi->Add(graphs.second);
  multi->Draw("ap");
  f->Draw("same");
  /*
  int Npoints = 5;
  TGraph* gZanalVsTime = new TGraph(Npoints);
  SetGraphStyle(gZanalVsTime, 8, 2, kBlack);
  TGraph* gZothVsTime = new TGraph(Npoints);
  SetGraphStyle(gZothVsTime, 21, 2, kRed);
  TGraph* gZanalVsTime_act1 = new TGraph(Npoints);
  SetGraphStyle(gZanalVsTime_act1, 8, 1, kGreen+2);
  TGraph* gZothVsTime_act1 = new TGraph(Npoints);
  SetGraphStyle(gZothVsTime_act1, 21, 1, kMagenta);
  for(int i = 0; i < Npoints; i++) {
    double timeFactor = (i+1)/1000.; 
    res->RescaleTime(timeFactor);
    res->Print();
    zAnal = res->CalcZanal();
    res->WriteOTHFile("OTHinput/inputYieldRescaled.dat");
    zOTH = res->CalcZoth("OTHinput/inputYieldRescaled.dat");
    cout << "i, Zanal, Zoth = " << i << " " << zAnal << " " << zOTH << endl;
    gZanalVsTime->SetPoint(i, res->m_time*60, zAnal);
    gZothVsTime->SetPoint(i, res->m_time*60, zOTH);
    
    res->RescaleActivityAndTime(0.8, timeFactor);
    res->Print();
    zAnal = res->CalcZanal();
    res->WriteOTHFile("OTHinput/inputYieldRescaled.dat");
    zOTH = res->CalcZoth("OTHinput/inputYieldRescaled.dat");
    cout << "i, Zanal, Zoth = " << i << " " << zAnal << " " << zOTH << endl;
    gZanalVsTime_act1->SetPoint(i, res->m_time*60, zAnal);
    gZothVsTime_act1->SetPoint(i, res->m_time*60, zOTH);
  }
  
  TMultiGraph* gMulti = new TMultiGraph();
  gMulti->Add(gZanalVsTime);
  gMulti->Add(gZothVsTime);
  gMulti->Add(gZanalVsTime_act1);
  gMulti->Add(gZothVsTime_act1);
  
  TCanvas* cMulti = new TCanvas();
  gMulti->Draw("ap");
  */
  /*
  res->RescaleActivity(0.05);
  
  res->Print();
  zAnal = res->CalcZanal();
  res->WriteOTHFile("OTHinput/inputYieldRescaled.dat");
  zOTH = res->CalcZoth("OTHinput/inputYieldRescaled.dat");
  cout << "zAnal, zOTH = " << zAnal << " " << zOTH << endl;
  TCanvas* c2 = new TCanvas();
  model->Plot(data, res, false);
  */
}

void FitLYSOPlusSig()
{
  FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "~/godaq_rootfiles/analysis_v2.10.0/run78.root", false);
}



// OTH stuff
//   
//   TH1F* hLYSO_gen = (TH1F*) m_LYSO.m_dh->createHistogram("E");
//   hLYSO_gen->Sumw2();
//   
//   //hLYSO_gen->Scale(lyso_yield->getVal()/hLYSO_gen->Integral());
//   hLYSO_gen->Scale(N1prime/hLYSO_gen->Integral());
//   
//   cout << "ERROR=" << hLYSO_gen->GetEntries() << "  " << hLYSO_gen->Integral() << "  " << hLYSO_gen->GetBinContent(hLYSO_gen->GetXaxis()->FindBin(400)) << " " << 
// hLYSO_gen->GetBinError(hLYSO_gen->GetXaxis()->FindBin(400)) << endl;
//   
//   //RooDataHist* dhSig_gen = sig_gaussian->generateBinned(*E, sig_yield->getVal(), ExpectedData());
//   RooDataHist* dhSig_gen = sig_gaussian->generateBinned(*E, N2prime, ExpectedData());
//   
//   TH1F* hSig_gen = (TH1F*) dhSig_gen->createHistogram("E");
//   cout << "hSig_gen: " << hSig_gen->Integral() << endl;
//     
//   TH1F* hLYSO_gen_range = MakeHistoInRange(hLYSO_gen, signalWindow_min, signalWindow_max, true);
//   TH1F* hSig_gen_range = MakeHistoInRange(hSig_gen, signalWindow_min, signalWindow_max, true);
//   TCanvas* c2 = new TCanvas("c2","c2",1200,400);
//   c2->Divide(3,1);
//   c2->cd(1);
//   hLYSO_gen->Draw("e");
//   hLYSO_gen_range->Draw("histsame");
//   hLYSO_gen_range->Draw("same");
//   c2->cd(2);
//   hSig_gen->Draw("e");
//   hSig_gen_range->Draw("histsame");
//   hSig_gen_range->Draw("same");
//   c2->cd(3);
//   TH1F* hData_gen_range = (TH1F*) hLYSO_gen_range->Clone("hData_gen_range");
//   hData_gen_range->Add(hSig_gen_range);
//   hData_gen_range->GetYaxis()->SetRangeUser(0, hData_gen_range->GetMaximum()*1.3);
//   hData_gen_range->Draw("e");
//   hLYSO_gen_range->Draw("same");
//   hSig_gen_range->Draw("same");
  
//   TFile* f = new TFile("OTHinput/histos.root", "recreate");
//   hLYSO_gen_range->Write();
//   hSig_gen_range->Write();
//   hData_gen_range->Write();
//   f->Write();
//   f->Close();