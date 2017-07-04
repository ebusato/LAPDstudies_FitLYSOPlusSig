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
// double eff_signal = 0.85;
// double eff_lyso = 0.0034;

// double eff_signal = 1;
// double eff_lyso = 1;

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

// returns run duration in seconds
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
  return timeEnd - timeBeg;
}


Model::Model(RooRealVar* E, Data* lyso, RooAbsPdf* sig_gaussian) : m_E(E), 
								   m_LYSO(lyso)
{
  m_lyso_yield = new RooRealVar("lyso_yield", "yield of lyso", 100000, 0, 1000000);
  m_sig_yield = new RooRealVar("sig_yield", "yield signal peak", 1000, 0, 10000000);
  
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
  
//   RooMsgService::instance().setSilentMode(true);


  m_sig_yield->setConstant();
  m_model->fitTo(*(data->m_dh), Extended(),Range("range_650_Max"), PrintEvalErrors(-1));

//  m_model->fitTo(*(data->m_dh), Extended(),Range("range_700_Max"), PrintEvalErrors(-1));

 Result* res = new Result(data->RunDuration(), data->m_activity, GetSigYield(data), m_lyso_yield->getVal(), m_LYSO->m_dh->sum(kFALSE));
 return res;
}

void Model::Plot(Data* data, Result* res, bool plotData)
{m_sig_yield->setConstant(0);
  m_sig_yield->setVal(peakEff*res->m_Nsig);
  m_lyso_yield->setVal(res->m_Nlyso);
  
  TGraph* data_graph, *sig_graph, *lyso_graph;
  
  RooPlot* frame = m_E->frame(Bins(100), Range("range_200_1000"));
  frame->GetYaxis()->SetTitle("Events");
  frame->GetXaxis()->SetTitle("E [keV]");
  if(plotData) {
    data->m_dh->plotOn(frame); //, DrawOption("PX"));
    data_graph = (TGraph*)frame->getObject( frame->numItems() - 1  );
  }
  
  m_model->plotOn(frame, Range("range_250_Max"));
  m_model->plotOn(frame, Range("range_250_Max"), Components("sig_gaussian"),LineColor(kRed));
  sig_graph = (TGraph*)frame->getObject( frame->numItems() - 1  );
  m_model->plotOn(frame, Range("range_200_Max"), Components("histpdf_LYSO"),LineColor(kGreen+2));
  lyso_graph = (TGraph*)frame->getObject( frame->numItems() - 1  );
  if(plotData) {
    data->m_dh->plotOn(frame); //, DrawOption("PX"));
  }
  frame->Draw();
  	frame->GetXaxis()->SetTitleSize(0.05);
	frame->GetYaxis()->SetTitleSize(0.05);
	frame->GetXaxis()->SetTitleOffset(1.25);
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->GetXaxis()->SetLabelSize(0.05);
	frame->GetYaxis()->SetLabelSize(0.05);
  double xText = 0.55;
  double yShift = 0.07;
  PutText(xText, 0.81, kBlack, "LAPD");
 // PutText(xText, 0.85-yShift, kBlack, "LPC");
  PutText(xText, 0.81-1*yShift, kBlack, "^{22}Na (14.4 kBq)");
  stringstream ss;
  ss.precision(3);
  ss << "Run duration: " << res->m_time << " min";
  PutText(xText, 0.81-2*yShift, kBlack, ss.str().c_str());
  
  TLegend* leg = new TLegend(0.6,0.42,0.75,0.62);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	if(data_graph) leg->AddEntry(data_graph, "Data", "lp");
	if(sig_graph) leg->AddEntry(sig_graph, "Signal photopeak", "l");
	if(lyso_graph) leg->AddEntry(lyso_graph, "LYSO background", "l");
// 	leg->AddEntry(h2.first, "z = 13.5 cm", "l");
// 	leg->AddEntry(h3.first, "z = 14 cm", "l");
	leg->Draw();
}

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
	cout << "  -> Time = " << m_time << " sec" << endl;
	cout << "  -> Nsig = " << m_Nsig << endl;
	cout << "  -> Nlyso = " << m_Nlyso << " +- " << m_NlysoErr << endl;
	cout << "  -> NlysoOrig = " << m_NlysoOrig << endl;
	//double rateSig = m_Nsig / m_time;
	//double rateLyso = m_Nlyso / m_time;
	cout << "  -> rate sig = " << SigRate() << endl;
	cout << "  -> rate lyso = " << LysoRate() << endl;
	cout << "  -> rate tot = " << SigRate()+LysoRate() << endl; 
}

void Result::ApplyEff(double effSig, double effLyso)
{
  m_Nsig *= effSig;
  m_Nlyso *= effLyso;
}

// time is in seconds
void Result::PrintYields(double activity, double time) {
  double mLYSOOrig = m_NlysoOrig/m_timeOrig;
  double mSigOrig = m_NsigOrig/m_timeOrig;

  double mLYSOprime = mLYSOOrig / (1+mSigOrig*m_deadTime*(activity/m_activityOrig-1));
  double mSigprime = activity/m_activityOrig*mSigOrig / (1+mSigOrig*m_deadTime*(activity/m_activityOrig-1));

  cout << "PrintYields result: " << endl;
  cout << "  -> activity=" << activity << endl;
  cout << "  -> time=" << time << endl;
  cout << "  -> mLYSOprime=" << mLYSOprime << endl;
  cout << "  -> mSigprime=" << mSigprime << endl;
  cout << "  -> NLYSOprime=" << mLYSOprime*time << endl;
  cout << "  -> NSigprime=" << mSigprime*time << endl;
}

void Result::RescaleActivity(double factor) {
  // compute rates in sec-1
  double mLYSO = m_Nlyso / m_time;
  double mSig = m_Nsig / m_time;
  
  // compute rescaled quantities
  double mLYSOprime = mLYSO / (1+mSig*m_deadTime*(factor-1));
  double mSigprime = factor*mSig / (1+mSig*m_deadTime*(factor-1));
  
  m_activity = m_activity * factor;
  m_Nlyso = mLYSOprime*m_time;
  m_Nsig = mSigprime*m_time;
}

void Result::RescaleTime(double factor) {
  m_Nlyso = m_Nlyso * factor;
  m_Nsig = m_Nsig * factor;
  m_time = m_time * factor;
}

void Result::SetTime(double time) {
  double factor = time/m_time;
  RescaleTime(factor);
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

TF1* Result::MakeFuncSignifTheoFromGamma(TString name, int color, int style, double eff_signal, double eff_lyso, double time)
{
  double max=10000.;
  if(time > 4)
    max=3000.;
  if(time > 10)
    //max=500.;
    max=1000.;
  if(time >50)
    //max=250.;
    max=400.;
  if(time > 60*8)
    //max=63.; 
    max=80;

  TF1* f = new TF1(name.Data(), "ROOT::Math::normal_quantile(1-ROOT::Math::gamma_cdf([1]/(1+[2]*[3]*(x/[4] - 1))*[5]*[6],[1]/(1+[2]*[3]*(x/[4] - 1))*[5]*[6] + x/[4]*[2]/(1+[2]*[3]*(x/[4] - 1))*[5]*[7], 1), [0])", 1, max);

  // [0] -> with of normal distribution used for quantile calculation
  // [1] -> mb_0
  // [2] -> ms_0
  // [3] -> dead time
  // [4] -> initial activity
  // [5] -> acquisition time
  // [6] -> eff lyso
  // [7] -> eff sig
  f->SetParameters(1, m_NlysoOrig/m_timeOrig, m_NsigOrig/m_timeOrig, m_deadTime, m_activityOrig, time, eff_lyso, eff_signal);
  f->SetLineColor(color);
  f->SetLineStyle(style);
  f->GetXaxis()->SetTitleSize(0.05);
  f->GetYaxis()->SetTitleSize(0.05);
  f->GetXaxis()->SetTitleOffset(1.25);
  f->GetYaxis()->SetTitleOffset(1.2);
  f->GetXaxis()->SetLabelSize(0.05);
  f->GetYaxis()->SetLabelSize(0.05);
  //f->SetLineWidth(3);
  //f->SetNpx(1e4);
  return f;
}


TF1* Result::MakeFuncSignifTheoFromGammaFromTrueRates(TString name, int color, int style, double eff_signal, double eff_lyso, double time, double deadtime)
{
  double max=10000.;
  if(time > 4)
    max=3000.;
  if(time > 10)
    //max=500.;
    max=1000.;
  if(time >50)
    //max=250.;
    max=400.;
  if(time > 60*8)
    //max=63.; 
    max=80;

  TF1* f = new TF1(name.Data(), "ROOT::Math::normal_quantile(1-ROOT::Math::gamma_cdf([1]/(1+[3]*([1] + [2]*x/[4]))*[5]*[6],[1]/(1+[3]*([1] + [2]*x/[4]))*[5]*[6] + x/[4]*[2]/(1+[3]*([1] + [2]*x/[4]))*[5]*[7], 1), [0])", 1, max);

  // [0] -> with of normal distribution used for quantile calculation
  // [1] -> rb
  // [2] -> rs_0
  // [3] -> dead time
  // [4] -> initial activity
  // [5] -> acquisition time
  // [6] -> eff lyso
  // [7] -> eff sig

  double deadT = 0;
  if(deadtime==-1) {
    deadT = m_deadTime;
  } else {
    deadT = deadtime;
  }

  f->SetParameters(1, BkgTrueRateOrig(), SigTrueRate0Orig(), deadT, m_activityOrig, time, eff_lyso, eff_signal);
  f->SetLineColor(color);
  f->SetLineStyle(style);
  f->GetXaxis()->SetTitleSize(0.05);
  f->GetYaxis()->SetTitleSize(0.05);
  f->GetXaxis()->SetTitleOffset(1.25);
  f->GetYaxis()->SetTitleOffset(1.2);
  f->GetXaxis()->SetLabelSize(0.05);
  f->GetYaxis()->SetLabelSize(0.05);
  //f->SetLineWidth(3);
  //f->SetNpx(1e4);
  return f;
}


double FindDetectionLimit(TF1* f) 
{
  double a=1.;
  double Z = f->Eval(a);
  while(fabs(Z-3) > 0.001) {
    if(Z < 3) {
      a*=1.2;
    } else {
      a/=1.1;
    }
    Z = f->Eval(a);
    //cout << "Testing a=" << a << " -> Z=" << Z << endl;
  }
  return a;
}

std::pair<double, double> ComputeEfficienciesInSignalRegion(Data* dataNa22PlusLYSO, Data* dataLYSO, Result* res)
{
	TCut signalRegion("NoPulses == 2 && E[0] > 425 && E[0] < 595 && E[1] > 425 && E[1] < 595 && fabs(T30[0] - T30[1]) < 3.6");
	
	// LYSO
	int nLYSO_init = dataLYSO->m_tree->GetEntries();
	int nLYSO_afterCut = dataLYSO->m_tree->GetEntries(signalRegion);
	cout << nLYSO_init << "  "  << nLYSO_afterCut << endl;
	double effLYSO = nLYSO_afterCut/float(nLYSO_init);
	
	// Signal
	// The efficiency for the signal is given by the following calculations:
	// nNa22PlusLYSO_init = Ns + Nb
	// nNa22PlusLYSO_afterCut = eff_signal*Ns + eff_lyso*Nb
	//  => eff_signal = (nNa22PlusLYSO_afterCut - eff_lyso*Nb)/Ns
	// where Nb and Ns are given by the fit whose result is stored in res
	int nNa22PlusLYSO_init = dataNa22PlusLYSO->m_tree->GetEntries(dataNa22PlusLYSO->m_cut);
	int nNa22PlusLYSO_afterCut = dataNa22PlusLYSO->m_tree->GetEntries(dataNa22PlusLYSO->m_cut && signalRegion);
	cout << nNa22PlusLYSO_init << "  "  << nNa22PlusLYSO_afterCut << endl;
	
	double Nb = res->m_Nlyso;
	double Ns = res->m_Nsig;
	cout << "Nb, Ns=" << Nb << "  " << Ns << endl;
	
	double effSignal = (nNa22PlusLYSO_afterCut - effLYSO*Nb)/Ns;
	cout << effSignal << endl;
	
	return make_pair(effSignal, effLYSO);
}

void FitLYSOPlusSig(string dataFile, string lysoFile, bool na22FromSimu=false)
{
  RooRealVar* E = new RooRealVar("E", "Energy", 0, 1200, "keV");
  E->setBins(150);
  E->setRange("range_whole", E->getMin(), E->getMax());
  E->setRange("range_200_1000", 200, 1000);
  E->setRange("range_250_Max", 250, E->getMax());
  E->setRange("range_300_Max", 300, E->getMax());
  E->setRange("range_400_Max", 400, E->getMax());
  E->setRange("range_450_Max", 450, E->getMax());
  E->setRange("range_650_Max", 650, E->getMax());
  E->setRange("range_700_Max", 700, E->getMax());
  E->setRange("range_750_Max", 750, E->getMax());
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
  //data->m_cut = "Evt < 1600 && Sat[0] == 0";
  
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
  cout << "HERE HERE" << endl;
  res->Print();
  
  TCanvas* c1 = new TCanvas("c1", "c1");
  model->Plot(data, res);
  c1->SaveAs("FitLYSOPlusSig_c1.png");
  
  TCanvas* c10 = new TCanvas("c10", "c10");
  c10->SetGridx(1);
  c10->SetGridy(1);
  TF1* fSigRateVsAlpha = new TF1("fSigRateVsAlpha", "x*[0]/(1+[0]*[1]*(x-1))",0,10);
 fSigRateVsAlpha->SetParameters(res->SigRate(), res->m_deadTime);
 TF1* fBkgRateVsAlpha = new TF1("fBkgRateVsAlpha", "[2]/(1+[0]*[1]*(x-1))",0,10);
 fBkgRateVsAlpha->SetParameters(res->SigRate(), res->m_deadTime, res->LysoRate());
 fSigRateVsAlpha->SetLineColor(kRed);
 fBkgRateVsAlpha->SetLineColor(kGreen+2);
 TF1* fSum = new TF1("fSum", "x*[0]/(1+[0]*[1]*(x-1)) + [2]/(1+[0]*[1]*(x-1))", 0, 10);
 fSum->SetParameters(res->SigRate(), res->m_deadTime, res->LysoRate());
 fSum->SetLineColor(kBlack);

 fSum->GetXaxis()->SetTitleSize(0.05);
 fSum->GetYaxis()->SetTitleSize(0.05);
 fSum->GetXaxis()->SetTitleOffset(1.25);
 fSum->GetYaxis()->SetTitleOffset(1.2);
 fSum->GetXaxis()->SetLabelSize(0.05);
 fSum->GetYaxis()->SetLabelSize(0.05);
 fSum->GetYaxis()->SetRangeUser(0,30);
 fSum->GetXaxis()->SetTitle("#alpha");
 fSum->GetYaxis()->SetTitle("Expected measured rates [Hz]");

 fSum->SetLineWidth(4);
 fSigRateVsAlpha->SetLineWidth(3);
 fSigRateVsAlpha->SetLineStyle(9);
 fBkgRateVsAlpha->SetLineWidth(3);
 fBkgRateVsAlpha->SetLineStyle(2);

 fSum->Draw();
 fSigRateVsAlpha->Draw("same");
 fBkgRateVsAlpha->Draw("same");

 PutText(0.2, 0.81, kBlack, "LAPD");
 PutText(0.2, 0.75, kBlack, "^{22}Na source");

  TLegend* leg = new TLegend(0.6,0.35,0.88,0.55);
  leg->SetBorderSize(1);
  leg->SetTextSize(0.04);	
  leg->AddEntry(fSum, "Total rate", "l");
  leg->AddEntry(fSigRateVsAlpha, "Signal rate", "l");
  leg->AddEntry(fBkgRateVsAlpha, "Background rate", "l");
  leg->Draw();
  c10->SaveAs("FitLYSOPlusSig_c10.png");

 
  cout << "--------------------------------------" << endl;
  cout << "Computing efficiencies" << endl;
  std::pair<double, double> effs = ComputeEfficienciesInSignalRegion(data, dataLYSO, res);
  cout << effs.first << "  " << effs.second << endl;
  
  double eff_signal = effs.first;
  double eff_lyso = effs.second;
 
//   double alpha = res->SolveForAlpha(3*sqrt(eff_lyso)/eff_signal);
  
  
  TCanvas* c2 = new TCanvas("c2", "c2");
  // draw 3 sigma line
  TF1* line3sigmas = new TF1("line3sigmas","3",10, 10000);
  line3sigmas->Draw();

  c2->SetLogx();
  //std::vector<double> detectionLimits = MakeZVsAlphaPlot(res, model, data, eff_signal, eff_lyso);

  gPad->SetGridx(1);
  gPad->SetGridy(1);
  
  std::vector<TF1*> funcs;
  std::vector<TF1*> funcsdT1;
  std::vector<TF1*> funcsdT2;

  std::vector<std::pair<int, int> > styles; // first element of pair: color, second element of pair: style
  styles.push_back(make_pair(kBlack, 9));
  styles.push_back(make_pair(kRed, 7));
  styles.push_back(make_pair(kGreen+2, 2));
  styles.push_back(make_pair(kBlue, 10));
  styles.push_back(make_pair(kMagenta, 10));
  
  std::vector<double> Times;
  Times.push_back(res->m_time);
  Times.push_back(2*60);
  Times.push_back(20.);
  Times.push_back(5.);
  Times.push_back(1.);
  
  std::vector<double> detectionLimits;

  for(int i=0; i<5; i++) {
	cout << "--------------------" << endl;
	cout << "RescaleTime " << i << endl;
	
	res->SetTime(Times[i]); res->Print();
	TString name("f");
	name+=i;
	TF1* f = res->MakeFuncSignifTheoFromGammaFromTrueRates(name, styles[i].first, 1, eff_signal, eff_lyso, Times[i]);
	funcs.push_back(f);

	TF1* f1 = res->MakeFuncSignifTheoFromGammaFromTrueRates(name, styles[i].first, 1, eff_signal, eff_lyso, Times[i]);
	funcsdT1.push_back(f1);

	double a0 = FindDetectionLimit(funcs[i]);
	cout << " Detection limit = " << a0 << endl;
	res->PrintYields(a0, Times[i]);
	detectionLimits.push_back(a0);
  }
  
  for(int i=0; i<5; i++) {
    funcs[i]->Draw("same");
  }
  
  TLatex l;
  l.SetTextColor(kBlack);
  l.SetTextSize(0.045);
  PutText(0.2, 0.81, kBlack, "LAPD");
  PutText(0.2, 0.81-0.071, kBlack, "^{22}Na source at (0,0,0)");
  l.SetTextAngle(50); 
  l.DrawLatex(24, 1.5, Form("%.1f min", Times[0]/60.));
  l.SetTextColor(kRed); l.DrawLatex(60, 1.5, Form("%.1f min", Times[1]/60.));
  l.SetTextColor(kGreen+2); l.DrawLatex(180, 1.5, Form("%.1f sec", Times[2]));
  l.SetTextColor(kBlue); l.DrawLatex(455, 1.5, Form("%.1f sec", Times[3]));
  l.SetTextColor(kMagenta); l.DrawLatex(1500, 1.5, Form("%.1f sec", Times[4]));
  
  //  line3sigmas->GetYaxis()->SetRangeUser(0.5,4.2);
  line3sigmas->GetXaxis()->SetRangeUser(10,8000);
  line3sigmas->SetLineWidth(3);
  line3sigmas->SetLineColor(kBlack);
  line3sigmas->GetXaxis()->SetTitle("activity [Bq]");
  line3sigmas->GetYaxis()->SetTitle("significance (expected average)");
  line3sigmas->GetXaxis()->SetTitleSize(0.05);
  line3sigmas->GetYaxis()->SetTitleSize(0.05);
  line3sigmas->GetXaxis()->SetTitleOffset(1.25);
  line3sigmas->GetYaxis()->SetTitleOffset(1.2);
  line3sigmas->GetXaxis()->SetLabelSize(0.05);
  line3sigmas->GetYaxis()->SetLabelSize(0.05);
  c2->Update();
  c2->SaveAs("FitLYSOPlusSig_c2.png");

  TCanvas* c3 = new TCanvas("c3", "c3");
  TGraph* g = new TGraph(Times.size());
  for(int i=0; i<Times.size(); i++) {
    g->SetPoint(i, Times[i], detectionLimits[i]);
  }
  
  g->Draw("apl");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  g->GetYaxis()->SetTitle("Detection limit [Bq]");
  g->GetXaxis()->SetTitle("Time [s]");
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetTitleOffset(1.25);
  g->GetYaxis()->SetTitleOffset(1.4);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetRangeUser(0.5, 2000);
  g->GetYaxis()->SetRangeUser(10, 5000);
  gPad->Update();
}

void FitLYSOPlusSigNoOTH()
{
//   FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "~/godaq_rootfiles/analysis_v2.10.0/run78.root", false);
FitLYSOPlusSig("analysis_v2.18-calibC1/run67.root", "analysis_v2.18-calibC1/run78.root", false);
// 	FitLYSOPlusSig("analysis_v2.18-calibC1/run34.root", "analysis_v2.18-calibC1/run78.root", false);
}


