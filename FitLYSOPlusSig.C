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
double eff_signal = 0.85;
double eff_lyso = 0.0034;

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
  
  WriteOTHFile(fileName.Data());
  OpTHyLiC oth(OTH::SystPolyexpo,OTH::StatLogN);
  oth.addChannel("ch1",fileName.Data());
  const int Nexp=5e5;
  std::pair<double, double> s = oth.significance(OTH::SignifExpectedMed,Nexp);
//    std::pair<double, double> s = oth.significance(OTH::SignifObserved,Nexp);
  const double p=s.first;
  const double z=s.second;
  
  /*
  TH1F* hLLRb=(TH1F*)oth.getHistoLLRb();
  hLLRb->GetXaxis()->SetTitleSize(0.05);
  hLLRb->GetXaxis()->SetTitleOffset(0.88);
  hLLRb->GetXaxis()->SetTitle("q_{ #mu}");
  hLLRb->GetYaxis()->SetTitle("");
  hLLRb->SetLineColor(kRed);
  TH1F* hLLRsb=(TH1F*)oth.getHistoLLRsb();
	 hLLRsb->SetLineColor(kBlue);
	 hLLRb->SetLineWidth(4);
  
  TCanvas *c1 = new TCanvas();
//   c1->SetLogy();
  hLLRb->Draw();
  hLLRsb->Draw("same");*/
  
  /*
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
  */
  
//   return std::make_pair(m_Nsig/sqrt(m_Nlyso), z);
  double zCowan = sqrt(2*((m_Nsig+m_Nlyso)*TMath::Log(1+m_Nsig/m_Nlyso) - m_Nsig));
  return std::make_pair(zCowan, z);
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
  m_Nlyso = m_Nlyso * factor;
  m_Nsig = m_Nsig * factor;
  m_time = m_time * factor;
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

TF1* Result::MakeFuncZvsActivity(int color, int style)
{
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Old stuff, could eventually be removed
  // [0]: sig yield
  // [1]: lyso yield
  // [2]: sig rate
  // [3]: dead time
  // [4]: activity ref
//   TF1* f = new TF1("f", "x/[4]*[0]/sqrt([1])*1/sqrt(1+[2]*[3]*(x/[4]-1))", 0, 1000);
//   f->SetParameters(m_Nsig*eff_signal, m_Nlyso*eff_lyso, SigRate(), m_deadTime, m_activity);
  //////////////////////////////////////////////////////////////////////////////////////////////////
	
  // [0]: sig rate
  // [1]: lyso rate
  // [2]: dead time
  // [3]: activity ref
  // [4]: time
  // [5]: eff signal
  // [6]: eff background
//   TF1* f = new TF1("f", "x/[3]*[0]*1/sqrt([1]*(1+[0]*[2]*(x/[3]-1))) * sqrt([4]) * [5]/sqrt([6])", 0, 1000);
  // Cowan's formula
  TF1* f = new TF1("f", "sqrt(2*([4]/(1+[0]*[2]*(x/[3]-1))*([1]*[6] + x/[3]*[0]*[5])*TMath::Log(1+([5]*x/[3]*[0])/([6]*[1])) - [5]*x/[3]*[0]*[4]/(1+[0]*[2]*(x/[3]-1))))", 0, 2000);
  f->SetParameters(SigRate(), LysoRate(), m_deadTime, m_activity, m_time*60., eff_signal, eff_lyso);
  f->SetLineColor(color);
  f->SetLineStyle(style);
  f->SetNpx(1e4);
  return f;
}

TGraph* GetOTHTGraphAroundAlpha(Result* res, double alpha, int color, Model* model, Data* data, bool removeFirstPoint = false)
{ 

  TGraph* g = new TGraph(4);
  for(int i=0; i < 4; i++) {
   Result* resClone = res->Clone();
   double factor;
   if(i==0) factor = 1/2.;
   else if(i==1) factor = 1/1.5;
   else if(i==2) factor = 1/1.1;
   else if(i==3) factor = 1.4;
   else if(i==4) factor = 1.4;
   resClone->RescaleActivity(alpha*factor);
   
   resClone->Print();
//    TCanvas* cc = new TCanvas();
//   model->Plot(data, resClone);
   resClone->ApplyEff(eff_signal, eff_lyso);
   std::pair<double, double> zs = resClone->CalcZ("OTHinput/inputYield.dat");
   resClone->Print();
   cout << "z(anal) = " << zs.first << " ; z(oth) = " << zs.second << endl;
   g->SetPoint(i, res->m_activity*alpha*factor, zs.second);
  
  }
  if(removeFirstPoint) g->RemovePoint(0);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  return g;
}

void MakeZVsAlphaPlot(Result* res, Model* model, Data* data)
{
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  
  // raw
  double time1 = res->m_time;
  TF1* f1 = res->MakeFuncZvsActivity(kBlack, 9);
  f1->GetYaxis()->SetRangeUser(0, 5);
  
  
  TGraph* g1 = GetOTHTGraphAroundAlpha(res, 35/14400., kBlack, model, data);
  
  // rescale
  res->RescaleTime(1/3.); res->Print();
  double time2 = res->m_time;
  TF1* f2 = res->MakeFuncZvsActivity(kRed, 7);
  
   TGraph* g2 = GetOTHTGraphAroundAlpha(res, 67/14400., kRed, model, data);
  
  // rescale
  res->RescaleTime(1/3.); res->Print();
  double time3 = res->m_time;
  TF1* f3 = res->MakeFuncZvsActivity(kGreen+2, 2);
  
   TGraph* g3 = GetOTHTGraphAroundAlpha(res, 130/14400., kGreen+2, model, data);
  
  // rescale
  res->RescaleTime(1/3.); res->Print();
  double time4 = res->m_time;
  TF1* f4 = res->MakeFuncZvsActivity(kBlue, 10);
  
   TGraph* g4 = GetOTHTGraphAroundAlpha(res, 260/14400., kBlue, model, data);
 
    // rescale
  res->RescaleTime(1/2.); res->Print();
  double time5 = res->m_time;
  TF1* f5 = res->MakeFuncZvsActivity(kMagenta, 10);
  
   TGraph* g5 = GetOTHTGraphAroundAlpha(res, 400/14400., kMagenta, model, data, true);
 
   
  TMultiGraph* g = new TMultiGraph();
  g->Add(g1);
  g->Add(g2);
  g->Add(g3);
  g->Add(g4);
  g->Add(g5);
  
  g->Draw("apl");
 	g->GetXaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(1.25);
	g->GetYaxis()->SetTitleOffset(1.2);
	g->GetXaxis()->SetLabelSize(0.05);
	g->GetYaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetTitle("activity [Bq]");
  g->GetYaxis()->SetTitle("significance (expected median)");
  g->GetYaxis()->SetRangeUser(0.5,4.2);
  g->GetXaxis()->SetRangeUser(0.,600);
  TLatex l;
  l.SetTextColor(kBlack);
  l.SetTextSize(0.045);
    PutText(0.55, 0.31, kBlack, "LAPD");
  PutText(0.55, 0.31-0.071, kBlack, "^{22}Na source at (0,0,0)");
  double x, y;
  g1->GetPoint(0, x, y); l.DrawLatex(x, y-0.3, Form("time = %.1f min", time1));
  g2->GetPoint(0, x, y); l.SetTextColor(kRed); l.DrawLatex(x+10, y-0.05, Form("time = %.1f min", time2));
  g3->GetPoint(0, x, y); l.SetTextColor(kGreen+2); l.DrawLatex(x+10, y-0.1, Form("time = %.1f min", time3));
  g4->GetPoint(0, x, y); l.SetTextColor(kBlue); l.DrawLatex(x+20, y, Form("time = %.1f sec", time4*60));
  g5->GetPoint(0, x, y); l.SetTextColor(kMagenta); l.DrawLatex(x+40, y, Form("time = %.1f sec", time5*60));
  
//   f1->Draw("same");
//   f2->Draw("same");
//   f3->Draw("same");
//   f4->Draw("same");
//   f5->Draw("same");
  
  
  // draw 3 sigma line
  TLine* line3sigmas = new TLine(0, 3, 590, 3);
  line3sigmas->SetLineWidth(3);
  line3sigmas->Draw("same");
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
  TCanvas* c1 = new TCanvas("c1", "c1");
  model->Plot(data, res);
 c1->SaveAs("FitLYSOPlusSig_c1.png");
  
//   double alpha = res->SolveForAlpha(3*sqrt(eff_lyso)/eff_signal);
  
  TCanvas* c2 = new TCanvas("c2", "c2");
  MakeZVsAlphaPlot(res, model, data);
  c2->SaveAs("FitLYSOPlusSig_c2.png");
  
}

void FitLYSOPlusSig()
{
//   FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "~/godaq_rootfiles/analysis_v2.10.0/run78.root", false);
FitLYSOPlusSig("~/godaq_rootfiles/analysis_v2.10.0/run67.root", "analysis_v2.18-calibC1/run78.root", false);
// 	FitLYSOPlusSig("analysis_v2.18-calibC1/run34.root", "analysis_v2.18-calibC1/run78.root", false);
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