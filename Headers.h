
class Data {
public:
	Data(TTree* t, RooRealVar* E);
	RooDataHist* GetDataHistFromTH1(string TH1Name, string histName);
	double RunDuration();
// 	void RescaleActivity(double factor);
	
	TTree* m_tree;
	TCut m_cut;
	
	RooDataHist* m_dh;
	RooRealVar* m_E;
	
	double m_activity;
};

class Result {
public:
	Result(double time, double activity, double Nsig, double Nlyso, double NlysoOrig);
	void WriteOTHFile(TString fileName);
	double CalcZoth(TString fileName);
	double CalcZanal();
	void RescaleActivity(double factor);
	void RescaleTime(double factor);
	void MakeCalculationsSensitivity(Data* data);
	void Print();
	
	double m_time;
	double m_activity;
	double m_Nsig;
	double m_Nlyso;
	double m_NlysoOrigSample;
	double m_NlysoErr;
	
	double m_timeOrig;
	double m_activityOrig;
	double m_NsigOrig;
	double m_NlysoOrig;
	double m_NlysoErrOrig;
	
	double m_deadTime;
};

class Model {
public:
  Model(RooRealVar* E, Data* lyso, RooAbsPdf* sig_gaussian);
  void MakeCalculationsSensitivity(Data* data);
  Result* Fit(Data* data);
  void Plot(Data* data, Result* res, bool plotData = true);
  double GetSigYield(Data* data);
  
  RooRealVar* m_lyso_yield;
  RooRealVar* m_sig_yield;
  
  RooAddPdf* m_model;
  
  RooRealVar* m_E;
  
  Data* m_LYSO;
};
