#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

using namespace std;

void analyzeKpKmCorr()
{

	double ptHatMin = 20.0;
	double ptHatMax = 50.0;

	TString inFileName = "data/pythiaJets_";
	inFileName += ptHatMin;
	inFileName += "_";
	inFileName += ptHatMax;
	inFileName += ".root";
	
	TString outFileName = "data/BjetHistos_";
	outFileName += ptHatMin;
	outFileName += "_";
	outFileName += ptHatMax;
	outFileName += ".root";	
	cout << "input file name : " << inFileName << endl;
  TFile *infile = TFile::Open(inFileName);
  TFile *outfile = new TFile(outFileName,"RECREATE");


  // Grab jet tree and relevant branches from the input file //
  int    eventNum, jetNum, jetNumInEvent, nConstituents;
  bool kpkmJet, kpOrkmJet;
  double jPx, jPy, jPz, jE, jPt, jEta, jPhi;
  TTree *jettree = (TTree*)infile->Get("jets");
  jettree->SetBranchAddress("eventNum", &eventNum);
  jettree->SetBranchAddress("jetNum", &jetNum);
  jettree->SetBranchAddress("jetNumInEvent", &jetNumInEvent);
  jettree->SetBranchAddress("nConstituents", &nConstituents);
  jettree->SetBranchAddress("jPx", &jPx);
  jettree->SetBranchAddress("jPy", &jPy);
  jettree->SetBranchAddress("jPz", &jPz);
  jettree->SetBranchAddress("jE" , &jE);
  
  // Make sure the naming conventions match your pythiaJet.C macro!! // 
  double bPx, bPy, bPz, bE;
  TTree *btree = (TTree*)infile->Get("bs");
  btree->SetBranchAddress("eventNum", &eventNum);
  btree->SetBranchAddress("bPx", &bPx);
  btree->SetBranchAddress("bPy", &bPy);
  btree->SetBranchAddress("bPz", &bPz);
  btree->SetBranchAddress("bE" , &bE);

  
  for (int i=0; i<btree->GetEntries(); ++i)
  {
    btree->GetEntry(i);
    TLorentzVector b4vec = TLorentzVector(bPx, bPy, bPz, bE);
    int bEvent = eventNum;
    for (int j=0; j<jettree->GetEntries(); ++j)
    {
      jettree->GetEntry(j);
      TLorentzVector jet4Vec = TLorentzVector(jPx, jPy, jPz, jE); 
      int jetEvent = eventNum;
      if (bEvent != jetEvent) {continue;}
      // Geometric matching goes here... // 
      
    }
	}
  outfile->Write();
  outfile->Close();
}
