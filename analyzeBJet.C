#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

using namespace std;

void analyzeBJet()
{

	double ptHatMin = 50.0;
	double ptHatMax = -1;

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
  double bbPx, bbPy, bbPz, bbE; 
  int bbevent, bbPID;
  TTree *bbtree = (TTree*)infile->Get("bbtree");
  bbtree->SetBranchAddress("bbevent", &bbevent);
  bbtree->SetBranchAddress("bbPx", &bbPx);
  bbtree->SetBranchAddress("bbPy", &bbPy);
  bbtree->SetBranchAddress("bbPz", &bbPz);
  bbtree->SetBranchAddress("bbE" , &bbE);
  bbtree->SetBranchAddress("bbPID",&bbPID);


  TH1D *rdiffbhist = new TH1D("rdiffbhist",";r;counts",100,0.,1.);
  TH1D *rbhist = new TH1D("rbhist",";r;counts",100,0.,10.);
  TH1D *zbhist = new TH1D("zbhist",";z;counts",100,0.,1.);
  TH1D *jtbhist = new TH1D("jtbhist",";j_{T} [GeV];counts",100,0.,10.);

  
  for (int i=0; i<bbtree->GetEntries(); ++i)
  {
    bbtree->GetEntry(i);
    TLorentzVector b4vec = TLorentzVector(bbPx, bbPy, bbPz, bbE);
	if (b4vec.Eta() < 2. || b4vec.Eta() > 5.)  {continue;}
	TVector3 b3vec(b4vec.Px(),b4vec.Py(),b4vec.Pz());
    for (int j=0; j<jettree->GetEntries(); ++j)
    {
      jettree->GetEntry(j);
      TLorentzVector jet4Vec = TLorentzVector(jPx, jPy, jPz, jE); 
	  if (jet4Vec.Eta() < 2. || jet4Vec.Eta() > 5.)  {continue;}
	  TVector3 jet3vec(jet4Vec.Px(),jet4Vec.Py(),jet4Vec.Pz());

      if (bbevent != eventNum) {continue;}
	  double r = sqrt(pow(jet4Vec.Eta() - b4vec.Eta(), 2));
	  rbhist->Fill(r);
      // Geometric matching goes here... // 
      double rdiffb = sqrt(pow(jet4Vec.Eta() - b4vec.Eta(), 2) + pow(jet4Vec.Phi() - b4vec.Phi(), 2) );
      rdiffbhist->Fill(rdiffb);
      if (rdiffb > 0.5){continue;}
      
        double z = (b3vec.Dot(jet3vec))/(jet3vec.Mag2());
        double jt = ((b3vec.Cross(jet3vec)).Mag())/ (jet3vec.Mag());
		zbhist->Fill(z);
		jtbhist->Fill(jt);

    }
	}
  outfile->Write();
  outfile->Close();
}

