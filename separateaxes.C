#include "TFile.h"
#include "TTree.h"
#include "TH1.h"



using namespace std;

void separateaxes()
{

	double ptHatMin = 50.0;
	double ptHatMax = -1;

	TString inFileName = "data/pythiaJets_";
	inFileName += ptHatMin;
	inFileName += "_";
	inFileName += ptHatMax;
	inFileName += ".root";
	
	TString outFileName = "data/jetHistos_";
	outFileName += ptHatMin;
	outFileName += "_";
	outFileName += ptHatMax;
	outFileName += "_new.root";	
	
  TFile *infile = TFile::Open(inFileName);
  TFile *outfile = new TFile(outFileName,"RECREATE");
  // z histos //
  TH1F *zH = new TH1F("z", ";z;", 100, 0., 1.);

  // jT histos //
  TH1F *jtH = new TH1F("jt", ";jt;", 120, 0., 60.);

  // r histos //
  TH1F *rH = new TH1F("r", ";r;", 100, 0., 1.);

  // z jt correlations //
  TH2D *zjt = new TH2D("zjt", ";z;j_{T} [GeV]", 100, 0., 1., 100, 0., 2. );
  TH2D *zjt_kpkmJet = new TH2D("zjt_kpkmJet", ";z;j_{T} [GeV]", 100, 0., 1., 100, 0., 2. );
  TH2D *zjt_kaons = new TH2D("zjt_kaons", ";z;j_{T} [GeV]", 100, 0., 1., 100, 0., 2. );
  TH2D *zjt_kaons_kpkmJet = new TH2D("zjt_kaons_kpkmJet", ";z;j_{T} [GeV]", 100, 0., 1., 100, 0., 2. );

  // z r correlations //
  TH2D *zr = new TH2D("zr", ";z;r", 100, 0., 1., 100, 0., 1. );
  TH2D *zr_kpkmJet = new TH2D("zr_kpkmJet", ";z;r", 100, 0., 1., 100, 0., 1. );
  TH2D *zr_kaons = new TH2D("zr_kaons", ";z;r", 100, 0., 1., 100, 0., 1. );
  TH2D *zr_kaons_kpkmJet = new TH2D("zr_kaons_kpkmJet", ";z;r", 100, 0., 1., 100, 0., 1. );

  // Declare histograms differential in jet pT //
  const int nPtBins = 3;
  double ptBins[nPtBins+1] = {15.0, 20.0, 25.0, 60.0};
  TH2D *zr_ptbinned[nPtBins];
  TH2D *zjt_ptbinned[nPtBins];
  for (int i=0; i<nPtBins; ++i)
  {
    TString histName = "zr";
    histName += ptBins[i]; histName+="_"; histName+=ptBins[i+1];
    zr_ptbinned[i] = new TH2D(histName, ";z;r", 100, 0., 1., 100, 0., 1.);
  }
  
  //kpkm correlations //
  TH2D *kpkmz = new TH2D("kpkmz",";zk; zk-", 100, 0., 1., 100, 0., 1.); 
  TH2D *kpkmjt = new TH2D("kpkmjt",";jtk+; jtk-", 100, 0., 2., 100, 0., 2.);
  TH2D *kpkmr = new TH2D("kpkmr",";rk+; rk-", 100, 0., 1., 100, 0., 1.);


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
  jettree->SetBranchAddress("jPt", &jPt);
  jettree->SetBranchAddress("jEta", &jEta);
  jettree->SetBranchAddress("jPhi", &jPhi);
  jettree->SetBranchAddress("kpkmJet", &kpkmJet);
  jettree->SetBranchAddress("kpOrkmJet", &kpOrkmJet);

  // Grab constituent tree and relevant branches from the input file //
  int    cPid;
  double cPx, cPy, cPz, cE, cPt, cEta, cPhi, cM;
  double z, jt, r;
  TTree *contree = (TTree*)infile->Get("cons");
  contree->SetBranchAddress("eventNum", &eventNum);
  contree->SetBranchAddress("jetNum", &jetNum);
  contree->SetBranchAddress("cPid", &cPid);
  contree->SetBranchAddress("cPx", &cPx);
  contree->SetBranchAddress("cPy", &cPy);
  contree->SetBranchAddress("cPz", &cPz);
  contree->SetBranchAddress("cE" , &cE);
  contree->SetBranchAddress("cPt", &cPt);
  contree->SetBranchAddress("cEta", &cEta);
  contree->SetBranchAddress("cPhi", &cPhi);
  contree->SetBranchAddress("z", &z);
  contree->SetBranchAddress("jt", &jt);
  contree->SetBranchAddress("r", &r);
  contree->SetBranchAddress("jPt", &jPt);
  
  int conEntry1 = 0;
  int conEntry2 = 0;
  for (int i=0; i<jettree->GetEntries(); ++i)
  {
    jettree->GetEntry(i);
	TLorentzVector jet = TLorentzVector(jPx, jPy, jPz, jE);
    // K+ K- Correlation Studies // 
    if (kpkmJet)
    {
      for (int j=conEntry1; j<(conEntry1 + nConstituents); ++j)
      {
        contree->GetEntry(j);
        if (jetNum != i) {cout << "something funky 1" << endl;}
        if (abs(cPid)!=321) {continue;} 
		

		if (cPid !=-321) {continue;}
		double kplusr = r;
		double kplusz = z;
		double kplusjt = jt;
	
        int kaon1Pid = cPid;
        TLorentzVector kaon1 = TLorentzVector(cPx, cPy, cPz, cE);
        for (int k=conEntry2; k<(conEntry2 + nConstituents); ++k)
        {
          if (k>=j) {continue;}	
          contree->GetEntry(k);
          if (abs(cPid) != 321) {continue;}
          int kaon2Pid = cPid;
          if (jetNum != i) {cout << "something funky 2" << endl;}
          if (kaon1Pid == kaon2Pid) {continue;}

		   //my work
		  if (cPid != 321) {continue;}
		  double kminusr = r;
		  double kminusz = z;
		  double kminusjt = jt;

          //TLorentzVector kaon2 = TLorentzVector(cPx, cPy, cPz, cE);
          //TLorentzVector kpair = kaon1 + kaon2;
          // Fill histos for K+ K- correlations here // 
		kpkmz->Fill(kplusz, kminusz);
		kpkmr->Fill(kplusr, kminusr);
		kpkmjt->Fill(kplusjt, kminusjt);
		
        } // end con loop 2
      } // end con loop 1
    } // end kpkmJet conditional
    
    conEntry1 += nConstituents;
    conEntry2 += nConstituents;
	}
  outfile->Write();
  outfile->Close();
}
