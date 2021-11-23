/// Jet hadronization studies with Pythia8 + FASTJET for K+ K- in jet correlation studies ///
/// Author(s): Dillon Fitzgerald, Nicole Kuchta, Julia Marchese ///
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>

using namespace fastjet;
using namespace std;

int verbosity = 0;
void pythia8Jets(Int_t nev  = 100, Int_t ndeb = 1)
{
	double ptHatMin = 20.0;
	double ptHatMax = 50.0;
	TString ptHatMin_str = "PhaseSpace:pTHatMin = ";
	TString ptHatMax_str = "PhaseSpace:pTHatMax = ";
	ptHatMin_str += ptHatMin;
	ptHatMax_str += ptHatMax;

	TString outFileName = "data/pythiaJets_";
	outFileName += ptHatMin;
	outFileName += "_";
	outFileName += ptHatMax;
	outFileName += ".root";

	TFile *outFile = new TFile(outFileName,"RECREATE");
	// Load libraries //
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

	// Setup histograms // 
	TH1D *ptHat_hist = new TH1D("ptHat", ";p_{T} [GeV/c];", 100, 0.0, 100.0);

	// Set up jet tree //
  int    eventNum, jetNum, jetNumInEvent, nConstituents, partonPid;
  bool kpkmJet, kpOrkmJet;
  double jPx, jPy, jPz, jE, jPt, jEta, jPhi, partondEta, partondPhi, rDiff;

  TTree* jettree;
  jettree = new TTree("jets","A tree with jet info");
  jettree->Branch("eventNum", &eventNum, "eventNum/I");
  jettree->Branch("jetNum", &jetNum, "jetNum/I");
  jettree->Branch("jetNumInEvent", &jetNumInEvent, "jetNumInEvent/I");
  jettree->Branch("nConstituents", &nConstituents, "nConstituents/I"); 
  jettree->Branch("jPx", &jPx, "jPx/D");
  jettree->Branch("jPy", &jPy, "jPy/D");
  jettree->Branch("jPz", &jPz, "jPz/D");
  jettree->Branch("jE" , &jE,  "jE/D" );
  jettree->Branch("jPt" , &jPt,  "jPt/D" );
  jettree->Branch("jEta" , &jEta,  "jEta/D" );
  jettree->Branch("jPhi", &jPhi, "jPhi/D"); 
  jettree->Branch("kpkmJet", &kpkmJet, "kpkmJet/O");
  jettree->Branch("kpOrkmJet", &kpOrkmJet, "kpOrkmJet/O");
	jettree->Branch("partondEta", &partondEta, "partondEta/D");
	jettree->Branch("partondPhi", &partondPhi, "partondPhi/D");
	jettree->Branch("partonPid", &partonPid, "partonPid/I");
	jettree->Branch("rDiff", &rDiff, "rDiff/D");

	// Set up jet constituent tree //
  int    cPid, cQ;
  double cPx, cPy, cPz, cE, cPt, cEta, cPhi, cM;
  double z, jt, r;
  TTree* contree;
  contree = new TTree("cons", "A tree with jet constituent information");
  contree->Branch("eventNum", &eventNum, "eventNum/I");
  contree->Branch("jetNum", &jetNum, "jetNum/I");
  contree->Branch("cPid", &cPid, "cPid/I");
  contree->Branch("cQ",  &cQ,   "cQ/I");
  contree->Branch("cPx", &cPx, "cPx/D");
  contree->Branch("cPy", &cPy, "cPy/D");
  contree->Branch("cPz", &cPz, "cPz/D");
  contree->Branch("cE" , &cE,  "cE/D" );
  contree->Branch("cPt" , &cPt,  "cPt/D" );
  contree->Branch("cEta" , &cEta,  "cEta/D" );
  contree->Branch("cPhi", &cPhi, "cPhi/D"); 	
  contree->Branch("cM", &cM, "cM/D");
  contree->Branch("z", &z, "z/D"); 	
  contree->Branch("jt", &jt, "jt/D"); 	
  contree->Branch("r", &r, "r/D"); 	 
  contree->Branch("jPt" , &jPt,  "jPt/D" );	
  contree->Branch("jEta" , &jEta,  "jEta/D" );	
  contree->Branch("jPhi" , &jPhi,  "jPhi/D" );	
	
//new tree october 5//
	int bbcount,bbevent,bbPID;
	double bbPx, bbPy, bbPz, bbE;
	TTree* bbtree;
	bbtree = new TTree ("bbtree","a tree with b/b-bar info");
	bbtree->Branch("bbcount",&bbcount,"bbcount/I");
	bbtree->Branch("bbPx",&bbPx,"bbPx/D");
	bbtree->Branch("bbPy",&bbPy,"bbPy/D");
	bbtree->Branch("bbPz",&bbPz,"bbPz/D");
	bbtree->Branch("bbE",&bbE,"bbE/D");
	bbtree->Branch("bbevent",&bbevent,"bbevent/I");
	bbtree->Branch("bbPID",&bbPID,"bbPID/I");


	int muplusPDG, kaonsPDG, muminusPDG;
	double muminusPx, muplusPx, kaonsPx, muminusPy, muplusPy, kaonsPy,muminusPz, muplusPz, kaonsPz, muminusE, muplusE, kaonsE;
	TTree *daughtertree;
	daughtertree = new TTree ("daughtertree","a tree with decay daughter info");
	daughtertree->Branch("muminusPx", &muminusPx, "muminusPx/D");
	daughtertree->Branch("muminusPy", &muminusPy, "muminusPy/D");
	daughtertree->Branch("muminusPz", &muminusPz, "muminusPz/D");
	daughtertree->Branch("muplusPx", &muplusPx, "muplusPx/D");
	daughtertree->Branch("muplusPy", &muplusPy, "muplusPy/D");
	daughtertree->Branch("muplusPz", &muplusPx, "muplusPz/D");
	daughtertree->Branch("kaonsPx", &kaonsPx, "kaonsPx/D");
	daughtertree->Branch("kaonsPy", &kaonsPy, "kaonsPy/D");
	daughtertree->Branch("kaonsPz", &kaonsPz, "kaonsPz/D");
	daughtertree->Branch("kaonsE", &kaonsE, "kaonsE/D");
	daughtertree->Branch("muplusE", &muplusE, "muplusE/D");
	daughtertree->Branch("muminusE", &muminusE, "muminusE/D");
	daughtertree->Branch("muplusPDG", &muplusPDG, "muplusPDG/I");
	daughtertree->Branch("muminusPDG", &muminusPDG, "muminusPDG/I");
	daughtertree->Branch("kaonsPDG", &kaonsPDG, "kaonsPDG/I");
	

  // choose a jet definition //
  double R = 0.5;
  JetDefinition jet_def(antikt_algorithm, R);
  vector<PseudoJet> parts;

	// Array of particles //
  TClonesArray* particles = new TClonesArray("TParticle", 1000000);
	// Create pythia8 object //
  TPythia8* pythia8 = new TPythia8();

  #if PYTHIA_VERSION_INTEGER == 8235
    // Pythia 8.235 is known to cause crashes: //
    printf("ABORTING PYTHIA8 TUTORIAL!\n");
   	printf("The version of Pythia you use is known to case crashes due to memory errors.\n");
   	printf("They have been reported to the authors; the Pythia versions 8.1... are known to work.\n");
   	return;
  #endif

  // Configure
  pythia8->ReadString("HardQCD:hardbbbar = on");
  pythia8->ReadString("Random:setSeed = on");
  // use a reproducible seed: always the same results for the tutorial.
  pythia8->ReadString("Random:seed = 42"); //make this random per event
  // Here is the pT hat cut... //
  pythia8->ReadString(ptHatMin_str);
  pythia8->ReadString(ptHatMax_str);
  pythia8->ReadString("521:oneChannel = 1 0.0010600 0 443 321");   
  pythia8->ReadString("443:oneChannel = 1 0.0593000 0 13 -13");





	// Initialize
	// RHIC energy //
  //pythia8->Initialize(2212 /* p */, 2212 /* p */, 200. /* GeV */);

	// LHC energy //
  pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* GeV */);
  
	// EIC beams and energies //
  //pythia8->Initialize(11 /* e ^-*/, 2212 /* p */, 18. /* GeV */, 275. /* GeV */);
	
  double parton1eta = 0;
  double parton2eta = 0;
  double parton1phi = 0;
  double parton2phi = 0;
  int parton1Pid = 0;
  int parton2Pid = 0;
  jetNum = 0;
  // Event loop
  for (Int_t iev = 0; iev < nev; iev++) 
	{
    pythia8->GenerateEvent();
    //cout << "pythia hard process... " << pythia8->Pythia8()->info.name() << " " << pythia8->Pythia8()->info.code() << endl;
    if (verbosity > 0)
      pythia8->EventListing();

    pythia8->ImportParticles(particles,"All");
    Int_t np = particles->GetEntriesFast();
    parts.clear();
	bool bplusflag = false;
	bool jpsifromb = false;

		// Particle loop (1)
    for (Int_t ip = 0; ip < np; ip++) 
    {
      TParticle* part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
      if (ip==4)
      {
        TLorentzVector parton1 = TLorentzVector(part->Px(), part->Py(), part->Pz(), part->Energy());
        parton1eta = parton1.Eta();
        parton1phi = parton1.Phi();
        parton1Pid = part->GetPdgCode();
        ptHat_hist->Fill(parton1.Pt());
        // store pT of parton in histogram -- add matched pT to jet tree 
      }
      if (ip==5)
      {
        TLorentzVector parton2 = TLorentzVector(part->Px(), part->Py(), part->Pz(), part->Energy());
        parton2eta = parton2.Eta();
        parton2phi = parton2.Phi();
        parton2Pid = part->GetPdgCode();
        ptHat_hist->Fill(parton2.Pt());
        // store pT of parton in histogram -- add matched pT to jet tree
      }
	  int bdecaydaughter1, bdecaydaughter2, jpsidecaydaughter1, jpsidecaydaughter2;
	  int bbcount = 0; 
	  if (part->Eta() < 2. || part->Eta() > 5.)  {continue;}
	  if (part->GetPdgCode() == abs(521)) 
	  {
		bplusflag = true;
	  	bbcount++;  
	  	bbPx=part->Px();
	  	bbPy=part->Py();
	  	bbPz=part->Pz();
	  	bbE=part->Energy();
	  	bbPID=part->GetPdgCode();
	  	bbevent=iev;
	  	bbtree->Fill();
	  	bdecaydaughter1=part->GetDaughter(0);
	  	bdecaydaughter2=part->GetDaughter(1);
	  }
	  if ((ip == bdecaydaughter1 || ip == bdecaydaughter2)&& bplusflag &&abs(part->GetPdgCode())==443)
	  { 
	  	//cout<<"decay daughter 1 code "<<part->GetPdgCode()<<endl;
		jpsifromb = true;
		jpsidecaydaughter1=part->GetDaughter(0);
		jpsidecaydaughter2=part->GetDaughter(1);
	  }
	  //kaon from b decay (and later muons)
		if ((ip == bdecaydaughter1 || ip == bdecaydaughter2)&& bplusflag &&abs(part->GetPdgCode())==321)
	  { 
	 kaonsPx=part->Px();
	 kaonsPy=part->Py();
	 kaonsPz=part->Pz();
	 kaonsE=part->Energy();
	 kaonsPDG=part->GetPdgCode();
	 cout<<"kaons"<<endl;
	  }
	if ((ip == jpsidecaydaughter1 || ip == jpsidecaydaughter2)&& jpsifromb && part->GetPdgCode()==13)
	  { 
	 muplusPx=part->Px();
	 muplusPy=part->Px();
	 muplusPz=part->Py();
	 muplusE=part->Energy();
	 muplusPDG=part->GetPdgCode();
	 cout<<"muplus"<<endl;
	 cout<<"decaydaughters " <<ip<<endl;
	 cout<<"jpsidecaydaughter1 " <<jpsidecaydaughter1<<endl;
	 cout<<"jpsidecaydaughter2 " <<jpsidecaydaughter2<<endl;
	  }
	if ((ip == jpsidecaydaughter1 || ip == jpsidecaydaughter2)&& jpsifromb && part->GetPdgCode()==-13)
	  { 
	 muminusPx=part->Px();
	 muminusPy=part->Py();
	 muminusPz=part->Pz();
	 muminusE=part->Energy();
	 muminusPDG=part->GetPdgCode();
	 cout<<"muminus"<<endl;
	 cout<<"decaydaughters " <<ip<<endl;
	 cout<<"jpsidecaydaughter1 " <<jpsidecaydaughter1<<endl;
	 cout<<"jpsidecaydaughter2 " <<jpsidecaydaughter2<<endl;
	 daughtertree->Fill();
	  }
//create trees and make sure amount of entries are the same !! show 100,000 show 3 to four percent we expect and me';ve moved onto investigating decay daughters show how to fix decays from chat
      // Positive codes are final particles.
      if (ist <= 0) continue;
      Int_t pdg = part->GetPdgCode();
      Float_t eta = part->Eta();
      Float_t pt  = part->Pt();
      if (pt < 0.2) continue;
      if (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) continue;
      parts.push_back(PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy() ));
      
    } // end particle loop (1)

    // Cluster the jets //
    ClusterSequence cs(parts, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    // Jet loop //
    for (unsigned i = 0; i < jets.size(); i++)
    {
      partonPid = 0;
      if (jets[i].pt() < (ptHatMin - 5.0)) {continue;}
	  if (jets[i].eta() < 2 + R || jets[i].eta() > 5 - R) {continue;}
      kpkmJet = false;
      kpOrkmJet = false;
      jetNumInEvent = i;
      
      // Have Nicole and Julia change this pT cut to 15 GeV! //
      if (verbosity == 1)
        cout << "jet num : " <<  jetNum << endl;
        
      PseudoJet jet = jets[i];
		  jPx = jet.px(); 
		  jPy = jet.py(); 
		  jPz = jet.pz(); 
		  jE = jet.E();
		  jPt = jet.pt(); 
		  jEta = jet.eta(); 
		  jPhi = jet.phi_std();
		  TVector3 jet3(jet.px(), jet.py(), jet.pz());
		  vector<PseudoJet> constituents = jet.constituents();
		  nConstituents = constituents.size();
		  
		  if (nConstituents < 3) {continue;}
		  eventNum = iev;
		  int kpcount = 0; int kmcount = 0;
		  TLorentzVector jet4Vec = TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.e());
		  double jetTagr1 = sqrt(pow(jet4Vec.Eta() - parton1eta, 2) + pow(jet4Vec.Phi() - parton1phi, 2) );
		  double jetTagr2 = sqrt(pow(jet4Vec.Eta() - parton2eta, 2) + pow(jet4Vec.Phi() - parton2phi, 2) );
		  rDiff = fabs(jetTagr1 - jetTagr2);
		  if (jetTagr1 < jetTagr2)
		  {
		    partonPid = parton1Pid;
		  }
      else if (jetTagr2 < jetTagr1)
      {
        partonPid = parton2Pid;
      }

      // Constituent loop //
      for (unsigned j = 0; j < constituents.size(); j++)
      {
        PseudoJet con = constituents[j];
        cPx= con.px();
        cPy = con.py();
        cPz = con.pz();
        cE = con.e();
        cPt = con.pt();
        cEta = con.eta();
        cPhi = con.phi_std();
        cM = con.m();
        
        TVector3 con3(con.px(), con.py(), con.pz());
        z = (con3.Dot(jet3))/(jet3.Mag2());
        jt = ((con3.Cross(jet3)).Mag())/ (jet3.Mag());
        
        // implement fix to get rid  of bumnp around 2pi in r distribution //
        double ucPhi = cPhi + 6.28;
        double lcPhi = cPhi - 6.28; 
        std::vector<double> phiDiffSqVec{pow(ucPhi-jPhi,2), pow(cPhi-jPhi,2), pow(lcPhi-jPhi,2)};
        std::vector<double>::iterator rPhiDiffSq = std::min_element(phiDiffSqVec.begin(), phiDiffSqVec.end());
        r = sqrt(pow(jet.eta() - con.eta(), 2) + phiDiffSqVec.at(std::distance(phiDiffSqVec.begin(), rPhiDiffSq)));
        
        // PID particle loop //
        for (int k = 0; k < np; k++) 
        {
          TParticle* part = (TParticle*) particles->At(k);
          int ist = part->GetStatusCode();
          // Positive codes are final particles.
          if (ist <= 0) continue;
          double partPx = part->Px();
          double partPy = part->Py();
          double partPz = part->Pz();
          double partE  = part->Energy();
          if (cPx == partPx && cPy == partPy && cPz == partPz && cE == partE)
          {
            cPid = part->GetPdgCode();
            cQ = TDatabasePDG::Instance()->GetParticle(cPid)->Charge();
            if (verbosity == 1)
              cout << "constituent PID : " << cPid << endl;
          }
        }
        
       // end PID particle loop //
        if (cPid ==  321) {kpcount++;}
        if (cPid == -321) {kmcount++;}
        contree->Fill();
      }
      // end constituent loop //
      
      if (verbosity == 1)
        cout << "k+ count : " << kpcount << " k- count " << kmcount << endl;
      
      if (kpcount == kmcount &&  kpcount > 0 && kmcount > 0 ) {kpkmJet = true;}
      if (kpcount != kmcount && (kpcount > 0 || kmcount > 0)) {kpOrkmJet = true;}
      jettree->Fill();
      jetNum += 1;

    }
    // end jet loop //
  }
  // end event loop //
  pythia8->PrintStatistics();
  outFile->Write();
  outFile->Close();
}
