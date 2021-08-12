/// \Adapted to cluster Jets using FastJet and Study K+ K- in jet correlations

/// \author Andreas Morsch
/// \adapted by Dillon Fitzgerald
/// \edited by Julia Marchese 5/21/21

#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <cmath> 

using namespace fastjet;
using namespace std;

bool verbosity = 0;
void pythia8Jets(Int_t nev  = 5000, Int_t ndeb = 1)
{
	// Outfiles for various bins
	double ptHatMin = 20.0;
	double ptHatMax = 50.0;
	TString ptHatMin_str = "PhaseSpace:pTHatMin = ";
	TString ptHatMax_str = "PhaseSpace:pTHatMax = ";
	ptHatMin_str += ptHatMin;
	ptHatMax_str += ptHatMax;

	TString outFileName = "pythiaJets_";
	outFileName += ptHatMin;
	outFileName += "_";
	outFileName += ptHatMax;
	outFileName += ".root";

	// Outfile
	TFile *outFile = new TFile(outFileName,"RECREATE");
	
	// Load libraries
	gSystem->Load("libEG");
	gSystem->Load("libEGPythia8");

	// ptHat_hist
	TH1D *ptHat_hist = new TH1D("ptHat", ";p_{T} [GeV/c];", 100, 0.0, 100.0);

 	// Jet tree
    int eventNum, jetNum, nConstituents, jettag, jetNumInEvent;
    bool kpkmJet, kpOrkmJet;
    double jPx, jPy, jPz, jE, jPt, jEta, jPhi, rDiff;

    vector<pair<TLorentzVector, vector<pair<TLorentzVector, int>>>> pythiaJets;
    vector<pair<TLorentzVector,int>> constituentList;

    TTree* jettree;
    jettree = new TTree("jet", "A tree with jet info");
    jettree->Branch("eventNum", &eventNum, "eventNum/I");
    jettree->Branch("jetNum", &jetNum, "jetNum/I");
    jettree->Branch("nConstituents", &nConstituents, "nConstituents/I");
    jettree->Branch("jPx", &jPx, "jPx/D");
    jettree->Branch("jPy", &jPy, "jPy/D");
    jettree->Branch("jPz", &jPz, "jPz/D");
    jettree->Branch("jE", &jE, "jE/D");
    jettree->Branch("jPt", &jPt, "jPt/D");
    jettree->Branch("jEta", &jEta, "jEta/D");
    jettree->Branch("jPhi", &jPhi, "jPhi/D");
    jettree->Branch("kpkmJet", &kpkmJet, "kpkmJet/O");
    jettree->Branch("kpOrkmJet", &kpOrkmJet, "kpOrkmJet/O");
    jettree->Branch("jettag", &jettag, "jettag/I");
	jettree->Branch("jetNumInEvent", &jetNumInEvent, "jetNumInEvent/I");
	jettree->Branch("rDiff", &rDiff, "rDiff/D");

    // Constituent tree
    int cPID;
    double r, jt, z, cPx, cPy, cPz, cE, cPt, cEta, cPhi, cM, cCharge;
    TTree* contree;
    contree = new TTree("con", "A tree with constituent info");
    contree->Branch("eventNum", &eventNum, "eventNum/I");
    contree->Branch("jetNum", &jetNum, "jetNum/I");
    contree->Branch("cPID", &cPID, "cPID/I");
    contree->Branch("r", &r, "r/D");
    contree->Branch("jt", &jt, "jt/D");
    contree->Branch("z", &z, "z/D");
    contree->Branch("cPx", &cPx, "cPx/D");
    contree->Branch("cPy", &cPy, "cPy/D");
    contree->Branch("cPz", &cPz, "cPz/D");
    contree->Branch("cE", &cE, "cE/D");
    contree->Branch("cPt" , &cPt,  "cPt/D" );
   	contree->Branch("cEta" , &cEta,  "cEta/D" );
   	contree->Branch("cPhi", &cPhi, "cPhi/D"); 	
   	contree->Branch("cM", &cM, "cM/D");
	contree->Branch("jPt", &jPt, "jPt/D");
	contree->Branch("jPhi", &jPhi, "jPhi/D");
	contree->Branch("jEta", &jEta, "jEta/D");
	contree->Branch("cCharge", &cCharge, "cCharge/D");

    // Choose a jet definition
    double R = 0.5;
    JetDefinition jet_def(antikt_algorithm, R);
    vector<PseudoJet> parts;

	// Create array of particles and pythia8 object
    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    TPythia8* pythia8 = new TPythia8();

	#if PYTHIA_VERSION_INTEGER == 8235
	//  Pythia 8.235 is known to cause crashes:
    printf("ABORTING PYTHIA8 TUTORIAL!\n");
    printf("The version of Pythia you use is known to case crashes due to memory errors.\n");
    printf("They have been reported to the authors; the Pythia versions 8.1... are known to work.\n");
    return;
	#endif

	// Configure and use a reproducible seed.
    pythia8->ReadString("HardQCD:all = on");
    pythia8->ReadString("Random:setSeed = on");
    pythia8->ReadString("Random:seed = 42");

	// pThat cut -- will do this in bins
	pythia8->ReadString(ptHatMin_str);
   	pythia8->ReadString(ptHatMax_str);

	// Initialize pythia8
    pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* GeV */);

	// Initialize variables
	double parton1eta = 0;
	double parton2eta = 0;
	double parton1phi = 0;
	double parton2phi = 0;
	int parton1Pid = 0;
	int parton2Pid = 0;
    jetNum = 0;

	// Initialize parton-type counters
	int quark = 0; // all quarks including strange quarks
	int gluon = 0; // all gluons
	int strange = 0; // strange quarks
	int total = 0; // all initial partons
	bool ptypeverb = 0; // couts for particle type

	// Event loop
    for (Int_t iev = 0; iev < nev; iev++) 
	{ 	
		pythia8->GenerateEvent();
    	
		if (verbosity >  0)
 	  		pythia8->EventListing(); 

    	pythia8->ImportParticles(particles,"All");
  		Int_t np = particles->GetEntriesFast();
    	parts.clear();

  	 	// Particle loop
    	for (Int_t ip = 0; ip < np; ip++)
		{
        	TParticle* partPID = (TParticle*) particles->At(ip);
        	Int_t ist = partPID->GetStatusCode();
			if (ip==4)
			{
				TLorentzVector parton1 = TLorentzVector(partPID->Px(), partPID->Py(), partPID->Pz(), partPID->Energy());
				parton1eta = parton1.Eta();
				parton1phi = parton1.Phi();
				parton1Pid = partPID->GetPdgCode();
				ptHat_hist->Fill(parton1.Pt());
			}
			if (ip==5)
			{
				TLorentzVector parton2 = TLorentzVector(partPID->Px(), partPID->Py(), partPID->Pz(), partPID->Energy());
				parton2eta = parton2.Eta();
				parton2phi = parton2.Phi();
				parton2Pid = partPID->GetPdgCode();
				ptHat_hist->Fill(parton2.Pt());
			}

			// Positive codes are final particles.
        	if (ist <= 0) continue;
        	Int_t pdg = partPID->GetPdgCode();
          	Float_t eta = partPID->Eta();
         	Float_t pt  = partPID->Pt();
  		    if (pt < 0.2) continue;
			if (abs(pdg) != 12 && abs(pdg) != 14 && abs(pdg) != 16)
			{
	    	parts.push_back(PseudoJet(partPID->Px(), partPID->Py(), partPID->Pz(), partPID->Energy() ));
			}

		} // End particle loop

		// cout << " parts size : " << parts.size() << endl;
		
		ClusterSequence cs(parts, jet_def);
    	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
		//cout << " num jets : " << jets.size() << endl;

	  	// Jet loop
   		for (unsigned i = 0; i < jets.size(); i++)
   		{
			jettag = 0;
      		if (jets[i].pt() < 15) continue;
     		kpkmJet = false;
			kpOrkmJet = false;
     		jetNumInEvent = i;

			PseudoJet jet = jets[i];
      		jPx = jet.px();
      		jPy = jet.py();
     		jPz = jet.pz();
    	 	jE = jet.E();
      		jPt = jet.pt();
      		jEta = jet.eta();
      		jPhi = jet.phi();

      		TVector3 jet3(jet.px(), jet.py(), jet.pz());
      		vector<PseudoJet> constituents = jet.constituents();
			nConstituents = constituents.size();

      		if (nConstituents < 3) {continue;}

      		eventNum = iev;
			int kpluscount = 0;
			int kminuscount = 0;
      		
			TLorentzVector jet4Vec = TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.e());
			double jetEta = jet4Vec.Eta();
			double jetPhi = jet4Vec.Phi();

			// New

			if (jEta<2 || jEta>5) {continue;}

//			cout << "jet eta : " << jet4Vec.Eta() << " jet phi : " << jet4Vec.Phi() << endl;
//			cout << "parton1 eta : " << parton1eta << " parton1 phi : " << parton1phi << endl;
//			cout << "parton2 eta : " << parton2eta << " parton2 phi : " << parton2phi << endl;
			
			double jettagr1 = sqrt(pow((jetPhi - parton1phi),2.0)+pow((jetEta-parton1eta), 2.0));
			double jettagr2 = sqrt(pow((jetPhi - parton2phi),2.0)+pow((jetEta-parton2eta), 2.0));
			rDiff = fabs(jettagr1 - jettagr2);

			// Jettag is now defined as the initializing Pid
			if (jettagr1 < jettagr2)
			{	jettag = parton1Pid;}
			else if (jettagr1 > jettagr2)
			{	jettag = parton2Pid;}

			// Counters and if-statements for Pids
			total++;
			if (abs(jettag) == 1 || abs(jettag) == 2 || abs(jettag) == 3 || abs(jettag) == 4 || abs(jettag) == 6 || abs(jettag) == 6)
			{	quark++;}
			if (jettag == 21)
			{	gluon++;}
			if (abs(jettag) == 3)
			{	strange++;}
			
			// Prints out each instance of strange/antistrange quarks -- Sanity check
			if (ptypeverb == 1)
			{	if (jettag == 3)
				{	cout << "initial particle: strange" << endl;}
				if (jettag == -3)
				{	cout << "initial particle: antistrange" << endl;}
			}

			if (ptypeverb == 1)
			{
				{	cout << "total:" << total << endl;}
				{	cout << "quark:" << quark << endl;}
				{	cout << "gluon:" << gluon << endl;}
				{	cout << "strange:" << strange << endl;}

			}

 			// Constituent loop within jet loop
          	for (unsigned j = 0; j < constituents.size(); j++)
			{
				PseudoJet con = constituents[j];
          		TLorentzVector con4Vec = TLorentzVector(con.px(), con.py(), con.pz(), con.e());
          		cPx= con.px();
          		cPy = con.py();
          		cPz = con.pz();
          		cE = con.e();
          		cPt = con.pt();
          		cEta = con.eta();
          		cPhi = con.phi();
          		cM = con.m();

         		TVector3 con3(con.px(), con.py(), con.pz());
		    	z = (con3.Dot(jet3))/(jet3.Mag2());       		
				jt = (con3.Cross(jet3).Mag()) / (jet3.Mag());

				// New calculation of r -- to fix bump around 2pi
				double ucPhi = cPhi + 6.28;
				double lcPhi = cPhi - 6.28; 
				std::vector<double> phiDiffSqVec{pow(ucPhi-jPhi,2), pow(cPhi-jPhi,2), pow(lcPhi-jPhi,2)};
				std::vector<double>::iterator rPhiDiffSq = std::min_element(phiDiffSqVec.begin(), phiDiffSqVec.end());
                r = sqrt(pow(jet.eta() - con.eta(), 2) + phiDiffSqVec.at(std::distance(phiDiffSqVec.begin(), rPhiDiffSq)));
		
				// Particle ID loop 
				for (int k=0;k < np; k++)
				{
					TParticle* partPID = (TParticle*) particles->At(k);
          			int ist = partPID->GetStatusCode();
					double partPx = partPID->Px(); // Parenthesis bc method from class
					double partPy = partPID->Py();
				 	double partPz = partPID->Pz();
				  	double partE = partPID->Energy();
					if (cPx == partPx && cPy == partPy && cPz == partPz && cE == partE)
		    		{
          				cPID = partPID->GetPdgCode();
						Float_t charge = TDatabasePDG::Instance()->GetParticle(cPID)->Charge();
						cCharge = charge;
		   			}
				} // End PID loop				

       			if (cPID == 321) {kpluscount++;}
      			if (cPID == -321) {kminuscount++;}
        		contree->Fill();
        		pair<TLorentzVector, int> conpair = make_pair(con4Vec, cPID);
				constituentList.push_back(conpair);

			} // End constituent loop
			
      		if (kpluscount == kminuscount &&  kpluscount > 0 && kminuscount > 0 ) {kpkmJet = true;}
      		if (kpluscount != kminuscount && (kpluscount > 0 || kminuscount > 0)) {kpOrkmJet = true;}
      		pair<TLorentzVector, vector<pair<TLorentzVector, int>>> pythiaJet = make_pair(jet4Vec, constituentList);
			pythiaJets.push_back(pythiaJet);
		  	jettree->Fill();
			jetNum += 1;
			constituentList.clear();

		} // End jet loop

 	} // End event loop

//  pythia8->PrintStatistics();
	outFile->Write();
	outFile->Close();

} // End macro
