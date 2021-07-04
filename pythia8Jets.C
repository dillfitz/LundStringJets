/// \author Andreas Morsch  -- Adapted to cluster Jets using FastJet and Study K+ K- in jet correlations by Dillon Fitzgerald //
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
void pythia8Jets(Int_t nev  = 1000, Int_t ndeb = 1)
{
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

	TFile *outFile = new TFile(outFileName,"RECREATE");
	// Load libraries //
    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");
	//ptHatMin = 20.0;
	//ptHatMax = 30.0;

	// Setup histograms // 
	TH1D *ptHat_hist = new TH1D("ptHat", ";p_{T} [GeV/c];", 100, 0.0, 100.0);

	// Set up jet tree //
   	int    eventNum, jetNum, jetNumInEvent, nConstituents, partonPid;
    bool kpkmJet, kpOrkmJet;
  	double jPx, jPy, jPz, jE, jPt, jEta, jPhi, partondEta, partondPhi, rDiff;
	vector<pair<TLorentzVector, vector<pair<TLorentzVector, int>>>> pythiaJets;
	vector<pair<TLorentzVector,int>> constituentList;

	TTree* eventtree;
	eventtree = new TTree("events", "A tree with event info");
	eventtree->Branch("pythiaJets", &pythiaJets);

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

  	// choose a jet definition //
  	double R = 0.5;
  	JetDefinition jet_def(antikt_algorithm, R);
  	vector<PseudoJet> parts;

	// Array of particles //
   	TClonesArray* particles = new TClonesArray("TParticle", 1000);
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
   	pythia8->ReadString("HardQCD:all = on");
   	pythia8->ReadString("Random:setSeed = on");
   	// use a reproducible seed: always the same results for the tutorial.
   	pythia8->ReadString("Random:seed = 42");
	// Here is the pT hat cut... //
   	pythia8->ReadString(ptHatMin_str);
   	pythia8->ReadString(ptHatMax_str);


	// Initialize
	// RHIC energy //
   	pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* GeV */);

	// LHC energy //
//   	pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* GeV */);
	// EIC beams and energies //
//   	pythia8->Initialize(11 /* e ^-*/, 2212 /* p */, 18. /* GeV */, 275. /* GeV */);
	
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
      	//if (iev < ndeb) pythia8->EventListing();

        if (verbosity > 0)
      	    pythia8->EventListing();

      	pythia8->ImportParticles(particles,"All");
      	Int_t np = particles->GetEntriesFast();
      	parts.clear();
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
			
         	// Positive codes are final particles.
         	if (ist <= 0) continue;
         	Int_t pdg = part->GetPdgCode();
         	//Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
         	//if (charge == 0.) continue;
         	Float_t eta = part->Eta();
         	Float_t pt  = part->Pt();
	 		if (pt < 0.2) continue;
			if (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) continue;
	 		parts.push_back(PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy() ));
      	} // end particle loop (1)

      	//cout << " parts size : " << parts.size() << endl;
		//cout << " event num : " << iev << endl;
		//cout << " num particles into jet reco : " << parts.size() << endl;
      	ClusterSequence cs(parts, jet_def);
      	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      	//cout << " num jets : " << jets.size() << endl;

      	// Jet loop //
      	for (unsigned i = 0; i < jets.size(); i++)
        {
			partonPid = 0;
            if (jets[i].pt() < (ptHatMin - 5.0)) {continue;}
            kpkmJet = false;
            kpOrkmJet = false;
		    jetNumInEvent = i;
// Apparently this would work if we were using full blown Pythia8 and not TPythia8 //
//            cout << "pythia hard process... " << pythia8->info()->name() << " " << pythia8->info()->code() << endl;
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
            
            // Have Nicole and Julia add this as well //
            if (nConstituents < 3) {continue;}

		    eventNum = iev;
            int kpcount = 0; int kmcount = 0;
			TLorentzVector jet4Vec = TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.e());
			//if (verbosity == 2)
			//{
				//cout << "jet eta : " << jet4Vec.Eta() << " jet phi : " << jet4Vec.Phi() << endl;
				//cout << "parton1 eta : " << parton1eta << " parton1 phi : " << parton1phi << endl;
				//cout << "parton2 eta : " << parton2eta << " parton2 phi : " << parton2phi << endl;
			//}
			double parton1dEta = jet4Vec.Eta() - parton1eta; 
			double parton2dEta = jet4Vec.Eta() - parton2eta; 
			double parton1dPhi = jet4Vec.Phi() - parton1phi; 
			double parton2dPhi = jet4Vec.Phi() - parton2phi; 
			//cout << " p1dEta : " << fabs(parton1dEta) << " p2dEta : " << fabs(parton2dEta) << endl;
			//cout << " p1dPhi : " << fabs(parton1dPhi) << " p2dPhi : " << fabs(parton2dPhi) << endl;
			if ( fabs(parton1dEta) < fabs(parton2dEta) && fabs(parton1dPhi) < fabs(parton2dPhi))
			{
				partondEta = parton1dEta;
				partondPhi = parton1dPhi;
				partonPid = parton1Pid;
			}
			else if ( fabs(parton1dEta) > fabs(parton2dEta) && fabs(parton1dPhi) > fabs(parton2dPhi))
			{
				partondEta = parton2dEta;
				partondPhi = parton2dPhi;
				partonPid = parton2Pid;
			}
			//else
			//{ 
				//cout << "unsuccesfull tag..." << endl;
				//partonPid = 0;
			//}
            // Constituent loop //
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
				pair<TLorentzVector, int> conpair = make_pair(con4Vec, cPid);
				constituentList.push_back(conpair);
	        }
            // end constituent loop //
            if (verbosity == 1)
                cout << "k+ count : " << kpcount << " k- count " << kmcount << endl;

            if (kpcount == kmcount &&  kpcount > 0 && kmcount > 0 ) {kpkmJet = true;}
            if (kpcount != kmcount && (kpcount > 0 || kmcount > 0)) {kpOrkmJet = true;}
			pair<TLorentzVector, vector<pair<TLorentzVector, int>>> pythiaJet = make_pair(jet4Vec, constituentList);
			pythiaJets.push_back(pythiaJet);
		    jettree->Fill();
			jetNum += 1;
			constituentList.clear();
        }
        // end jet loop //
		eventtree->Fill();
		pythiaJets.clear();
    }
    // end event loop //
    pythia8->PrintStatistics();
    outFile->Write();
    outFile->Close();
}
