//Read and store all the variables
//

// c++ -o completereading completereading.cpp `root-config --cflags --glibs`
#include "LHEF.h"
#include<iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include <fstream>
#include <algorithm>
#include "TLorentzVector.h"

using namespace std;

int main(int argc, char** argv){
    if (argc < 3){
        cout << "Usage: " << argv[0] << " output.root file1.lhe.." << endl;
        return 1;
    }

    char* rootfile = argv[1];     
    TFile output( rootfile, "RECREATE");
    TTree* tree = new TTree("tree", "Background events");
    
    double lep_eta, lep_pt, lep_E, lep_Et, lep_phi;
    double n_eta, n_pt;
    double j1_eta, j1_pt, j1_phi, j2_eta, j2_pt, j2_phi;
	double b1_eta, b1_pt, b1_phi, b2_eta, b2_pt, b2_phi;
    double mww, mjj, mbb;
	double ww_pt, j_pt, bb_pt;
	double j_Et, bb_Et, ww_Et;
	double Ht, Htnu, Ptm;
	double deltar_ljj, deltar_bbljj, deltar_bb;
	double deltaphi_ljj, deltaphi_bbljj, deltaphi_bb, deltaphi_bwbw;
    
    TH1F *q1 = new TH1F("q1", "PGID of j1", 20, -10, 10);
	TH1F *q2 = new TH1F("q2", "PGID of j2", 20, -10, 10);   
    TH1F* l1 = new TH1F("l1", "PGID of lepton", 40, -20, 20);
   
    tree->Branch("lep_pt", &lep_pt);
	tree->Branch("lep_eta", &lep_eta);
	tree->Branch("lep_E", &lep_E);
	tree->Branch("lep_Et", &lep_Et);
    tree->Branch("n_pt", &n_pt);
    tree->Branch("j1_pt", &j1_pt);
	tree->Branch("j1_eta", &j1_eta);
    tree->Branch("j2_pt", &j2_pt);	
    tree->Branch("b1_pt", &b1_pt);
	tree->Branch("b1_eta", &b1_eta);
    tree->Branch("b2_pt", &b2_pt);

	tree->Branch("mjj", &mjj);
	tree->Branch("j_pt", &j_pt);
	tree->Branch("j_Et", &j_Et);
	tree->Branch("mww", &mww);
	tree->Branch("ww_pt", &ww_pt);
	tree->Branch("ww_Et", &ww_Et);
	tree->Branch("mbb", &mbb);
	tree->Branch("bb_pt", &bb_pt);
	tree->Branch("bb_Et", &bb_Et);
	
	tree->Branch("Ht", &Ht);
	tree->Branch("Htnu", &Htnu);
	tree->Branch("Ptm", &Ptm);

	tree->Branch("deltaphi_bb", &deltaphi_bb);
	tree->Branch("deltaphi_ljj", &deltaphi_ljj);
	tree->Branch("deltaphi_bbljj", &deltaphi_bbljj);
	tree->Branch("deltaphi_bwbw", &deltaphi_bwbw);

	tree->Branch("deltar_bb", &deltar_bb);
	tree->Branch("deltar_ljj", &deltar_ljj);
	tree->Branch("deltar_bbljj", &deltar_bbljj);
	
    long iEv = 0;
	int tau = 0;
    for (int i = 2; i < argc; i++){

        //Open a stream connected to an event file:
        ifstream ifs(argv[i]);
        
        //Create the Reader object:
        LHEF::Reader reader(ifs);
        
        //Now loop over all the events     
        while ( reader.readEvent() ){
            iEv++;
//            cout << "Event " << iEv << endl;
            
            vector<int> charged_leptons; 
            vector<int> neutrinos;
            vector<int> leptons;
            vector<int> quarks_jet; 
			vector<int> quarks_b; 
            vector<int> jets;
            
            //TlorentVector of j and bosons w
            TLorentzVector jet_j1_mom, jet_j2_mom, w1_mom;
            //TlorentVector of lepton and neutrino
            TLorentzVector lep_mom, nu_mom; 
			//TlorentVector of quarks beauty and antibeauty
            TLorentzVector b1_mom, b2_mom, bb_mom, bw1_mom, bw2_mom; 
            
            //Save quadrimomentum of all particles mapped with position in the event
            map<int, TLorentzVector> momenta;
            
            //Save eta and pt of particles            
            //PG loop over particles in the event
            for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); iPart++){
//                cout    << "\t part type [" << iPart << "] " << reader.hepeup.IDUP.at (iPart)
//                         << "\t status " << reader.hepeup.ISTUP.at (iPart)
//                         << endl;
                        
                //Saving info of final state particles only
                if (reader.hepeup.ISTUP.at (iPart) == 1){
                    int ID = reader.hepeup.IDUP.at(iPart);

					//PUP --> quadrimomentum                   
 					TLorentzVector momentum(
                        reader.hepeup.PUP.at (iPart).at (0), //PG px
                        reader.hepeup.PUP.at (iPart).at (1), //PG py
                        reader.hepeup.PUP.at (iPart).at (2), //PG pz
                        reader.hepeup.PUP.at (iPart).at (3) //PG E
                    );
      
                    // Save momentum
                    momenta[iPart] =  momentum;
					// eletrons, positrons, muon, antimuon 
					// Events with a Tau or an antitau are not considered
                    if( abs(ID) == 11 || abs(ID) == 13 ){  
                        charged_leptons.push_back(iPart);
                        leptons.push_back(iPart);
                        lep_mom = momentum;
                        lep_eta = momentum.Eta();
                        lep_pt = momentum.Pt();
						lep_E = momentum.E();
						lep_Et = momentum.Et();
						lep_phi = momentum.Phi();
                        l1->Fill(ID);
                    }
					if (abs(ID) == 15){ tau = 1;}

                    // neutrinos
                    if ( abs(ID) == 12 || abs(ID) == 14 ){  
                        neutrinos.push_back(iPart);
                        leptons.push_back(iPart);
                        nu_mom = momentum;
                        n_eta = momentum.Eta();
                        n_pt = momentum.Pt();
                    }
               		// quark beauty 
                    if ( ID == 5 ){  
                        quarks_b.push_back(iPart);
                        b1_mom = momentum;
                        b1_eta = momentum.Eta();
                        b1_pt = momentum.Pt();
						b1_phi = momentum.Phi();
                    }
					// quark antibeauty 
                    if ( ID == -5 ){  
                        quarks_b.push_back(iPart);
                        b2_mom = momentum;
                        b2_eta = momentum.Eta();
                        b2_pt = momentum.Pt();
						b2_phi = momentum.Phi();
                    }
                    // Other quarks in final state are from the jet 
                    if ( abs(ID) < 6 && abs(ID) != 5 ){  
                        quarks_jet.push_back(iPart);
                        jets.push_back(iPart);
                    }   
                }// info of final state particles
            }// end of loop over particles	
			
			// mass and pt of the two beauty           
			mbb = ( b1_mom + b2_mom ).M();
			bb_pt = ( b1_mom + b2_mom ).Pt();
			bb_Et = ( b1_mom + b2_mom ).Et();    			
        
            // The j jet are saved looking at the element of the jet vector.
			// // Events with a DOUBLE JETS or NO JET are not saved.   
			if(jets.size() > 0 && jets.size() < 3 && tau == 0){      
            	for (int i=0; i<jets.size(); i += 2){
                       	jet_j1_mom = momenta[jets[i]];
                      	jet_j2_mom = momenta[jets[i+1]];
						j1_pt = jet_j1_mom.Pt();
            			j1_eta = jet_j1_mom.Eta();
						j1_phi = jet_j1_mom.Phi();
           				j2_pt = jet_j2_mom.Pt();
            			j2_eta = jet_j2_mom.Eta();
						j2_phi = jet_j2_mom.Phi();						
						mjj = (jet_j1_mom + jet_j2_mom).M();
           				j_pt = (jet_j1_mom + jet_j2_mom).Pt();
						j_Et = (jet_j1_mom + jet_j2_mom).Et();

						//bosons mass and momentum, the neutrinos are not measured  
            			mww = ( jet_j1_mom + jet_j2_mom + lep_mom ).M();
						ww_pt = ( jet_j1_mom + jet_j2_mom + lep_mom ).Pt();

						//flavour of jets are saved in an histogram
						q1->Fill(reader.hepeup.IDUP.at(jets[i]));
            			q2->Fill(reader.hepeup.IDUP.at(jets[i+1]));
						
						// saving Ht, Htnu, Ptm
						Ht = b1_pt + b2_pt + j1_pt + j2_pt + lep_pt;
						Htnu = b1_pt + b2_pt + j1_pt + j2_pt + lep_pt + n_pt;  
						Ptm = ( b1_mom + b2_mom + jet_j1_mom + jet_j2_mom + lep_mom ).Pt();

						//Angolar variables 
						deltaphi_bb = b1_mom.DeltaPhi( b2_mom );
						deltar_bb = b1_mom.DeltaR( b2_mom );
						deltaphi_ljj = lep_mom.DeltaPhi( jet_j1_mom + jet_j2_mom );
						deltar_ljj = lep_mom.DeltaR( jet_j1_mom + jet_j2_mom );

						bb_mom = b1_mom + b2_mom;
						w1_mom = jet_j1_mom + jet_j2_mom; 
						deltaphi_bbljj = bb_mom.DeltaPhi( lep_mom + w1_mom);
						deltar_bbljj = bb_mom.DeltaR( lep_mom + w1_mom);

						bw1_mom = b1_mom + lep_mom; 
						bw2_mom = b2_mom + w1_mom;
						deltaphi_bwbw = bw1_mom.DeltaPhi( bw2_mom );
						tree->Fill();
                 }//end of jet construction 
			}else{
//			cout << "Event " << iEv << " is not saved. " << endl;
			tau = 0; 
			}     
        }//end of loop over events        
        ifs.close();
    }//end of loop over files
    l1->Write();
    q1->Write();
    q2->Write();
    tree->Write();
    output.Close();    
    return 0;
}
