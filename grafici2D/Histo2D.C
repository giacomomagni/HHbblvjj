
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

using namespace std;

void histo_cutmass( TTree * tree_s, TTree * tree_b, TH2D * hs,  TH2D * hb, string var1, string var2 ){	


	//Dal TTree del background vengono selezionati gli eventi con:  
	// -- mww < 125 [GeV]/c^2
	// -- 120 < mbb < 125 [GeV]/c^2
	//Gli eventi sono poi divisi in due categorie: mjj <= 50 e mjj > 50 [GeV]/c^2
	// -- DeltaR_ljj < 0.9
	//Vengono studiati GLI EVENTI A MASSA mjj ALTA ( entambi i W on-Shell ) 
	//Viene poi fatto un istogramma della pdf di ogni variabile scelta
 
	TTreeReader reader_s( tree_s );
	TTreeReader reader_b( tree_b );
	TTreeReaderValue<double> mbb(reader_b, "mww"); 
	TTreeReaderValue<double> mww(reader_b, "mbb"); 
	TTreeReaderValue<double> mjj_b(reader_b, "mjj");
	
	TTreeReaderValue<double> mjj_s(reader_s, "mjj");
	
	TTreeReaderValue<double>  variable1_b( reader_b, var1.c_str());
	TTreeReaderValue<double>  variable1_s( reader_s, var1.c_str());

	TTreeReaderValue<double>  variable2_b( reader_b, var2.c_str());
	TTreeReaderValue<double>  variable2_s( reader_s, var2.c_str());


	//Selezione eventi backgrond e scrittura degli istogrammi 
	while ( reader_b.Next() ) {
		if ( * mww <= 125. && * mbb >= 120. && * mbb <= 130. && * mjj_b >= 50.){
				hb->Fill(* variable1_b, * variable2_b);
		}
	}	
	//Scrittura degli istogrammi del segnale
	while ( reader_s.Next() ) {
		if ( * mjj_s >= 50.){
				hs->Fill(* variable1_s, * variable2_s);			
		}
	}
	return;
}

void Histo2D( ){	

	int N = 2;
	string * var1 = new string [N]; 
	string * var2 = new string [N]; 

	//AGGIUNGERE LE VARIABILI DA PLOTTARE
	var1[0] = "deltaphi_bbljj";
	var2[0] = "deltaphi_ljj";
	
	var1[1] = "deltar_bbljj";
	var2[1] = "deltar_ljj";

	//LETTURA DEI TTree
    TFile * sinput = TFile::Open( "HH.root" );
    TFile * binput = TFile::Open( "ttbar.root" );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	

	TCanvas ** c = new TCanvas * [N];
	TH2D *** h = new  TH2D ** [N];

	for(int i = 0; i < N; i++){

		//INSERIRE RANGE DEGLI ISTOGRAMMI
		double min1, max1, min2, max2;
		int nbin1, nbin2;
		cout << "Inserire range della grandezza 1 da visualizzare: " << var1[i] <<endl;
		cout << "Min: "; cin >> min1;
		cout << "Max: "; cin >> max1;
		cout << "Nbin: "; cin >> nbin1; 
		cout << "Inserire range della grandezza 2 da visualizzare: " << var2[i] <<endl;
		cout << "Min: "; cin >> min2;
		cout << "Max: "; cin >> max2;
		cout << "Nbin: "; cin >> nbin2; 

//		min1=0; max1=300; nbin1=100;
//		min2=0; max2=300; nbin2=100;

		string title = "Graph " + var1[i] + " vs " + var2[i] ;
		h[i] = new TH2D * [2];
		c[i] = new TCanvas(Form("c%d",i), title.c_str() );	

		//ISTOGRAMMI DELLE DISTIBUZIONI
		string title1 = "Background " + var1[i] + " vs " + var2[i];
		string title2 = "Signal " + var1[i] + " vs " + var2[i];	
 		h[i][1] = new TH2D( title1.c_str(), "Background", nbin1,min1,max1, nbin2,min2,max2);
		h[i][2] = new TH2D( title2.c_str(), "Signal", nbin1,min1,max1, nbin2,min2,max2); 
   		histo_cutmass(signal, background, h[i][2], h[i][1], var1[i], var2[i]);

   		//ISTOGRAMMA pdf
		double N_b = h[i][1]->Integral("width");
		double N_s = h[i][2]->Integral("width");
    	for(int j = 1; j <= nbin1; j++){
			for(int k = 1; k <= nbin2; k++){
    			h[i][1]->SetBinContent(j, k, h[i][1]->GetBinContent(j,k)/N_b);
				h[i][2]->SetBinContent(j, k, h[i][2]->GetBinContent(j,k)/N_s);
			}
		}
		cout<<"Total probability signal in this range : "<< h[i][2]->Integral("width") <<"\n"; 
		cout<<"Total probability background in this range: "<< h[i][1]->Integral("width") <<"\n"; 
	
		c[i]->Divide(1,2);
		c[i]->cd(1);
		gPad->SetGrid(1,1);	
		h[i][2]->SetTitle(title2.c_str());
		h[i][2]->GetXaxis()->SetTitle(var1[i].c_str());
		h[i][2]->GetYaxis()->SetTitle(var2[i].c_str());
		h[i][2]->Draw("COLZ");
		c[i]->cd(2);
		gPad->SetGrid(1,1);
		h[i][1]->SetTitle(title1.c_str());
		h[i][1]->GetXaxis()->SetTitle(var1[i].c_str());
		h[i][1]->GetYaxis()->SetTitle(var2[i].c_str());
		h[i][1]->Draw("COLZ");	
	}
	return;
}
