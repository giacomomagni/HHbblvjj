//Programma per individuare le vaiabili più significative con il taglio sulle masse.
// -- Selezione delgli eventi con tagli sulle masse invarianti
// -- Test di Kolmogrov
// -- Plot della pdf(x) e fit con distribuzione gamma
// -- Plot di tre fattori di merito 

// c++ -o mcut mcut.cpp `root-config --cflags --glibs`
//./tcut q1 q2...

#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <THStack.h>
#include <TApplication.h>
#include <fstream>

using namespace std;

void histo_cutmass( TTree * tree_s, TTree * tree_b, TH1D * hs,  TH1D * hb, string var ){	

	//Dal TTree del background vengono selezionati gli eventi con:  
	// -- mww < 125 [GeV]/c^2
	// -- 120 < mbb < 125 [GeV]/c^2
	//Gli eventi sono poi divisi in due categorie: mvbs <= 50 e mvbs > 50 [GeV]/c^2
	// -- DeltaR_ljj < 0.9
	//Vengono studiati GLI EVENTI A MASSA MVBS ALTA ( entambi i W on-Shell ) 
	//Viene poi fatto un istogramma della pdf di ogni variabile scelta
 
	TTreeReader reader_s( tree_s );
	TTreeReader reader_b( tree_b );
	TTreeReaderValue<double> mbb(reader_b, "mww"); 
	TTreeReaderValue<double> mww(reader_b, "mbb"); 
	TTreeReaderValue<double> mvbs_b(reader_b, "mvbs");
	TTreeReaderValue<double>  b1_phi_b( reader_b, "b1_phi");
	TTreeReaderValue<double>  b2_phi_b( reader_b, "b2_phi");
	TTreeReaderValue<double>  vbs1_phi_b( reader_b, "vbs1_phi");
	TTreeReaderValue<double>  vbs2_phi_b( reader_b, "vbs2_phi");
	TTreeReaderValue<double>  lep_phi_b( reader_b, "lep_phi");
	TTreeReaderValue<double>  b1_eta_b( reader_b, "b1_eta");
	TTreeReaderValue<double>  b2_eta_b( reader_b, "b2_eta");
	TTreeReaderValue<double>  vbs1_eta_b( reader_b, "vbs1_eta");
	TTreeReaderValue<double>  vbs2_eta_b( reader_b, "vbs2_eta");
	TTreeReaderValue<double>  lep_eta_b( reader_b, "lep_eta");

	TTreeReaderValue<double> mvbs_s(reader_s, "mvbs");
	TTreeReaderValue<double>  b1_phi_s( reader_s, "b1_phi");
	TTreeReaderValue<double>  b2_phi_s( reader_s, "b2_phi");
	TTreeReaderValue<double>  vbs1_phi_s( reader_s, "vbs1_phi");
	TTreeReaderValue<double>  vbs2_phi_s( reader_s, "vbs2_phi");
	TTreeReaderValue<double>  lep_phi_s( reader_s, "lep_phi");
	TTreeReaderValue<double>  b1_eta_s( reader_s, "b1_eta");
	TTreeReaderValue<double>  b2_eta_s( reader_s, "b2_eta");
	TTreeReaderValue<double>  vbs1_eta_s( reader_s, "vbs1_eta");
	TTreeReaderValue<double>  vbs2_eta_s( reader_s, "vbs2_eta");
	TTreeReaderValue<double>  lep_eta_s( reader_s, "lep_eta");
	TTreeReaderValue<double>  variable_b( reader_b, var.c_str());
	TTreeReaderValue<double>  variable_s( reader_s, var.c_str());

	//Selezione eventi backgrond e scrittura degli istogrammi 
	while ( reader_b.Next() ) {
		if ( * mww <= 125. && * mbb >= 120. && * mbb <= 130. && * mvbs_b >= 50.){
				hb->Fill(* variable_b);
		}
	}	
	//Scrittura degli istogrammi del segnale
	while ( reader_s.Next() ) {
		if ( * mvbs_s >= 50.){
				hs->Fill(* variable_s);			
		}
	}
	return;
}

double gammadist(double * x, double * p ){
	return pow(x[0]/p[1], p[0]-1)*exp(-1.*x[0]/p[1])/(p[1]*TMath::Gamma(p[0]));
}

void setstack(TH1D * hs, TH1D * hb, THStack * stack ){	
	hb->SetFillStyle(3003);
	hs->SetFillStyle(3004);
    hs->SetFillColor(2);
	hb->SetFillColor(4);
	hs->SetLineColor(2);
	hb->SetLineColor(4);	
	stack->Add(hs, "S");
	stack->Add(hb, "S");
    stack->Draw("nostack");
	gPad->SetGrid(1,1);		
	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
    legend1->AddEntry(hs,"Signal", "f");
    legend1->AddEntry(hb,"Background", "f");
    legend1->Draw("SAME");
	return;
}

void DrawHisto(TTree * signal, TTree * background, string * var, TCanvas ** c, TH1D *** h, int N ){	
	for(int i = 0; i < N; i++){

		//INSERIRE RANGE DEGLI ISTOGRAMMI
		double min, max;
		int nbin;
		cout << "Inserire range della grandezza da visualizzare: " << var[i] <<endl;
		cout << "Min: "; cin >> min;
		cout << "Max: "; cin >> max;
		cout << "Nbin: "; cin >> nbin; 
//		min=0; max=300; nbin=100;

		string title = "Graph_" + var[i];
		h[i] = new TH1D * [8];
		c[i] = new TCanvas(Form("c%d",i), title.c_str() );	
		c[i]->Divide(2,2);

		//ISTOGRAMMI DELLE DISTIBUZIONI
//		c[i]->cd(1);
		string title1 = "Background_" + var[i];
		string title2 = "Signal_" + var[i];
//		THStack * dist = new THStack("dist", "Distribution histo");
   		h[i][1] = new TH1D( title2.c_str(), "Signal", nbin,min,max);
 		h[i][2] = new TH1D( title1.c_str(), "Background", nbin,min,max); 
   		histo_cutmass(signal, background, h[i][1], h[i][2], var[i]);
//		setstack( h[i][1], h[i][2], dist );
//    	dist->SetTitle("Distribution;  ; Events");

   		//ISTOGRAMMA pdf
		c[i]->cd(1);
		THStack * pdf = new THStack("pdf", "Pdf histo");
   		h[i][3] = new TH1D(Form("ps%d", i),"Signal Probability distribution", nbin,min,max);
		h[i][4] = new TH1D(Form("pb%d", i),"Background Probability distribution", nbin,min,max);
    	for(int j = 1; j <= nbin; j++){
    		h[i][3]->SetBinContent(j, h[i][1]->GetBinContent(j)/h[i][1]->Integral("width"));
			h[i][4]->SetBinContent(j, h[i][2]->GetBinContent(j)/h[i][2]->Integral("width"));
		}
//		cout<<"Total probability signal in this range : "<< h[i][3]->Integral("width") <<"\n"; 
//		cout<<"Total probability background in this range: "<< h[i][4]->Integral("width") <<"\n"; 
	
		//TEST DI KOLMOGROF 
		//poichè i dati sarebbero unbinned questo test è significativo solo se la larghezza dei bin è sufficientemente più piccola 
		//della precisione con cui si vuole ottenere in risultato 
		TH1D * hs = new TH1D(Form("cps%d", i),"Signal Cumulative Probability distribution", nbin,min,max);
		TH1D * hb = new TH1D(Form("cpb%d", i),"Background Cumulative Probability distribution", nbin,min,max);
    	for(int j = 1; j <= nbin; j++){
    		hs->SetBinContent(j, h[i][3]->Integral(1,j,""));
			hb->SetBinContent(j, h[i][4]->Integral(1,j,""));
		}
		cout<< "KOLMOGROF TEST: probability that the two distribuions comes from the same phenomena:"  << hs->KolmogorovTest( hb,"" ) << endl;

		//FIT DELLA PDF viene usata una distribuzione Gamma
		TF1 * pdf_funz = new TF1("pdf_funz", gammadist, min, max, 2);
		pdf_funz->SetParameter(0, 2.);
		pdf_funz->SetParameter(1, 10.);
		gStyle->SetOptFit(1111);
		pdf_funz->SetParName(0, "k");
		pdf_funz->SetParName(1, "theta");
		pdf_funz->SetLineColor(2);
		h[i][3]->Fit("pdf_funz", "Q");
		pdf_funz->SetLineColor(4);
		h[i][4]->Fit("pdf_funz", "Q");
		setstack( h[i][3], h[i][4], pdf );
		gPad->Update();
		TPaveStats * sb1 = (TPaveStats * )(h[i][3]->GetListOfFunctions()->FindObject("stats"));
		TPaveStats * sb2 = (TPaveStats * )(h[i][4]->GetListOfFunctions()->FindObject("stats"));
		sb1->SetX1NDC(.7);
		sb1->SetX2NDC(0.99);
		sb1->SetY1NDC(.4);
		sb1->SetY2NDC(.7);
		sb1->SetTextColor(2);
		sb2->SetX1NDC(.7);
		sb2->SetX2NDC(0.99);
		sb2->SetY1NDC(.7);
		sb2->SetY2NDC(1.0);
		sb2->SetTextColor(4);
		gPad->Modified();
		pdf->SetTitle("Probability distribution; ; pdf ");

		//FATTORI DI MERITO
		//-- h[5] Product between significance level (alpha) and power (beta) that must be minimized
		//-- h[6] Fattore di merito S/sqrt(B)
		//-- h[7] Signal purity: purity = (1-alpha)*N_s / ( (1-alpha)*N_s + beta*N_b )
    	h[i][5] = new TH1D(Form("S*B %d", i)," S_sx*B_dx", nbin, min, max);
		h[i][6] = new TH1D(Form("S/sqrt(B) %d", i),"T cut S_sx/sqrt(B_sx)", nbin, min, max);
		h[i][7] = new TH1D(Form("Purity %d", i),"Purity ", nbin, min, max);
    	for(int j = 1; j <= nbin; j++){
			double s_sx = h[i][3]->Integral(1,j, "width");
			double s_dx = h[i][3]->Integral(j, nbin, "width");
			double b_sx = h[i][4]->Integral(1,j, "width");
			double b_dx = h[i][4]->Integral(j, nbin, "width");
			double purity_s = s_dx*h[i][1]->Integral() / (s_dx*h[i][1]->Integral()+b_sx*h[i][2]->Integral());
	   		h[i][5]->SetBinContent(j, s_sx * b_dx);
			if ( b_sx != 0) h[i][7]->SetBinContent(j, s_sx / sqrt(b_sx) );
			h[i][6]->SetBinContent(j, purity_s);
		}
		for (int j = 0; j < 3; j++){
			c[i]->cd(2+j);
			gPad->SetGrid(1,1);
			gStyle->SetOptStat(0000);	
			h[i][5+j]->SetLineWidth(3);
			h[i][5+j]->SetLineColor(2);
			h[i][5+j]->Draw("C");
			if (j == 0 ) h[i][5]->SetTitle(" Significance level (S_sx)* Power (B_dx); ;  ");
    		if (j == 1 ) h[i][7]->SetTitle(" T cut S_sx/sqrt(B_sx); ; T cut ");
			if (j == 2 ) h[i][6]->SetTitle(" Signal Purity; ; Purity ");	
		}

	}
	return;
}

//Le grandezze di cui fare i grafici vanno passate sulla riga di comando.
//Per vedere i nomi con cui sono salvate nei TTree, consultare il file di lettura reading.cpp
//Il programma chiede di inserire il range di ogni variabile. Le grandezze hanno  circa i seguenti range:
//	var_phi: -pi, pi
//	var_eta: -4, 4
//	var_Et: 0, 350
//	var_pt:	0, 350
//	Ht, Htnu: 120, 800
//	Ptm: 0, 150

int main(int argc, char** argv){
    if (argc < 1){
        cout << "Usage: " << argv[0] << " q1 q2 q3 ... " << endl;
        return 1;
    }

	TApplication * Grafica = new TApplication("App", 0, 0);
	//LETTURA DEI TTree
    TFile * sinput = TFile::Open( "HH.root" );
    TFile * binput = TFile::Open( "ttbar.root" );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	
	int N = argc-1;
	TCanvas ** c = new TCanvas * [N];
	TH1D *** h = new  TH1D ** [N];
	string * var = new string [N]; 
	for( int i = 0; i < N; i++) var[i] = argv[i+1];
	DrawHisto(signal, background, var, c, h, N);
    Grafica->Run();
	return 0;
}
