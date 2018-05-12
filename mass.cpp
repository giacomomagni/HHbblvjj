//plot of invatiant masses.
//
//
// c++ -o mass mass.cpp `root-config --cflags --glibs`
//./mass HH.root ttbar.root

#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <THStack.h>
#include <TApplication.h>
#include <fstream>

using namespace std;

void setstack(TH1F * hs, TH1F * hb, THStack * stack ){	
	hb->SetFillStyle(3003);
	hs->SetFillStyle(3004);
    hs->SetFillColor(2);
	hb->SetFillColor(4);
	hs->SetLineColor(2);
	hb->SetLineColor(4);	
	stack->Add(hs, "S");
	stack->Add(hb, "S");
    stack->Draw("nostack");
	gStyle->SetOptStat(1111);	
	gPad->SetGrid(1,1);
	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
    legend1->AddEntry(hs,"Signal", "f");
    legend1->AddEntry(hb,"Background", "f");
    legend1->Draw("SAME");
	return;
}

void DrawHisto	(TTree * signal, TTree * background, string * var, int nbin, float min, float max, int N ){	
		TH1F *** h = new TH1F ** [N];		
		TCanvas * c1 = new TCanvas("c", "Garph_mass" );	
		c1->Divide(N, 1);
	
	for(int i = 0; i < N; i++){
		h[i] = new TH1F * [2];
		
		//SCRITTURA ISTOGRAMMA DISTRIBUZIONE
		c1->cd(i+1);
		string title1 = "Background_" + var[i];
		string title2 = "Signal_" + var[i];
		string var1 = var[i] + ">>" + title1;
		string var2 = var[i] + ">>" + title2;
		THStack * dist = new THStack("dist", "Distribution histo");
   		h[i][1] = new TH1F( title2.c_str(), "Signal", nbin,min,max);
 		h[i][2] = new TH1F( title1.c_str(), "Background", nbin,min,max); 
   		signal->Draw( var2.c_str(), "" );
    	background->Draw( var1.c_str(), "" );	
		for( int j=1; j<=nbin; j++ ){
			h[i][1]->SetBinContent(j, h[i][1]->GetBinContent(j) / (h[i][1]->GetEntries()*(max-min)/nbin) );
			h[i][2]->SetBinContent(j, h[i][2]->GetBinContent(j) / (h[i][2]->GetEntries()*(max-min)/nbin) );		
		}
		setstack( h[i][1], h[i][2], dist );
		if ( i == 0 )	dist->SetTitle(" b bbar invariant mass; Mass [GeV/c^2]; dN/dm");
		if ( i == 1 )	dist->SetTitle(" jet components invariant mass;  Mass [GeV/c^2]; dN/dm");
		if ( i == 2 )	dist->SetTitle(" Visible WW invariant mass ;  Mass [GeV/c^2]; dN/dm");
	}
	return;
}

int main(int argc, char** argv){
    if (argc < 3){
        cout << "Usage: " << argv[0] << " HH.root ttbar.root " << endl;
        return 1;
    } 

	TApplication * Grafica = new TApplication("App", 0, 0);
	//LETTURA DEI TTree
    TFile * sinput = TFile::Open( argv[1] );
    TFile * binput = TFile::Open( argv[2] );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	
	//Var contiene le variabili lette
	//N = numero di grandezze prese in cosiderazione
	int N = 3;
	string * var = new string [N] {  
				"mbb",  "mvbs", "mww"
	}; 

	//MASS  
	int nbin = 100;
	double min = 0, max = 300;
	DrawHisto(signal, background, var, nbin, min, max, N);
	
    Grafica->Run();
	return 0;
}
