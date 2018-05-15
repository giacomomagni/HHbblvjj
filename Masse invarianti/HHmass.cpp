//HH invariant mass.
//
//
// c++ -o HHmass HHmass.cpp `root-config --cflags --glibs`
//./HHmass HH.root ttbar.root

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
#include <TLegend.h>
#include <THStack.h>
#include <TApplication.h>
#include <fstream>

using namespace std;

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
	gStyle->SetOptStat(1111);	
	gPad->SetGrid(1,1);
	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
    legend1->AddEntry(hs,"Signal", "f");
    legend1->AddEntry(hb,"Background", "f");
    legend1->Draw("SAME");
	return;
}

void DrawHisto	(TTree * signal, TTree * background, string var){	
	
	TTreeReader reader_s( signal );
	TTreeReader reader_b( background );

	TTreeReaderValue<double> mbb_b(reader_b, "mww"); 
	TTreeReaderValue<double> mww_b(reader_b, "mbb"); 
	TTreeReaderValue<double> mbb_s(reader_s, "mww"); 
	TTreeReaderValue<double> mww_s(reader_s, "mbb"); 

	int nbin = 100;
	double min = 0, max = 800;
	string title2 = "Background_" + var;
	string title1 = "Signal_" + var;
	TCanvas *can = new TCanvas("can", "HH invariant mass");
	THStack * dist = new THStack("dist", "Distribution histo");
	TH1D * hs = new TH1D( title1.c_str(), "Signal", nbin,min,max);
 	TH1D *hb = new TH1D( title2.c_str(), "Background", nbin,min,max); 
	
	//SCRITTURA ISTOGRAMMA DISTRIBUZIONE
	while ( reader_s.Next() ) {
		hs->Fill( *mbb_s + *mww_s );
	}
	while ( reader_b.Next() ) {
		hb->Fill( *mbb_b + *mww_b );
	}
	
	for( int j=1; j<=nbin; j++ ){
		hb->SetBinContent(j, hb->GetBinContent(j) / (hb->GetEntries()*(max-min)/nbin) );
		hs->SetBinContent(j, hs->GetBinContent(j) / (hs->GetEntries()*(max-min)/nbin) );		
	}
	setstack( hs, hb, dist );
	dist->SetTitle(" HH invariant mass; Mass [GeV/c^2]; dN/dm");
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
	
	DrawHisto(signal, background, "mHH");
	
    Grafica->Run();
	return 0;
}
