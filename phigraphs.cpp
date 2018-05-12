//Grafici di bb_phi , jj_lep_phi, bb_ jj+lep_deltaphi
//
//
// c++ -o phigraphs phigraphs.cpp `root-config --cflags --glibs`
//./phigraphs

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
	gPad->SetGrid(1,1);		
	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
    legend1->AddEntry(hs,"Signal", "f");
    legend1->AddEntry(hb,"Background", "f");
    legend1->Draw("SAME");
	return;
}

void DrawHisto ( TH1D * hs, TH1D * hb, THStack * stack, TCanvas * c, string var, int N_s, int N_b){
	c->Divide(2,2);
	c->cd(1);
	setstack( hs , hb , stack );
	int nbin = 50;
	string title0 = "hs_cdf_" + var;
	string title1 = "hb_cdf_" + var;
	string title2 = "S*B_" + var;
	string title3 = "S/sqrt(B)_" + var;
	string title4 = "Purity_" + var;
	//TEST DI KOLMOGROF 
	//poichè i dati sarebbero unbinned questo test è significativo solo se la larghezza dei bin è sufficientemente più piccola 
	//della precisione con cui si vuole ottenere in risultato 
	TH1D * h0 = new TH1D( title0.c_str(),"hs_cdf", nbin, 0, TMath::Pi() ); 
	TH1D * h1 = new TH1D( title1.c_str(),"hb_cdf", nbin, 0, TMath::Pi() ); 
    for(int j = 1; j <= nbin; j++){
    	h0->SetBinContent(j, hs->Integral(1,j,""));
		h1->SetBinContent(j, hb->Integral(1,j,""));
	}
	cout<< "KOLMOGROF TEST: probability that the two distribuions comes from the same phenomena:"  << h0->KolmogorovTest( h1,"" ) << endl;

	//FATTORI DI MERITO
	//-- h2 Product between significance level (alpha) and power (beta) that must be 0imized
	//-- h3 Fattore di merito S/sqrt(B)
	//-- h4 Signal purity: purity = (1-alpha)*N_s / ( (1-alpha)*N_s + beta*N_b )
	TH1D * h2 = new TH1D( title2.c_str()," S_sx*B_dx", nbin, 0, TMath::Pi());
	TH1D * h3 = new TH1D( title3.c_str(),"T cut S_sx/sqrt(B_sx)", nbin, 0, TMath::Pi());
	TH1D * h4 = new TH1D( title4.c_str(), "Purity ", nbin, 0, TMath::Pi());
	for(int j = 1; j <= nbin; j++){
		double s_sx = hs->Integral(1,j, "width");
		double s_dx = hs->Integral(j, nbin, "width");
		double b_sx = hb->Integral(1,j, "width");
		double b_dx = hb->Integral(j, nbin, "width");
		double purity_s = (double) s_dx*N_s / (s_dx*N_s+b_sx*N_b);
		h2->SetBinContent(j, s_sx * b_dx);
		if ( sqrt(b_sx) != 0 ) h3->SetBinContent(j, s_sx / sqrt(b_sx) );
		h4->SetBinContent(j, purity_s);
	}
	c->cd(2);
	gPad->SetGrid(1,1);
	gStyle->SetOptStat(0000);	
	h2->SetLineWidth(3);
	h2->SetLineColor(2);
	h2->Draw("C");
	h2->SetTitle(" Significance level (S_sx)* Power (B_dx); ;  ");
	c->cd(3);
	gPad->SetGrid(1,1);
	gStyle->SetOptStat(0000);	
	h4->SetLineWidth(3);
	h4->SetLineColor(2);
	h4->Draw("C");
	h4->SetTitle(" Signal Purity; ; Purity ");
	c->cd(4);
	gPad->SetGrid(1,1);
	gStyle->SetOptStat(0000);	
	h3->SetLineWidth(3);
	h3->SetLineColor(2);
	h3->Draw("C");
    h3->SetTitle(" T cut S_sx/sqrt(B_sx); ; T cut ");
	return;
}

int main(int argc, char** argv){
    if (argc < 1){
        cout << "Usage: " << argv[0] << " ./phigraphs " << endl;
        return 1;
    }

	TApplication * Grafica = new TApplication("App", 0, 0);
	//LETTURA DEI TTree
    TFile * sinput = TFile::Open( "HH.root" );
    TFile * binput = TFile::Open( "ttbar.root" );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	int N_s=0, N_b=0;

	TCanvas ** c = new TCanvas * [2];
	THStack * bbphi = new THStack("bbphi", "Delta phi b-bbar" );
	THStack * ljphi = new THStack("ljphi", "Delta phi lep-jj" );
	THStack * wbphi = new THStack("wbphi", "Delta phi ljj-bb" );
	int nbin = 50;
	TH1D * hs_bb = new TH1D("hs_bb","hs_bb", nbin, 0, TMath::Pi() ); 
	TH1D * hb_bb = new TH1D("hb_bb","hs_bb", nbin, 0, TMath::Pi() ); 
	TH1D * hs_lj = new TH1D("hs_lj","hs_lj", nbin, 0, TMath::Pi() ); 
	TH1D * hb_lj = new TH1D("hb_lj","hs_lj", nbin, 0, TMath::Pi() ); 
	TH1D * hs_wb = new TH1D("hs_wb","hs_wb", nbin, 0, TMath::Pi() ); 
	TH1D * hb_wb = new TH1D("hb_wb","hs_wb", nbin, 0, TMath::Pi() );
 
	//Dal TTree del background vengono selezionati gli eventi con:  
	// -- mww < 125 [GeV]/c^2
	// -- 120 < mbb < 125 [GeV]/c^2
	//Gli eventi sono poi divisi in due categorie: mvbs <= 50 e mvbs > 50 [GeV]/c^2
	//Vengono studiati gli eventi a mvbs alta.
	//Viene poi fatto un istogramma della pdf di ogni variabile scelta
 
	TTreeReader reader_s( signal );
	TTreeReader reader_b( background );

	TTreeReaderValue<double> mbb(reader_b, "mww"); 
	TTreeReaderValue<double> mww(reader_b, "mbb"); 
	TTreeReaderValue<double> mvbs_b(reader_b, "mvbs");
	TTreeReaderValue<double>  b1_phi_b( reader_b, "b1_phi");
	TTreeReaderValue<double>  b2_phi_b( reader_b, "b2_phi");
	TTreeReaderValue<double>  vbs1_phi_b( reader_b, "vbs1_phi");
	TTreeReaderValue<double>  vbs2_phi_b( reader_b, "vbs2_phi");
	TTreeReaderValue<double>  lep_phi_b( reader_b, "lep_phi");

	TTreeReaderValue<double> mvbs_s(reader_s, "mvbs");
	TTreeReaderValue<double>  b1_phi_s( reader_s, "b1_phi");
	TTreeReaderValue<double>  b2_phi_s( reader_s, "b2_phi");
	TTreeReaderValue<double>  vbs1_phi_s( reader_s, "vbs1_phi");
	TTreeReaderValue<double>  vbs2_phi_s( reader_s, "vbs2_phi");
	TTreeReaderValue<double>  lep_phi_s( reader_s, "lep_phi");
	

	//Selezione eventi backgrond e scrittura degli istogrammi 
	while ( reader_b.Next() ) {
		if ( * mww <= 125. && * mbb >= 120. && * mbb <= 130. && * mvbs_b >= 50.){				
				double deltabb = abs( * b1_phi_b - *  b2_phi_b );
				double deltalj = abs( ( * vbs1_phi_b + *  vbs2_phi_b ) * 0.5 - * lep_phi_b );
				double deltawb = abs( (( * vbs1_phi_b + *  vbs2_phi_b ) * 0.5 + * lep_phi_b )*0.5 - (* b1_phi_b + *  b2_phi_b)*0.5 );
				if (deltabb > TMath::Pi()) deltabb = TMath::TwoPi()-deltabb;
				if (deltalj > TMath::Pi()) deltalj = TMath::TwoPi()-deltalj;
				if (deltawb > TMath::Pi()) deltawb = TMath::TwoPi()-deltawb;
				hb_bb->Fill( deltabb );
				hb_lj->Fill( deltalj );
				hb_wb->Fill( deltawb );
				N_b++;
		}
	}	
	//Scrittura degli istogrammi del segnale
	while ( reader_s.Next() ) {
		if ( * mvbs_s >= 50.){
			double deltabb = abs( * b1_phi_s - *  b2_phi_s );
			double deltalj = abs( ( * vbs1_phi_s + *  vbs2_phi_s ) * 0.5 - * lep_phi_s );
			double deltawb = abs( (( * vbs1_phi_s + *  vbs2_phi_s ) * 0.5 + * lep_phi_s )*0.5 - (* b1_phi_s + *  b2_phi_s)*0.5 );
			if (deltabb > TMath::Pi()) deltabb = TMath::TwoPi()-deltabb;
			if (deltalj > TMath::Pi()) deltalj = TMath::TwoPi()-deltalj;
			if (deltawb > TMath::Pi()) deltawb = TMath::TwoPi()-deltawb;
			hs_bb->Fill( deltabb );
			hs_lj->Fill( deltalj );	
			hs_wb->Fill( deltawb );
			N_s++;		
		}
	}
	for(int j = 1; j <= nbin; j++){
   		hb_bb->SetBinContent(j, hb_bb->GetBinContent(j) / (N_b*(TMath::Pi()/nbin)) );	
		hs_bb->SetBinContent(j, hs_bb->GetBinContent(j) / (N_s*(TMath::Pi()/nbin)) );
		hb_lj->SetBinContent(j, hb_lj->GetBinContent(j) / (N_b*(TMath::Pi()/nbin)) );	
		hs_lj->SetBinContent(j, hs_lj->GetBinContent(j) / (N_s*(TMath::Pi()/nbin)) );
		hb_wb->SetBinContent(j, hb_wb->GetBinContent(j) / (N_b*(TMath::Pi()/nbin)) );	
		hs_wb->SetBinContent(j, hs_wb->GetBinContent(j) / (N_s*(TMath::Pi()/nbin)) );
	}
	
	//ANALISI di b-bbar delta phi 
	c[0] = new TCanvas("c0", "bb deltaPhi");
	DrawHisto( hs_bb, hb_bb, bbphi, c[0], "deltaPhi_bb", N_s, N_b);
	bbphi->SetTitle("Delta Phi b-bbar; Phi [rad]; dN/dphi");
		
	//ANALISI di lep-vbs delta phi 	
	c[1] = new TCanvas("c1", "lep-jj deltaPhi");
	DrawHisto( hs_lj, hb_lj, ljphi, c[1], "deltaPhi_ljj", N_s, N_b);
	ljphi->SetTitle("Delta Phi lep-jj; Phi [rad]; dN/dphi");

	//ANALISI di bb-lepjj delta phi 
	c[2] = new TCanvas("c2", "bb-ljj deltaPhi");
	DrawHisto( hs_wb, hb_wb, wbphi, c[2], "deltaPhi_bb-ljj", N_s, N_b);
	wbphi->SetTitle("Delta Phi bb-ljj; Phi [rad]; dN/dphi");

    Grafica->Run();
	return 0;
}
