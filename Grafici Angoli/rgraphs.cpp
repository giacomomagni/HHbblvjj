//Grafici di bb_deltar,  jj_lep_deltar, bb_ jj+lep_deltar,
//
//
// c++ -o rgraphs rgraphs.cpp `root-config --cflags --glibs`
//./rgraphs

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
#define nbin 50

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
        cout << "Usage: " << argv[0] << " ./rgraphs " << endl;
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
	THStack * bbr = new THStack("bbr", "Delta R b-bbar" );
	THStack * ljr = new THStack("ljr", "Delta R lep-jj" );
	THStack * wbr = new THStack("wbr", "Delta R ljj-bb" );
	TH1D * hs_bb = new TH1D("hs_bb","hs_bb", nbin, 0, 2*TMath::Pi() ); 
	TH1D * hb_bb = new TH1D("hb_bb","hs_bb", nbin, 0, 2*TMath::Pi() ); 
	TH1D * hs_lj = new TH1D("hs_lj","hs_lj", nbin, 0, 2*TMath::Pi() ); 
	TH1D * hb_lj = new TH1D("hb_lj","hs_lj", nbin, 0, 2*TMath::Pi() ); 
	TH1D * hs_wb = new TH1D("hs_wb","hs_wb", nbin, 0, 2*TMath::Pi() ); 
	TH1D * hb_wb = new TH1D("hb_wb","hs_wb", nbin, 0, 2*TMath::Pi() ); 

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
	
	

	//Selezione eventi backgrond e scrittura degli istogrammi 
	while ( reader_b.Next() ) {
		if ( * mww <= 125. && * mbb >= 120. && * mbb <= 130. && * mvbs_b >= 50.){				
			double deltaphi_bb = abs( * b1_phi_b - *  b2_phi_b );
			double deltaphi_lj = abs( ( * vbs1_phi_b + *  vbs2_phi_b ) * 0.5 - * lep_phi_b );
			double deltaphi_wb = abs( (( * vbs1_phi_b + *  vbs2_phi_b ) * 0.5 + * lep_phi_b )*0.5 - (* b1_phi_b + *  b2_phi_b)*0.5 );
			double deltaeta_bb = abs( * b1_eta_b - *  b2_eta_b );
			double deltaeta_lj = abs( ( * vbs1_eta_b + *  vbs2_eta_b ) * 0.5 - * lep_eta_b );
			double deltaeta_wb = abs( (( * vbs1_eta_b + *  vbs2_eta_b ) * 0.5 + * lep_eta_b )*0.5 - (* b1_eta_b + *  b2_eta_b)*0.5 );
			if (deltaphi_bb > TMath::Pi()) deltaphi_bb = TMath::TwoPi()-deltaphi_bb;
			if (deltaphi_lj > TMath::Pi()) deltaphi_lj = TMath::TwoPi()-deltaphi_lj;
			if (deltaphi_wb > TMath::Pi()) deltaphi_wb = TMath::TwoPi()-deltaphi_wb;
			double deltar_bb = sqrt( pow(deltaphi_bb, 2) + pow(deltaeta_bb, 2) );
			double deltar_lj = sqrt( pow(deltaphi_lj, 2) + pow(deltaeta_lj, 2) );
			double deltar_wb = sqrt( pow(deltaphi_wb, 2) + pow(deltaeta_wb, 2) );
			hb_bb->Fill( deltar_bb );
			hb_lj->Fill( deltar_lj );
			hb_wb->Fill( deltar_wb );
			N_b++;		
		}
	}	
	//Scrittura degli istogrammi del segnale
	while ( reader_s.Next() ) {
		if ( * mvbs_s >= 50.){
			double deltaphi_bb = abs( * b1_phi_s - *  b2_phi_s );
			double deltaphi_lj = abs( ( * vbs1_phi_s + *  vbs2_phi_s ) * 0.5 - * lep_phi_s );
			double deltaphi_wb = abs( (( * vbs1_phi_s + *  vbs2_phi_s ) * 0.5 + * lep_phi_s )*0.5 - (* b1_phi_s + *  b2_phi_s)*0.5 );
			double deltaeta_bb = abs( * b1_eta_s - *  b2_eta_s );
			double deltaeta_lj = abs( ( * vbs1_eta_s + *  vbs2_eta_s ) * 0.5 - * lep_eta_s );
			double deltaeta_wb = abs( (( * vbs1_eta_s + *  vbs2_eta_s ) * 0.5 + * lep_eta_s )*0.5 - (* b1_eta_s + *  b2_eta_s)*0.5 );
			if (deltaphi_bb >= TMath::Pi()) deltaphi_bb = TMath::TwoPi()-deltaphi_bb;
			if (deltaphi_lj >= TMath::Pi()) deltaphi_lj = TMath::TwoPi()-deltaphi_lj;
			if (deltaphi_wb >= TMath::Pi()) deltaphi_wb = TMath::TwoPi()-deltaphi_wb;
			double deltar_bb = sqrt( pow(deltaphi_bb, 2) + pow(deltaeta_bb, 2) );
			double deltar_lj = sqrt( pow(deltaphi_lj, 2) + pow(deltaeta_lj, 2) );
			double deltar_wb = sqrt( pow(deltaphi_wb, 2) + pow(deltaeta_wb, 2) );
			hs_bb->Fill( deltar_bb );
			hs_lj->Fill( deltar_lj );
			hs_wb->Fill( deltar_wb );
			N_s++;	
		}
	}
	for(int j = 1; j <= nbin; j++){
   		hb_bb->SetBinContent(j, hb_bb->GetBinContent(j) / (N_b*(2*TMath::Pi()/nbin)) );	
		hs_bb->SetBinContent(j, hs_bb->GetBinContent(j) / (N_s*(2*TMath::Pi()/nbin)) );
		hb_lj->SetBinContent(j, hb_lj->GetBinContent(j) / (N_b*(2*TMath::Pi()/nbin)) );	
		hs_lj->SetBinContent(j, hs_lj->GetBinContent(j) / (N_s*(2*TMath::Pi()/nbin)) );
		hb_wb->SetBinContent(j, hb_wb->GetBinContent(j) / (N_b*(2*TMath::Pi()/nbin)) );	
		hs_wb->SetBinContent(j, hs_wb->GetBinContent(j) / (N_s*(2*TMath::Pi()/nbin)) );
	}
	
	//ANALISI di b-bbar delta R 
	c[0] = new TCanvas("c0", "b-bbar deltaR");
	DrawHisto( hs_bb, hb_bb, bbr, c[0], "deltaR_bb", N_s, N_b);
	bbr->SetTitle("Delta R b-bbar; R [rad]; dN/dR");
	
	//ANALISI di lep-jj delta R 	
	c[1] = new TCanvas("c1", "lep-jj deltaR");
	DrawHisto( hs_lj, hb_lj, ljr, c[1], "deltaR_lj", N_s, N_b);
	ljr->SetTitle("Delta R lep-jj; R [rad]; dN/dR");
		
	//ANALISI di bb-lepjj delta R 	
	c[2] = new TCanvas("c2", "bb-ljj deltaR");
	DrawHisto( hs_wb, hb_wb, wbr, c[2], "deltaR_bbljj", N_s, N_b);
	wbr->SetTitle("Delta R bb-lepjj; R [rad]; dN/dR");

    Grafica->Run();
	return 0;
}
