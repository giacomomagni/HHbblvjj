
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TApplication.h"

using namespace std;
//cross sections
#define cr_s 33.5
#define cr_b 831760
//branching ratios
#define br_s 0.0717
#define br_b 0.238

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
	TLegend *legend = new TLegend(0.8,0.2,0.98,0.38);
    legend->AddEntry(hs,"Signal", "f");
    legend->AddEntry(hb,"Background", "f");
    legend->Draw("SAME");
	gPad->Update();
	TPaveStats * sb1 = (TPaveStats * )(hs->GetListOfFunctions()->FindObject("stats"));
	TPaveStats * sb2 = (TPaveStats * )(hb->GetListOfFunctions()->FindObject("stats"));
	sb1->SetX1NDC(.7);	sb1->SetX2NDC(0.99);	sb1->SetY1NDC(.4);	sb1->SetY2NDC(.7);
	sb2->SetX1NDC(.7);	sb2->SetX2NDC(0.99);	sb2->SetY1NDC(.7);	sb2->SetY2NDC(1.0);
	sb1->SetTextColor(2);	
	sb2->SetTextColor(4);
	gPad->Modified();
	return;
}

void setstack1(TH1D * hs, TH1D * hb, THStack * stack ){	
	gStyle->SetOptStat(0000);
	gStyle->SetOptFit(1111);
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
    legend1->AddEntry(hb,"Background + Signal", "f");
    legend1->Draw("SAME");
	gPad->Update();
	TPaveStats * sb = (TPaveStats * )(hb->GetListOfFunctions()->FindObject("stats"));
	sb->SetX1NDC(.55);	sb->SetX2NDC(0.98);	sb->SetY1NDC(.4);	sb->SetY2NDC(.9);
	sb->SetTextColor(1);
	gPad->Modified();
	return;
}

void toymc( string var  ){	

	gRandom->SetSeed(0);

	//INSERIRE RANGE in CUI GENERARE EVENTI 
	double min, max;
	int nbin, dof;
	cout << "Inserire range in cui generare i nuovi eventi di: " << var <<endl;
	cout << "Min: "; cin >> min;
	cout << "Max: "; cin >> max;
	cout << "Nbin: "; cin >> nbin; 
	cout << "Polinomial degree: "; cin >> dof;


	//ISTOGRAMMI DELLE DISTIBUZIONI	
    TFile * sinput = TFile::Open( "HH.root" );
    TFile * binput = TFile::Open( "ttbar.root" );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	string title1 = "Background_" + var;
	string title2 = "Signal_" + var;

	TCanvas * can = new TCanvas( "can" , "pdf lhe files");
	THStack * dist = new THStack("dist", var.c_str() );
	TH1D * hs = new TH1D( title2.c_str(), "Signal", nbin,min,max);
	TH1D * hb = new TH1D( title1.c_str(), "Background", nbin,min,max); 

	//Dal TTree del background vengono selezionati gli eventi con:  
	// -- mww < 125 [GeV]/c^2
	// -- 120 < mbb < 125 [GeV]/c^2
	//Gli eventi sono poi divisi in due categorie: mjj <= 50 e mjj > 50 [GeV]/c^2
	//Vengono studiati GLI EVENTI A MASSA mjj ALTA ( entambi i W on-Shell ) 
	//Viene poi fatto un istogramma della pdf di ogni variabile scelta
 
	TTreeReader reader_s( signal );
	TTreeReader reader_b( background );
	TTreeReaderValue<double> mbb(reader_b, "mww"); 
	TTreeReaderValue<double> mww(reader_b, "mbb"); 
	TTreeReaderValue<double> mjj_b(reader_b, "mjj");
	TTreeReaderValue<double> mjj_s(reader_s, "mjj");
	TTreeReaderValue<double>  variable_b( reader_b, var.c_str());
	TTreeReaderValue<double>  variable_s( reader_s, var.c_str());

	//contatore degli eventi salvati
	int tot_ev=0, s_ev=0, b_ev=0;
	//Selezione eventi backgrond e scrittura degli istogrammi 
	while ( reader_b.Next() ) {
		tot_ev++;
		if ( * mww <= 125. && * mbb >= 120. && * mbb <= 130. && * mjj_b >= 50.){
			hb->Fill( log(* variable_b) );
			b_ev++;
		}
	}
	double evb_saved = (double) b_ev/tot_ev;	
	tot_ev=0;
	//Scrittura degli istogrammi del segnale
	while ( reader_s.Next() ) {
		tot_ev += 1 ;	
		if ( * mjj_s >= 50.){
			hs->Fill( log(* variable_s) );
			s_ev ++;						
		}
	}
	double evs_saved = (double) s_ev/tot_ev;

	//NORMALIZZAZIONE PDF A 1
	double Ns = hs->Integral("width");
	double Nb = hb->Integral("width");
   	for(int j = 1; j <= nbin; j++){
   		hs->SetBinContent(j, hs->GetBinContent(j)/Ns);
		hb->SetBinContent(j, hb->GetBinContent(j)/Nb);
	}

	//FIT DELLA PDF viene usata una distribuzione polinomiale
	TF1 * pdf_funzs = new TF1("pdf_funzs", Form("abs(pol%i)", dof), min, max);
	TF1 * pdf_funzb = new TF1("pdf_funzb", Form("abs(pol%i)", dof), min, max );
	pdf_funzs->SetLineColor(2);
	pdf_funzb->SetLineColor(4);
	hs->Fit("pdf_funzs", "RQEMLI");
	hb->Fit("pdf_funzb", "RQEMLI");
	setstack( hs, hb, dist);

	//GENERAZIONE DEI NUOVI VALORI CON METODO MC
	// PLOT DEL NUOVO HISTO S+B al variariare della luminosità int
	// FIT con funzione S+B
	//GRAFICO ERRORE RELATIVO N_s(LUMINOSITY)

	double lum [4] = { 375., 750., 1500., 3000. };	
//	double lum [4] = { 100., 150., 300., 3000. };	

	TH1D ** bs = new TH1D * [4];
	TH1D ** s = new TH1D * [4];
	TCanvas ** c = new TCanvas * [4];
	TGraphErrors * gr = new TGraphErrors();
	TGraphErrors * gr_exp = new TGraphErrors();
	TMultiGraph *mg = new TMultiGraph();

	int dof1;
	cout << "degree S+B: "; cin>> dof1;

	for(int i=0; i< 4; i++){
		cout << "Luminosity: " << lum[i] << " fb^-1" << endl;

		string title = "Background + Signal at luminosity " + to_string(lum[i]) + " fb^-1";

		THStack * dist = new THStack("dist","Signal+Background");
		c[i] = new TCanvas( Form("c%i", i) , title.c_str() );
		bs[i] = new TH1D( title.c_str(), title.c_str() , nbin, min, max);
		s[i] = new TH1D( Form("signal, %i",i), Form("signal, %i",i) , nbin, min, max);

		double evs_exp = cr_s*br_s*lum[i]*evs_saved;
		double evb_exp = cr_b*br_b*lum[i]*evb_saved;

		bs[i]->FillRandom( "pdf_funzb", evb_exp );
		bs[i]->FillRandom( "pdf_funzs", evs_exp );
		s[i]->FillRandom( "pdf_funzs",  evs_exp );
		
		
		string funz = Form("[0]*abs(pol%i(2))+[1]*abs(pol%i(%i))", dof1, dof1, dof1+3);
		TF1 * funz_bs = new TF1("funz_bs", funz.c_str() , min, max );
		funz_bs->SetParameter(0, evs_exp);
		funz_bs->SetParameter(1, evb_exp);
		for (int j=2; j<dof1+3; j++) 			funz_bs->SetParameter(j, pdf_funzs->GetParameter(j-2));
		for (int j=dof1+3; j<2*dof1+4; j++) 	funz_bs->SetParameter(j, pdf_funzb->GetParameter(j-(dof1+3)));
		funz_bs->SetParName(0, "N_s");
		funz_bs->SetParName(1, "N_b");
		funz_bs->SetLineColor(4);
		bs[i]->Fit("funz_bs", "QRLM");
	
		setstack1( s[i], bs[i], dist );	
		c[i]->SetLogy();
		dist->SetTitle( title.c_str() );
		string xtitle = "log( " + var + " )";
		dist->GetXaxis()->SetTitle( xtitle.c_str() );
		dist->GetYaxis()->SetTitle("dN/dx");	

		gr->SetPoint(i, lum[i], funz_bs->GetParameter(0) );
		gr->SetPointError(i, 0.0, funz_bs->GetParError(0) );
		gr_exp->SetPoint(i, lum[i],  evs_exp );
		gr_exp->SetPointError(i, 0.0, 0.01*evs_exp );
		
		cout<< "Relative uncertainty on N_signal: " << funz_bs->GetParError(0)/funz_bs->GetParameter(0)*100 <<" %" <<endl;
	}

	//Disegno del grafico N_s vs luminosity
	TCanvas * c5 = new TCanvas("c5", "N_s vs luminosity");
	mg->SetTitle("N_signal vs Luminosity; Luminosity [fb^-1]; N_s");
	gPad->SetGrid(1,1);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.7);
	gr->Fit("pol1", "Q");
	gr_exp->SetMarkerStyle(2);
	gr_exp->SetMarkerSize(2);
	gr_exp->SetMarkerColor(4);
	mg->Add( gr_exp );
	mg->Add( gr );
	mg->Draw("AP");
	return;
}
