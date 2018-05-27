
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TApplication.h"
#include "TRandom3.h"

using namespace std;
//cross sections
#define cr_s 33.5
#define cr_b 831760
//branching ratios
#define br_s 0.0717
#define br_b 0.238

void setstack(TH1D * hs, TH1D * hb, THStack * stack ){
	gStyle->SetOptFit(0000);
	gStyle->SetOptStat(0000);	
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
	TLegend *legend = new TLegend(0.8,0.7,0.98,0.88);
    legend->AddEntry(hs,"Signal", "f");
    legend->AddEntry(hb,"Background", "f");
    legend->Draw("SAME");
	gPad->Update();
//	TPaveStats * sb1 = (TPaveStats * )(hs->GetListOfFunctions()->FindObject("stats"));
//	TPaveStats * sb2 = (TPaveStats * )(hb->GetListOfFunctions()->FindObject("stats"));
//	sb1->SetX1NDC(.7);	sb1->SetX2NDC(0.99);	sb1->SetY1NDC(.4);	sb1->SetY2NDC(.7);
//	sb2->SetX1NDC(.7);	sb2->SetX2NDC(0.99);	sb2->SetY1NDC(.7);	sb2->SetY2NDC(1.0);
//	sb1->SetTextColor(2);	
//	sb2->SetTextColor(4);
	gPad->Modified();
	return;
}

void toymc( string var  ){	
	
	if(gRandom) delete gRandom;
	gRandom = new TRandom3(0);

	int nbin, dof;
	double min, max;
	bool correct = 1;
	bool logscale = 1;

	//INSERIRE RANGE in CUI GENERARE EVENTI 
	cout << "------------------------------------------------------------" << endl
		 << "Type the MC simulation range of the variable: " << var <<endl;
	cout << "Log scale? [Yes type 1, No type 0]: "; cin >> logscale; 
	cout << "Min: "; cin >> min;
	cout << "Max: "; cin >> max;
	cout << "Nbin: "; cin >> nbin; 

	//ISTOGRAMMI DELLE DISTIBUZIONI	
    TFile * sinput = TFile::Open( "HH.root" );
    TFile * binput = TFile::Open( "newttbar.root" );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	string title0 = "Pdf distribution of " + var;
	string title1 = "Background_" + var;
	string title2 = "Signal_" + var;

	TCanvas * can = new TCanvas( "can" , title0.c_str() );
	THStack * dist = new THStack("dist", title0.c_str() );	
	TH1D * hb = new TH1D( title1.c_str(), "Background", nbin,min,max);
	TH1D * hs = new TH1D( title2.c_str(), "Signal", nbin,min,max); 

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
	TTreeReaderValue<double> deltaphi_bbljj_b(reader_b, "deltaphi_bbljj");
	TTreeReaderValue<double> deltaphi_bbljj_s(reader_s, "deltaphi_bbljj");
	TTreeReaderValue<double>  variable_b( reader_b, var.c_str());
	TTreeReaderValue<double>  variable_s( reader_s, var.c_str());

	//contatore degli eventi salvati
	int tot_ev=0, s_ev=0, b_ev=0;
	//Selezione eventi backgrond e scrittura degli istogrammi 
	while ( reader_b.Next() ) {
		tot_ev++;
		if ( * mww <= 125. && * mbb >= 124. && * mbb <= 126. /*&& abs(* deltaphi_bbljj_b) >= 2.*/){
			if(logscale == 1){ 
				hb->Fill( log(* variable_b) );
			}else{
				hb->Fill( * variable_b );
			}
			b_ev++;
		}
	}
	double evb_saved = (double) b_ev/tot_ev;	
	tot_ev=0;
	//Scrittura degli istogrammi del segnale
	while ( reader_s.Next() ) {
		tot_ev += 1 ;	
		//if ( abs(* deltaphi_bbljj_s)>=2 ){
			if(logscale == 1){ 
				hs->Fill( log(* variable_s) );
			}else{
				hs->Fill( * variable_s );
			}
			s_ev ++;						
	//	}
	}
	double evs_saved = (double) s_ev/tot_ev;

	//NORMALIZZAZIONE PDF A 1
	//hs->Smooth(100);
	//hb->Smooth(100);
	double Ns = hs->Integral( );
	double Nb = hb->Integral( );
   	for(int j = 1; j <= nbin; j++){
   		hs->SetBinContent(j, hs->GetBinContent(j)/Ns);
		hb->SetBinContent(j, hb->GetBinContent(j)/Nb);
	}


	//FIT DELLA PDF viene usata una distribuzione polinomiale
	cout << "Fit polinomial degree: "; cin >> dof;
	TF1 * pdf_funzs = new TF1("pdf_funzs", Form("pol%i", dof), min, max);
	TF1 * pdf_funzb = new TF1("pdf_funzb", Form("pol%i", dof), min, max );
	pdf_funzs->SetLineColor(2);
	pdf_funzb->SetLineColor(4);
	hs->Fit("pdf_funzs", "LRME");
	hb->Fit("pdf_funzb", "LRME");
	setstack( hs, hb, dist);
	if(logscale == 1){
		string xtitle = "log( " + var + " )";
		dist->GetXaxis()->SetTitle( xtitle.c_str() );
	}
	if(logscale == 0)		dist->GetXaxis()->SetTitle( var.c_str() );
	dist->GetYaxis()->SetTitle("dN/dx");	

	string title3 = var + ".png";
	can->Print( title3.c_str() );

	cout<< "Cut on Background:" << evb_saved << endl;
	cout<< "Cut on Signal: " << evs_saved << endl;

	cout<< "Is the fit correct? [type 0 to exit, 1 to continue]: ";	cin >> correct;
	if( correct == 0) return;		

	sinput->Close();
	binput->Close();

	//GENERAZIONE DEI NUOVI VALORI CON METODO MC
	//GRAFICO ERRORE RELATIVO N_s vs LUMINOSITY

	double lum [4] = { 375., 750., 1500., 3000. };	
//	double lum [4] = { 100., 150., 300., 3000. };	

	string fileoutname = var + " simulation.root";
	TFile * fileout = new TFile( fileoutname.c_str() ,"recreate");
	TCanvas * can2 = new TCanvas("can2", "can");

	// FIT con funzione S+B PER UN MUMERO DI Nsim VOLTE PER CIASCUNA LUMINOSITA'
	for(int i=0; i<4; i++){

		int evs_exp = (int) cr_s*br_s*lum[i]*evs_saved;
		int evb_exp = (int) cr_b*br_b*lum[i]*evb_saved;		
		//int evb_exp = 100*evs_exp;		
		int Nsim=100;
		cout<< "------------------------------------------------------------" << endl
			<< "Luminosity: " << lum[i] << " fb^-1" << endl;
//			<< "Type the number of MC simulation: "; cin >> Nsim;

		string nametree = "tree_" + to_string(lum[i]);
		TTree * treesim = new TTree( nametree.c_str(), nametree.c_str());

		double ns;
		double nb;
		double nbs;
		
		treesim->Branch("ns",&ns);
		treesim->Branch("nb",&nb);
		treesim->Branch("nbs",&nbs);
	
		for(int k=1; k <= Nsim; k++ ){
			cout<<"I'm woking on simulation number: " << k <<endl;
			string title4 = "Background + Signal at luminosity " + to_string(lum[i]) + " fb^-1" + to_string(k);	
			TH1D * bs = new TH1D(title4.c_str(), title4.c_str(), nbin, min, max );

			gRandom->SetSeed(0);
			//cout << evs_exp <<endl;
			//cout<< evb_exp << endl;
			TRandom3 r;
			bs->FillRandom( "pdf_funzs", r.Poisson(evs_exp));	
			bs->FillRandom( "pdf_funzb", r.Poisson(evb_exp));

			string funz = Form("[0]*pol%i(2)+[1]*pol%i(%i)", dof, dof, dof+3);
			TF1 * funz_bs = new TF1("funz_bs", funz.c_str(), min, max );
			funz_bs->SetParameter(0, evb_exp);
			funz_bs->SetParameter(1, evs_exp);


			funz_bs->SetParLimits(0,  evb_exp-5*sqrt(evb_exp), evb_exp+5*sqrt(evb_exp) );
			//funz_bs->SetParLimits(1,  evs_exp-10*sqrt(evs_exp), evs_exp+10*sqrt(evs_exp) );
			//funz_bs->SetParLimits(1, 0, 2*evs_exp);
			//funz_bs->SetParLimits(0, 0, 2*evb_exp);

			for (int j=2; j<dof+3; j++) 		funz_bs->FixParameter(j, pdf_funzb->GetParameter(j-2));
			for (int j=dof+3; j<2*dof+4; j++) 	funz_bs->FixParameter(j, pdf_funzs->GetParameter(j-(dof+3)));
			
			bs->Fit("funz_bs", "LQMERB");
			gPad->Update();
			gStyle->SetOptFit(1111);
			bs->Draw("");
			gPad->Modified();
			//cout<< funz_bs->GetParameter(1) <<endl;
			//cout<< funz_bs->GetParameter(0) <<endl;
			ns = funz_bs->GetParameter(1);
			nb = funz_bs->GetParameter(0);
			nbs = ns + nb;
			treesim->Fill();	
		}
		treesim->Write();
	}
	fileout->Close();
	return;
}
