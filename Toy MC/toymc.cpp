
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

	double binsize = abs(max-min)/nbin;

//	min=0; max=4.5; nbin=67; dof=8;
	//ISTOGRAMMI DELLE DISTIBUZIONI	
    TFile * sinput = TFile::Open( "HH.root" );
    TFile * binput = TFile::Open( "ttbar.root" );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	string title0 = "Pdf distribution of " + var;
	string title1 = "Background_" + var;
	string title2 = "Signal_" + var;

	TCanvas * can = new TCanvas( "can" , title0.c_str() );
	THStack * dist = new THStack("dist", title0.c_str() );
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
			if(logscale == 1){ 
				hb->Fill( log(* variable_b) );
			}else{
				hb->Fill( * variable_b);
			}
			b_ev++;
		}
	}
	double evb_saved = (double) b_ev/tot_ev;	
	tot_ev=0;
	//Scrittura degli istogrammi del segnale
	while ( reader_s.Next() ) {
		tot_ev += 1 ;	
		if ( * mjj_s >= 50.){
			if(logscale == 1){ 
				hs->Fill( log(* variable_s) );
			}else{
				hs->Fill( * variable_s);
			}
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
//	cout<<hs->Integral() << endl;
	//FIT DELLA PDF viene usata una distribuzione polinomiale
	cout << "Fit polinomial degree: "; cin >> dof;
	TF1 * pdf_funzs = new TF1("pdf_funzs", Form("abs(pol%i)", dof), min, max);
	TF1 * pdf_funzb = new TF1("pdf_funzb", Form("abs(pol%i)", dof), min, max );
	pdf_funzs->SetLineColor(2);
	pdf_funzb->SetLineColor(4);
	hs->Fit("pdf_funzs", "RQLEM");
	hb->Fit("pdf_funzb", "RQLEM");
	setstack( hs, hb, dist);
	dist->GetXaxis()->SetTitle( var.c_str() );
	dist->GetYaxis()->SetTitle("dN/dx");	
	can->Modified();
	can->Update();
	string name = var + ".png";
//	can->Print( name.c_str() );

//	cout<< "Is the fit correct? [type 0 to exit, 1 to continue]: ";	cin >> correct;
//	if( correct == 0) return;		

	//GENERAZIONE DEI NUOVI VALORI CON METODO MC
	//GRAFICO ERRORE RELATIVO N_s vs LUMINOSITY

	double lum [4] = { 375., 750., 1500., 3000. };	
//	double lum [4] = { 100., 150., 300., 3000. };	

	TCanvas * can1 = new TCanvas("can1", "N_s simulated distribution" );
	can1->Divide(1,4);
	TGraphErrors * gr = new TGraphErrors();
	TGraphErrors * gr_b = new TGraphErrors();
	TGraphErrors * gr_exp = new TGraphErrors();
	TGraphErrors * gr_exp_b = new TGraphErrors();
	TMultiGraph *mg = new TMultiGraph();
	TMultiGraph *mg_b = new TMultiGraph();

	// FIT con funzione S+B PER UN MUMERO DI Nsim VOLTE PER CIASCUNA LUMINOSITA'
	for(int i=0; i<4; i++){

		int evs_exp = (int) cr_s*br_s*lum[i]*evs_saved;
		int evb_exp = (int) cr_b*br_b*lum[i]*evb_saved;
		int Nsim;
		cout<< "------------------------------------------------------------" << endl
			<< "Luminosity: " << lum[i] << " fb^-1" << endl
			<< "Type the number of MC simulation: "; cin >> Nsim;

		TVectorT<double> ns(Nsim);
		TVectorT<double> nb(Nsim);
		TF1 * f = new TF1("f", "gaus", 0, 3*evs_exp);
		
		for(int k=1; k <= Nsim; k++ ){
			cout<<"I'm woking on simulation number: " << k <<endl;
			string title5 = "Background + Signal at luminosity " + to_string(lum[i]) + " fb^-1" + to_string(k);	
			TH1D * bs = new TH1D(title5.c_str(), title5.c_str() , nbin, min, max );
	
			gRandom->SetSeed(0);
			//cout << evs_exp <<endl;
			//cout<< evb_exp << endl;
			bs->FillRandom( "pdf_funzs", evs_exp );	
			bs->FillRandom( "pdf_funzb", evb_exp );

			string funz = Form("[0]*abs(pol%i(2))+[1]*abs(pol%i(%i))", dof, dof, dof+3);
			TF1 * funz_bs = new TF1("funz_bs", funz.c_str() , min, max );
			funz_bs->SetParameter(0, evb_exp);
			funz_bs->SetParameter(1, evs_exp);

			for (int j=2; j<dof+3; j++) 		funz_bs->FixParameter(j, pdf_funzb->GetParameter(j-2));
			for (int j=dof+3; j<2*dof+4; j++) 	funz_bs->FixParameter(j, pdf_funzs->GetParameter(j-(dof+3)));
			
			bs->Fit("funz_bs", "QRLM");
			//cout<< funz_bs->GetParameter(1)/binsize <<endl;
			//cout<< funz_bs->GetParameter(0)/binsize <<endl;
			ns[k-1] = funz_bs->GetParameter(1)/binsize;
			nb[k-1] = funz_bs->GetParameter(0)/binsize;
		}

		
		string title4 = "Signal events simulated at luminosity " + to_string(lum[i]) + " fb^-1";
		string title6 = "Background events simulated at luminosity " + to_string(lum[i]) + " fb^-1";
		TH1D * histo_Ns = new TH1D( title4.c_str(), title4.c_str(), Nsim, ns.Min(), ns.Max());
		TH1D * histo_Nb = new TH1D( title6.c_str(), title6.c_str(), Nsim, nb.Min(), nb.Max());
		for(int k=0; k< Nsim; k++){
			histo_Ns->Fill( ns[k] );
			histo_Nb->Fill( nb[k] );		
		}
		histo_Ns->GetXaxis()->SetTitle( "Simulated signal events" );
		histo_Ns->GetYaxis()->SetTitle(	"Frequency [Events]");
		histo_Ns->Fit("f", "Q");
		double Ns_avg = f->GetParameter(1);
		double Ns_avg_err = f->GetParameter(2);
		
		can1->cd(i+1);	 
		histo_Ns->Draw();
		can1->Modified();
		can1->Update();

		histo_Nb->Fit("f", "QN");
		double Nb_avg = f->GetParameter(1);
		double Nb_avg_err = f->GetParameter(2);
		
		cout << "------------------------------------------------------------" << endl
			 << "Expected signal events: " << evs_exp << endl
			 << "Expected background events: " << evb_exp << endl
			 << "Simulated signal events on average: " << Ns_avg <<  endl
			 << "Relative error: " << Ns_avg_err/Ns_avg*100 << "%" << endl
			 << "Simulated bg events on average: " << Nb_avg <<  endl;

		gr->SetPoint(i, lum[i], Ns_avg );
		gr->SetPointError(i, 0.0, Ns_avg_err );
		gr_b->SetPoint(i, lum[i], Nb_avg );
		gr_b->SetPointError(i, 0.0, Nb_avg_err );
		gr_exp->SetPoint(i, lum[i],  evs_exp );
		gr_exp->SetPointError(i, 0.0, 0.01*evs_exp );
		gr_exp_b->SetPoint(i, lum[i],  evb_exp );
		gr_exp_b->SetPointError(i, 0.0, 0.01*evb_exp );
		
		//cout<< "Relative uncertainty on N_signal: " << funz_bs->GetParError(0)/funz_bs->GetParameter(0)*100 <<" %" <<endl;
	}
	string title6 = var + " Simulated Ns_frequency.png";
//	can1->Print( title6.c_str() );

	//Disegno del grafico N_s vs luminosity
	TCanvas * can2 = new TCanvas("can2", "N_s vs luminosity");
	mg->SetTitle("N_signal vs Luminosity; Luminosity [fb^-1]; N_s");
	gStyle->SetOptFit(0000);
	gPad->SetGrid(1,1);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.7);
	gr->Fit("pol1", "Q");
	gr_exp->SetMarkerStyle(2);
	gr_exp->SetMarkerSize(3);
	gr_exp->SetMarkerColor(4);
	mg->Add( gr_exp );
	mg->Add( gr );
	mg->Draw("AP");
	string title7 = var + " N_signal_luminostiy.png";
	can2->Print( title7.c_str() );

	TCanvas * can3 = new TCanvas("can3", "N_b vs luminosity");
	mg_b->SetTitle("N_backgruond vs Luminosity; Luminosity [fb^-1]; N_b");
	gStyle->SetOptFit(0000);
	gPad->SetGrid(1,1);
	gr_b->SetMarkerStyle(20);
	gr_b->SetMarkerSize(0.7);
	gr_b->Fit("pol1", "Q");
	gr_exp_b->SetMarkerStyle(2);
	gr_exp_b->SetMarkerSize(3);
	gr_exp_b->SetMarkerColor(4);
	mg_b->Add( gr_exp );
	mg_b->Add( gr );
	mg_b->Draw("AP");
	string title10 = var + " N_bg_luminostiy.png";
	can2->Print( title10.c_str() );
	return;
}
