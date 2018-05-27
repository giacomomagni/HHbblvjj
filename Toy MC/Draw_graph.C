
using namespace std;
//cross sections
#define cr_s 33.5
#define cr_b 831760
//branching ratios
#define br_s 0.0717
#define br_b 0.238
//cut fraction
#define evs_saved 1
#define evb_saved 0.006176

void Draw_graph( string var ){

	string name = var + "/" + var + " simulation.root";
	TFile * input = TFile::Open( name.c_str() );

	double lum [4] = { 375., 750., 1500., 3000. };	
//	double lum [4] = { 100., 150., 300., 3000. };	


	//mg_s 		Grafico segnale
	//mg_b  	Grafico background
	//gr_err_b	Grafici errore background
	//gr_err_s	Grafici errore segnale
	//gr_rel_b	Grafici errore relativo background
	//gr_rel_s	Grafici errore relativo segnale

	TGraphErrors * gr_s = new TGraphErrors();
	TGraphErrors * gr_b = new TGraphErrors();
	TGraphErrors * gr_err_s = new TGraphErrors();
	TGraphErrors * gr_err_b = new TGraphErrors();
	TGraphErrors * gr_rel_s = new TGraphErrors();
	TGraphErrors * gr_rel_b = new TGraphErrors();
	TGraphErrors * gr_exp_s = new TGraphErrors();
	TGraphErrors * gr_exp_b = new TGraphErrors();
	TMultiGraph *mg_s = new TMultiGraph();
	TMultiGraph *mg_b = new TMultiGraph();
	mg_s->SetTitle(" SM signal Events; Luminosity [fb^-1]; Ns");
	mg_b->SetTitle(" ttbar background; Luminosity [fb^-1]; Nb");
	gr_err_s->SetTitle("SM signal error; Luminosity [fb^-1]; sigma Ns");
	gr_err_b->SetTitle("ttbar error; Luminosity [fb^-1]; sigma Nb");
	gr_rel_s->SetTitle("SM signal relative error; Luminosity [fb^-1]; sigma Ns/Ns");
	gr_rel_b->SetTitle("ttbar relative error; Luminosity [fb^-1]; sigma Nb/Nb");

	TF1 * f0 = new TF1("f0", "gaus");		
	TF1 * f1 = new TF1("f1", "gaus");	
	
	TCanvas ** c = new TCanvas * [4];
	c[0] = new TCanvas("c0", "SM signal frequency"); 
	c[0]->Divide(4,1);

	for(int i=0; i<4;i++){

		c[0]->cd(i+1);
		cout << "--------------------------------------------------------------------------" << endl
			 << "                    *** " <<"LUMINOSITY " << lum[i] << " fb^-1 ***" << endl;

		double evs_exp = cr_s*br_s*lum[i]*evs_saved;
		double evb_exp = cr_b*br_b*lum[i]*evb_saved;

		string nametree = "tree_" + to_string(lum[i]);
   		TTree * tree = (TTree*) input->Get( nametree.c_str() );
		string title0 = "Signal events simulated at luminosity " + to_string(lum[i]) + " fb^-1";
		string title1 = "Background events simulated at luminosity " + to_string(lum[i]) + " fb^-1";
		string title2 = "Background + Signal events simulated at luminosity " + to_string(lum[i]) + " fb^-1";

		string comand = "nb>>" + title1; 
		tree->Draw( comand.c_str() );
		TH1F * histo_Nb = (TH1F*) gDirectory->Get( title1.c_str() );
		histo_Nb->GetXaxis()->SetTitle( "Simulated background events" );
		histo_Nb->GetYaxis()->SetTitle(	"Frequency [Events]");
		histo_Nb->Fit("f1", "NQ");
		double Nb_avg = f1->GetParameter(1);
		double Nb_avg_err = f1->GetParameter(2);
		
		comand = "ns>>" + title0; 
		tree->Draw( comand.c_str() );
		TH1F * histo_Ns = (TH1F*) gDirectory->Get( title0.c_str() );
		histo_Ns->GetXaxis()->SetTitle( "Simulated signal events" );
		histo_Ns->GetYaxis()->SetTitle(	"Frequency [Events]");
		histo_Ns->Fit("f0", "Q");
		double Ns_avg = f0->GetParameter(1);
		double Ns_avg_err = f0->GetParameter(2);

		gr_s->SetPoint(i, lum[i], Ns_avg );
		gr_s->SetPointError(i, 0.0, Ns_avg_err );
		gr_exp_s->SetPoint(i, lum[i],  evs_exp );
		gr_exp_s->SetPointError(i, 0.0, 0.09*evs_exp );

		gr_b->SetPoint(i, lum[i], Nb_avg );
		gr_b->SetPointError(i, 0.0, Nb_avg_err );
		gr_exp_b->SetPoint(i, lum[i],  evb_exp );
		gr_exp_b->SetPointError(i, 0.0, 0.09*evb_exp );

		gr_err_s->SetPoint(i, lum[i], Ns_avg_err);
		gr_err_b->SetPoint(i,lum[i], Nb_avg_err);
		gr_err_s->SetPointError(i, 0.0, f0->GetParError(2) );
		gr_err_b->SetPointError(i, 0.0, f1->GetParError(2) );

		gr_rel_s->SetPoint(i, lum[i], Ns_avg_err/Ns_avg);
		gr_rel_b->SetPoint(i,lum[i], Nb_avg_err/Nb_avg);
		
		cout << "---------------" << endl
			 << "Expected signal events: " << evs_exp << endl			 
			 << "Simulated signal events: " << Ns_avg <<  endl
			 << "Relative error: " << Ns_avg_err/Ns_avg*100 << "%" << endl
			 << "---------------" << endl
			 << "Expected background events: " << evb_exp << endl
			 << "Simulated background events: " << Nb_avg <<  endl
			 << "Relative error: " << Nb_avg_err/Nb_avg*100 << "%" << endl
			 << "----------------" << endl;
	}

//	string title3 = var + " Simulated Nb Ns frequency.png";
//	c[0]->Print( title3.c_str() );

	//Disegno del grafico Ns vs luminosity
	c[1] = new TCanvas("c1", "Ns vs luminosity");	
	gStyle->SetOptFit(0000);
	gPad->SetGrid(1,1);
	gr_s->SetMarkerStyle(20);
	gr_s->SetMarkerSize(2);
	gr_s->Fit("pol1", "Q");
	gr_exp_s->SetMarkerStyle(34);
	gr_exp_s->SetMarkerSize(3);
	gr_exp_s->SetMarkerColor(4);
	mg_s->Add( gr_exp_s );
	mg_s->Add( gr_s );
	mg_s->Draw("AP");
	
	TLegend *legend1 = new TLegend(0.6,0.2,0.98,0.38);
    legend1->AddEntry(gr_exp_s,"Signal expected", "lep");
    legend1->AddEntry(gr_s," Signal simulated", "lep");
    legend1->Draw("SAME");

//	string title4 = var + " Ns luminostiy.png";
//	c[1]->Print( title4.c_str() );

	//Disegno del grafico Nb vs luminosity	
	c[2] = new TCanvas("c2", "Nb vs luminosity");	
	gStyle->SetOptFit(0000);
	gPad->SetGrid(1,1);
	gr_b->SetMarkerStyle(20);
	gr_b->SetMarkerSize(2);
	gr_b->Fit("pol1", "Q");
	gr_exp_b->SetMarkerStyle(34);
	gr_exp_b->SetMarkerSize(3);
	gr_exp_b->SetMarkerColor(4);
	mg_b->Add( gr_exp_b );
	mg_b->Add( gr_b );
	mg_b->Draw("AP");

	TLegend *legend2 = new TLegend(0.6,0.2,0.98,0.38);
    legend2->AddEntry(gr_exp_b,"Background expected", "lep");
    legend2->AddEntry(gr_b,"Background expected", "lep");
    legend2->Draw("SAME");

//	string title5 = var + " Nb luminostiy.png";
//	c[2]->Print( title5.c_str() );

	c[3] = new TCanvas("c3", "Signal errors vs luminosity");	
	gStyle->SetOptFit(0000);
	gr_err_s->SetMarkerStyle(20);
	gr_err_s->SetMarkerSize(2);
	gr_err_s->SetMarkerColor(4);
	gr_err_s->SetLineColor(4);
	gr_err_s->SetLineWidth(3);
	gr_rel_s->SetMarkerStyle(20);
	gr_rel_s->SetMarkerSize(2);
	gr_rel_s->SetMarkerColor(4);
	gr_rel_s->SetLineColor(4);
	gr_rel_s->SetLineWidth(3);

	c[3]->Divide(2,1);
	c[3]->cd(1);
	gPad->SetGrid(1,1);
	gr_err_s->Draw("APL");
	c[3]->cd(2);
	gPad->SetGrid(1,1);
	gr_rel_s->Draw("APL");

//	string title6 = var + " relative error.png";
//	c[3]->Print( title6.c_str() );
	
	c[4] = new TCanvas("c4", "Background errors vs luminosity");	
	gStyle->SetOptFit(0000);
	gr_err_b->SetMarkerStyle(20);
	gr_err_b->SetMarkerSize(2);
	gr_err_b->SetMarkerColor(2);
	gr_err_b->SetLineColor(2);
	gr_err_b->SetLineWidth(3);
	gr_rel_b->SetMarkerStyle(20);
	gr_rel_b->SetMarkerSize(2);
	gr_rel_b->SetMarkerColor(2);
	gr_rel_b->SetLineColor(2);
	gr_rel_b->SetLineWidth(3);
	
	c[4]->Divide(2,1);
	c[4]->cd(1);
	gPad->SetGrid(1,1);
	gr_err_b->Draw("APL");
	c[4]->cd(2);
	gPad->SetGrid(1,1);
	gr_rel_b->Draw("APL");

//	string title7 = var + " relative error.png";
//	c[4]->Print( title7.c_str() )
	return;
}
