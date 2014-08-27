#define radmonanb_cxx
#include "radmonanb.h"

#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TMultiGraph.h>

// Analysis of PHENIX radiation monitors in Run 13
// John Haggerty, BNL, 2013.02.22

// Modified for Run 14
// John Haggerty, BNL, 2014.08.16

// 5 channels of sensors each of which have 4 subchannels:
// N negative temperature coefficient thermistor
// K 1 kohm resistor
// R radfet (read voltage with 160 microamp current)
// Si sensor (read voltage with 25 mA current)

// This works on radmon.root which you create with
// ./dumpradmon

Double_t temperature( Double_t v, Double_t i, Double_t beta = 3530.0 ) {

//T_NTC = 1/(1/297+1/BETA*(LN((VOLTAGE/2.5)/10000)))-273.15
  return ( 1.0/(1.0/297.0 + 1.0/beta*(TMath::Log((v/i)/10000.0)))-273.15 );

}

TMultiGraph *time_series( TGraph *g[], TString legend[], TString ptitle, TString ytitle, Int_t nsensor ) {
  TMultiGraph *mg = new TMultiGraph( "mg", "mg" );
  TLegend *leg = new TLegend(0.8680651,0.7870257,0.9982345,0.9987,NULL,"brNDC");
  
  Int_t icolor = 1;
  
  for ( Int_t i = 0; i < nsensor; i++ ) {
    g[i]->SetMarkerColor( icolor );
    g[i]->SetLineColor( icolor++ );
    leg->AddEntry( g[i], legend[i], "L" );
    mg->Add( g[i], "p" );
  }
  
  mg->Draw("a");
   
  // Make changes to axis after drawing, otherwise they don't exist
  mg->SetTitle( ptitle );
  mg->GetYaxis()->SetTitle( ytitle );

  mg->GetXaxis()->SetTimeDisplay(1);
  mg->GetXaxis()->SetNdivisions(-504); 
  mg->GetXaxis()->SetTitleOffset( 0.4 );
  mg->GetXaxis()->SetLabelOffset( 0.05 );
  mg->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}");
  mg->GetXaxis()->SetTimeOffset(0); 
  leg->SetFillColor(0);
  leg->Draw();
  
  return mg;
}

void radmonanb::Loop( Int_t idfirst, Int_t idlast, TString plotfile ) {

  //   In a ROOT session, you can do:
  //      Root > .L radmonanb.C
  //      Root > radmonanb t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  const Int_t nsensor = 7;

   // N, K, S, R
   std::vector<double> V_s[nsensor];
   std::vector<double> V_r[nsensor];
   std::vector<double> R_n[nsensor];
   std::vector<double> T_n[nsensor];
   std::vector<double> dV_r[nsensor];
   std::vector<double> dR_n[nsensor];
   std::vector<double> R_k[nsensor];
   std::vector<double> R_r[nsensor];
   std::vector<double> time[nsensor];
   std::vector<double> record[nsensor];
   std::vector<double> V_r_corrected[nsensor];
   std::vector<double> V_s_corrected[nsensor];

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if ( id < idfirst || id > idlast ) continue;

      // cuts to eliminate data pathologies

      if ( TMath::Abs( (i_r/i_r_set) - 1.0 ) < 0.01 && 
	   TMath::Abs( (i_s/i_s_set) - 1.0 ) < 0.01 && 
	   TMath::Abs( (i_n/i_n_set) - 1.0 ) < 0.01 &&
	   0.75 < v_n && v_n < 5.0 && 
	   v_r < 10.0 && v_s < 10.0 ) {
	// fill vectors for plotting
	V_s[channel].push_back(v_s);
	V_r[channel].push_back(v_r);
	dV_r[channel].push_back( 0.00005 );
	R_n[channel].push_back(v_n/i_n);
	T_n[channel].push_back( temperature(v_n,i_n) );
	dR_n[channel].push_back( 5.0 );
	R_k[channel].push_back(v_k/i_k);
	R_r[channel].push_back(v_r/i_r);
	time[channel].push_back(t);
	record[channel].push_back(id);
      }
   }

   gSystem->Setenv("TZ","EDT");
   gStyle->SetOptFit(111);

   // zero point and temperature correction for radfet 
   // Run 14 after 373780 May 28, 2014
   Double_t a0[nsensor] = { 3.923, 7.513, 5.395, 6.038, 4.140, 3.306, 3.263 };
   Double_t b0[nsensor] = { 0.0 };
   Double_t c0[nsensor] = { 0.0 };

   // zero point correction for Si detector
   // Run 14 after 373780 May 28, 2014
   Double_t s0[nsensor] = { 4.290, 7.220, 6.250, 2.694, 2.050, 1.054, 1.051 };
   
   TString legend[nsensor] = { "Channel 0: 150", "Channel 1: 151", "Channel 2: 152", "Channel 3: 148", "Channel 4: 149", "Channel 5", "Channel 6" };
   
   TGraphErrors *gr[nsensor];
   TGraph *grna[nsensor];

   TGraph *gr_vr[nsensor];
   TGraph *gr_vrc[nsensor];
   TGraph *gr_vs[nsensor];
   TGraph *gr_vsc[nsensor];
   TGraph *gr_rn[nsensor];
   TGraph *gr_tn[nsensor];
   TGraph *gr_rk[nsensor];
   TGraph *gr_vs_tn[nsensor];
   TGraph *gr_vr_tn[nsensor];
   
   TMultiGraph *mg_vs_tn = new TMultiGraph("mg_vs_tn","mg_vs_tn");;
   TMultiGraph *mg_vr_tn = new TMultiGraph("mg_vr_tn","mg_vr_tn");;

   Int_t icolor = 1;
   Double_t corrected;
   for ( Int_t i = 0; i < nsensor; i++ ) {

     std::cout << i << ": " << a0[i] << "," << b0[i] << "," << c0[i] << std::endl;
     for ( Int_t j = 0; j < (int) R_n[i].size(); j++ ) {
       corrected = V_r[i][j] - 
       	 ( a0[i] + b0[i]*R_n[i][j] + c0[i]*TMath::Power( R_n[i][j],2.0 ) );
       V_r_corrected[i].push_back( corrected );
       V_s_corrected[i].push_back( V_s[i][j] - s0[i] );
     }

     gr[i] = new TGraphErrors( R_n[i].size(), &R_n[i][0], &V_r[i][0], &dR_n[i][0], &dV_r[i][0] );
     grna[i] = new TGraph( R_n[i].size(), &(record[i][0]), &(R_n[i][0]) );

     gr_vs_tn[i] = new TGraph( T_n[i].size(), &(T_n[i][0]), &(V_s[i][0]) );
     gr_vs_tn[i]->SetMarkerColor( icolor );
     gr_vs_tn[i]->SetLineColor( icolor );
     mg_vs_tn->Add( gr_vs_tn[i], "p" );
    
     gr_vr_tn[i] = new TGraph( T_n[i].size(), &(T_n[i][0]), &(V_r[i][0]) );
     gr_vr_tn[i]->SetMarkerColor( icolor );
     gr_vr_tn[i]->SetLineColor( icolor );
     mg_vr_tn->Add( gr_vr_tn[i], "p" );

     // time series graphs

     gr_vr[i] = new TGraph( R_n[i].size(), &(time[i][0]), &(V_r[i][0]) );
     gr_vrc[i] = new TGraph( R_n[i].size(), &(time[i][0]), &(V_r_corrected[i][0]) );
     gr_vs[i] = new TGraph( V_s[i].size(), &(time[i][0]), &(V_s[i][0]) );
     gr_vsc[i] = new TGraph( V_s[i].size(), &(time[i][0]), &(V_s_corrected[i][0]) );
     gr_rn[i] = new TGraph( R_n[i].size(), &(time[i][0]), &(R_n[i][0]) );
     gr_tn[i] = new TGraph( T_n[i].size(), &(time[i][0]), &(T_n[i][0]) );
     gr_rk[i] = new TGraph( R_n[i].size(), &(time[i][0]), &(R_k[i][0]) );

     icolor++;
   }

   TCanvas *correct = new TCanvas("correct","correct",960,480);
   correct->Divide(2,1);

   correct->cd(1);
   mg_vs_tn->Draw("a");
   mg_vs_tn->GetXaxis()->SetTitle("T (#circC)");
   mg_vs_tn->GetYaxis()->SetTitle("V (V)");

   correct->cd(2);
   mg_vr_tn->Draw("a");
   mg_vr_tn->GetXaxis()->SetTitle("T (#circC)");
   mg_vr_tn->GetYaxis()->SetTitle("V (V)");

   TCanvas *cradmon = new TCanvas("cradmon","radmon",640,480);
   cradmon->Divide(2,3);

   // Draw the time series summary

   cradmon->cd(1);
   TMultiGraph *mg_rk;
   mg_rk = time_series( gr_rk, legend, "PHENIX Dosimeter Resistor", "R ( #Omega )", nsensor ); 
   mg_rk->GetYaxis()->SetRangeUser( 995, 1015.0 );

   cradmon->cd(2);
   TMultiGraph *mg_tn;
   mg_tn = time_series( gr_tn, legend, "PHENIX Dosimeter Temperature", "T (#circC)", nsensor ); 

   cradmon->cd(3);
   TMultiGraph *mg_vr;
   mg_vr = time_series( gr_vr, legend, "PHENIX Dosimeter RadFET", "V@160 #mu A (V)", nsensor ); 

   cradmon->cd(4);
   TMultiGraph *mg_vs;
   mg_vs = time_series( gr_vs, legend, "PHENIX Dosimeter p-i-n Diode", "V@25mA (V)", nsensor ); 

   cradmon->cd(5);
   TMultiGraph *mg_vrc;
   mg_vrc = time_series( gr_vrc, legend, "PHENIX Corrected RadFET", "#DeltaV (V)", nsensor ); 

   cradmon->cd(6);
   TMultiGraph *mg_vsc;
   mg_vsc = time_series( gr_vsc, legend, "PHENIX Corrected p-i-n Diode", "#DeltaV (V)", nsensor ); 

   cradmon->Print(plotfile);

   // save graphs
   
   TFile *rootfile = new TFile( "radmonanb.root", "recreate");
   Char_t graphname[20];
   for ( Int_t i = 0; i < nsensor; i++ ) {
     sprintf( graphname, "gr_vr%d", i );
     gr_vr[i]->Write( graphname );
     sprintf( graphname, "gr_vs%d", i );
     gr_vs[i]->Write( graphname );
   }
   rootfile->Close();

}
