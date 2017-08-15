#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TString.h"
#include "TH2.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include <TSystem.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSpline.h>

#include <iostream>
#include <bitset>

using namespace std;

const Int_t nsensor = 12;

//Bool_t plotme[nsensor] = { true, true, true, false, false, false, false, false, false, false, false, false };
// all of them
//Bool_t plotme[nsensor] = { true, true, true, true, true, true, true, true, true, true, true, true };
// central?
Bool_t plotme[nsensor] = { true, true, true, true, true, false, false, true, true, false, false, false };

TMultiGraph *mg_tn;
TMultiGraph *mg_rk;
TMultiGraph *mg_vr;
TMultiGraph *mg_vs;
TMultiGraph *mg_vrc;
TMultiGraph *mg_vsc;
TMultiGraph *mg_dose_r;
TMultiGraph *mg_dose_s;

// Functions which extract the data from the CERN dosimeters installed in PHENIX
// directly from the database and plot time series

// Example: plotradmon("2014-05-29 00:00:00", "2014-07-01 00:00:00")

// John Haggerty, BNL, 2014-08-25

// N, K, S, R
std::vector<double> V_s[nsensor];
std::vector<double> V_r[nsensor];
std::vector<double> R_n[nsensor];
std::vector<double> T_n[nsensor];
std::vector<double> R_k[nsensor];
std::vector<double> R_r[nsensor];
std::vector<double> read_time[nsensor];
std::vector<double> V_r_corrected[nsensor];
std::vector<double> dose_r[nsensor];
std::vector<double> V_s_corrected[nsensor];
std::vector<double> dose_s[nsensor];
std::vector<double> rs_ratio[nsensor];

TGraph *gr_vr[nsensor];
TGraph *gr_vrc[nsensor];
TGraph *gr_dose_r[nsensor];
TGraph *gr_vs[nsensor];
TGraph *gr_vsc[nsensor];
TGraph *gr_dose_s[nsensor];
TGraph *gr_rn[nsensor];
TGraph *gr_tn[nsensor];
TGraph *gr_rk[nsensor];

const Int_t nscaler = 2;

std::vector<double> sread_time[nscaler];
std::vector<double> cum_scaler_db[nscaler];
std::vector<double> cum_scaler_sum[nscaler];

TGraph *gr_scaler[nscaler];
TGraph *gr_scaler_sum[nscaler];

// These are ratios of the corrected voltages to the two scalers

std::vector<double> vsc_0[nsensor];
std::vector<double> vsc_1[nsensor];
std::vector<double> vrc_0[nsensor];
std::vector<double> vrc_1[nsensor];

TGraph *gr_sc0[nsensor];
TGraph *gr_sc1[nsensor];
TGraph *gr_rc0[nsensor];
TGraph *gr_rc1[nsensor];

// The SiPM leakage current data

const Int_t nsipm = 1;

std::vector<double> read_time_sipm[nsipm];
std::vector<double> sipm_current[nsipm];

TGraph *gr_sipm[nsipm];

TGraph *grxtime( TGraph *grin, TString title = "Title;x;y" )
{
  
  // generic help function for making time series TGraphs

  grin->SetTitle( title );

  grin->GetXaxis()->SetTimeDisplay(1);
  grin->GetXaxis()->SetNdivisions(-504); 
  grin->GetXaxis()->SetLabelOffset( 0.05 );
  grin->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}");
  //  grin->GetXaxis()->SetTimeOffset(-4*3600,"gmt"); 
  //  grin->GetYaxis()->SetTitleOffset( 1.4 );

  return grin;

}

Double_t radfet_dose( Double_t *deltav, Double_t *par )
{

  if ( deltav[0] < 0.0 ) return 0;

  // thesis-2007-013.pdf p. 96 eq. 6.7
  Double_t alow = 0.01854; // V/Gy
  Double_t blow = 0.91072;
  
  Double_t ahigh = 0.02921; // V/Gy
  Double_t bhigh = 0.78778;
  
  Double_t deltav_high_low = 0.5343;

  Double_t dose;
  if ( deltav[0] < deltav_high_low ) {
    dose = TMath::Exp( TMath::Log(deltav[0]/alow)/blow );
  } else {
    dose = TMath::Exp( TMath::Log(deltav[0]/ahigh)/bhigh );
  };
  
  return dose;
}

Double_t si_dose( Double_t *deltav, Double_t *par )
{
  // thesis-2007-013.pdf p. 135, S-1
  // 1/c = 1.6E8 cm^-2/mV, note mV
  Double_t one_over_c = 1.6E11; // 

  // return in units of 10^12 cm^-2  
  return one_over_c*deltav[0]/1.0E12;
}

void plotradfet()
{
  TCanvas *c1 = new TCanvas();
  TF1 *d = new TF1("d",radfet_dose,0.001,100.0,0);
  d->SetTitle("REM-501C Type K 250 nm RadFET;#DeltaV [V];Dose [Gy]");
  d->Draw();
  c1->SetLogx();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  c1->Print("radfet_dose.pdf","pdf portrait");
}

void plotsi()
{
  TCanvas *c1 = new TCanvas();
  TF1 *d = new TF1("d",si_dose,0.001,10.0,0);
  d->SetTitle("Si-1 LBSD diode;#DeltaV [V];#Phi [10^{12} cm^{-2}]");
  d->Draw();
  c1->SetLogx();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  c1->Print("si_dose.pdf","pdf portrait");
}


Double_t temperature( Double_t v, Double_t i, Double_t beta = 3530.0 ) {

  //T_NTC = 1/(1/297+1/BETA*(LN((VOLTAGE/2.5)/10000)))-273.15
  return ( 1.0/(1.0/297.0 + 1.0/beta*(TMath::Log((v/i)/10000.0)))-273.15 );

}

TMultiGraph *time_series( TGraph *g[], TString legend[], TString ptitle, TString ytitle, Int_t nsensor ) {

  TMultiGraph *mg = new TMultiGraph( "mg", "mg" );
  
  TLegend *leg = new TLegend(0.8680651,0.7870257,0.9982345,0.9987,NULL,"brNDC");

  Int_t icolor = 1;
  
  for ( Int_t i = 0; i < nsensor; i++ ) {
    if ( g[i] != 0 ) {
      g[i]->SetMarkerColor( icolor );
      g[i]->SetLineColor( icolor++ );
      g[i]->SetMarkerStyle( kFullCircle );
      g[i]->SetMarkerSize( 0.35 );
      g[i]->SetTitle(legend[i]);
      g[i]->SetFillColor(0);
      g[i]->SetFillStyle(0);
      leg->AddEntry( g[i], legend[i], "L" );
      mg->Add( g[i], "P" );
      //    mg->Add( g[i], "L" );
    }
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
  //  mg->GetXaxis()->SetTimeOffset(-4*3600,"gmt"); 
  leg->SetFillColor(0);
  leg->Draw();
  
  return mg;
}

Int_t getradmon( Char_t *start_datetime, Char_t *end_datetime )
{

  TSQLServer *serv = TSQLServer::Connect("pgsql://phnxdb0.phenix.bnl.gov/daq", "phnxrc", "");
  
  // Create the sql query
  
  TString columns = "id, EXTRACT(EPOCH FROM read_datetime)::INT AS read_timestamp, channel, i_n_set, i_n, v_n, i_k_set, i_k, v_k, i_s_set, i_s, v_s, i_r_set, i_r, v_r";
  
  TString sql = "SELECT ";
  sql += columns;
  sql += " FROM radmon WHERE read_datetime>=\'";
  sql += start_datetime;
  sql += "\' AND read_datetime<=\'";
  sql += end_datetime;
  sql += "\'";
  sql += " AND ABS( (i_r/i_r_set) - 1.0 ) < 0.01"; 
  sql += " AND ABS( (i_s/i_s_set) - 1.0 ) < 0.01";
  sql += " AND ABS( (i_n/i_n_set) - 1.0 ) < 0.01";
  sql += " AND v_n > 0.4 AND v_n < 2.0";
  // key
  sql += " AND v_r < 20.0 AND v_s < 20.0";
  //  sql += " AND v_r < 50.0 AND v_s < 50.0";
  sql += ";";
  cout << "sql query: " << sql << endl;

  TSQLResult *res;
  res = serv->Query(sql);

  // Extract the result of the query into vectors

  Int_t nrows = res->GetRowCount();
  Int_t nfields = res->GetFieldCount();
  cout << "rows: " << nrows << " columns: " << nfields << endl;

  TString fieldname;
  TString field;
  TSQLRow *row;

  Int_t channel = 0;
  Double_t read_timestamp;
  Double_t v_k = 0.0, v_n = 0.0, v_s = 0.0, v_r = 0.0;
  Double_t i_k = 0.0, i_n = 0.0, i_s = 0.0, i_r = 0.0;

  // zero point and temperature correction for radfet 
  // Run 14 after 373780 May 28, 2014
  //  Double_t a0[nsensor] = { 3.923, 7.513, 5.395, 6.038, 4.140, 3.306, 3.263 };
  // Beginning of Run 15 Jan 14, 2015
  //  Double_t a0[nsensor] = { 3.971, 7.871, 5.645, 7.566, 3.322, 3.108, 3.214 };
  // after adding CAN 1 and CAN 2 2015.04.15 
  //  Double_t a0[nsensor] = { 3.971, 7.871, 5.645, 7.566, 3.322, 3.108, 3.214, 3.274, 3.643 };
  Double_t a0[nsensor] = { 3.36, 3.36, 3.36, 3.36, 3.36, 3.36, 3.36, 3.36, 3.36 };

  // Here's at the beginning of Run 14 in January
  //  Double_t a0[nsensor] = { 3.78, 6.74, 4.94, 3.40, 3.45, 3.306, 3.263 };


  // zero point correction for Si detector
  // Run 14 after 373780 May 28, 2014
  //  Double_t s0[nsensor] = { 4.290, 7.220, 6.250, 2.694, 2.050, 1.054, 1.051 };
  // Begining of Run 15 Jan 14, 2015  
  //  Double_t s0[nsensor] = { 4.326, 7.078, 6.208, 3.349, 1.108, 1.048, 1.058 };
  // after adding CAN 1 and CAN 2 2015.04.15 
  //  Double_t s0[nsensor] = { 4.326, 7.078, 6.208, 3.349, 1.108, 1.048, 1.058, 1.046, 1.048 };
  // Here's the beginning of Run 14 in January
  //  Double_t s0[nsensor] = { 3.60, 6.20, 5.26, 1.16, 1.18, 1.054, 1.051 };
  Double_t s0[nsensor] = { 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05 };

  for ( channel = 0; channel < nsensor; channel++ ) {
    V_s[channel].clear();
    V_r[channel].clear();
    R_n[channel].clear();
    T_n[channel].clear();
    R_k[channel].clear();
    R_r[channel].clear();
    read_time[channel].clear();
    V_r_corrected[channel].clear();
    V_s_corrected[channel].clear();
    dose_r[channel].clear();
    dose_s[channel].clear();
    rs_ratio[channel].clear();
  }

  Double_t v_r_c = 0.0;
  Double_t v_s_c = 0.0;

  for (Int_t i = 0; i < nrows; i++) {
    row = res->Next();
    for (Int_t j = 0; j < nfields; j++) {
      fieldname = TString( res->GetFieldName(j) );
      field = TString( row->GetField(j) );
      // Extract all columns of each row
      // std::cout << "fieldname: " << fieldname << " field: " << field << std::endl;
      if ( fieldname == "read_timestamp" ) read_timestamp = field.Atof();
      if ( fieldname == "channel" ) channel = field.Atoi();
      if ( fieldname == "v_k" ) v_k = field.Atof();
      if ( fieldname == "v_n" ) v_n = field.Atof();
      if ( fieldname == "v_s" ) v_s = field.Atof();
      if ( fieldname == "v_r" ) v_r = field.Atof();
      if ( fieldname == "i_k" ) i_k = field.Atof();
      if ( fieldname == "i_n" ) i_n = field.Atof();
      if ( fieldname == "i_s" ) i_s = field.Atof();
      if ( fieldname == "i_r" ) i_r = field.Atof();
    }
    // Save all the columns in this row in vectors for plotting
      V_s[channel].push_back(v_s);
      V_r[channel].push_back(v_r);
      R_n[channel].push_back(v_n/i_n);
      T_n[channel].push_back( temperature(v_n,i_n) );
      R_k[channel].push_back(v_k/i_k);
      R_r[channel].push_back(v_r/i_r);
      read_time[channel].push_back(read_timestamp);
      v_r_c = v_r - a0[channel];
      V_r_corrected[channel].push_back( v_r_c );
      //      std::cout << channel << ": " << v_r_c << " " << radfet_dose( &v_r_c, 0 ) << std::endl; 
      dose_r[channel].push_back( radfet_dose( &v_r_c, 0 ) );
      v_s_c = v_s - s0[channel];
      V_s_corrected[channel].push_back( v_s_c );
      dose_s[channel].push_back( si_dose( &v_s_c, 0 ) );
      if ( TMath::Abs( v_s_c ) > 1E-9 ) {
	rs_ratio[channel].push_back( (v_r - a0[channel])/v_s_c );
      } else {
	rs_ratio[channel].push_back( 0.0 );
	};
  }
  
  return nrows;
  
}

Int_t getscaler( Int_t scalervector, Int_t scalerchannel, Char_t *start_datetime, Char_t *end_datetime )
{

  TSQLServer *serv = TSQLServer::Connect("mysql://phnxdb1.phenix.bnl.gov/scalers", "phoncs", "phenix7815");
  
  // Create the sql query
  
  Int_t which_table, which_rate;

  which_table = (scalerchannel/16) + 1;
  which_rate = ( scalerchannel%16 ) + 1;
  
  TString scaler_field_name = "rate";
  scaler_field_name += which_rate;

  TString columns = "UNIX_TIMESTAMP(rs.read_datetime) AS read_timestamp,rs.";
  columns += scaler_field_name;
  columns += ", (@csum:=@csum+60.0*rs.";
  columns += scaler_field_name;
  columns += ") as cum";
  
  TString sql = "SELECT ";
  sql += columns;
  sql += " FROM rhicscaler";
  sql += which_table;
  sql += " AS rs WHERE";
  sql += " rs.read_datetime>=\"";
  sql += start_datetime;
  sql += "\" AND rs.read_datetime<=\"";
  sql += end_datetime;
  sql += "\";";

  cout << "sql query: " << sql << endl;

  TSQLResult *res;
  res = serv->Query("SET @csum=0.0;");
  res = serv->Query(sql);

  // Extract the result of the query into vectors

  Int_t nrows = res->GetRowCount();
  Int_t nfields = res->GetFieldCount();
  cout << "rows: " << nrows << " columns: " << nfields << endl;

  TString fieldname;
  TString field;
  TSQLRow *row;

  Double_t read_timestamp;
  Double_t cum;
  Double_t running_sum = 0.0;
  Double_t summand = 0.0;

  sread_time[scalervector].clear();
  cum_scaler_db[scalervector].clear();
  cum_scaler_sum[scalervector].clear();

  for (Int_t i = 0; i < nrows; i++) {
    row = res->Next();
    for (Int_t j = 0; j < nfields; j++) {
      fieldname = TString( res->GetFieldName(j) );
      field = TString( row->GetField(j) );
      // Extract all columns of each row
      // std::cout << "fieldname: " << fieldname << " field: " << field << std::endl;
      if ( fieldname == "read_timestamp" ) read_timestamp = field.Atof();
      if ( fieldname == "cum" ) cum = field.Atof();
      if ( fieldname == scaler_field_name ) summand = field.Atof();
    }
    // Save all the columns in this row in vectors for plotting
    sread_time[scalervector].push_back(read_timestamp);
    cum_scaler_db[scalervector].push_back(cum);
    running_sum += 60.0*summand;
    cum_scaler_sum[scalervector].push_back(running_sum);
  }
  
  return nrows;
  
}

Int_t getsipm( Char_t *start_datetime, Char_t *end_datetime )
{

  TSQLServer *serv = TSQLServer::Connect("pgsql://phnxdb1.phenix.bnl.gov/daq", "phnxrc", "");
  
  // Create the sql query
  TString sql = "SELECT EXTRACT(EPOCH FROM read_datetime)::INT AS read_timestamp,current FROM sipm WHERE read_datetime >= \'";
  sql += start_datetime;
  sql += "\' AND read_datetime < \'";
  sql += end_datetime;
  sql += "\'";
  std::cout << "sql query: " << sql << std::endl;

  TSQLResult *res;
  res = serv->Query(sql);
  
  // Extract the result of the query into vectors
  
  Int_t nrows = res->GetRowCount();
  Int_t nfields = res->GetFieldCount();
  std::cout << "rows: " << nrows << " columns: " << nfields << std::endl;
  
  TString fieldname;
  TString field;
  TSQLRow *row;
  //  std::vector<double> read_timestamp;
  //  std::vector<double> read_current;

 for (Int_t i = 0; i < nrows; i++) {
    row = res->Next();
    for (Int_t j = 0; j < nfields; j++) {
      fieldname = TString( res->GetFieldName(j) );
      field = TString( row->GetField(j) );
      //      std::cout << "fieldname: " << fieldname << " field: " << field << std::endl;
      if ( fieldname == "read_timestamp" ) read_time_sipm[0].push_back( field.Atof() );
      if ( fieldname.Contains( "current" ) ) sipm_current[0].push_back( field.Atof() );
    }
  }

 return nrows;

}

Int_t plotradmon( Char_t *start_datetime, Char_t *end_datetime ) {
  
  gROOT->SetStyle("Modern");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.4);
  gStyle->SetTimeOffset(6*3600);


  gROOT->SetStyle("Modern");

  TString legend[nsensor] = { "Channel 0: 150", "Channel 1: 151", "Channel 2: 152", "Channel 3: 148", "Channel 4: 149", "Channel 5", "Channel 6", "CAN 1", "CAN 2", "X-1", "X-2", "X-3" };

  Int_t nlines;
  nlines = getradmon( start_datetime, end_datetime );

  /*
  TGraph *gr_vr[nsensor];
  TGraph *gr_vrc[nsensor];
  TGraph *gr_vs[nsensor];
  TGraph *gr_vsc[nsensor];
  TGraph *gr_rn[nsensor];
  TGraph *gr_tn[nsensor];
  TGraph *gr_rk[nsensor];
  */

  //  TGraph *gr_vs_tn[nsensor];
  //  TGraph *gr_vr_tn[nsensor];

  //  TMultiGraph *mg_vs_tn = new TMultiGraph("mg_vs_tn","mg_vs_tn");;
  //  TMultiGraph *mg_vr_tn = new TMultiGraph("mg_vr_tn","mg_vr_tn");;
  
  for ( Int_t i = 0; i < nsensor; i++ ) {
    if ( plotme[i] ) {
      if ( R_n[i].size() > 0 ) {
	gr_vr[i] = new TGraph( R_n[i].size(), &(read_time[i][0]), &(V_r[i][0]) );
	gr_vrc[i] = new TGraph( R_n[i].size(), &(read_time[i][0]), &(V_r_corrected[i][0]) );
	gr_dose_r[i] = new TGraph( R_n[i].size(), &(read_time[i][0]), &(dose_r[i][0]) );
	gr_rn[i] = new TGraph( R_n[i].size(), &(read_time[i][0]), &(R_n[i][0]) );
	gr_rk[i] = new TGraph( R_n[i].size(), &(read_time[i][0]), &(R_k[i][0]) );
      } else {
	std::cout << "no R_n entries " << i << std::endl;
      }
      if (  V_s[i].size() > 0 ) {
	gr_vs[i] = new TGraph( V_s[i].size(), &(read_time[i][0]), &(V_s[i][0]) );
	gr_vsc[i] = new TGraph( V_s[i].size(), &(read_time[i][0]), &(V_s_corrected[i][0]) );
	gr_dose_s[i] = new TGraph( V_s[i].size(), &(read_time[i][0]), &(dose_s[i][0]) );
      } else {
	std::cout << "no R_n entries " << i << std::endl;
      }
      if (  T_n[i].size() > 0 ) {
	gr_tn[i] = new TGraph( T_n[i].size(), &(read_time[i][0]), &(T_n[i][0]) );
      } else { 
	std::cout << "no T_n entries " << i << std::endl;
      }
    }
  };

  TCanvas *cradmon = new TCanvas("cradmon","radmon",640,480);
  cradmon->Divide(2,4);
  Int_t ipad = 0;

// Draw the time series summary

  cradmon->cd(++ipad);
  mg_rk = time_series( gr_rk, legend, "PHENIX Dosimeter Resistor", "R ( #Omega )", nsensor ); 
  mg_rk->GetYaxis()->SetRangeUser( 995, 1015.0 );
  
  cradmon->cd(++ipad);
  mg_tn = time_series( gr_tn, legend, "PHENIX Dosimeter Temperature", "T (#circC)", nsensor ); 
  
  cradmon->cd(++ipad);
  mg_vr = time_series( gr_vr, legend, "PHENIX Dosimeter RadFET", "V@160 #mu A (V)", nsensor ); 
  
  cradmon->cd(++ipad);
  mg_vs = time_series( gr_vs, legend, "PHENIX Dosimeter p-i-n Diode", "V@25mA (V)", nsensor ); 

  cradmon->cd(++ipad);
  mg_vrc = time_series( gr_vrc, legend, "PHENIX Subtracted RadFET", "#DeltaV (V)", nsensor ); 
  
  cradmon->cd(++ipad);
  mg_vsc = time_series( gr_vsc, legend, "PHENIX Subtracted p-i-n Diode", "#DeltaV (V)", nsensor ); 

  cradmon->cd(++ipad);
  mg_dose_r = time_series( gr_dose_r, legend, "RadFET Dose", "Dose (Gy)", nsensor ); 
  
  cradmon->cd(++ipad);
  mg_dose_s = time_series( gr_dose_s, legend, "Fluence from PIN Diode", "Fluence (10^{12}/cm^{2})", nsensor ); 

  TCanvas *cradmon1 = new TCanvas("cradmon1","radmon",640,480);
  cradmon1->Divide(2,2);

  ipad = 0;
  cradmon1->cd(++ipad);
  mg_vrc = time_series( gr_vrc, legend, "PHENIX Subtracted RadFET", "#DeltaV (V)", nsensor ); 
  
  cradmon1->cd(++ipad);
  mg_vsc = time_series( gr_vsc, legend, "PHENIX Subtracted p-i-n Diode", "#DeltaV (V)", nsensor ); 

  cradmon1->cd(++ipad);
  mg_dose_r = time_series( gr_dose_r, legend, "RadFET Dose", "Dose (Gy)", nsensor ); 
  
  cradmon1->cd(++ipad);
  mg_dose_s = time_series( gr_dose_s, legend, "Fluence from LBSD Si-1 PIN Diode (k=1)", "Fluence (10^{12}/cm^{2})", nsensor ); 


  TCanvas *cradmon2 = new TCanvas("cradmon2","radmon",640,480);
  mg_dose_r = time_series( gr_dose_r, legend, "250 nm Thin Film RadFET Dose (160 #muA)", "Dose (Gy)", nsensor ); 
  
  TCanvas *cradmon3 = new TCanvas("cradmon3","radmon",640,480);
  mg_dose_s = time_series( gr_dose_s, legend, "LBSD Si-1 PIN Diode Fluence (25 mA)", "Fluence (10^{12}/cm^{2})", nsensor ); 

  return 0;
}

Int_t plotscaler( Char_t *start_datetime, Char_t *end_datetime ) {
  
  gROOT->SetStyle("Modern");

  Int_t nlines;
  // ZDCNS
  nlines = getscaler( 0, 0, start_datetime, end_datetime );
  // He3
  nlines = getscaler(1, 22, start_datetime, end_datetime );
  
  TCanvas *cscaler = new TCanvas("cscaler","cscaler",640,480);
  cscaler->Divide(1,2);
  Int_t ipad = 0;

  //  TGraph *gr_scaler[nscaler];
  for ( Int_t i = 0; i < nscaler; i++ ) {
    gr_scaler[i] = new TGraph( sread_time[i].size(), &(sread_time[i][0]), &(cum_scaler_db[i][0]) );
  }

  //  TGraph *gr_scaler_sum[nscaler];
  for ( Int_t i = 0; i < nscaler; i++ ) {
    gr_scaler_sum[i] = new TGraph( sread_time[i].size(), &(sread_time[i][0]), &(cum_scaler_sum[i][0]) );
  }

// Draw the time series summary

  TString scaler_legend[nscaler] = { "ZDC", "^{3}He tube" };

  cscaler->cd(++ipad);
  TMultiGraph *mg_sc;
  mg_sc = time_series( gr_scaler, scaler_legend, "PHENIX Scalers", "Cumulative counts", nscaler ); 

  //  TString scaler_legend[nscaler] = { "ZDC", "GE 3He-CO2 Neutron Counter" };

  cscaler->cd(++ipad);
  TMultiGraph *mg_sc_sum;
  mg_sc_sum = time_series( gr_scaler_sum, scaler_legend, "PHENIX Scalers--summed", "Cumulative counts", nscaler ); 

  return 0;
}

Int_t plotsipm( Char_t *start_datetime, Char_t *end_datetime ) {
  
  gROOT->SetStyle("Modern");

  Int_t nlines;
  nlines = getsipm( start_datetime, end_datetime );
  
  gr_sipm[0] = new TGraph( read_time_sipm[0].size(), &read_time_sipm[0][0], &sipm_current[0][0] );
  gr_sipm[0]->GetYaxis()->SetRangeUser(0.0,3500.0);
  
  TCanvas *csipm = new TCanvas("csipm","sipm",640,480);
  csipm->Divide(1,1);
  gr_sipm[0]->Draw("ap");
   
  // Make changes to axis after drawing, otherwise they don't exist
  gr_sipm[0]->SetTitle( "PHENIX SiPM Test" );
  gr_sipm[0]->GetYaxis()->SetTitle( "Leakage current (nA)" );

  gr_sipm[0]->GetXaxis()->SetTimeDisplay(1);
  gr_sipm[0]->GetXaxis()->SetNdivisions(-504); 
  gr_sipm[0]->GetXaxis()->SetTitleOffset( 0.4 );
  gr_sipm[0]->GetXaxis()->SetLabelOffset( 0.05 );
  gr_sipm[0]->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}");
  //  gr_sipm[0]->GetXaxis()->SetTimeOffset(-4*3600,"gmt"); 
  //  gr_sipm[0]->GetYaxis()->SetTitleOffset(1.5);
  
  return 0;
}

Int_t plotratio( void )
{
  
  TString legend[nsensor] = { "Channel 0: 150", "Channel 1: 151", "Channel 2: 152", "Channel 3: 148", "Channel 4: 149", "Channel 5", "Channel 6" };
  
  TSpline3 *sp0 = new TSpline3("sp0",gr_scaler[0]);
  TSpline3 *sp1 = new TSpline3("sp1",gr_scaler[1]);
  Double_t s0;
  Double_t s1;

  for ( Int_t i = 0; i < nsensor; i++ ) {
    for (UInt_t j = 0; j < read_time[i].size(); j++) {
      s0 = sp0->Eval( read_time[i][j] )*1E-9;
      s1 = sp1->Eval( read_time[i][j] )*1E-9;
      vsc_0[i].push_back( V_s_corrected[i][j]/s0 );
      vrc_0[i].push_back( V_r_corrected[i][j]/s0 );
      vsc_1[i].push_back( V_s_corrected[i][j]/s1 );
      vrc_1[i].push_back( V_r_corrected[i][j]/s1 );
    }
    gr_sc0[i] = new TGraph( read_time[i].size(), &(read_time[i][0]), &(vsc_0[i][0]) );
    gr_sc1[i] = new TGraph( read_time[i].size(), &(read_time[i][0]), &(vsc_1[i][0]) );
    gr_rc0[i] = new TGraph( read_time[i].size(), &(read_time[i][0]), &(vrc_0[i][0]) );
    gr_rc1[i] = new TGraph( read_time[i].size(), &(read_time[i][0]), &(vrc_1[i][0]) );
  }

  TCanvas *srats = new TCanvas("srats","sratmon",640,480);
  srats->Divide(2,2);
  Int_t ipad = 0;

  Bool_t logy = false;
  Double_t ymin = 0.0;
  Double_t ymax = 0.1;

  srats->cd(++ipad);
  TMultiGraph *mgs_0;
  mgs_0 = time_series( gr_sc0, legend, "PHENIX p-i-n diode", "#Delta V_{s}/ZDNS", nsensor ); 
  mgs_0->GetYaxis()->SetRangeUser(ymin,ymax);
  if ( logy ) gPad->SetLogy();

  srats->cd(++ipad);
  TMultiGraph *mgs_1;
  mgs_1 = time_series( gr_sc1, legend, "PHENIX p-i-n diode", "#Delta V_{s}/^{3}He counter", nsensor ); 
  mgs_1->GetYaxis()->SetRangeUser(ymin,ymax);
   if ( logy ) gPad->SetLogy();

  srats->cd(++ipad);
  TMultiGraph *mgr_0;
  mgr_0 = time_series( gr_rc0, legend, "PHENIX RadFET", "#Delta V_{r}/ZDCNS", nsensor ); 
  mgr_0->GetYaxis()->SetRangeUser(ymin,ymax);
   if ( logy ) gPad->SetLogy();

  srats->cd(++ipad);
  TMultiGraph *mgr_1;
  mgr_1 = time_series( gr_rc1, legend, "PHENIX RadFET", "#Delta V_{r}/^{3}He counter", nsensor ); 
  mgr_1->GetYaxis()->SetRangeUser(ymin,ymax);
   if ( logy ) gPad->SetLogy();

  return 0;
}


