

////
//  combined.c
//
//


#include <stdio.h>
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TText.h"
#include <TString.h>
#include "TMath.h"







void combined(){
    
    TFile *f = TFile::Open("combined_out.root");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    //gPad->SetTicks(1,1);
    //////////////////////////Barrel//////////////////////////////
    
    //////////////////////////Denominator////////////////////////
    TH1F *Wjets_den_Barrel = (TH1F*)f->Get("Wjets_den_pt_Barrel");
    Wjets_den_Barrel->Sumw2();
    TH1F *Zjets_den_Barrel = (TH1F*)f->Get("Zjets_den_pt_Barrel");
    Zjets_den_Barrel->Sumw2();
    TH1F *ttbar_den_Barrel = (TH1F*)f->Get("ttbar_den_pt_Barrel");
    ttbar_den_Barrel->Sumw2();
    TH1F *WW_den_Barrel = (TH1F*)f->Get("WW_den_pt_Barrel");
    WW_den_Barrel->Sumw2();
    TH1F *WZ_den_Barrel = (TH1F*)f->Get("WZ_den_pt_Barrel");
    WZ_den_Barrel->Sumw2();
    TH1F *ZZ_den_Barrel = (TH1F*)f->Get("ZZ_den_pt_Barrel");
    ZZ_den_Barrel->Sumw2();
    TH1F *tW_den_Barrel = (TH1F*)f->Get("tW_den_pt_Barrel");
    tW_den_Barrel->Sumw2();
    TH1F *Wantitop_den_Barrel = (TH1F*)f->Get("Wantitop_den_pt_Barrel");
    Wantitop_den_Barrel->Sumw2();
    TH1F *dy_den_Barrel = (TH1F*)f->Get("dy_den_pt_Barrel");
    dy_den_Barrel->Sumw2();
    TH1F *Data_den_Barrel = (TH1F*)f->Get("Data_den_pt_Barrel");
    Data_den_Barrel->Sumw2();
    
    TH1F *Wjets_den_Barrel_fx = (TH1F*)f->Get("Wjets_den_pt_Barrel_fx");
    Wjets_den_Barrel_fx->Sumw2();
    TH1F *Zjets_den_Barrel_fx = (TH1F*)f->Get("Zjets_den_pt_Barrel_fx");
    Zjets_den_Barrel_fx->Sumw2();
    TH1F *ttbar_den_Barrel_fx = (TH1F*)f->Get("ttbar_den_pt_Barrel_fx");
    ttbar_den_Barrel_fx->Sumw2();
    TH1F *WW_den_Barrel_fx = (TH1F*)f->Get("WW_den_pt_Barrel_fx");
    WW_den_Barrel_fx->Sumw2();
    TH1F *WZ_den_Barrel_fx = (TH1F*)f->Get("WZ_den_pt_Barrel_fx");
    WZ_den_Barrel_fx->Sumw2();
    TH1F *ZZ_den_Barrel_fx = (TH1F*)f->Get("ZZ_den_pt_Barrel_fx");
    ZZ_den_Barrel_fx->Sumw2();
    TH1F *tW_den_Barrel_fx = (TH1F*)f->Get("tW_den_pt_Barrel_fx");
    tW_den_Barrel_fx->Sumw2();
    TH1F *Wantitop_den_Barrel_fx = (TH1F*)f->Get("Wantitop_den_pt_Barrel_fx");
    Wantitop_den_Barrel_fx->Sumw2();
    TH1F *dy_den_Barrel_fx = (TH1F*)f->Get("dy_den_pt_Barrel_fx");
    dy_den_Barrel_fx->Sumw2();
    TH1F *Data_den_Barrel_fx = (TH1F*)f->Get("Data_den_pt_Barrel_fx");
    Data_den_Barrel_fx->Sumw2();
    //////////////////////////Denominator////////////////////////
    
    //////////////////////////Numerator////////////////////////
    TH1F *Wjets_num_Barrel = (TH1F*)f->Get("Wjets_num_pt_Barrel");
    Wjets_num_Barrel->Sumw2();
    TH1F *Zjets_num_Barrel = (TH1F*)f->Get("Zjets_num_pt_Barrel");
    Zjets_num_Barrel->Sumw2();
    TH1F *ttbar_num_Barrel = (TH1F*)f->Get("ttbar_num_pt_Barrel");
    ttbar_num_Barrel->Sumw2();
    TH1F *WW_num_Barrel = (TH1F*)f->Get("WW_num_pt_Barrel");
    WW_num_Barrel->Sumw2();
    TH1F *WZ_num_Barrel = (TH1F*)f->Get("WZ_num_pt_Barrel");
    WZ_num_Barrel->Sumw2();
    TH1F *ZZ_num_Barrel = (TH1F*)f->Get("ZZ_num_pt_Barrel");
    ZZ_num_Barrel->Sumw2();
    TH1F *tW_num_Barrel = (TH1F*)f->Get("tW_num_pt_Barrel");
    tW_num_Barrel->Sumw2();
    TH1F *Wantitop_num_Barrel = (TH1F*)f->Get("Wantitop_num_pt_Barrel");
    Wantitop_num_Barrel->Sumw2();
    TH1F *dy_num_Barrel = (TH1F*)f->Get("dy_num_pt_Barrel");
    dy_num_Barrel->Sumw2();
    TH1F *Data_num_Barrel = (TH1F*)f->Get("Data_num_pt_Barrel");
    Data_num_Barrel->Sumw2();
    
    TH1F *Wjets_num_Barrel_fx = (TH1F*)f->Get("Wjets_num_pt_Barrel_fx");
    Wjets_num_Barrel_fx->Sumw2();
    TH1F *Zjets_num_Barrel_fx = (TH1F*)f->Get("Zjets_num_pt_Barrel_fx");
    Zjets_num_Barrel_fx->Sumw2();
    TH1F *ttbar_num_Barrel_fx = (TH1F*)f->Get("ttbar_num_pt_Barrel_fx");
    ttbar_num_Barrel_fx->Sumw2();
    TH1F *WW_num_Barrel_fx = (TH1F*)f->Get("WW_num_pt_Barrel_fx");
    WW_num_Barrel_fx->Sumw2();
    TH1F *WZ_num_Barrel_fx = (TH1F*)f->Get("WZ_num_pt_Barrel_fx");
    WZ_num_Barrel_fx->Sumw2();
    TH1F *ZZ_num_Barrel_fx = (TH1F*)f->Get("ZZ_num_pt_Barrel_fx");
    ZZ_num_Barrel_fx->Sumw2();
    TH1F *tW_num_Barrel_fx = (TH1F*)f->Get("tW_num_pt_Barrel_fx");
    tW_num_Barrel_fx->Sumw2();
    TH1F *Wantitop_num_Barrel_fx = (TH1F*)f->Get("Wantitop_num_pt_Barrel_fx");
    Wantitop_num_Barrel_fx->Sumw2();
    TH1F *dy_num_Barrel_fx = (TH1F*)f->Get("dy_num_pt_Barrel_fx");
    dy_num_Barrel_fx->Sumw2();
    TH1F *Data_num_Barrel_fx = (TH1F*)f->Get("Data_num_pt_Barrel_fx");
    Data_num_Barrel_fx->Sumw2();
    //////////////////////////Numerator/////////////////////////
    
    //////////////////////////Barrel//////////////////////////////
    
    
    
    //////////////////////////Endcap//////////////////////////////
    
    //////////////////////////Denominator////////////////////////
    TH1F *Wjets_den_Endcap = (TH1F*)f->Get("Wjets_den_pt_Endcap");
    Wjets_den_Endcap->Sumw2();
    TH1F *Zjets_den_Endcap = (TH1F*)f->Get("Zjets_den_pt_Endcap");
    Zjets_den_Endcap->Sumw2();
    TH1F *ttbar_den_Endcap = (TH1F*)f->Get("ttbar_den_pt_Endcap");
    ttbar_den_Endcap->Sumw2();
    TH1F *WW_den_Endcap = (TH1F*)f->Get("WW_den_pt_Endcap");
    WW_den_Endcap->Sumw2();
    TH1F *WZ_den_Endcap = (TH1F*)f->Get("WZ_den_pt_Endcap");
    WZ_den_Endcap->Sumw2();
    TH1F *ZZ_den_Endcap = (TH1F*)f->Get("ZZ_den_pt_Endcap");
    ZZ_den_Endcap->Sumw2();
    TH1F *tW_den_Endcap = (TH1F*)f->Get("tW_den_pt_Endcap");
    tW_den_Endcap->Sumw2();
    TH1F *Wantitop_den_Endcap = (TH1F*)f->Get("Wantitop_den_pt_Endcap");
    Wantitop_den_Endcap->Sumw2();
    TH1F *dy_den_Endcap = (TH1F*)f->Get("dy_den_pt_Endcap");
    dy_den_Endcap->Sumw2();
    TH1F *Data_den_Endcap = (TH1F*)f->Get("Data_den_pt_Endcap");
    Data_den_Endcap->Sumw2();
    
    TH1F *Wjets_den_Endcap_fx = (TH1F*)f->Get("Wjets_den_pt_Endcap_fx");
    Wjets_den_Endcap_fx->Sumw2();
    TH1F *Zjets_den_Endcap_fx = (TH1F*)f->Get("Zjets_den_pt_Endcap_fx");
    Zjets_den_Endcap_fx->Sumw2();
    TH1F *ttbar_den_Endcap_fx = (TH1F*)f->Get("ttbar_den_pt_Endcap_fx");
    ttbar_den_Endcap_fx->Sumw2();
    TH1F *WW_den_Endcap_fx = (TH1F*)f->Get("WW_den_pt_Endcap_fx");
    WW_den_Endcap_fx->Sumw2();
    TH1F *WZ_den_Endcap_fx = (TH1F*)f->Get("WZ_den_pt_Endcap_fx");
    WZ_den_Endcap_fx->Sumw2();
    TH1F *ZZ_den_Endcap_fx = (TH1F*)f->Get("ZZ_den_pt_Endcap_fx");
    ZZ_den_Endcap_fx->Sumw2();
    TH1F *tW_den_Endcap_fx = (TH1F*)f->Get("tW_den_pt_Endcap_fx");
    tW_den_Endcap_fx->Sumw2();
    TH1F *Wantitop_den_Endcap_fx = (TH1F*)f->Get("Wantitop_den_pt_Endcap_fx");
    Wantitop_den_Endcap_fx->Sumw2();
    TH1F *dy_den_Endcap_fx = (TH1F*)f->Get("dy_den_pt_Endcap_fx");
    dy_den_Endcap_fx->Sumw2();
    TH1F *Data_den_Endcap_fx = (TH1F*)f->Get("Data_den_pt_Endcap_fx");
    Data_den_Endcap_fx->Sumw2();
    //////////////////////////Denominator////////////////////////
    
    //////////////////////////Numerator////////////////////////
    TH1F *Wjets_num_Endcap = (TH1F*)f->Get("Wjets_num_pt_Endcap");
    Wjets_num_Endcap->Sumw2();
    TH1F *Zjets_num_Endcap = (TH1F*)f->Get("Zjets_num_pt_Endcap");
    Zjets_num_Endcap->Sumw2();
    TH1F *ttbar_num_Endcap = (TH1F*)f->Get("ttbar_num_pt_Endcap");
    ttbar_num_Endcap->Sumw2();
    TH1F *WW_num_Endcap = (TH1F*)f->Get("WW_num_pt_Endcap");
    WW_num_Endcap->Sumw2();
    TH1F *WZ_num_Endcap = (TH1F*)f->Get("WZ_num_pt_Endcap");
    WZ_num_Endcap->Sumw2();
    TH1F *ZZ_num_Endcap = (TH1F*)f->Get("ZZ_num_pt_Endcap");
    ZZ_num_Endcap->Sumw2();
    TH1F *tW_num_Endcap = (TH1F*)f->Get("tW_num_pt_Endcap");
    tW_num_Endcap->Sumw2();
    TH1F *Wantitop_num_Endcap = (TH1F*)f->Get("Wantitop_num_pt_Endcap");
    Wantitop_num_Endcap->Sumw2();
    TH1F *dy_num_Endcap = (TH1F*)f->Get("dy_num_pt_Endcap");
    dy_num_Endcap->Sumw2();
    TH1F *Data_num_Endcap = (TH1F*)f->Get("Data_num_pt_Endcap");
    Data_num_Endcap->Sumw2();
    
    TH1F *Wjets_num_Endcap_fx = (TH1F*)f->Get("Wjets_num_pt_Endcap_fx");
    Wjets_num_Endcap_fx->Sumw2();
    TH1F *Zjets_num_Endcap_fx = (TH1F*)f->Get("Zjets_num_pt_Endcap_fx");
    Zjets_num_Endcap_fx->Sumw2();
    TH1F *ttbar_num_Endcap_fx = (TH1F*)f->Get("ttbar_num_pt_Endcap_fx");
    ttbar_num_Endcap_fx->Sumw2();
    TH1F *WW_num_Endcap_fx = (TH1F*)f->Get("WW_num_pt_Endcap_fx");
    WW_num_Endcap_fx->Sumw2();
    TH1F *WZ_num_Endcap_fx = (TH1F*)f->Get("WZ_num_pt_Endcap_fx");
    WZ_num_Endcap_fx->Sumw2();
    TH1F *ZZ_num_Endcap_fx = (TH1F*)f->Get("ZZ_num_pt_Endcap_fx");
    ZZ_num_Endcap_fx->Sumw2();
    TH1F *tW_num_Endcap_fx = (TH1F*)f->Get("tW_num_pt_Endcap_fx");
    tW_num_Endcap_fx->Sumw2();
    TH1F *Wantitop_num_Endcap_fx = (TH1F*)f->Get("Wantitop_num_pt_Endcap_fx");
    Wantitop_num_Endcap_fx->Sumw2();
    TH1F *dy_num_Endcap_fx = (TH1F*)f->Get("dy_num_pt_Endcap_fx");
    dy_num_Endcap_fx->Sumw2();
    TH1F *Data_num_Endcap_fx = (TH1F*)f->Get("Data_num_pt_Endcap_fx");
    Data_num_Endcap_fx->Sumw2();
    //////////////////////////Numerator/////////////////////////
    
    //////////////////////////Endcap//////////////////////////////
    
    
    
    
    //////////////////////////Barrel//////////////////////////////
    
    TCanvas *c1 = new TCanvas("rate", "rate", 1000, 800);
    c1->cd();
    
    
    
    TH1F* total_den_Barrel = (TH1F*) Data_den_Barrel->Clone();
    total_den_Barrel->Sumw2();
    
    /////////////////////For average/////////////////////////////
    TH1F* total_den_Barrel_fx = (TH1F*) Data_den_Barrel_fx->Clone();
    total_den_Barrel_fx->Sumw2();
    TH1F * den_Barrel_avg = (TH1F*) total_den_Barrel_fx->Clone();
    den_Barrel_avg->Sumw2();
    den_Barrel_avg->Divide(total_den_Barrel);
    
    //const int Barrel_bins = den_Barrel_avg->GetNbinsX();
    const int Barrel_bins = 4;
    std::cout<<"Number of bins : "<<Barrel_bins<<std::endl;
    double avg_pt_Barrel[Barrel_bins] = {0.0};
    double FR_Barrel[Barrel_bins] = {0.0};
    double FR_low_Barrel[Barrel_bins] = {0.0};
    double FR_high_Barrel[Barrel_bins] = {0.0};
    double avg_pt_low_Barrel[Barrel_bins] = {0.0};
    double avg_pt_high_Barrel[Barrel_bins] = {0.0};
    
    
    for(int i=1; i<=Barrel_bins; i++ ){
        avg_pt_Barrel[i-1] = den_Barrel_avg->GetBinContent(i);
        std::cout<<"bin : "<<i<<"   avg_pt : "<<avg_pt_Barrel[i-1]<<std::endl;
    }
    /////////////////////For average/////////////////////////////
    
    
    
    TH1F* MC_den_Barrel = (TH1F*) Wjets_den_Barrel->Clone();
    MC_den_Barrel->Sumw2();
    MC_den_Barrel->Add(dy_den_Barrel);
    MC_den_Barrel->Add(Zjets_den_Barrel);
    MC_den_Barrel->Add(ttbar_den_Barrel);
    MC_den_Barrel->Add(WW_den_Barrel);
    MC_den_Barrel->Add(WZ_den_Barrel);
    MC_den_Barrel->Add(ZZ_den_Barrel);
    MC_den_Barrel->Add(Wantitop_den_Barrel);
    MC_den_Barrel->Add(tW_den_Barrel);
    total_den_Barrel->Add(MC_den_Barrel,-1);
    
    
    TH1F* total_num_Barrel = (TH1F*) Data_num_Barrel->Clone();
    total_num_Barrel->Sumw2();
    TH1F* MC_num_Barrel = (TH1F*) Wjets_num_Barrel->Clone();
    MC_num_Barrel->Sumw2();
    MC_num_Barrel->Add(dy_num_Barrel);
    MC_num_Barrel->Add(Zjets_num_Barrel);
    MC_num_Barrel->Add(ttbar_num_Barrel);
    MC_num_Barrel->Add(WW_num_Barrel);
    MC_num_Barrel->Add(WZ_num_Barrel);
    MC_num_Barrel->Add(ZZ_num_Barrel);
    MC_num_Barrel->Add(Wantitop_num_Barrel);
    MC_num_Barrel->Add(tW_num_Barrel);
    total_num_Barrel->Add(MC_num_Barrel,-1);
    
    
    TH1F* MC_den_Barrel_fx = (TH1F*) Wjets_den_Barrel_fx->Clone();
    MC_den_Barrel_fx->Sumw2();
    MC_den_Barrel_fx->Add(dy_den_Barrel_fx);
    MC_den_Barrel_fx->Add(Zjets_den_Barrel_fx);
    MC_den_Barrel_fx->Add(ttbar_den_Barrel_fx);
    MC_den_Barrel_fx->Add(WW_den_Barrel_fx);
    MC_den_Barrel_fx->Add(WZ_den_Barrel_fx);
    MC_den_Barrel_fx->Add(ZZ_den_Barrel_fx);
    MC_den_Barrel_fx->Add(Wantitop_den_Barrel_fx);
    MC_den_Barrel_fx->Add(tW_den_Barrel_fx);
    total_den_Barrel_fx->Add(MC_den_Barrel_fx,-1);
    
    
    
    
    
    
    
    /*
     TH1F* total_num_Barrel_fx = (TH1F*) Data_num_Barrel_fx->Clone();
     total_num_Barrel_fx->Sumw2();
     TH1F* MC_num_Barrel_fx = (TH1F*) Wjets_num_Barrel_fx->Clone();
     MC_num_Barrel_fx->Sumw2();
     MC_num_Barrel_fx->Add(dy_num_Barrel_fx);
     MC_num_Barrel_fx->Add(Zjets_num_Barrel_fx);
     MC_num_Barrel_fx->Add(ttbar_num_Barrel_fx);
     MC_num_Barrel_fx->Add(WW_num_Barrel_fx);
     MC_num_Barrel_fx->Add(WZ_num_Barrel_fx);
     MC_num_Barrel_fx->Add(ZZ_num_Barrel_fx);
     MC_num_Barrel_fx->Add(Wantitop_num_Barrel_fx);
     MC_num_Barrel_fx->Add(tW_num_Barrel_fx);
     total_num_Barrel_fx->Add(MC_num_Barrel_fx,-1);
     
     TH1F * num_Barrel_avg = (TH1F*) total_num_Barrel_fx->Clone();
     num_Barrel_avg->Sumw2();
     num_Barrel_avg->Divide(total_num_Barrel);
     */
    
    
    
    TH1F* FR_pt_Barrel= (TH1F*) total_num_Barrel->Clone();
    FR_pt_Barrel->Sumw2();
    
    std::cout<<"FR info"<<std::endl;
    
    double dn_B=0.0, mn_B=0.0, dd_B=0.0, md_B=0.0 ;
    double num_B=0.0;
    double num_err_B = 0.0;
    double den_B=0.0;
    double den_err_B = 0.0;
    double r_B =0.0;
    double err_B = 0.0;
    
    for(int i=1; i<=Barrel_bins; i++){
        
        dn_B = Data_num_Barrel->GetBinContent(i);
        dd_B = Data_den_Barrel->GetBinContent(i);
        mn_B = MC_num_Barrel->GetBinError(i)*MC_num_Barrel->GetBinError(i);
        md_B = MC_den_Barrel->GetBinError(i)*MC_den_Barrel->GetBinError(i);
        num_B = abs(total_num_Barrel->GetBinContent(i));
        den_B = abs(total_den_Barrel->GetBinContent(i));
        num_err_B = total_num_Barrel->GetBinError(i);
        den_err_B = total_den_Barrel->GetBinError(i);
        r_B = num_B/den_B;
        //double err = r * sqrt(((Data_num->GetBinContent(i)+MC_num->GetBinContent(i))/(total_num->GetBinContent(i)*total_num->GetBinContent(i))) + ((Data_den->GetBinContent(i)+MC_den->GetBinContent(i))/(total_den->GetBinContent(i)*total_den->GetBinContent(i))));
        //double err = r * sqrt((num_err*num_err)/(num*num) + (den_err*den_err)/(den*den));
        
        err_B = r_B * sqrt(((dn_B+mn_B)/(num_B*num_B)) + ((dd_B+md_B)/(den_B*den_B)) - 2*(dn_B+mn_B)/(num_B*den_B));
        
        std::cout<<i<<" "<<num_B<<"   "<<total_num_Barrel->GetBinError(i)<<"   "<<den_B<<"   "<<total_den_Barrel->GetBinError(i)<<"   "<<dd_B<<"    "<<Data_den_Barrel->GetBinError(i)<<"  "<<md_B<<"    "<<MC_den_Barrel->GetBinError(i)<<"    "<<r_B<<" "<<err_B<<std::endl;
        FR_pt_Barrel->SetBinContent(i,r_B);
        FR_pt_Barrel->SetBinError(i,err_B);
        
        FR_Barrel[i-1] = r_B;
        FR_low_Barrel[i-1] = err_B;
        FR_high_Barrel[i-1] = err_B;
        avg_pt_low_Barrel[i-1] = avg_pt_Barrel[i-1] - (FR_pt_Barrel->GetXaxis()->GetBinLowEdge(i));
        avg_pt_high_Barrel[i-1] = ((FR_pt_Barrel->GetXaxis()->GetBinLowEdge(i)) + (FR_pt_Barrel->GetXaxis()->GetBinWidth(i))) - avg_pt_Barrel[i-1];
        if(i==4) avg_pt_high_Barrel[i-1]= 5000 - avg_pt_Barrel[i-1];
        
        
    }
    
    
    
    
    
    
    
    std::cout<<"Bin"<<" "<<"Avg pt"<<"   "<<"FR"<<" "<<"FR error"<<"    "<<"pt low"<<"  "<<"pt high"<<std::endl;
    for(int i=1; i<=Barrel_bins; i++) {
        std::cout<<i<<" "<<avg_pt_Barrel[i-1]<<"   "<<FR_Barrel[i-1]<<" "<<FR_low_Barrel[i-1]<<"    "<<avg_pt_low_Barrel[i-1]<<"    "<<avg_pt_high_Barrel[i-1]<<std::endl;
    }
    
    FR_pt_Barrel->SetMarkerStyle(20);
    FR_pt_Barrel->SetMarkerColor(kBlack);
    FR_pt_Barrel->SetMarkerSize(0.95);
    
    
    FR_pt_Barrel->GetYaxis()->SetTitle("Fake Rate");
    FR_pt_Barrel->GetXaxis()->SetTitle("pT(#mu)[GeV]");
    FR_pt_Barrel->SetTitle("Fake_Rate_Barrel");
    FR_pt_Barrel->GetYaxis()->SetRangeUser(0, 1.5);
    c1->SetGrid();
    
    FR_pt_Barrel->Draw("PE1");
    
    c1->Update();
    
    TPaveText* tText1 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText1->SetBorderSize(0);
    tText1->SetFillColor(0);
    tText1->SetFillStyle(0);
    TText *t1 = tText1->AddText("CMS 2017,2018 Data 103.4 fb^{-1} (13TeV)");
    tText1->SetTextSize(0.035);
    tText1->Draw();
    c1->Update();
    
    
    // FR_pt_Barrel->SaveAs("Fake_Rate_Barrel.root","root");
    //c1->SaveAs("Fake_Rate_Barrel.pdf","pdf");
    //c1->SaveAs("Fake_Rate_Barrel.png","png");
    
    //////////////////////////Barrel//////////////////////////////
    
    
    
    
    
    //////////////////////////Endcap//////////////////////////////
    
    TCanvas *c2 = new TCanvas("rate", "rate", 1000, 800);
    c2->cd();
    
    
    
    TH1F* total_den_Endcap = (TH1F*) Data_den_Endcap->Clone();
    total_den_Endcap->Sumw2();
    
    /////////////////////For average/////////////////////////////
    
    TH1F* total_den_Endcap_fx = (TH1F*) Data_den_Endcap_fx->Clone();
    total_den_Endcap_fx->Sumw2();
    TH1F * den_Endcap_avg = (TH1F*) total_den_Endcap_fx->Clone();
    den_Endcap_avg->Sumw2();
    den_Endcap_avg->Divide(total_den_Endcap);
    
    //const int Endcap_bins = den_Endcap_avg->GetNbinsX();
    const int Endcap_bins = 4;
    std::cout<<"Number of bins : "<<Endcap_bins<<std::endl;
    double avg_pt_Endcap[Endcap_bins] = {0.0};
    double FR_Endcap[Endcap_bins] = {0.0};
    double FR_low_Endcap[Endcap_bins] = {0.0};
    double FR_high_Endcap[Endcap_bins] = {0.0};
    double avg_pt_low_Endcap[Endcap_bins] = {0.0};
    double avg_pt_high_Endcap[Endcap_bins] = {0.0};
    
    
    for(int i=1; i<=Endcap_bins; i++ ){
        avg_pt_Endcap[i-1] = den_Endcap_avg->GetBinContent(i);
        std::cout<<"bin : "<<i<<"   avg_pt : "<<avg_pt_Endcap[i-1]<<std::endl;
    }
    /////////////////////For average/////////////////////////////
    
    TH1F* MC_den_Endcap = (TH1F*) Wjets_den_Endcap->Clone();
    MC_den_Endcap->Sumw2();
    MC_den_Endcap->Add(dy_den_Endcap);
    MC_den_Endcap->Add(Zjets_den_Endcap);
    MC_den_Endcap->Add(ttbar_den_Endcap);
    MC_den_Endcap->Add(WW_den_Endcap);
    MC_den_Endcap->Add(WZ_den_Endcap);
    MC_den_Endcap->Add(ZZ_den_Endcap);
    MC_den_Endcap->Add(Wantitop_den_Endcap);
    MC_den_Endcap->Add(tW_den_Endcap);
    total_den_Endcap->Add(MC_den_Endcap,-1);
    
    
    
    TH1F* total_num_Endcap = (TH1F*) Data_num_Endcap->Clone();
    total_num_Endcap->Sumw2();
    TH1F* MC_num_Endcap = (TH1F*) Wjets_num_Endcap->Clone();
    MC_num_Endcap->Sumw2();
    MC_num_Endcap->Add(dy_num_Endcap);
    MC_num_Endcap->Add(Zjets_num_Endcap);
    MC_num_Endcap->Add(ttbar_num_Endcap);
    MC_num_Endcap->Add(WW_num_Endcap);
    MC_num_Endcap->Add(WZ_num_Endcap);
    MC_num_Endcap->Add(ZZ_num_Endcap);
    MC_num_Endcap->Add(Wantitop_num_Endcap);
    MC_num_Endcap->Add(tW_num_Endcap);
    total_num_Endcap->Add(MC_num_Endcap,-1);
    
    
    TH1F* MC_den_Endcap_fx = (TH1F*) Wjets_den_Endcap_fx->Clone();
    MC_den_Endcap_fx->Sumw2();
    MC_den_Endcap_fx->Add(dy_den_Endcap_fx);
    MC_den_Endcap_fx->Add(Zjets_den_Endcap_fx);
    MC_den_Endcap_fx->Add(ttbar_den_Endcap_fx);
    MC_den_Endcap_fx->Add(WW_den_Endcap_fx);
    MC_den_Endcap_fx->Add(WZ_den_Endcap_fx);
    MC_den_Endcap_fx->Add(ZZ_den_Endcap_fx);
    MC_den_Endcap_fx->Add(Wantitop_den_Endcap_fx);
    MC_den_Endcap_fx->Add(tW_den_Endcap_fx);
    total_den_Endcap_fx->Add(MC_den_Endcap_fx,-1);
    
    
    
    
    
    
    
    /*
     TH1F* total_num_Endcap_fx = (TH1F*) Data_num_Endcap_fx->Clone();
     total_num_Endcap_fx->Sumw2();
     TH1F* MC_num_Endcap_fx = (TH1F*) Wjets_num_Endcap_fx->Clone();
     MC_num_Endcap_fx->Sumw2();
     MC_num_Endcap_fx->Add(dy_num_Endcap_fx);
     MC_num_Endcap_fx->Add(Zjets_num_Endcap_fx);
     MC_num_Endcap_fx->Add(ttbar_num_Endcap_fx);
     MC_num_Endcap_fx->Add(WW_num_Endcap_fx);
     MC_num_Endcap_fx->Add(WZ_num_Endcap_fx);
     MC_num_Endcap_fx->Add(ZZ_num_Endcap_fx);
     MC_num_Endcap_fx->Add(Wantitop_num_Endcap_fx);
     MC_num_Endcap_fx->Add(tW_num_Endcap_fx);
     total_num_Endcap_fx->Add(MC_num_Endcap_fx,-1);
     
     TH1F * num_Endcap_avg = (TH1F*) total_num_Endcap_fx->Clone();
     num_Endcap_avg->Sumw2();
     num_Endcap_avg->Divide(total_num_Endcap);
     */
    
    
    TH1F* FR_pt_Endcap= (TH1F*) total_num_Endcap->Clone();
    FR_pt_Endcap->Sumw2();
    
    std::cout<<"FR info"<<std::endl;
    
    double dn_E=0.0, mn_E=0.0, dd_E=0.0, md_E=0.0 ;
    double num_E=0.0;
    double num_err_E = 0.0;
    double den_E=0.0;
    double den_err_E = 0.0;
    double r_E =0.0;
    double err_E = 0.0;
    
    for(int i=1; i<=Endcap_bins; i++){
        
        dn_E = Data_num_Endcap->GetBinContent(i);
        dd_E = Data_den_Endcap->GetBinContent(i);
        mn_E = MC_num_Endcap->GetBinError(i)*MC_num_Endcap->GetBinError(i);
        md_E = MC_den_Endcap->GetBinError(i)*MC_den_Endcap->GetBinError(i);
        num_E = abs(total_num_Endcap->GetBinContent(i));
        den_E = abs(total_den_Endcap->GetBinContent(i));
        num_err_E = total_num_Endcap->GetBinError(i);
        den_err_E = total_den_Endcap->GetBinError(i);
        r_E = num_E/den_E;
        //double err = r * sqrt(((Data_num->GetBinContent(i)+MC_num->GetBinContent(i))/(total_num->GetBinContent(i)*total_num->GetBinContent(i))) + ((Data_den->GetBinContent(i)+MC_den->GetBinContent(i))/(total_den->GetBinContent(i)*total_den->GetBinContent(i))));
        //double err = r * sqrt((num_err*num_err)/(num*num) + (den_err*den_err)/(den*den));
        
        err_E = r_E * sqrt(((dn_E+mn_E)/(num_E*num_E)) + ((dd_E+md_E)/(den_E*den_E)) - 2*(dn_E+mn_E)/(num_E*den_E));
        
        std::cout<<i<<" "<<num_E<<"   "<<total_num_Endcap->GetBinError(i)<<"   "<<den_E<<"   "<<total_den_Endcap->GetBinError(i)<<"   "<<dd_E<<"    "<<Data_den_Endcap->GetBinError(i)<<"  "<<md_E<<"    "<<MC_den_Endcap->GetBinError(i)<<"    "<<r_E<<" "<<err_E<<std::endl;
        
        FR_pt_Endcap->SetBinContent(i,r_E);
        FR_pt_Endcap->SetBinError(i,err_E);
        
        FR_Endcap[i-1] = r_E;
        FR_low_Endcap[i-1] = err_E;
        FR_high_Endcap[i-1] = err_E;
        avg_pt_low_Endcap[i-1] = avg_pt_Endcap[i-1] - (FR_pt_Endcap->GetXaxis()->GetBinLowEdge(i));
        avg_pt_high_Endcap[i-1] = ((FR_pt_Endcap->GetXaxis()->GetBinLowEdge(i)) + (FR_pt_Endcap->GetXaxis()->GetBinWidth(i))) - avg_pt_Endcap[i-1];
        if(i==4) avg_pt_high_Endcap[i-1]= 5000 - avg_pt_Endcap[i-1];
        
        
    }
    
    
    
    
    
    
    std::cout<<"Bin"<<" "<<"Avg pt"<<"   "<<"FR"<<" "<<"FR error"<<"    "<<"pt low"<<"  "<<"pt high"<<std::endl;
    for(int i=1; i<=Endcap_bins; i++) {
        std::cout<<i<<" "<<avg_pt_Endcap[i-1]<<"   "<<FR_Endcap[i-1]<<" "<<FR_low_Endcap[i-1]<<"    "<<avg_pt_low_Endcap[i-1]<<"    "<<avg_pt_high_Endcap[i-1]<<std::endl;
    }
    
    FR_pt_Endcap->SetMarkerStyle(20);
    FR_pt_Endcap->SetMarkerColor(kBlack);
    FR_pt_Endcap->SetMarkerSize(0.95);
    
    
    FR_pt_Endcap->GetYaxis()->SetTitle("Fake Rate");
    FR_pt_Endcap->GetXaxis()->SetTitle("pT(#mu)[GeV]");
    FR_pt_Endcap->SetTitle("Fake_Rate_Endcap");
    FR_pt_Endcap->GetYaxis()->SetRangeUser(0, 1.5);
    c2->SetGrid();
    
    FR_pt_Endcap->Draw("PE1");
    c2->Update();
    
    TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText2->SetBorderSize(0);
    tText2->SetFillColor(0);
    tText2->SetFillStyle(0);
    TText *t2 = tText2->AddText("CMS 2017,2018 Data 103.4 fb^{-1} (13TeV)");
    tText2->SetTextSize(0.035);
    tText2->Draw();
    c2->Update();
    
    
    //FR_pt_Endcap->SaveAs("Fake_Rate_Endcap.root","root");
    // c2->SaveAs("Fake_Rate_Endcap.pdf","pdf");
    //c2->SaveAs("Fake_Rate_Endcap.png","png");
    
    //////////////////////////Endcap//////////////////////////////
    
    
    
    
    
    /////////////////////////////////////TGraps/////////////////////////////////////////////////////////////
    
    
    /////////////////////////////////////Fitting/////////////////////////////////////////////////
    
    TF1 *brl = new TF1("brl","[0]/(1+[1]*TMath::Exp(-x*[2]))",0,5000);
    brl->SetParLimits(0,0.545,0.621);
    brl->SetParLimits(1,9.102,17.434);
    brl->SetParLimits(2,0.0056,0.0068);
   
    TF1 *endcp = new TF1("endcp","landau",0,1200);
    TF1 *constant  = new TF1("constant","0.2",1200,5000);
    
    
    
    TCanvas *c3 = new TCanvas("c3","Fake_Rate_Barrel",1000,800);
    c3->SetGrid();
    
    TPad* pad1 = new TPad("pad1", "pad1",  0, 0, 1, 1);
    pad1->Draw();
    pad1->cd();
    pad1->SetTicks();
    pad1->SetGrid();
    
    // create TGraphAsymmErrors with the arrays
    TGraphAsymmErrors* gr_Barrel = new TGraphAsymmErrors(Barrel_bins, avg_pt_Barrel ,FR_Barrel, avg_pt_low_Barrel, avg_pt_high_Barrel, FR_low_Barrel, FR_high_Barrel);
    gr_Barrel->SetTitle(" ");
    gr_Barrel->SetMarkerColor(4);
    gr_Barrel->SetMarkerStyle(21);
    gr_Barrel->GetYaxis()->SetTitle("Fake Rate");
    gr_Barrel->GetXaxis()->SetTitle("p_{T}(#mu)[GeV]");
    gr_Barrel->Draw("AP");
    gr_Barrel->Fit("brl", "RW");
    brl->Draw("SAME");
    gr_Barrel->GetYaxis()->SetRangeUser(0.0,2.0);
    
    c3->Update();
    
    TLegend *l1 = new TLegend(0.2,0.6,0.6,0.8);
    l1->SetBorderSize(0);
    l1->AddEntry(gr_Barrel, "Muon Barrel", "lep");
    l1->AddEntry(brl, "Fit function : #frac{p0}{1 + p1*exp(-p2*p_{T})}", "l");
    
    
    l1->Draw();
    c3->Update();
    
    TPaveText* tText3 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText3->SetBorderSize(0);
    tText3->SetFillColor(0);
    tText3->SetFillStyle(0);
    TText *t3 = tText3->AddText("CMS 2017,2018 Data 103.4 fb^{-1} (13TeV)");
    tText3->SetTextSize(0.035);
    tText3->Draw();
    pad1->Update();
    c3->Update();
    
    
    c3->Update();
    
    //gr_Barrel->SaveAs("gr_Barrel.root","root");
    
    
    //c3->SaveAs("gr_Barrel.pdf","pdf");
    // c3->SaveAs("gr_Barrel.png","png");
    
    
    /////////////////////////////////////////////////////////////////////////////////////////
    
    TCanvas *c4 = new TCanvas("c4","Fake_Rate_Endcap",1000,800);
    c4->SetGrid();
    
    TPad* pad2 = new TPad("pad1", "pad1",  0, 0, 1, 1);
    pad2->Draw();
    pad2->cd();
    pad2->SetTicks();
    pad2->SetGrid();
    
    // create TGraphAsymmErrors with the arrays
    TGraphAsymmErrors* gr_Endcap = new TGraphAsymmErrors(Endcap_bins, avg_pt_Endcap ,FR_Endcap, avg_pt_low_Endcap, avg_pt_high_Endcap, FR_low_Endcap, FR_high_Endcap);
    gr_Endcap->SetTitle(" ");
    gr_Endcap->SetMarkerColor(4);
    gr_Endcap->SetMarkerStyle(21);
    gr_Endcap->GetYaxis()->SetTitle("Fake Rate");
    gr_Endcap->GetXaxis()->SetTitle("p_{T}(#mu)[GeV]");
    gr_Endcap->Draw("AP");
    gr_Endcap->Fit("endcp", "RW");
    endcp->Draw("SAME");
    gr_Endcap->Fit("constant", "RW");
    constant->Draw("SAME");
    gr_Endcap->GetYaxis()->SetRangeUser(0.0,2.0);
    
    c4->Update();
    
    TLegend *l2 = new TLegend(0.2,0.6,0.6,0.8);
    l2->SetBorderSize(0);
    l2->AddEntry(gr_Endcap, "Muon Endcap", "lep");
    l2->AddEntry(brl, "Fit function : landau", "l");
    
    
    l2->Draw();
    c4->Update();
    
    TPaveText* tText4 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText4->SetBorderSize(0);
    tText4->SetFillColor(0);
    tText4->SetFillStyle(0);
    TText *t4 = tText4->AddText("CMS 2017,2018 Data 103.4 fb^{-1} (13TeV)");
    tText4->SetTextSize(0.035);
    tText4->Draw();
    
    pad2->Update();
    c4->Update();
    
    
    
    
    pad2->Update();
    c4->Update();
    
    
    gr_Endcap->SaveAs("gr_Endcap.root","root");
    c4->SaveAs("gr_Endcap.pdf","pdf");
    c4->SaveAs("gr_Endcap.png","png");
    
}









