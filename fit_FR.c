

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







void fit_FR(){
    
    
    /////////////////////////////////////Fitting Barrel/////////////////////////////////////////////////
    
    TF1 *brl = new TF1("brl","[0]/(1+[1]*TMath::Exp(-x*[2]))",0,5000);
    brl->SetParLimits(0,0.545,0.621);
    brl->SetParLimits(1,9.102,17.434);
    brl->SetParLimits(2,0.0056,0.0068);
    
    
    
    
    TCanvas *c3 = new TCanvas("c3","Fake_Rate_Barrel",1000,800);
    c3->SetGrid();
    
    TPad* pad1 = new TPad("pad1", "pad1",  0, 0, 1, 1);
    pad1->Draw();
    pad1->cd();
    pad1->SetTicks();
    pad1->SetGrid();
    
    // create TGraphAsymmErrors with the arrays
    
    const int Barrel_bins = 4;
    double avg_pt_Barrel[Barrel_bins] = {68.6058, 252.985, 481.872, 915.982};
    double FR_Barrel[Barrel_bins] = {0.0466836, 0.369449, 0.419311, 0.676781};
    double FR_low_Barrel[Barrel_bins] = {0.00920205, 0.041077, 0.143134, 0.161144};
    double FR_high_Barrel[Barrel_bins] = {0.00920205, 0.041077, 0.143134, 0.161144};
    double avg_pt_low_Barrel[Barrel_bins] = {18.6058, 52.9846, 81.8719, 215.982};
    double avg_pt_high_Barrel[Barrel_bins] = {131.394, 147.015, 218.128, 4084.02};
    
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
    
    gr_Barrel->SaveAs("gr_Barrel.root","root");
    
    
    c3->SaveAs("gr_Barrel.pdf","pdf");
    c3->SaveAs("gr_Barrel.png","png");
    
    
    //////////////////////////////////////////Fitting Endcap///////////////////////////////////////////////
    
    TF1 *endcp = new TF1("endcp","landau",0,1200);
    TF1 *constant  = new TF1("constant","0.2",1200,5000);
    
    
    TCanvas *c4 = new TCanvas("c4","Fake_Rate_Endcap",1000,800);
    c4->SetGrid();
    
    TPad* pad2 = new TPad("pad1", "pad1",  0, 0, 1, 1);
    pad2->Draw();
    pad2->cd();
    pad2->SetTicks();
    pad2->SetGrid();
    
    // create TGraphAsymmErrors with the arrays
    
    const int Endcap_bins = 4;
    std::cout<<"Number of bins : "<<Endcap_bins<<std::endl;
    double avg_pt_Endcap[Endcap_bins] = {67.0077, 251.203, 476.057, 1121.94};
    double FR_Endcap[Endcap_bins] = {0.194999, 0.527867, 0.477867, 0.22148};
    double FR_low_Endcap[Endcap_bins] = {0.0142154, 0.0557531, 0.091, 0.13};
    double FR_high_Endcap[Endcap_bins] = {0.0142154, 0.0557531, 0.091, 0.13};
    double avg_pt_low_Endcap[Endcap_bins] = {17.0077, 51.2029, 76.0565, 421.939};
    double avg_pt_high_Endcap[Endcap_bins] = {132.992, 148.797, 223.943, 3878.06};
    
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










