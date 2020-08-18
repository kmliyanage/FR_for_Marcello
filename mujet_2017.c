//
//  mujet.c
//  To estimate the muon+jet background to the SR coming from Fake Muons. Here we consider the dimon events(not single muons as in FR estimation calculations)
//
//  Created by Kalpanie Liyanage on 8/6/20.
//

#include <stdio.h>
#include "TFile.h"
#include <iostream>
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
#include "TPad.h"
#include <TString.h>
#include "TEfficiency.h"






void mujet(){
    
    
    TFile *f = new TFile("out.root", "RECREATE");
    
    //////////////////////FR fitting functions - 2017 & 2018 combined/////////////////////////
    TF1 *f_barrel = new TF1("f_barrel","0.621/(1 + (9.102*TMath::Exp(-0.0068*x)))",0,7000);
    TF1 *f_endcap = new TF1("f_endcap","3.011*TMath::Landau(x,411.8,197.3,false)",0,7000);
    
    float FR = 0.0;
    float W = 0.0;
    
    
    //////////////////////FR fitting functions - 2017 & 2018 combined/////////////////////////
    
    const int mc = 33;
    const int data = 5;
    const double PI = 3.141592654;
    
    
    
    //k-factor applied for DY samples only to accamodate NNPDF corrections and PI bkgs
    float gM = 1.0;
    float kFactor = 1.0;
    float kFactor_BB = 1.0;
    float kFactor_BE = 1.0;
    
    TString MC_samples[mc] =  {
        "qcd15to30",        //0
        "qcd30to50",        //1
        "qcd50to80",        //2
        "qcd80to120",       //3
        "qcd120to170",      //4
        "qcd170to300",      //5
        "qcd300to470",      //6
        "qcd470to600",      //7
        "qcd600to800",      //8
        "qcd800to1000",     //9
        "qcd1000to1400",    //10
        "qcd1400to1800",    //11
        "qcd1800to2400",    //12
        "qcd2400to3200",    //13
        "qcd3200toInf",     //14
        "Wjets",            //15
        "Zjets_filtered",   //16
        "ttbar",            //17
        "WW",               //18
        "WZ",               //19
        "ZZ",               //20
        "tW",               //21
        "Wantitop",         //22
        "dy50to120",    //23
        "dy120to200",   //24
        "dy200to400",   //25
        "dy400to800",   //26
        "dy800to1400",  //27
        "dy1400to2300", //28
        "dy2300to3500", //29
        "dy3500to4500", //30
        "dy4500to6000", //31
        "dy6000toInf",  //32
    };
    
    double events[mc] = {
        19986190,   //0 qcd15to30
        19873760,   //1 qcd30to50
        19176400,   //2 qcd50to80
        28430940,   //3 qcd80to120
        29854280,   //4 qcd120to170
        29829920,   //5 qcd170to300
        53798780,   //6 qcd300to470
        27881030,   //7 qcd470to600
        66134960,   //8 qcd600to800
        39529010,   //9 qcd800to1000
        19631810,   //10    qcd1000to1400
        5685270,    //11    qcd1400to1800
        2923941,    //12    qcd1800to2400
        1910526,    //13    qcd2400to3200
        757837,     //14    qcd3200toInf  151000
        30008250,   //15    Wjets
        15296830,  //16    Zjets
        960752,   //17    ttbar
        7765828,    //18    WW
        3928630,    //19    WZ
        1949766,    //20    ZZ
        4974435,    //21    tW
        5635539,    //22    Wantitop
        2863000,    //23    dy50to120
        100000,     //24    dy120to200
        100000,     //25    dy200to400
        100000,     //26    dy400to800
        100000,     //27    dy800to1400
        100000,     //28    dy1400to2300
        100000,     //29    dy2300to3500
        100000,     //30    dy3500to4500
        100000,     //31    dy4500to6000
        100000,     //32    dy6000toInf
        
    };
    
    //Using cross sections in XSDB
    double sigma[mc] = {
        1246000000.0,   //0 qcd15to30
        106900000.0,    //1 qcd30to50
        15710000.0,     //2 qcd50to80
        2336000.0,      //3 qcd80to120
        407300.0,       //4 qcd120to170
        103500.0,       //5 qcd170to300
        6830.0,         //6 qcd300to470
        552.1,          //7 qcd470to600
        156.5,          //8 qcd600to800
        26.28,          //9 qcd800to1000
        7.477,           //10    qcd1000to1400
        0.6484,         //11    qcd1400to1800
        0.0875,        //12    qcd1800to2400
        0.005236,       //13    qcd2400to3200
        0.0001357,      //14    qcd3200toInf
        61526.7,        //15    Wjets
        6077.22,        //16    Zjets
        88.29,          //17    ttbar
        118.7,           //18    WW
        50.2,           //19    WZ
        16.523,          //20    ZZ
        35.6,          //21    tW
        35.6,          //22    Wantitop
        2112.904,         //23    dy50to120
        20.553,          //24    dy120to200
        2.886,          //25    dy200to400
        0.2517,         //26    dy400to800
        0.01707,        //27    dy800to1400
        1.366E-3,       //28    dy1400to2300
        8.178E-5,       //29    dy2300to3500
        3.191E-6,       //30    dy3500to4500
        2.787E-7,       //31    dy4500to6000
        9.569E-9,       //32    dy6000toInf
    };
    
    //Using cross sections in XSDB and AN_2018_011
    /*  double sigma[mc] = {
     1246000000.0,   //0 qcd15to30
     106900000.0,    //1 qcd30to50
     15710000.0,     //2 qcd50to80
     2336000.0,      //3 qcd80to120
     407300.0,       //4 qcd120to170
     103500.0,       //5 qcd170to300
     6830.0,         //6 qcd300to470
     552.1,          //7 qcd470to600
     156.5,          //8 qcd600to800
     26.28,          //9 qcd800to1000
     7.47,           //10    qcd1000to1400
     0.6484,         //11    qcd1400to1800
     0.08743,        //12    qcd1800to2400
     0.005236,       //13    qcd2400to3200
     0.0001357,      //14    qcd3200toInf
     52940.0,        //15    Wjets
     5765.4,        //16    Zjets
     87.31,          //17    ttbar
     118.7,           //18    WW
     47.13,           //19    WZ
     16.523,          //20    ZZ
     35.6,          //21    tW
     35.6,          //22    Wantitop
     2112.904,         //23    dy50to120
     20.553,          //24    dy120to200
     2.886,          //25    dy200to400
     0.2517,         //26    dy400to800
     0.01707,        //27    dy800to1400
     1.366E-3,       //28    dy1400to2300
     8.178E-5,       //29    dy2300to3500
     3.191E-6,       //30    dy3500to4500
     2.787E-7,       //31    dy4500to6000
     9.569E-9,       //32    dy6000toInf
     };
     */
    
    
    
    
    
    
    TString samp[data] =  {
        "Run2017MuonsOnly_SingleMuonRun2017B-31Mar2018-v1",
        "Run2017MuonsOnly_SingleMuonRun2017C-17Nov2017-v1",
        "Run2017MuonsOnly_SingleMuonRun2017D-17Nov2017-v1",
        "Run2017MuonsOnly_SingleMuonRun2017E-17Nov2017-v1",
        "Run2017MuonsOnly_SingleMuonRun2017F-17Nov2017-v1",
    };
    
    
    TString DATA_samples[data] =  {
        "Run2017B",
        "Run2017C",
        "Run2017D",
        "Run2017E",
        "Run2017F",
    };
    
    double LUMINOSITY = 42135.25562;
    
    
    
    float dil_mass;
    float dil_pt;
    float dil_eta;
    float dil_rap;
    float dil_phi;
    float cos_angle;
    float vertex_chi2;
    int dil_chosen;
    float met_pt;
    float met_phi;
    
    float mt = -1.0; // for transverse mass
    float leading_pt = -1.0; // leading pt of the dimuon
    float leading_phi = -1.0;
    float delta_phi = -1.0;
    Int_t c_event;
    Int_t n_event;
    
    UInt_t event;
    UInt_t run;
    unsigned lumi;
    Int_t prev_event = -88;
    
    float lep_pt[2];
    float lep_phi[2];
    int lep_id[2];
    float lep_pt_err[2];
    float lep_eta[2];
    float lep_tk_pt[2];
    float lep_glb_pt[2];
    float lep_picky_pt[2];
    float lep_tpfms_pt[2];
    float lep_dB[2];
    float lep_tk_dz[2];
    float lep_tuneP_pt[2];
    float lep_sumPt[2];
    float lep_triggerMatchPt[2];
    short lep_glb_numberOfValidTrackerLayers[2];
    short lep_glb_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidMuonHits[2];
    short lep_TuneP_numberOfValidMuonHits[2];
    short lep_numberOfMatchedStations[2];
    short lep_expectedNnumberOfMatchedStations[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];
    bool GoodVtx;
    short lep_numberOfMatchedRPCLayers[2];
    unsigned int lep_stationMask[2];
    float vertex_m;
    float gen_dil_mass;
    float genWeight;
    float gen_lep_qOverPt[2];
    float lep_qOverPt[2]; //for data
    float gen_lep_eta[2];
    float gen_lep_pt[2];
    
    //ratios between DATA and MC in the range 60 - 120 GeV (for normalizing the background events to z peak), considering also the luminosity sections applied for all backgrounds
    float Z_peak = 1.0282;
    float Z_peak_BB = 1.0286;
    float Z_peak_BE = 1.0278;
    
    
    float weight[mc] = {0.0};
    float weight_BB[mc] = {0.0};
    float weight_BE[mc] = {0.0};
    
    for(int i = 0; i<mc; i++){ //Normalizing MC to the luminosity of DATA
        weight[i] = LUMINOSITY * (float) sigma[i] / (float) events[i];
        
        weight_BB[i] = weight[i] * (float) Z_peak_BB ;
        weight_BE[i] = weight[i] * (float) Z_peak_BE;
        weight[i] *= Z_peak;
        std::cout<<MC_samples[i]<<" "<<weight[i]<<"    "<<weight_BB[i]<<"    "<<weight_BE[i]<<std::endl;
    }
    
    int mass_bins=200;
    
    const int    NMBINS = 51;
    const double MMIN = 70., MMAX =4000.;
    double logMbins[NMBINS+1];
    
    for (int ibin = 0; ibin <= NMBINS; ibin++){
        logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
        //       if(ibin != 0 && logMbins[ibin] <= logMbins[ibin-1]) std::cout<<"PROBLEMA"<<std::endl;
        //ss  std::cout<<logMbins[ibin]<<",";
    }
    
    
    
    
    ///////////////Histograms/////////////////////////////////
    
    ////////////////////qcd///////////////////////////////////
    
    
    TH1F* qcd_mass_linear = new TH1F("qcd_mass_linear", "qcd_mass_linear;qcd_mass_linear;#Events", mass_bins, 70, 4000);
    qcd_mass_linear->Sumw2();
    TH1F* qcd_mass_linear_BB = new TH1F("qcd_mass_linear_BB", "qcd_mass_linear_BB;qcd_mass_linear_BB;#Events", mass_bins, 70, 4000);
    qcd_mass_linear_BB->Sumw2();
    TH1F* qcd_mass_linear_BE = new TH1F("qcd_mass_linear_BE", "qcd_mass_linear_BE;qcd_mass_linear_BE;#Events", mass_bins, 70, 4000);
    qcd_mass_linear_BE->Sumw2();
    
    TH1F* qcd_mass_log = new TH1F("qcd_mass_log", "qcd_mass_log;qcd_mass_log;#Events", NMBINS, logMbins);
    qcd_mass_log->Sumw2();
    TH1F* qcd_mass_log_BB = new TH1F("qcd_mass_log_BB", "qcd_mass_log_BB;qcd_mass_log_BB;#Events", NMBINS, logMbins);
    qcd_mass_log_BB->Sumw2();
    TH1F* qcd_mass_log_BE = new TH1F("qcd_mass_log_BE", "qcd_mass_log_BE;qcd_mass_log_BE;#Events", NMBINS, logMbins);
    qcd_mass_log_BE->Sumw2();
    
    TH1F* qcd_mass_cumulative = new TH1F("qcd_mass_cumulative", "qcd_mass_cumulative;qcd_mass_cumulative;#Events", mass_bins, 70, 4000);
    qcd_mass_cumulative->Sumw2();
    TH1F* qcd_mass_cumulative_BB = new TH1F("qcd_mass_cumulative_BB", "qcd_mass_cumulative_BB;qcd_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    qcd_mass_cumulative_BB->Sumw2();
    TH1F* qcd_mass_cumulative_BE = new TH1F("qcd_mass_cumulative_BE", "qcd_mass_cumulative_BE;qcd_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    qcd_mass_cumulative_BE->Sumw2();
    
    TH1F* qcd_mass_log_cumulative = new TH1F("qcd_mass_log_cumulative", "qcd_mass_log_cumulative;qcd_mass_log_cumulative;#Events", NMBINS, logMbins);
    qcd_mass_log_cumulative->Sumw2();
    TH1F* qcd_mass_log_cumulative_BB = new TH1F("qcd_mass_log_cumulative_BB", "qcd_mass_log_cumulative_BB;qcd_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    qcd_mass_log_cumulative_BB->Sumw2();
    TH1F* qcd_mass_log_cumulative_BE = new TH1F("qcd_mass_log_cumulative_BE", "qcd_mass_log_cumulative_BE;qcd_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    qcd_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////Wjets////////////////////////////////////////////
    TH1F* Wjets_mass_linear = new TH1F("Wjets_mass_linear", "Wjets_mass_linear;Wjets_mass_linear;#Events", mass_bins, 70, 4000);
    Wjets_mass_linear->Sumw2();
    TH1F* Wjets_mass_linear_BB = new TH1F("Wjets_mass_linear_BB", "Wjets_mass_linear_BB;Wjets_mass_linear_BB;#Events", mass_bins, 70, 4000);
    Wjets_mass_linear_BB->Sumw2();
    TH1F* Wjets_mass_linear_BE = new TH1F("Wjets_mass_linear_BE", "Wjets_mass_linear_BE;Wjets_mass_linear_BE;#Events", mass_bins, 70, 4000);
    Wjets_mass_linear_BE->Sumw2();
    
    TH1F* Wjets_mass_log = new TH1F("Wjets_mass_log", "Wjets_mass_log;Wjets_mass_log;#Events", NMBINS, logMbins);
    Wjets_mass_log->Sumw2();
    TH1F* Wjets_mass_log_BB = new TH1F("Wjets_mass_log_BB", "Wjets_mass_log_BB;Wjets_mass_log_BB;#Events", NMBINS, logMbins);
    Wjets_mass_log_BB->Sumw2();
    TH1F* Wjets_mass_log_BE = new TH1F("Wjets_mass_log_BE", "Wjets_mass_log_BE;Wjets_mass_log_BE;#Events", NMBINS, logMbins);
    Wjets_mass_log_BE->Sumw2();
    
    TH1F* Wjets_mass_cumulative = new TH1F("Wjets_mass_cumulative", "Wjets_mass_cumulative;Wjets_mass_cumulative;#Events", mass_bins, 70, 4000);
    Wjets_mass_cumulative->Sumw2();
    TH1F* Wjets_mass_cumulative_BB = new TH1F("Wjets_mass_cumulative_BB", "Wjets_mass_cumulative_BB;Wjets_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    Wjets_mass_cumulative_BB->Sumw2();
    TH1F* Wjets_mass_cumulative_BE = new TH1F("Wjets_mass_cumulative_BE", "Wjets_mass_cumulative_BE;Wjets_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    Wjets_mass_cumulative_BE->Sumw2();
    
    TH1F* Wjets_mass_log_cumulative = new TH1F("Wjets_mass_log_cumulative", "Wjets_mass_log_cumulative;Wjets_mass_log_cumulative;#Events", NMBINS, logMbins);
    Wjets_mass_log_cumulative->Sumw2();
    TH1F* Wjets_mass_log_cumulative_BB = new TH1F("Wjets_mass_log_cumulative_BB", "Wjets_mass_log_cumulative_BB;Wjets_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    Wjets_mass_log_cumulative_BB->Sumw2();
    TH1F* Wjets_mass_log_cumulative_BE = new TH1F("Wjets_mass_log_cumulative_BE", "Wjets_mass_log_cumulative_BE;Wjets_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    Wjets_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////Zjets////////////////////////////////////////////
    TH1F* Zjets_mass_linear = new TH1F("Zjets_mass_linear", "Zjets_mass_linear;Zjets_mass_linear;#Events", mass_bins, 70, 4000);
    Zjets_mass_linear->Sumw2();
    TH1F* Zjets_mass_linear_BB = new TH1F("Zjets_mass_linear_BB", "Zjets_mass_linear_BB;Zjets_mass_linear_BB;#Events", mass_bins, 70, 4000);
    Zjets_mass_linear_BB->Sumw2();
    TH1F* Zjets_mass_linear_BE = new TH1F("Zjets_mass_linear_BE", "Zjets_mass_linear_BE;Zjets_mass_linear_BE;#Events", mass_bins, 70, 4000);
    Zjets_mass_linear_BE->Sumw2();
    
    TH1F* Zjets_mass_log = new TH1F("Zjets_mass_log", "Zjets_mass_log;Zjets_mass_log;#Events", NMBINS, logMbins);
    Zjets_mass_log->Sumw2();
    TH1F* Zjets_mass_log_BB = new TH1F("Zjets_mass_log_BB", "Zjets_mass_log_BB;Zjets_mass_log_BB;#Events", NMBINS, logMbins);
    Zjets_mass_log_BB->Sumw2();
    TH1F* Zjets_mass_log_BE = new TH1F("Zjets_mass_log_BE", "Zjets_mass_log_BE;Zjets_mass_log_BE;#Events", NMBINS, logMbins);
    Zjets_mass_log_BE->Sumw2();
    
    TH1F* Zjets_mass_cumulative = new TH1F("Zjets_mass_cumulative", "Zjets_mass_cumulative;Zjets_mass_cumulative;#Events", mass_bins, 70, 4000);
    Zjets_mass_cumulative->Sumw2();
    TH1F* Zjets_mass_cumulative_BB = new TH1F("Zjets_mass_cumulative_BB", "Zjets_mass_cumulative_BB;Zjets_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    Zjets_mass_cumulative_BB->Sumw2();
    TH1F* Zjets_mass_cumulative_BE = new TH1F("Zjets_mass_cumulative_BE", "Zjets_mass_cumulative_BE;Zjets_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    Zjets_mass_cumulative_BE->Sumw2();
    
    TH1F* Zjets_mass_log_cumulative = new TH1F("Zjets_mass_log_cumulative", "Zjets_mass_log_cumulative;Zjets_mass_log_cumulative;#Events", NMBINS, logMbins);
    Zjets_mass_log_cumulative->Sumw2();
    TH1F* Zjets_mass_log_cumulative_BB = new TH1F("Zjets_mass_log_cumulative_BB", "Zjets_mass_log_cumulative_BB;Zjets_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    Zjets_mass_log_cumulative_BB->Sumw2();
    TH1F* Zjets_mass_log_cumulative_BE = new TH1F("Zjets_mass_log_cumulative_BE", "Zjets_mass_log_cumulative_BE;Zjets_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    Zjets_mass_log_cumulative_BE->Sumw2();
    
    
    
    ////////////////////////////////ttbar////////////////////////////////////////////
    TH1F* ttbar_mass_linear = new TH1F("ttbar_mass_linear", "ttbar_mass_linear;ttbar_mass_linear;#Events", mass_bins, 70, 4000);
    ttbar_mass_linear->Sumw2();
    TH1F* ttbar_mass_linear_BB = new TH1F("ttbar_mass_linear_BB", "ttbar_mass_linear_BB;ttbar_mass_linear_BB;#Events", mass_bins, 70, 4000);
    ttbar_mass_linear_BB->Sumw2();
    TH1F* ttbar_mass_linear_BE = new TH1F("ttbar_mass_linear_BE", "ttbar_mass_linear_BE;ttbar_mass_linear_BE;#Events", mass_bins, 70, 4000);
    ttbar_mass_linear_BE->Sumw2();
    
    TH1F* ttbar_mass_log = new TH1F("ttbar_mass_log", "ttbar_mass_log;ttbar_mass_log;#Events", NMBINS, logMbins);
    ttbar_mass_log->Sumw2();
    TH1F* ttbar_mass_log_BB = new TH1F("ttbar_mass_log_BB", "ttbar_mass_log_BB;ttbar_mass_log_BB;#Events", NMBINS, logMbins);
    ttbar_mass_log_BB->Sumw2();
    TH1F* ttbar_mass_log_BE = new TH1F("ttbar_mass_log_BE", "ttbar_mass_log_BE;ttbar_mass_log_BE;#Events", NMBINS, logMbins);
    ttbar_mass_log_BE->Sumw2();
    
    TH1F* ttbar_mass_cumulative = new TH1F("ttbar_mass_cumulative", "ttbar_mass_cumulative;ttbar_mass_cumulative;#Events", mass_bins, 70, 4000);
    ttbar_mass_cumulative->Sumw2();
    TH1F* ttbar_mass_cumulative_BB = new TH1F("ttbar_mass_cumulative_BB", "ttbar_mass_cumulative_BB;ttbar_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    ttbar_mass_cumulative_BB->Sumw2();
    TH1F* ttbar_mass_cumulative_BE = new TH1F("ttbar_mass_cumulative_BE", "ttbar_mass_cumulative_BE;ttbar_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    ttbar_mass_cumulative_BE->Sumw2();
    
    TH1F* ttbar_mass_log_cumulative = new TH1F("ttbar_mass_log_cumulative", "ttbar_mass_log_cumulative;ttbar_mass_log_cumulative;#Events", NMBINS, logMbins);
    ttbar_mass_log_cumulative->Sumw2();
    TH1F* ttbar_mass_log_cumulative_BB = new TH1F("ttbar_mass_log_cumulative_BB", "ttbar_mass_log_cumulative_BB;ttbar_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    ttbar_mass_log_cumulative_BB->Sumw2();
    TH1F* ttbar_mass_log_cumulative_BE = new TH1F("ttbar_mass_log_cumulative_BE", "ttbar_mass_log_cumulative_BE;ttbar_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    ttbar_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////WW////////////////////////////////////////////
    TH1F* WW_mass_linear = new TH1F("WW_mass_linear", "WW_mass_linear;WW_mass_linear;#Events", mass_bins, 70, 4000);
    WW_mass_linear->Sumw2();
    TH1F* WW_mass_linear_BB = new TH1F("WW_mass_linear_BB", "WW_mass_linear_BB;WW_mass_linear_BB;#Events", mass_bins, 70, 4000);
    WW_mass_linear_BB->Sumw2();
    TH1F* WW_mass_linear_BE = new TH1F("WW_mass_linear_BE", "WW_mass_linear_BE;WW_mass_linear_BE;#Events", mass_bins, 70, 4000);
    WW_mass_linear_BE->Sumw2();
    
    TH1F* WW_mass_log = new TH1F("WW_mass_log", "WW_mass_log;WW_mass_log;#Events", NMBINS, logMbins);
    WW_mass_log->Sumw2();
    TH1F* WW_mass_log_BB = new TH1F("WW_mass_log_BB", "WW_mass_log_BB;WW_mass_log_BB;#Events", NMBINS, logMbins);
    WW_mass_log_BB->Sumw2();
    TH1F* WW_mass_log_BE = new TH1F("WW_mass_log_BE", "WW_mass_log_BE;WW_mass_log_BE;#Events", NMBINS, logMbins);
    WW_mass_log_BE->Sumw2();
    
    TH1F* WW_mass_cumulative = new TH1F("WW_mass_cumulative", "WW_mass_cumulative;WW_mass_cumulative;#Events", mass_bins, 70, 4000);
    WW_mass_cumulative->Sumw2();
    TH1F* WW_mass_cumulative_BB = new TH1F("WW_mass_cumulative_BB", "WW_mass_cumulative_BB;WW_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    WW_mass_cumulative_BB->Sumw2();
    TH1F* WW_mass_cumulative_BE = new TH1F("WW_mass_cumulative_BE", "WW_mass_cumulative_BE;WW_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    WW_mass_cumulative_BE->Sumw2();
    
    TH1F* WW_mass_log_cumulative = new TH1F("WW_mass_log_cumulative", "WW_mass_log_cumulative;WW_mass_log_cumulative;#Events", NMBINS, logMbins);
    WW_mass_log_cumulative->Sumw2();
    TH1F* WW_mass_log_cumulative_BB = new TH1F("WW_mass_log_cumulative_BB", "WW_mass_log_cumulative_BB;WW_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    WW_mass_log_cumulative_BB->Sumw2();
    TH1F* WW_mass_log_cumulative_BE = new TH1F("WW_mass_log_cumulative_BE", "WW_mass_log_cumulative_BE;WW_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    WW_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////WZ////////////////////////////////////////////
    TH1F* WZ_mass_linear = new TH1F("WZ_mass_linear", "WZ_mass_linear;WZ_mass_linear;#Events", mass_bins, 70, 4000);
    WZ_mass_linear->Sumw2();
    TH1F* WZ_mass_linear_BB = new TH1F("WZ_mass_linear_BB", "WZ_mass_linear_BB;WZ_mass_linear_BB;#Events", mass_bins, 70, 4000);
    WZ_mass_linear_BB->Sumw2();
    TH1F* WZ_mass_linear_BE = new TH1F("WZ_mass_linear_BE", "WZ_mass_linear_BE;WZ_mass_linear_BE;#Events", mass_bins, 70, 4000);
    WZ_mass_linear_BE->Sumw2();
    
    TH1F* WZ_mass_log = new TH1F("WZ_mass_log", "WZ_mass_log;WZ_mass_log;#Events", NMBINS, logMbins);
    WZ_mass_log->Sumw2();
    TH1F* WZ_mass_log_BB = new TH1F("WZ_mass_log_BB", "WZ_mass_log_BB;WZ_mass_log_BB;#Events", NMBINS, logMbins);
    WZ_mass_log_BB->Sumw2();
    TH1F* WZ_mass_log_BE = new TH1F("WZ_mass_log_BE", "WZ_mass_log_BE;WZ_mass_log_BE;#Events", NMBINS, logMbins);
    WZ_mass_log_BE->Sumw2();
    
    TH1F* WZ_mass_cumulative = new TH1F("WZ_mass_cumulative", "WZ_mass_cumulative;WZ_mass_cumulative;#Events", mass_bins, 70, 4000);
    WZ_mass_cumulative->Sumw2();
    TH1F* WZ_mass_cumulative_BB = new TH1F("WZ_mass_cumulative_BB", "WZ_mass_cumulative_BB;WZ_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    WZ_mass_cumulative_BB->Sumw2();
    TH1F* WZ_mass_cumulative_BE = new TH1F("WZ_mass_cumulative_BE", "WZ_mass_cumulative_BE;WZ_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    WZ_mass_cumulative_BE->Sumw2();
    
    TH1F* WZ_mass_log_cumulative = new TH1F("WZ_mass_log_cumulative", "WZ_mass_log_cumulative;WZ_mass_log_cumulative;#Events", NMBINS, logMbins);
    WZ_mass_log_cumulative->Sumw2();
    TH1F* WZ_mass_log_cumulative_BB = new TH1F("WZ_mass_log_cumulative_BB", "WZ_mass_log_cumulative_BB;WZ_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    WZ_mass_log_cumulative_BB->Sumw2();
    TH1F* WZ_mass_log_cumulative_BE = new TH1F("WZ_mass_log_cumulative_BE", "WZ_mass_log_cumulative_BE;WZ_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    WZ_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////ZZ////////////////////////////////////////////
    TH1F* ZZ_mass_linear = new TH1F("ZZ_mass_linear", "ZZ_mass_linear;ZZ_mass_linear;#Events", mass_bins, 70, 4000);
    ZZ_mass_linear->Sumw2();
    TH1F* ZZ_mass_linear_BB = new TH1F("ZZ_mass_linear_BB", "ZZ_mass_linear_BB;ZZ_mass_linear_BB;#Events", mass_bins, 70, 4000);
    ZZ_mass_linear_BB->Sumw2();
    TH1F* ZZ_mass_linear_BE = new TH1F("ZZ_mass_linear_BE", "ZZ_mass_linear_BE;ZZ_mass_linear_BE;#Events", mass_bins, 70, 4000);
    ZZ_mass_linear_BE->Sumw2();
    
    TH1F* ZZ_mass_log = new TH1F("ZZ_mass_log", "ZZ_mass_log;ZZ_mass_log;#Events", NMBINS, logMbins);
    ZZ_mass_log->Sumw2();
    TH1F* ZZ_mass_log_BB = new TH1F("ZZ_mass_log_BB", "ZZ_mass_log_BB;ZZ_mass_log_BB;#Events", NMBINS, logMbins);
    ZZ_mass_log_BB->Sumw2();
    TH1F* ZZ_mass_log_BE = new TH1F("ZZ_mass_log_BE", "ZZ_mass_log_BE;ZZ_mass_log_BE;#Events", NMBINS, logMbins);
    ZZ_mass_log_BE->Sumw2();
    
    TH1F* ZZ_mass_cumulative = new TH1F("ZZ_mass_cumulative", "ZZ_mass_cumulative;ZZ_mass_cumulative;#Events", mass_bins, 70, 4000);
    ZZ_mass_cumulative->Sumw2();
    TH1F* ZZ_mass_cumulative_BB = new TH1F("ZZ_mass_cumulative_BB", "ZZ_mass_cumulative_BB;ZZ_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    ZZ_mass_cumulative_BB->Sumw2();
    TH1F* ZZ_mass_cumulative_BE = new TH1F("ZZ_mass_cumulative_BE", "ZZ_mass_cumulative_BE;ZZ_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    ZZ_mass_cumulative_BE->Sumw2();
    
    TH1F* ZZ_mass_log_cumulative = new TH1F("ZZ_mass_log_cumulative", "ZZ_mass_log_cumulative;ZZ_mass_log_cumulative;#Events", NMBINS, logMbins);
    ZZ_mass_log_cumulative->Sumw2();
    TH1F* ZZ_mass_log_cumulative_BB = new TH1F("ZZ_mass_log_cumulative_BB", "ZZ_mass_log_cumulative_BB;ZZ_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    ZZ_mass_log_cumulative_BB->Sumw2();
    TH1F* ZZ_mass_log_cumulative_BE = new TH1F("ZZ_mass_log_cumulative_BE", "ZZ_mass_log_cumulative_BE;ZZ_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    ZZ_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////tW////////////////////////////////////////////
    TH1F* tW_mass_linear = new TH1F("tW_mass_linear", "tW_mass_linear;tW_mass_linear;#Events", mass_bins, 70, 4000);
    tW_mass_linear->Sumw2();
    TH1F* tW_mass_linear_BB = new TH1F("tW_mass_linear_BB", "tW_mass_linear_BB;tW_mass_linear_BB;#Events", mass_bins, 70, 4000);
    tW_mass_linear_BB->Sumw2();
    TH1F* tW_mass_linear_BE = new TH1F("tW_mass_linear_BE", "tW_mass_linear_BE;tW_mass_linear_BE;#Events", mass_bins, 70, 4000);
    tW_mass_linear_BE->Sumw2();
    
    TH1F* tW_mass_log = new TH1F("tW_mass_log", "tW_mass_log;tW_mass_log;#Events", NMBINS, logMbins);
    tW_mass_log->Sumw2();
    TH1F* tW_mass_log_BB = new TH1F("tW_mass_log_BB", "tW_mass_log_BB;tW_mass_log_BB;#Events", NMBINS, logMbins);
    tW_mass_log_BB->Sumw2();
    TH1F* tW_mass_log_BE = new TH1F("tW_mass_log_BE", "tW_mass_log_BE;tW_mass_log_BE;#Events", NMBINS, logMbins);
    tW_mass_log_BE->Sumw2();
    
    TH1F* tW_mass_cumulative = new TH1F("tW_mass_cumulative", "tW_mass_cumulative;tW_mass_cumulative;#Events", mass_bins, 70, 4000);
    tW_mass_cumulative->Sumw2();
    TH1F* tW_mass_cumulative_BB = new TH1F("tW_mass_cumulative_BB", "tW_mass_cumulative_BB;tW_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    tW_mass_cumulative_BB->Sumw2();
    TH1F* tW_mass_cumulative_BE = new TH1F("tW_mass_cumulative_BE", "tW_mass_cumulative_BE;tW_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    tW_mass_cumulative_BE->Sumw2();
    
    TH1F* tW_mass_log_cumulative = new TH1F("tW_mass_log_cumulative", "tW_mass_log_cumulative;tW_mass_log_cumulative;#Events", NMBINS, logMbins);
    tW_mass_log_cumulative->Sumw2();
    TH1F* tW_mass_log_cumulative_BB = new TH1F("tW_mass_log_cumulative_BB", "tW_mass_log_cumulative_BB;tW_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    tW_mass_log_cumulative_BB->Sumw2();
    TH1F* tW_mass_log_cumulative_BE = new TH1F("tW_mass_log_cumulative_BE", "tW_mass_log_cumulative_BE;tW_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    tW_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////Wantitop////////////////////////////////////////////
    TH1F* Wantitop_mass_linear = new TH1F("Wantitop_mass_linear", "Wantitop_mass_linear;Wantitop_mass_linear;#Events", mass_bins, 70, 4000);
    Wantitop_mass_linear->Sumw2();
    TH1F* Wantitop_mass_linear_BB = new TH1F("Wantitop_mass_linear_BB", "Wantitop_mass_linear_BB;Wantitop_mass_linear_BB;#Events", mass_bins, 70, 4000);
    Wantitop_mass_linear_BB->Sumw2();
    TH1F* Wantitop_mass_linear_BE = new TH1F("Wantitop_mass_linear_BE", "Wantitop_mass_linear_BE;Wantitop_mass_linear_BE;#Events", mass_bins, 70, 4000);
    Wantitop_mass_linear_BE->Sumw2();
    
    TH1F* Wantitop_mass_log = new TH1F("Wantitop_mass_log", "Wantitop_mass_log;Wantitop_mass_log;#Events", NMBINS, logMbins);
    Wantitop_mass_log->Sumw2();
    TH1F* Wantitop_mass_log_BB = new TH1F("Wantitop_mass_log_BB", "Wantitop_mass_log_BB;Wantitop_mass_log_BB;#Events", NMBINS, logMbins);
    Wantitop_mass_log_BB->Sumw2();
    TH1F* Wantitop_mass_log_BE = new TH1F("Wantitop_mass_log_BE", "Wantitop_mass_log_BE;Wantitop_mass_log_BE;#Events", NMBINS, logMbins);
    Wantitop_mass_log_BE->Sumw2();
    
    TH1F* Wantitop_mass_cumulative = new TH1F("Wantitop_mass_cumulative", "Wantitop_mass_cumulative;Wantitop_mass_cumulative;#Events", mass_bins, 70, 4000);
    Wantitop_mass_cumulative->Sumw2();
    TH1F* Wantitop_mass_cumulative_BB = new TH1F("Wantitop_mass_cumulative_BB", "Wantitop_mass_cumulative_BB;Wantitop_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    Wantitop_mass_cumulative_BB->Sumw2();
    TH1F* Wantitop_mass_cumulative_BE = new TH1F("Wantitop_mass_cumulative_BE", "Wantitop_mass_cumulative_BE;Wantitop_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    Wantitop_mass_cumulative_BE->Sumw2();
    
    TH1F* Wantitop_mass_log_cumulative = new TH1F("Wantitop_mass_log_cumulative", "Wantitop_mass_log_cumulative;Wantitop_mass_log_cumulative;#Events", NMBINS, logMbins);
    Wantitop_mass_log_cumulative->Sumw2();
    TH1F* Wantitop_mass_log_cumulative_BB = new TH1F("Wantitop_mass_log_cumulative_BB", "Wantitop_mass_log_cumulative_BB;Wantitop_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    Wantitop_mass_log_cumulative_BB->Sumw2();
    TH1F* Wantitop_mass_log_cumulative_BE = new TH1F("Wantitop_mass_log_cumulative_BE", "Wantitop_mass_log_cumulative_BE;Wantitop_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    Wantitop_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////dy////////////////////////////////////////////
    TH1F* DY_mass_linear = new TH1F("DY_mass_linear", "DY_mass_linear;DY_mass_linear;#Events", mass_bins, 70, 4000);
    DY_mass_linear->Sumw2();
    TH1F* DY_mass_linear_BB = new TH1F("DY_mass_linear_BB", "DY_mass_linear_BB;DY_mass_linear_BB;#Events", mass_bins, 70, 4000);
    DY_mass_linear_BB->Sumw2();
    TH1F* DY_mass_linear_BE = new TH1F("DY_mass_linear_BE", "DY_mass_linear_BE;DY_mass_linear_BE;#Events", mass_bins, 70, 4000);
    DY_mass_linear_BE->Sumw2();
    
    TH1F* DY_mass_log = new TH1F("DY_mass_log", "DY_mass_log;DY_mass_log;#Events", NMBINS, logMbins);
    DY_mass_log->Sumw2();
    TH1F* DY_mass_log_BB = new TH1F("DY_mass_log_BB", "DY_mass_log_BB;DY_mass_log_BB;#Events", NMBINS, logMbins);
    DY_mass_log_BB->Sumw2();
    TH1F* DY_mass_log_BE = new TH1F("DY_mass_log_BE", "DY_mass_log_BE;DY_mass_log_BE;#Events", NMBINS, logMbins);
    DY_mass_log_BE->Sumw2();
    
    TH1F* DY_mass_cumulative = new TH1F("DY_mass_cumulative", "DY_mass_cumulative;DY_mass_cumulative;#Events", mass_bins, 70, 4000);
    DY_mass_cumulative->Sumw2();
    TH1F* DY_mass_cumulative_BB = new TH1F("DY_mass_cumulative_BB", "DY_mass_cumulative_BB;DY_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    DY_mass_cumulative_BB->Sumw2();
    TH1F* DY_mass_cumulative_BE = new TH1F("DY_mass_cumulative_BE", "DY_mass_cumulative_BE;DY_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    DY_mass_cumulative_BE->Sumw2();
    
    TH1F* DY_mass_log_cumulative = new TH1F("DY_mass_log_cumulative", "DY_mass_log_cumulative;DY_mass_log_cumulative;#Events", NMBINS, logMbins);
    DY_mass_log_cumulative->Sumw2();
    TH1F* DY_mass_log_cumulative_BB = new TH1F("DY_mass_log_cumulative_BB", "DY_mass_log_cumulative_BB;DY_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    DY_mass_log_cumulative_BB->Sumw2();
    TH1F* DY_mass_log_cumulative_BE = new TH1F("DY_mass_log_cumulative_BE", "DY_mass_log_cumulative_BE;DY_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    DY_mass_log_cumulative_BE->Sumw2();
    
    
    
    /////////////////////////////////Data////////////////////////////////////////////
    TH1F* DATA_mass_linear = new TH1F("DATA_mass_linear", "DATA_mass_linear;DATA_mass_linear;#Events", mass_bins, 70, 4000);
    DATA_mass_linear->Sumw2();
    TH1F* DATA_mass_linear_BB = new TH1F("DATA_mass_linear_BB", "DATA_mass_linear_BB;DATA_mass_linear_BB;#Events", mass_bins, 70, 4000);
    DATA_mass_linear_BB->Sumw2();
    TH1F* DATA_mass_linear_BE = new TH1F("DATA_mass_linear_BE", "DATA_mass_linear_BE;DATA_mass_linear_BE;#Events", mass_bins, 70, 4000);
    DATA_mass_linear_BE->Sumw2();
    
    TH1F* DATA_mass_log = new TH1F("DATA_mass_log", "DATA_mass_log;DATA_mass_log;#Events", NMBINS, logMbins);
    DATA_mass_log->Sumw2();
    TH1F* DATA_mass_log_BB = new TH1F("DATA_mass_log_BB", "DATA_mass_log_BB;DATA_mass_log_BB;#Events", NMBINS, logMbins);
    DATA_mass_log_BB->Sumw2();
    TH1F* DATA_mass_log_BE = new TH1F("DATA_mass_log_BE", "DATA_mass_log_BE;DATA_mass_log_BE;#Events", NMBINS, logMbins);
    DATA_mass_log_BE->Sumw2();
    
    TH1F* DATA_mass_cumulative = new TH1F("DATA_mass_cumulative", "DATA_mass_cumulative;DATA_mass_cumulative;#Events", mass_bins, 70, 4000);
    DATA_mass_cumulative->Sumw2();
    TH1F* DATA_mass_cumulative_BB = new TH1F("DATA_mass_cumulative_BB", "DATA_mass_cumulative_BB;DATA_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    DATA_mass_cumulative_BB->Sumw2();
    TH1F* DATA_mass_cumulative_BE = new TH1F("DATA_mass_cumulative_BE", "DATA_mass_cumulative_BE;DATA_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    DATA_mass_cumulative_BE->Sumw2();
    
    TH1F* DATA_mass_log_cumulative = new TH1F("DATA_mass_log_cumulative", "DATA_mass_log_cumulative;DATA_mass_log_cumulative;#Events", NMBINS, logMbins);
    DATA_mass_log_cumulative->Sumw2();
    TH1F* DATA_mass_log_cumulative_BB = new TH1F("DATA_mass_log_cumulative_BB", "DATA_mass_log_cumulative_BB;DATA_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    DATA_mass_log_cumulative_BB->Sumw2();
    TH1F* DATA_mass_log_cumulative_BE = new TH1F("DATA_mass_log_cumulative_BE", "DATA_mass_log_cumulative_BE;DATA_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    DATA_mass_log_cumulative_BE->Sumw2();
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    //strat of looping over MC samples
    for(int j=0; j<mc; j++){
        
        
        std::cout<<"opening.. "<<MC_samples[j]<<" --- "<<" No of evnets: "<<events[j]<<" Cross section: "<<sigma[j]<<std::endl;
        
        TChain *treeMC = new TChain("SimpleNtupler/t");
        
        treeMC->Add("/eos/user/k/kaliyana/2017_MC/FileBased/"+ MC_samples[j] +"/"+ MC_samples[j] +".root");
        
        treeMC->SetBranchAddress("genWeight",&genWeight);
        treeMC->SetBranchAddress("event",&event);
        treeMC->SetBranchAddress("run",&run);
        treeMC->SetBranchAddress("lumi",&lumi);
        treeMC->SetBranchAddress("dil_mass",&dil_mass);
        treeMC->SetBranchAddress("dil_pt",&dil_pt);
        treeMC->SetBranchAddress("dil_eta",&dil_eta);
        treeMC->SetBranchAddress("dil_rap",&dil_rap);
        treeMC->SetBranchAddress("dil_phi",&dil_phi);
        treeMC->SetBranchAddress("cos_angle",&cos_angle);
        treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
        treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
        treeMC->SetBranchAddress("lep_pt",lep_pt);
        treeMC->SetBranchAddress("lep_id",lep_id);
        treeMC->SetBranchAddress("lep_eta",lep_eta);
        treeMC->SetBranchAddress("lep_phi",lep_phi);
        treeMC->SetBranchAddress("lep_dB",lep_dB);
        treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
        treeMC->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
        treeMC->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
        treeMC->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
        treeMC->SetBranchAddress("lep_expectedNnumberOfMatchedStations",lep_expectedNnumberOfMatchedStations);
        treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
        treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
        treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
        treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
        treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
        treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
        treeMC->SetBranchAddress("lep_glb_pt",lep_glb_pt);
        treeMC->SetBranchAddress("lep_picky_pt",lep_picky_pt);
        treeMC->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
        treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
        treeMC->SetBranchAddress("vertex_m",&vertex_m);
        treeMC->SetBranchAddress("GoodVtx",&GoodVtx);
        treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
        treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
        treeMC->SetBranchAddress("gen_dil_mass", &gen_dil_mass);
        treeMC->SetBranchAddress("gen_lep_qOverPt", gen_lep_qOverPt);
        treeMC->SetBranchAddress("gen_lep_eta", gen_lep_eta);
        treeMC->SetBranchAddress("gen_lep_pt", gen_lep_pt);
        
        treeMC->SetBranchAddress("met_pt",&met_pt);
        treeMC->SetBranchAddress("met_phi",&met_phi);
        treeMC->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
        treeMC->SetBranchAddress("lep_tk_dz",lep_tk_dz);
        
        Long64_t ne = treeMC->GetEntries();
        
        
        std::cout<<"START: "<<ne<<" --For mujet background"<<std::endl;
        
        
        //start looping over entries
        for(int p=0; p<ne ;p++){
            
            if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            treeMC->GetEntry(p); //takes pth event
            
            if (!(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50)) continue; //dimuon evnts should pass the trigger
            
            // c_event=event;
            
            // if((j==15 && run==1 && lumi==20725 && event==51812270) || (j==15 && run==1 && lumi==17494 && event==43733785)) continue;
            
            
            
            
            
            //start mujet selection
            if(
               //both muons passing pre-selection
               ((lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1) &&
                (lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1) &&
                (fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2) &&
                (lep_glb_numberOfValidPixelHits[0] > 0.0 && lep_glb_numberOfValidPixelHits[1] > 0.0) &&
                (lep_glb_numberOfValidTrackerLayers[0] > 5.0 && lep_glb_numberOfValidTrackerLayers[1] > 5.0)
                //(fabs(lep_tk_dz[0]) < 1.0 && fabs(lep_tk_dz[1]) < 1.0)
                )
               &&
               //one muon passing and one muon failing high-pt ID selection
               //Muon0 failing high-pt ID selection
               ( ( !((lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
                     (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0)  &&
                     lep_glb_numberOfValidPixelHits[0]>0 &&
                     lep_glb_numberOfValidTrackerLayers[0]>5 &&
                     lep_sumPt[0]/lep_tk_pt[0]<0.10 &&
                     lep_pt_err[0]/lep_pt[0]<0.3 &&
                     lep_triggerMatchPt[0]>50)
                  &&
                  //Muon1 passing high-pt ID selection
                  ((lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
                   (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)  &&
                   lep_glb_numberOfValidPixelHits[1]>0 &&
                   lep_glb_numberOfValidTrackerLayers[1]>5 &&
                   lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
                   lep_pt_err[1]/lep_pt[1]<0.3 &&
                   lep_triggerMatchPt[1]>50) )
                ||
                //Muon0 passing high-pt ID selection
                ( ((lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
                   (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0)  &&
                   lep_glb_numberOfValidPixelHits[0]>0 &&
                   lep_glb_numberOfValidTrackerLayers[0]>5 &&
                   lep_sumPt[0]/lep_tk_pt[0]<0.10 &&
                   lep_pt_err[0]/lep_pt[0]<0.3 &&
                   lep_triggerMatchPt[0]>50)
                 &&
                 //Muon1 failing high-pt ID selection
                 !((lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
                   (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)  &&
                   lep_glb_numberOfValidPixelHits[1]>0 &&
                   lep_glb_numberOfValidTrackerLayers[1]>5 &&
                   lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
                   lep_pt_err[1]/lep_pt[1]<0.3 &&
                   lep_triggerMatchPt[1]>50) )    )
               
               ){
                
                if (!GoodVtx) continue;
                if (cos_angle<=-0.9998) continue;
                if (lep_id[0]*lep_id[1]>0) continue; //dimuon evnts with opposite-sign muons
                if (vertex_chi2 >= 20) continue;
                if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
                prev_event = event;
                
                
                
                
                
                
              
                double gen_lead_pt = 0.;
                gM = gen_dil_mass;
                gen_lead_pt = (gen_lep_pt[1]*(gen_lep_pt[1]>gen_lep_pt[0]) + gen_lep_pt[0]*(gen_lep_pt[0]>gen_lep_pt[1]));
                
               /*
                
                /////////////////mt cut to prevent Wjets//////////////////
                if(lep_pt[0]>lep_pt[1]){
                    leading_pt = lep_pt[0];
                    leading_phi = lep_phi[0];
                }
                else{
                    leading_pt = lep_pt[1];
                    leading_phi = lep_phi[1];
                }
                delta_phi = std::abs(met_phi-leading_phi);
                if(delta_phi>PI) delta_phi=float(2*PI)-delta_phi;
                mt = sqrt(2*leading_pt*met_pt*(1-cos(delta_phi)));
                
                if(mt>=35.0) continue;
                
                 /////////////////mt cut to prevent Wjets//////////////////
               */
                
                weight[j] *= genWeight;
                weight_BB[j] *= genWeight;
                weight_BE[j] *= genWeight;
                
                
               
                /////////////To check which muon passes the high-pt ID selection//////////////////////
                bool mu0 = false;
                bool mu1 = false;
                if(((lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
                    (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0)  &&
                    lep_glb_numberOfValidPixelHits[0]>0 &&
                    lep_glb_numberOfValidTrackerLayers[0]>5 &&
                    lep_sumPt[0]/lep_tk_pt[0]<0.10 &&
                    lep_pt_err[0]/lep_pt[0]<0.3 &&
                    lep_triggerMatchPt[0]>50)){
                    ////Muon0 passes the high-pt ID selection
                    mu0 = true;
                }
                
                if(((lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
                    (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)  &&
                    lep_glb_numberOfValidPixelHits[1]>0 &&
                    lep_glb_numberOfValidTrackerLayers[1]>5 &&
                    lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
                    lep_pt_err[1]/lep_pt[1]<0.3 &&
                    lep_triggerMatchPt[1]>50)){
                    ////Muon0 passes the high-pt ID selection
                    mu1 = true;
                }
                //std::cout<<m0<<"    "<<m1<<std::endl;
                
                
                if(j==15){ //Wjets
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        Wjets_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        Wjets_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        Wjets_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        Wjets_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
            
                
                else if(j==16){ //Zjets
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        Zjets_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        Zjets_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        Zjets_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        Zjets_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
                
                
                else if(j==17){//ttbar
                    
                    
                    /*
                     kFactor = 0.991403 + 3.05593e-05 * gM - 2.21967e-07 * pow(gM,2) + 6.63658e-11 * pow(gM,3);
                     kFactor_BB = 0.990973 + 6.17124e-06 * gM - 3.31244e-07 * pow(gM,2) + 1.2125e-10 * pow(gM,3);
                     kFactor_BE = 0.990038 + 5.44269e-05 * gM - 2.43311e-07 * pow(gM,2) + 5.9748e-11 * pow(gM,3);
                     */
                    
                    // New NNPDF  from Jan for 2017
                    
                    double NNPDFFac = 1.;
                    double NNPDFFac_bb = 1.;
                    double NNPDFFac_be = 1.;
                    
                    NNPDFFac = ((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
                    NNPDFFac_bb = ((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
                    NNPDFFac_be =((gM<120)*1. + (gM>120 && gM<3000)*((0.994078695151) + (2.64819793287e-05)*pow(gM,1) + (-3.73996461024e-08)*pow(gM,2) + (-1.11452866827e-11)*pow(gM,3)) + (gM>3000)*(0.436005));
                    
                    kFactor = 1*NNPDFFac;
                    kFactor_BB = 1*NNPDFFac_bb;
                    kFactor_BE = 1*NNPDFFac_be;
                    
                    weight[j] *= kFactor;
                    weight_BB[j] *= kFactor_BB;
                    weight_BE[j] *= kFactor_BE;
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        ttbar_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        ttbar_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        ttbar_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        ttbar_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                    
                    
                    weight[j] /=  (float) kFactor;
                    weight_BB[j] /=  (float) kFactor_BB;
                    weight_BE[j] /=  (float) kFactor_BE;
                }
                
                
                
                else if(j==18){ //WW
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        WW_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        WW_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        WW_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        WW_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
                
                
                else if(j==19){ //WZ
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        WZ_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        WZ_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        WZ_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        WZ_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
                
                
                else if(j==20){ //ZZ
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        ZZ_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        ZZ_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        ZZ_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        ZZ_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
                
                
                else if(j==21){ //tW
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        tW_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        tW_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        tW_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        tW_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
                
                
                else if(j==22){ //Wantitop
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        Wantitop_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        Wantitop_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        Wantitop_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        Wantitop_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
                
                
                else if(j==23 || j==24 || j==25 || j==26 || j==27 || j==28 || j==29 || j==30 || j==31 || j==32){ //DY
                    
                    
                    /*
                     gM = gen_dil_mass - 400.0;
                     kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
                     kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
                     kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
                     */
                    
                    //New K-factor and NNDF from Jan for DY
                    
                    double NNPDFFac = 1.;
                    double NNPDFFac_bb = 1.;
                    double NNPDFFac_be = 1.;
                    
                    
                    ////////////NNPDF//////////////////
                    if (gM < 120){
                        
                        NNPDFFac =((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.8245728118) + (-0.0537728412909)*pow(gen_lead_pt,1) + (0.000731365981935)*pow(gen_lead_pt,2) + (7.16669312495e-06)*pow(gen_lead_pt,3) + (-1.99723894101e-07)*pow(gen_lead_pt,4) + (1.0112316789e-09)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.01849023288));
                        NNPDFFac_bb = ((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.91383074609) + (-0.0596201865777)*pow(gen_lead_pt,1) + (0.000811074027001)*pow(gen_lead_pt,2) + (7.90677720686e-06)*pow(gen_lead_pt,3) + (-2.21489848717e-07)*pow(gen_lead_pt,4) + (1.12700571973e-09)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.00484010198));
                        NNPDFFac_be = ((gen_lead_pt<30)*0.9 + (gen_lead_pt>30 && gen_lead_pt<100)*((1.71913319508) + (-0.0481243962238)*pow(gen_lead_pt,1) + (0.000666286154366)*pow(gen_lead_pt,2) + (6.45776405133e-06)*pow(gen_lead_pt,3) + (-1.82202504311e-07)*pow(gen_lead_pt,4) + (9.24567381899e-10)*pow(gen_lead_pt,5)) + (gen_lead_pt>100)*(1.02790393101));
                        
                    }
                    else{
                        NNPDFFac = ((0.918129) + (6.92702e-05)*pow(gM,1) + (1.62175e-08)*pow(gM,2) + (-2.47833e-11)*pow(gM,3) + (8.75707e-15)*pow(gM,4) + (-7.53019e-19)*pow(gM,5));
                        NNPDFFac_bb = ((0.914053) + (7.91618e-05)*pow(gM,1) + (2.19722e-08)*pow(gM,2) + (-3.49212e-11)*pow(gM,3) + (1.22504e-14)*pow(gM,4) + (-1.07347e-18)*pow(gM,5));
                        NNPDFFac_be = ((0.933214) + (3.76813e-05)*pow(gM,1) + (1.95612e-08)*pow(gM,2) + (-1.2688e-11)*pow(gM,3) + (3.69867e-15)*pow(gM,4) + (-2.62212e-19)*pow(gM,5));
                    }
                    ////////////NNPDF//////////////////
                    
                    
                    ////////////////K-factor///////////////
                    if(gM < 150){
                        kFactor = 1*NNPDFFac;
                        kFactor_BB = 1*NNPDFFac_bb;
                        kFactor_BE = 1*NNPDFFac_be;
                    }
                    
                    if(gM > 150){
                        kFactor = (1.053 - 0.0001552 * gM + 5.661e-08 * pow(gM,2) - 8.382e-12 * pow(gM,3))*NNPDFFac;
                        kFactor_BB = (1.032 - 0.000138 * gM + 4.827e-08 * pow(gM,2) - 7.321e-12 * pow(gM,3))*NNPDFFac_bb;
                        kFactor_BE = (1.064 - 0.0001674 * gM + 6.599e-08 * pow(gM,2) - 9.657e-12 * pow(gM,3))*NNPDFFac_be;
                    }
                    
                    weight[j] *= kFactor;
                    weight_BB[j] *= kFactor_BB;
                    weight_BE[j] *= kFactor_BE;
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        DY_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        DY_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        DY_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        DY_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                    weight[j] /=  (float) kFactor;
                    weight_BB[j] /=  (float) kFactor_BB;
                    weight_BE[j] /=  (float) kFactor_BE;
                }
                
                
                else{ //qcd
                    
                    
                    //Barrel-Barrel
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        W = FR/(1-FR);
                        qcd_mass_linear_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        qcd_mass_log_BB->Fill(vertex_m,  weight_BB[j]*Z_peak_BB*W);
                        
                    }
                    
                    //Barrel-Endcap or Endcap-Endcap
                    else {
                        
                        //Endcap-Endcap
                        if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        //Barrel-Endcap
                        else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                            if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                            else FR = (float) f_barrel->Eval(lep_eta[0]);
                            
                        }
                        
                        //Endcap-Barrel
                        else{
                            if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                            else FR = (float) f_endcap->Eval(lep_eta[0]);
                        }
                        
                        W = FR/(1-FR);
                        qcd_mass_linear_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        qcd_mass_log_BE->Fill(vertex_m,  weight_BE[j]*Z_peak_BE*W);
                        
                    }
                    
                }
                
                
                
                
                
                weight[j] /= genWeight;
                weight_BB[j] /= genWeight;
                weight_BE[j] /= genWeight;
                
            } //end mujet selection
            
            
            
            
        }//end looping over entries
        
        
    }//end looping over MC samples
    
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = DY_mass_linear_BB->Integral(i, mass_bins);
        DY_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DY_mass_linear_BE->Integral(i, mass_bins);
        DY_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = ttbar_mass_linear_BB->Integral(i, mass_bins);
        ttbar_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ttbar_mass_linear_BE->Integral(i, mass_bins);
        ttbar_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = tW_mass_linear_BB->Integral(i, mass_bins);
        tW_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = tW_mass_log_BE->Integral(i, mass_bins);
        tW_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = Wantitop_mass_linear_BB->Integral(i, mass_bins);
        Wantitop_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wantitop_mass_log_BE->Integral(i, mass_bins);
        Wantitop_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = WW_mass_linear_BB->Integral(i, mass_bins);
        WW_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WW_mass_linear_BE->Integral(i, mass_bins);
        WW_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = WZ_mass_linear_BB->Integral(i, mass_bins);
        WZ_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WZ_mass_linear_BE->Integral(i, mass_bins);
        WZ_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = ZZ_mass_linear_BB->Integral(i, mass_bins);
        ZZ_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ZZ_mass_linear_BE->Integral(i, mass_bins);
        ZZ_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = Wjets_mass_linear_BB->Integral(i, mass_bins);
        Wjets_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wjets_mass_linear_BE->Integral(i, mass_bins);
        Wjets_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = Zjets_mass_linear_BB->Integral(i, mass_bins);
        Zjets_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zjets_mass_linear_BE->Integral(i, mass_bins);
        Zjets_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = qcd_mass_linear_BB->Integral(i, mass_bins);
        qcd_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = qcd_mass_linear_BE->Integral(i, mass_bins);
        qcd_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = DY_mass_log_BB->Integral(i, NMBINS);
        DY_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DY_mass_log_BE->Integral(i, NMBINS);
        DY_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = ttbar_mass_log_BB->Integral(i, NMBINS);
        ttbar_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ttbar_mass_log_BE->Integral(i, NMBINS);
        ttbar_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = tW_mass_log_BB->Integral(i, NMBINS);
        tW_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = tW_mass_log_BE->Integral(i, NMBINS);
        tW_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = Wantitop_mass_log_BB->Integral(i, NMBINS);
        Wantitop_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wantitop_mass_log_BE->Integral(i, NMBINS);
        Wantitop_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = WW_mass_log_BB->Integral(i, NMBINS);
        WW_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WW_mass_log_BE->Integral(i, NMBINS);
        WW_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = WZ_mass_log_BB->Integral(i, NMBINS);
        WZ_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WZ_mass_log_BE->Integral(i, NMBINS);
        WZ_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = ZZ_mass_log_BB->Integral(i, NMBINS);
        ZZ_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ZZ_mass_log_BE->Integral(i, NMBINS);
        ZZ_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = Wjets_mass_log_BB->Integral(i, NMBINS);
        Wjets_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wjets_mass_log_BE->Integral(i, NMBINS);
        Wjets_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = Zjets_mass_log_BB->Integral(i, NMBINS);
        Zjets_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zjets_mass_log_BE->Integral(i, NMBINS);
        Zjets_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = qcd_mass_log_BB->Integral(i, NMBINS);
        qcd_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = qcd_mass_log_BE->Integral(i, NMBINS);
        qcd_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    
    
    
    
    
    
    
    
    //strat of looping over Data
    for(int j=0; j<data; j++) {
        
        
        std::cout<<"opening.. "<<DATA_samples[j]<<std::endl;
        
        
        
        TChain *treeDATA = new TChain("SimpleNtupler/t");
        
        treeDATA->Add("/eos/user/k/kaliyana/2017_data/"+ samp[j] +"/"+ DATA_samples[j] +".root");
        
        
        //treeDATA->SetBranchAddress("genWeight",&genWeight);
        treeDATA->SetBranchAddress("event",&event);
        treeDATA->SetBranchAddress("run",&run);
        treeDATA->SetBranchAddress("lumi",&lumi);
        treeDATA->SetBranchAddress("dil_mass",&dil_mass);
        treeDATA->SetBranchAddress("dil_pt",&dil_pt);
        treeDATA->SetBranchAddress("dil_eta",&dil_eta);
        treeDATA->SetBranchAddress("dil_rap",&dil_rap);
        treeDATA->SetBranchAddress("dil_phi",&dil_phi);
        treeDATA->SetBranchAddress("cos_angle",&cos_angle);
        treeDATA->SetBranchAddress("vertex_chi2",&vertex_chi2);
        treeDATA->SetBranchAddress("dil_chosen",&dil_chosen);
        treeDATA->SetBranchAddress("lep_pt",lep_pt);
        treeDATA->SetBranchAddress("lep_id",lep_id);
        treeDATA->SetBranchAddress("lep_eta",lep_eta);
        treeDATA->SetBranchAddress("lep_phi",lep_phi);
        treeDATA->SetBranchAddress("lep_dB",lep_dB);
        treeDATA->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
        treeDATA->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
        treeDATA->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
        treeDATA->SetBranchAddress("lep_expectedNnumberOfMatchedStations",lep_expectedNnumberOfMatchedStations);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
        treeDATA->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
        treeDATA->SetBranchAddress("lep_sumPt",lep_sumPt);
        treeDATA->SetBranchAddress("lep_tk_pt",lep_tk_pt);
        treeDATA->SetBranchAddress("lep_glb_pt",lep_glb_pt);
        treeDATA->SetBranchAddress("lep_picky_pt",lep_picky_pt);
        treeDATA->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
        treeDATA->SetBranchAddress("lep_pt_err",lep_pt_err);
        treeDATA->SetBranchAddress("vertex_m",&vertex_m);
        treeDATA->SetBranchAddress("GoodVtx",&GoodVtx);
        treeDATA->SetBranchAddress("lep_stationMask",&lep_stationMask);
        treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
        //treeDATA->SetBranchAddress("gen_dil_mass", &gen_dil_mass);
        //treeDATA->SetBranchAddress("gen_lep_qOverPt", gen_lep_qOverPt);
        //treeDATA->SetBranchAddress("gen_lep_eta", gen_lep_eta);
        //treeDATA->SetBranchAddress("gen_lep_pt", gen_lep_pt);
        
        treeDATA->SetBranchAddress("met_pt",&met_pt);
        treeDATA->SetBranchAddress("met_phi",&met_phi);
        treeDATA->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
        treeDATA->SetBranchAddress("lep_tk_dz",lep_tk_dz);
        
        Long64_t ne = treeDATA->GetEntries();
        
        std::cout<<"START: "<<"<<ne<<"<<" -For mujet background"<<std::endl;
        
        
        
        //looping over entries
        for ( int p=0; p<ne ;p++) {
            if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            
            
            treeDATA->GetEntry(p); //takes pth event
            
            
            
            if (!(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50)) continue; //dimuon evnts should pass the trigger
            
            
            //  c_event=event;
            
            //start mujet selection
            if(
               //both muons passing pre-selection
               ((lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1) &&
                (lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1) &&
                (fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2) &&
                (lep_glb_numberOfValidPixelHits[0] > 0.0 && lep_glb_numberOfValidPixelHits[1] > 0.0) &&
                (lep_glb_numberOfValidTrackerLayers[0] > 5.0 && lep_glb_numberOfValidTrackerLayers[1] > 5.0)
                //(fabs(lep_tk_dz[0]) < 1.0 && fabs(lep_tk_dz[1]) < 1.0)
                )
               &&
               //one muon passing and one muon failing high-pt ID selection
               //Muon0 failing high-pt ID selection
               ( ( !((lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
                     (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0)  &&
                     lep_glb_numberOfValidPixelHits[0]>0 &&
                     lep_glb_numberOfValidTrackerLayers[0]>5 &&
                     lep_sumPt[0]/lep_tk_pt[0]<0.10 &&
                     lep_pt_err[0]/lep_pt[0]<0.3 &&
                     lep_triggerMatchPt[0]>50)
                  &&
                  //Muon1 passing high-pt ID selection
                  ((lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
                   (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)  &&
                   lep_glb_numberOfValidPixelHits[1]>0 &&
                   lep_glb_numberOfValidTrackerLayers[1]>5 &&
                   lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
                   lep_pt_err[1]/lep_pt[1]<0.3 &&
                   lep_triggerMatchPt[1]>50) )
                ||
                //Muon0 passing high-pt ID selection
                ( ((lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
                   (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0)  &&
                   lep_glb_numberOfValidPixelHits[0]>0 &&
                   lep_glb_numberOfValidTrackerLayers[0]>5 &&
                   lep_sumPt[0]/lep_tk_pt[0]<0.10 &&
                   lep_pt_err[0]/lep_pt[0]<0.3 &&
                   lep_triggerMatchPt[0]>50)
                 &&
                 //Muon1 failing high-pt ID selection
                 !((lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
                   (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)  &&
                   lep_glb_numberOfValidPixelHits[1]>0 &&
                   lep_glb_numberOfValidTrackerLayers[1]>5 &&
                   lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
                   lep_pt_err[1]/lep_pt[1]<0.3 &&
                   lep_triggerMatchPt[1]>50) )    )
               
               ){
                
                if (!GoodVtx) continue;
                if (cos_angle<=-0.9998) continue;
                if (lep_id[0]*lep_id[1]>0) continue; //dimuon evnts with opposite-sign muons
                if (vertex_chi2 >= 20) continue;
                if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
                prev_event = event;
                
               /*
                
                /////////////////mt cut to prevent Wjets//////////////////
                if(lep_pt[0]>lep_pt[1]){
                    leading_pt = lep_pt[0];
                    leading_phi = lep_phi[0];
                }
                else{
                    leading_pt = lep_pt[1];
                    leading_phi = lep_phi[1];
                }
                delta_phi = std::abs(met_phi-leading_phi);
                if(delta_phi>PI) delta_phi=float(2*PI)-delta_phi;
                mt = sqrt(2*leading_pt*met_pt*(1-cos(delta_phi)));
                
                if(mt>=35.0) continue;
                
                /////////////////mt cut to prevent Wjets//////////////////
                
               */
                
                
                /////////////To check which muon passes the high-pt ID selection//////////////////////
                bool mu0 = false;
                bool mu1 = false;
                if(((lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
                    (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0)  &&
                    lep_glb_numberOfValidPixelHits[0]>0 &&
                    lep_glb_numberOfValidTrackerLayers[0]>5 &&
                    lep_sumPt[0]/lep_tk_pt[0]<0.10 &&
                    lep_pt_err[0]/lep_pt[0]<0.3 &&
                    lep_triggerMatchPt[0]>50)){
                    ////Muon0 passes the high-pt ID selection
                    mu0 = true;
                }
                
                if(((lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
                    (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)  &&
                    lep_glb_numberOfValidPixelHits[1]>0 &&
                    lep_glb_numberOfValidTrackerLayers[1]>5 &&
                    lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
                    lep_pt_err[1]/lep_pt[1]<0.3 &&
                    lep_triggerMatchPt[1]>50)){
                    ////Muon0 passes the high-pt ID selection
                    mu1 = true;
                }
                //std::cout<<m0<<"    "<<m1<<std::endl;
                
                //Barrel-Barrel
                if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                    
                    if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                    else FR = (float) f_barrel->Eval(lep_eta[0]);
                    W = FR/(1-FR);
                    DATA_mass_linear_BB->Fill(vertex_m,  W);
                    DATA_mass_log_BB->Fill(vertex_m,  W);
                    
                }
                
                //Barrel-Endcap or Endcap-Endcap
                else {
                    
                    //Endcap-Endcap
                    if(fabs(lep_eta[0]) >= 1.2 && fabs(lep_eta[1]) >= 1.2){
                        if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                        else FR = (float) f_endcap->Eval(lep_eta[0]);
                    }
                    
                    //Barrel-Endcap
                    else if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) >= 1.2){
                        if(mu0) FR = (float) f_endcap->Eval(lep_eta[1]);
                        else FR = (float) f_barrel->Eval(lep_eta[0]);
                        
                    }
                    
                    //Endcap-Barrel
                    else{
                        if(mu0) FR = (float) f_barrel->Eval(lep_eta[1]);
                        else FR = (float) f_endcap->Eval(lep_eta[0]);
                    }
                    
                    W = FR/(1-FR);
                    DATA_mass_linear_BE->Fill(vertex_m,  W);
                    DATA_mass_log_BE->Fill(vertex_m,  W);
                    
                }
                
            } //end mujet selection
            
            
            
        }//end looping over entries
        
        
    }//end looping over Data samples
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = DATA_mass_linear_BB->Integral(i, mass_bins);
        DATA_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DATA_mass_linear_BE->Integral(i, mass_bins);
        DATA_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = DATA_mass_log_BB->Integral(i, NMBINS);
        DATA_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DATA_mass_log_BE->Integral(i, NMBINS);
        DATA_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    /*
    
     ///////////////////////////////////////////Data-driven mujet estimate////////////////////////////////////////////////////////
     
     double val=0.0;
     
     ///////////////////////////////////////////Barrel-Barrel/////////////////////////////////////
     
     TH1F* dijet_linear_BB = (TH1F*) DATA_mass_linear_BB->Clone();
     dijet_linear_BB->Sumw2();
     TH1F* MC_linear_BB = (TH1F*) Wjets_mass_linear_BB->Clone();
     //TH1F* MC_linear_BB = (TH1F*) Zjets_mass_linear_BB->Clone();
     MC_linear_BB->Sumw2();
     MC_linear_BB->Add(DY_mass_linear_BB);
     MC_linear_BB->Add(Zjets_mass_linear_BB);
     MC_linear_BB->Add(ttbar_mass_linear_BB);
     MC_linear_BB->Add(WW_mass_linear_BB);
     MC_linear_BB->Add(WZ_mass_linear_BB);
     MC_linear_BB->Add(ZZ_mass_linear_BB);
     MC_linear_BB->Add(Wantitop_mass_linear_BB);
     MC_linear_BB->Add(tW_mass_linear_BB);
     dijet_linear_BB->Add(MC_linear_BB,-1);
     
     for(int i=1; i<=mass_bins; i++){
         val=dijet_linear_BB->GetBinContent(i);
         if(val<0) dijet_linear_BB->SetBinContent(i,val*(-1));
         else dijet_linear_BB->SetBinContent(i,val);
     }
     
     dijet_linear_BB->SaveAs("dijet_linear_BB.root","root");
     
     TH1F* dijet_log_BB = (TH1F*) DATA_mass_log_BB->Clone();
     dijet_log_BB->Sumw2();
     TH1F* MC_log_BB = (TH1F*) Wjets_mass_log_BB->Clone();
     //TH1F* MC_log_BB = (TH1F*) Zjets_mass_log_BB->Clone();
     MC_log_BB->Sumw2();
     MC_log_BB->Add(DY_mass_log_BB);
     MC_log_BB->Add(Zjets_mass_log_BB);
     MC_log_BB->Add(ttbar_mass_log_BB);
     MC_log_BB->Add(WW_mass_log_BB);
     MC_log_BB->Add(WZ_mass_log_BB);
     MC_log_BB->Add(ZZ_mass_log_BB);
     MC_log_BB->Add(Wantitop_mass_log_BB);
     MC_log_BB->Add(tW_mass_log_BB);
     dijet_log_BB->Add(MC_log_BB,-1);
     
     for(int i=1; i<=NMBINS; i++){
         val=dijet_log_BB->GetBinContent(i);
         if(val<0) dijet_log_BB->SetBinContent(i,val*(-1));
         else dijet_log_BB->SetBinContent(i,val);
     }
     
     dijet_log_BB->SaveAs("dijet_log_BB.root","root");
     
     TH1F* dijet_cumulative_BB = (TH1F*) DATA_mass_cumulative_BB->Clone();
     dijet_cumulative_BB->Sumw2();
     TH1F* MC_cumulative_BB = (TH1F*) Wjets_mass_cumulative_BB->Clone();
     //TH1F* MC_cumulative_BB = (TH1F*) Zjets_mass_cumulative_BB->Clone();
     MC_cumulative_BB->Sumw2();
     MC_cumulative_BB->Add(DY_mass_cumulative_BB);
     MC_cumulative_BB->Add(Zjets_mass_cumulative_BB);
     MC_cumulative_BB->Add(ttbar_mass_cumulative_BB);
     MC_cumulative_BB->Add(WW_mass_cumulative_BB);
     MC_cumulative_BB->Add(WZ_mass_cumulative_BB);
     MC_cumulative_BB->Add(ZZ_mass_cumulative_BB);
     MC_cumulative_BB->Add(Wantitop_mass_cumulative_BB);
     MC_cumulative_BB->Add(tW_mass_cumulative_BB);
     dijet_cumulative_BB->Add(MC_cumulative_BB,-1);
     
     for(int i=1; i<=mass_bins; i++){
     val=dijet_cumulative_BB->GetBinContent(i);
     if(val<0) dijet_cumulative_BB->SetBinContent(i,val*(-1));
         else dijet_cumulative_BB->SetBinContent(i,val);
     }
     
     dijet_cumulative_BB->SaveAs("dijet_cumulative_BB.root","root");
     
     TH1F* dijet_log_cumulative_BB = (TH1F*) DATA_mass_log_cumulative_BB->Clone();
     dijet_log_cumulative_BB->Sumw2();
     TH1F* MC_log_cumulative_BB = (TH1F*) Wjets_mass_log_cumulative_BB->Clone();
     //TH1F* MC_log_cumulative_BB = (TH1F*) Zjets_mass_log_cumulative_BB->Clone();
     MC_log_cumulative_BB->Sumw2();
     MC_log_cumulative_BB->Add(DY_mass_log_cumulative_BB);
     MC_log_cumulative_BB->Add(Zjets_mass_log_cumulative_BB);
     MC_log_cumulative_BB->Add(ttbar_mass_log_cumulative_BB);
     MC_log_cumulative_BB->Add(WW_mass_log_cumulative_BB);
     MC_log_cumulative_BB->Add(WZ_mass_log_cumulative_BB);
     MC_log_cumulative_BB->Add(ZZ_mass_log_cumulative_BB);
     MC_log_cumulative_BB->Add(Wantitop_mass_log_cumulative_BB);
     MC_log_cumulative_BB->Add(tW_mass_log_cumulative_BB);
     dijet_log_cumulative_BB->Add(MC_log_cumulative_BB,-1);
     
     for(int i=1; i<=NMBINS; i++){
         val=dijet_log_cumulative_BB->GetBinContent(i);
         if(val<0) dijet_log_cumulative_BB->SetBinContent(i,val*(-1));
         else dijet_log_cumulative_BB->SetBinContent(i,val);
     }
     
     dijet_log_cumulative_BB->SaveAs("dijet_log_cumulative_BB.root","root");
     
     ///////////////////////////////////////////Barrel-Barrel/////////////////////////////////////
     
     
     
     ///////////////////////////////////////////Barrel-Endcap or Endcap-Endcap/////////////////////////////////////
     
     TH1F* dijet_linear_BE = (TH1F*) DATA_mass_linear_BE->Clone();
     dijet_linear_BE->Sumw2();
     TH1F* MC_linear_BE = (TH1F*) Wjets_mass_linear_BE->Clone();
     //TH1F* MC_linear_BE = (TH1F*) Zjets_mass_linear_BE->Clone();
     MC_linear_BE->Sumw2();
     MC_linear_BE->Add(DY_mass_linear_BE);
     MC_linear_BE->Add(Zjets_mass_linear_BE);
     MC_linear_BE->Add(ttbar_mass_linear_BE);
     MC_linear_BE->Add(WW_mass_linear_BE);
     MC_linear_BE->Add(WZ_mass_linear_BE);
     MC_linear_BE->Add(ZZ_mass_linear_BE);
     MC_linear_BE->Add(Wantitop_mass_linear_BE);
     MC_linear_BE->Add(tW_mass_linear_BE);
     dijet_linear_BE->Add(MC_linear_BE,-1);
     
     for(int i=1; i<=mass_bins; i++){
         val=dijet_linear_BE->GetBinContent(i);
         if(val<0) dijet_linear_BE->SetBinContent(i,val*(-1));
         else dijet_linear_BE->SetBinContent(i,val);
     }
     
     dijet_linear_BE->SaveAs("dijet_linear_BE.root","root");
     
     TH1F* dijet_log_BE = (TH1F*) DATA_mass_log_BE->Clone();
     dijet_log_BE->Sumw2();
     TH1F* MC_log_BE = (TH1F*) Wjets_mass_log_BE->Clone();
     //TH1F* MC_log_BE = (TH1F*) Zjets_mass_log_BE->Clone();
     MC_log_BE->Sumw2();
     MC_log_BE->Add(DY_mass_log_BE);
     MC_log_BE->Add(Zjets_mass_log_BE);
     MC_log_BE->Add(ttbar_mass_log_BE);
     MC_log_BE->Add(WW_mass_log_BE);
     MC_log_BE->Add(WZ_mass_log_BE);
     MC_log_BE->Add(ZZ_mass_log_BE);
     MC_log_BE->Add(Wantitop_mass_log_BE);
     MC_log_BE->Add(tW_mass_log_BE);
     dijet_log_BE->Add(MC_log_BE,-1);
     
     for(int i=1; i<=NMBINS; i++){
         val=dijet_log_BE->GetBinContent(i);
         if(val<0) dijet_log_BE->SetBinContent(i,val*(-1));
         else dijet_log_BE->SetBinContent(i,val);
     }
     
     dijet_log_BE->SaveAs("dijet_log_BE.root","root");
     
     TH1F* dijet_cumulative_BE = (TH1F*) DATA_mass_cumulative_BE->Clone();
     dijet_cumulative_BE->Sumw2();
     TH1F* MC_cumulative_BE = (TH1F*) Wjets_mass_cumulative_BE->Clone();
     //TH1F* MC_cumulative_BE = (TH1F*) Zjets_mass_cumulative_BE->Clone();
     MC_cumulative_BE->Sumw2();
     MC_cumulative_BE->Add(DY_mass_cumulative_BE);
     MC_cumulative_BE->Add(Zjets_mass_cumulative_BE);
     MC_cumulative_BE->Add(ttbar_mass_cumulative_BE);
     MC_cumulative_BE->Add(WW_mass_cumulative_BE);
     MC_cumulative_BE->Add(WZ_mass_cumulative_BE);
     MC_cumulative_BE->Add(ZZ_mass_cumulative_BE);
     MC_cumulative_BE->Add(Wantitop_mass_cumulative_BE);
     MC_cumulative_BE->Add(tW_mass_cumulative_BE);
     dijet_cumulative_BE->Add(MC_cumulative_BE,-1);
     
     for(int i=1; i<=mass_bins; i++){
         val=dijet_cumulative_BE->GetBinContent(i);
         if(val<0) dijet_cumulative_BE->SetBinContent(i,val*(-1));
         else dijet_cumulative_BE->SetBinContent(i,val);
     }
     
     dijet_cumulative_BE->SaveAs("dijet_cumulative_BE.root","root");
     
     TH1F* dijet_log_cumulative_BE = (TH1F*) DATA_mass_log_cumulative_BE->Clone();
     dijet_log_cumulative_BE->Sumw2();
     TH1F* MC_log_cumulative_BE = (TH1F*) Wjets_mass_log_cumulative_BE->Clone();
     //TH1F* MC_log_cumulative_BE = (TH1F*) Zjets_mass_log_cumulative_BE->Clone();
     MC_log_cumulative_BE->Sumw2();
     MC_log_cumulative_BE->Add(DY_mass_log_cumulative_BE);
     MC_log_cumulative_BE->Add(Zjets_mass_log_cumulative_BE);
     MC_log_cumulative_BE->Add(ttbar_mass_log_cumulative_BE);
     MC_log_cumulative_BE->Add(WW_mass_log_cumulative_BE);
     MC_log_cumulative_BE->Add(WZ_mass_log_cumulative_BE);
     MC_log_cumulative_BE->Add(ZZ_mass_log_cumulative_BE);
     MC_log_cumulative_BE->Add(Wantitop_mass_log_cumulative_BE);
     MC_log_cumulative_BE->Add(tW_mass_log_cumulative_BE);
     dijet_log_cumulative_BE->Add(MC_log_cumulative_BE,-1);
     
     for(int i=1; i<=NMBINS; i++){
         val=dijet_log_cumulative_BE->GetBinContent(i);
         if(val<0) dijet_log_cumulative_BE->SetBinContent(i,val*(-1));
         else dijet_log_cumulative_BE->SetBinContent(i,val);
     }
     
     dijet_log_cumulative_BE->SaveAs("dijet_log_cumulative_BE.root","root");
     
     ///////////////////////////////////////////Barrel-Endcap or Endcap-Endcap/////////////////////////////////////
     
     ///////////////////////////////////////Inclusive/////////////////////////
     
     TH1F* dijet_linear = (TH1F*) dijet_linear_BB->Clone();
     dijet_linear->Sumw2();
     dijet_linear->Add(dijet_linear_BE);
     dijet_linear->SaveAs("dijet_linear.root","root");
     
     TH1F* dijet_log = (TH1F*) dijet_log_BB->Clone();
     dijet_log->Sumw2();
     dijet_log->Add(dijet_log_BE);
     dijet_log->SaveAs("dijet_log.root","root");
     
     TH1F* dijet_cumulative = (TH1F*) dijet_cumulative_BB->Clone();
     dijet_cumulative->Sumw2();
     dijet_cumulative->Add(dijet_cumulative_BE);
     dijet_cumulative->SaveAs("dijet_cumulative.root","root");
     
     TH1F* dijet_log_cumulative = (TH1F*) dijet_log_cumulative_BB->Clone();
     dijet_log_cumulative->Sumw2();
     dijet_log_cumulative->Add(dijet_log_cumulative_BE);
     dijet_log_cumulative->SaveAs("dijet_log_cumulative.root","root");
     
    */
    
    
    f->Write();
    
    f->Close();
    
}
