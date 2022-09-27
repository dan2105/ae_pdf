#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TChain.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include <math.h>
#include <algorithm>
#include <list>
#include <stdbool.h>
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <cstdio>

/*
float MinDeltaPhi(float phi1, float phi2) {
dphi = phi1 - phi2;
while (dphi < -M_PI) dphi += 2*M_PI;
while (dphi >  M_PI) dphi -= 2*M_PI;
return dphi;
}
*/
void double_atau_analysis_v2(const double proc = 1., string collider ="", string pdf=""){

  TChain * chain = new TChain("T");

  if(proc == 1 && pdf=="luxlep")
  {
    cout << "Ok! Analyzing p p (a tau > a tau)"<<endl;
    chain->Add("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/pp_luxlep_13TeV_hepmc.root/T");
  }
  if(proc == 1 && pdf =="mrst")
  {
    cout << "Ok! Analyzing p p (a tau > a tau)"<<endl;
    chain->Add("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/pp_mrst_13TeV_hepmc.root/T");
  }
  if(proc == 1 && pdf =="nn23nlo" )
  {
    cout << "Ok! Analyzing p p (a tau > a tau)"<<endl;
    chain->Add("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/pp_nn23nlo_13TeV_hepmc.root/T");
  }
  if(proc == 2 && pdf=="luxlep")
  {
    cout << "Ok! Analyzing Pb p (a tau > a tau)"<<endl;
    chain->Add("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/Pbp_luxlep_tau_8TeV_hepmc.root/T");
  }
  if(proc == 2 && pdf =="mrst")
  {
    cout << "Ok! Analyzing Pb p (a tau > a tau)"<<endl;
    chain->Add("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/Pbp_mrst_tau_8TeV_hepmc.root/T");
  }
  if(proc == 2 && pdf =="nn23nlo" )
  {
    cout << "Ok! Analyzing Pb p (a tau > a tau)"<<endl;
    chain->Add("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/Pbp_nn23nlo_tau_8TeV_hepmc.root/T");
  }


  int const NTAU1MAX=100;
  int n_tau1=0;
  double tau1_pt[NTAU1MAX];
  double tau1_px[NTAU1MAX];
  double tau1_py[NTAU1MAX];
  double tau1_pz[NTAU1MAX];
  double tau1_energy[NTAU1MAX];
  double tau1_eta[NTAU1MAX];
  double tau1_phi[NTAU1MAX];
  double tau1_charge[NTAU1MAX];


  int const NAMAX=100;
  int n_a=0;
  double a_pt[NAMAX];
  double a_px[NAMAX];
  double a_py[NAMAX];
  double a_pz[NAMAX];
  double a_energy[NAMAX];
  double a_eta[NAMAX];
  double a_phi[NAMAX];
  double a_charge[NAMAX];


  // Branch of Eons
  chain->SetBranchAddress("n_tau1",&n_tau1);
  chain->SetBranchAddress("tau1_pt",&tau1_pt);
  chain->SetBranchAddress("tau1_px",&tau1_px);
  chain->SetBranchAddress("tau1_py",&tau1_py);
  chain->SetBranchAddress("tau1_pz",&tau1_pz);
  chain->SetBranchAddress("tau1_energy",&tau1_energy);
  chain->SetBranchAddress("tau1_eta",&tau1_eta);
  chain->SetBranchAddress("tau1_phi",&tau1_phi);
  chain->SetBranchAddress("tau1_charge",&tau1_charge);


  chain->SetBranchAddress("n_a",&n_a);
  chain->SetBranchAddress("a_pt",&a_pt);
  chain->SetBranchAddress("a_px",&a_px);
  chain->SetBranchAddress("a_py",&a_py);
  chain->SetBranchAddress("a_pz",&a_pz);
  chain->SetBranchAddress("a_energy",&a_energy);
  chain->SetBranchAddress("a_eta",&a_eta);
  chain->SetBranchAddress("a_phi",&a_phi);
  chain->SetBranchAddress("a_charge",&a_charge);


  //===============================================================
  TH1F* h_ntau1            = new TH1F("ntau1","ntau1",10,0.,10.);
  TH1F* h_tau1_pt          = new TH1F("tau1_pt","tau1_pt",100,0.,50.);
  TH1F* h_tau1_px          = new TH1F("tau1_px","tau1_px",100,-10.,10.);
  TH1F* h_tau1_py          = new TH1F("tau1_py","tau1_py",1000,-10.,10.);
  TH1F* h_tau1_pz          = new TH1F("tau1_pz","tau1_pz",1000,-3000.,3000.);
  TH1F* h_tau1_eta         = new TH1F("tau1_eta","tau1_eta",1000,-10,10);
  TH1F* h_tau1_phi         = new TH1F("tau1_phi","tau1_phi",1000,-M_PI,M_PI);
  TH1F* h_tau1_energy      = new TH1F("tau1_energy","tau1_energy",100,0.,10.);
  TH1F* h_tau1_charge = new TH1F("tau1_charge","tau1_charge",20,-10.,10.);

  TH1F* h_tau1_eta_frontal         = new TH1F("tau1_eta_frontal","tau1_eta",1000,-10,10);
  TH1F* h_tau1_eta_central         = new TH1F("tau1_eta_central","tau1_eta",1000,-10,10);

  TH1F* h_pttaua_eta_central         = new TH1F("pttaua_eta_central","pttaua_eta_central",1000,0,100);
  TH1F* h_pttaua_eta_frontal         = new TH1F("pttaua_eta_frontal","pttaua_eta_frontal",1000,0,100);

  TH1F* h_atau_eta_central         = new TH1F("atau_eta_central","atau_eta_central",1000,-10,10);
  TH1F* h_atau_eta_frontal         = new TH1F("atau_eta_frontal","atau_eta_frontal",1000,-10,10);

  TH1F* h_atau_m_inv          = new TH1F("atau_minv","atau_minv",1000,0,500);
  TH1F* h_atau_m_inv_central          = new TH1F("atau_minv_central","atau_minv_central",1000,0,500);
  TH1F* h_atau_m_inv_frontal          = new TH1F("atau_minv_frontal","atau_minv_frontal",1000,0,500);

  TH1F* h_na            = new TH1F("na","na",10,0.,10.);
  TH1F* h_a_pt          = new TH1F("a_pt","a_pt",100,0.,50.);
  TH1F* h_a_px          = new TH1F("a_px","a_px",100,-10.,10.);
  TH1F* h_a_py          = new TH1F("a_py","a_py",1000,-10.,10.);
  TH1F* h_a_pz          = new TH1F("a_pz","a_pz",1000,-3000.,3000.);
  TH1F* h_a_eta         = new TH1F("a_eta","a_eta",1000,-10,10);
  TH1F* h_a_phi         = new TH1F("a_phi","a_phi",1000,-M_PI,M_PI);
  TH1F* h_a_energy      = new TH1F("a_energy","a_energy",100,0.,10.);
  TH1F* h_a_charge = new TH1F("a_charge","a_charge",20,-10.,10.);

  TH1F* h_aco_taua         = new TH1F("aco_taua","aco_taua",1000,0.,1.);
  TH1F* h_aco_taua_central = new TH1F("aco_taua_cen","aco_taua",1000,0.,1.);
  TH1F* h_aco_taua_frontal = new TH1F("aco_taua_for","aco_taua",1000,0.,1.);

  TH1F* h_tau1_pt_central          = new TH1F("tau1_pt_central","tau1_pt_central",100,0.,50.);
  TH1F* h_tau1_pt_frontal          = new TH1F("tau1_pt_frontal","tau1_pt_frontal",100,0.,50.);


  TH1F* h_a_eta_central_efixed  = new TH1F("a_eta_central_efixed","a_eta_central_efixed",1000,-10,10);
  TH1F* h_tau_eta_central_afixed  = new TH1F("e_eta_central_afixed","e_eta_central_afixed",1000,-10,10);
  TH1F* h_a_eta_frontal_efixed  = new TH1F("a_eta_frontal_efixed","a_eta_frontal_efixed",1000,-10,10);
  TH1F* h_tau_eta_frontal_afixed  = new TH1F("e_eta_frontal_afixed","e_eta_frontal_afixed",1000,-10,10);
  TH1F* h_nevents_cuts         = new TH1F("nevents_cuts","nevents_cuts",10,0.,10.);

  TH2D* h_pte_pta_2D_correlation_central = new TH2D("pte_pta_central","pte_pta_central", 1000,0.0,50.0,1000,0.0,50.0);
  TH2D* hme_ma_2D_correlation_central = new TH2D("me_ma_central","me_ma_central", 1000,0.0,50,1000,0.0,50);

  double deltaphi = 0., dphi =0;
  //====================================================================
  //tree loop
  //====================================================================
  Int_t nentries = (Int_t)chain->GetEntries();
  for (Int_t i=0; i<nentries; i++) {
    chain->GetEntry(i);
    //===================================================================
    //======================================================================
    //electron
    //======================================================================
    TLorentzVector tau1p4;
    int itau1_maxpt = -1;
    double tau1_pt_max=0.;
    for(int itau1 = 0; itau1 < n_tau1; itau1++) {
      h_tau1_charge->Fill(tau1_charge[itau1]);
      //      if(i%2==0) e_pz[ie]=-e_pz[ie];
      if( tau1_pt[itau1] > tau1_pt_max ){
        itau1_maxpt = itau1 ;
        tau1_pt_max = tau1_pt[itau1];
      }
    }

    double tau1pt = tau1_pt[0];
    double tau1px = tau1_px[0];
    double tau1py = tau1_py[0];
    double tau1pz = tau1_pz[0];
    double tau1eta = tau1_eta[0];
    double tau1phi = tau1_phi[0];
    double tau1energy = tau1_energy[0];

    //======================================================================
    //Photon
    //======================================================================
    TLorentzVector ap4;
    int ia_maxpt = -1;
    double a_pt_max=0.;
    for(int ia = 0; ia < n_a; ia++) {
      h_a_charge->Fill(a_charge[ia]);
      //      if(i%2==0) e_pz[ie]=-e_pz[ie];
      if( a_pt[ia] > a_pt_max ){
        ia_maxpt = ia ;
        a_pt_max = a_pt[ia];
      }
    }

    double apt = a_pt[0];
    double apx = a_px[0];
    double apy = a_py[0];
    double apz = a_pz[0];
    double aeta = a_eta[0];
    double aphi = a_phi[0];
    double aenergy = a_energy[0];


    //general eletron vector
    //======================================================================
    h_ntau1->Fill(n_tau1);
    h_na->Fill(n_a);

    tau1p4.SetPxPyPzE(tau1px,tau1py,tau1pz,tau1energy);
    ap4.SetPtEtaPhiE(apx,apy,apz,aenergy);


    h_tau1_pt->Fill(tau1p4.Pt());
    h_tau1_px->Fill(tau1p4.Px());
    h_tau1_py->Fill(tau1p4.Py());
    h_tau1_pz->Fill(tau1p4.Pz());
    h_tau1_eta->Fill(tau1p4.Eta());
    h_tau1_phi->Fill(tau1p4.Phi());
    h_tau1_energy->Fill(tau1p4.E());


    h_a_pt->Fill(ap4.Pt());
    h_a_px->Fill(ap4.Px());
    h_a_py->Fill(ap4.Py());
    h_a_pz->Fill(ap4.Pz());
    h_a_eta->Fill(ap4.Eta());
    h_a_phi->Fill(ap4.Phi());
    h_a_energy->Fill(ap4.E());


    //===================================================================
    h_nevents_cuts->Fill(0);

    deltaphi =  tau1p4.Phi() - ap4.Phi();
    while (deltaphi < -M_PI) deltaphi += 2*M_PI;
    while (deltaphi >  M_PI) deltaphi -= 2*M_PI;
    double aco_taua = 1 -fabs(deltaphi/TMath::Pi());

    h_aco_taua->Fill(aco_taua);

    bool eta_tau_fixed = false;
    bool eta_a_fixed = false;
    bool eta_tau_fixed_frontal = false;
    bool eta_a_fixed_frontal = false;


    if(fabs(tau1p4.Eta()) < 2.5){
      eta_tau_fixed=true;
    }
    if(fabs(ap4.Eta()) < 2.5){
      eta_a_fixed=true;
    }
    if(tau1p4.Eta() > 2.0 && tau1p4.Eta() < 4.5){
      eta_tau_fixed_frontal=true;
    }
    if(ap4.Eta() > 2.0 && ap4.Eta() < 4.5){
      eta_a_fixed_frontal=true;
    }

    //=====================================================================
    //CMS cuts
    //=====================================================================
    // Invariant mass of the diphoton pair (> 5 GeV);
    //- Rapidity of the diphoton pair  (|eta|< 2.5);
    //- Transverse momentum of the pair ( 0 < pT < 1) GeV;
    //- Acoplanarity (< 0.01);

    h_atau_m_inv->Fill((ap4+tau1p4).M());

    if(proc == 1 && collider =="central"){

      if(tau1p4.Pt() < 2 && ap4.Pt() < 2) continue;
      h_nevents_cuts->Fill(1);
      if(eta_tau_fixed){
        h_a_eta_central_efixed->Fill(ap4.Eta());
      }
      if(eta_a_fixed){
        h_tau_eta_central_afixed->Fill(tau1p4.Eta());
      }
      //======================================================================
      //electron and photon in the CMS acceptance
      //======================================================================

      if(fabs(tau1p4.Eta()) > 2.5) continue;
      if(fabs(ap4.Eta()) > 2.5) continue;
      h_nevents_cuts->Fill(2);
      h_tau1_eta_central->Fill(tau1p4.Eta());
      h_pttaua_eta_central->Fill((tau1p4+ap4).Pt());
      h_aco_taua_central->Fill(aco_taua);
      h_atau_eta_central->Fill((tau1p4+ap4).Eta());
      h_atau_m_inv_central->Fill((ap4+tau1p4).M());
      h_tau1_pt_central->Fill(tau1p4.Pt());

      h_pte_pta_2D_correlation_central->Fill(tau1p4.Pt(),ap4.Pt());
      hme_ma_2D_correlation_central->Fill((tau1p4+ap4).Pt(),aco_taua);

      if((ap4+tau1p4).M() < 150) continue;
      h_nevents_cuts->Fill(3);


    }

    if(proc == 1 && collider =="frontal"){
      if(tau1p4.Pt() < 2 && ap4.Pt() < 2) continue;
      h_nevents_cuts->Fill(1);
      //======================================================================
      //electron and photon in the LHCb acceptance
      //======================================================================
      if(eta_tau_fixed_frontal){
        h_a_eta_frontal_efixed->Fill(ap4.Eta());
      }
      if(eta_a_fixed){
        h_tau_eta_frontal_afixed->Fill(tau1p4.Eta());
      }

      if(tau1p4.Eta() < 2.0 || tau1p4.Eta() > 4.5) continue;
      if(ap4.Eta() < 2.0 || ap4.Eta() > 4.5)continue;
      h_nevents_cuts->Fill(2);
      h_aco_taua_frontal->Fill(aco_taua);
      h_tau1_eta_frontal->Fill(tau1p4.Eta());
      h_pttaua_eta_frontal->Fill((tau1p4+ap4).Pt());
      h_atau_eta_frontal->Fill((tau1p4+ap4).Eta());
      h_atau_m_inv_frontal->Fill((ap4+tau1p4).M());
      h_tau1_pt_frontal->Fill(tau1p4.Pt());
    }

    if(proc == 2 && collider =="central"){
      if(tau1p4.Pt() < 2 && ap4.Pt() < 2) continue;
      h_nevents_cuts->Fill(1);
      //======================================================================
      //electron and photon in the CMS acceptance
      //======================================================================
      if(eta_tau_fixed){
        h_a_eta_central_efixed->Fill(ap4.Eta());
      }
      if(eta_a_fixed){
        h_tau_eta_central_afixed->Fill(tau1p4.Eta());
      }
      if(fabs(tau1p4.Eta()) > 2.5) continue;
      if(fabs(ap4.Eta()) > 2.5) continue;

      h_nevents_cuts->Fill(2);
      h_tau1_eta_central->Fill(tau1p4.Eta());
      h_aco_taua_central->Fill(aco_taua);
      h_pttaua_eta_central->Fill((tau1p4+ap4).Pt());
      h_atau_eta_central->Fill((tau1p4 + ap4).Eta());
      h_atau_m_inv_central->Fill((ap4+tau1p4).M());
      h_tau1_pt_central->Fill(tau1p4.Pt());

      if((ap4+tau1p4).M() < 150) continue;
      h_nevents_cuts->Fill(3);

    }

    if(proc == 2 && collider =="frontal"){
      if(tau1p4.Pt() < 2 && ap4.Pt() < 2) continue;
      h_nevents_cuts->Fill(1);

      if(eta_tau_fixed_frontal){
        h_a_eta_frontal_efixed->Fill(ap4.Eta());
      }
      if(eta_a_fixed){
        h_tau_eta_frontal_afixed->Fill(tau1p4.Eta());
      }
      //======================================================================
      //electron and photon in the LHCb acceptance
      //======================================================================
      if(tau1p4.Eta() < 2.0 || tau1p4.Eta() > 4.5 ) continue;
      if(ap4.Eta() < 2.0 || ap4.Eta() > 4.5 ) continue;

      h_nevents_cuts->Fill(2);
      h_aco_taua_frontal->Fill(aco_taua);
      h_tau1_eta_frontal->Fill(tau1p4.Eta());
      h_pttaua_eta_frontal->Fill((tau1p4+ap4).Pt());
      h_atau_eta_frontal->Fill((tau1p4+ap4).Eta());
      h_atau_m_inv_frontal->Fill((ap4+tau1p4).M());
      h_tau1_pt_frontal->Fill(tau1p4.Pt());



    }




  } //Fim do loop na tree

  if(proc == 1 && collider =="central" && pdf=="luxlep"  ){
    cout << "Ok! Analyzing p p (a tau > a tau) inside the central collider range with LUXLEP pdf"<<endl;
  }
  if(proc == 1 && collider =="frontal" && pdf=="luxlep" ){
    cout << "Ok! Analyzing p p (a tau > a tau) inside the frontal collider range with LUXLEP pdf"<<endl;
  }
  if(proc == 1 && collider =="central" && pdf=="nn23nlo"  ){
    cout << "Ok! Analyzing p p (a tau > a tau) inside the central collider range with NN23NLO pdf"<<endl;
  }
  if(proc == 1 && collider =="frontal" && pdf=="nn23nlo" ){
    cout << "Ok! Analyzing p p (a tau > a tau) inside the frontal collider range with NN23NLO pdf"<<endl;
  }
  if(proc == 1 && collider =="central" && pdf=="mrst"  ){
    cout << "Ok! Analyzing p p (a tau > a tau) inside the central collider range with MRST pdf"<<endl;
  }
  if(proc == 1 && collider =="frontal" && pdf=="mrst" ){
    cout << "Ok! Analyzing p p (a tau > a tau) inside the frontal collider range with MRST pdf"<<endl;
  }
//=========================================================================================================
  if(proc == 2 && collider =="central" && pdf=="luxlep"  ){
    cout << "Ok! Analyzing Pb p (a tau > a tau) inside the central collider range with LUXLEP pdf"<<endl;
  }
  if(proc == 2 && collider =="frontal" && pdf=="luxlep" ){
    cout << "Ok! Analyzing Pb p (a tau > a tau) inside the frontal collider range with LUXLEP pdf"<<endl;
  }
  if(proc == 2 && collider =="central" && pdf=="nn23nlo"  ){
    cout << "Ok! Analyzing Pb p (a tau > a tau) inside the central collider range with NN23NLO pdf"<<endl;
  }
  if(proc == 2 && collider =="frontal" && pdf=="nn23nlo" ){
    cout << "Ok! Analyzing Pb p (a tau > a tau) inside the frontal collider range with NN23NLO pdf"<<endl;
  }
  if(proc == 2 && collider =="central" && pdf=="mrst"  ){
    cout << "Ok! Analyzing Pb p (a tau > a tau) inside the central collider range with MRST pdf"<<endl;
  }
  if(proc == 2 && collider =="frontal" && pdf=="mrst" ){
    cout << "Ok! Analyzing Pb p (a tau > a tau) inside the frontal collider range with MRST pdf"<<endl;
  }




  double sigma =0;

  if(proc == 1 && pdf=="luxlep")
  {
    sigma = (1.91); // pp > atau
  }
  if(proc == 1 && pdf=="nn23nlo")
  {
    sigma = (0.280); // pp > atau
  }
  if(proc == 1 && pdf=="mrst")
  {
    sigma = (0.62); // pp > atau
  }

  if(proc == 2 && pdf=="luxlep")
  {
    sigma = (265.6); // Pbp > atau
  }
  if(proc == 2 && pdf=="nn23nlo")
  {
    sigma = (471.3); // Pbp > atau
  }
  if(proc == 2 && pdf=="mrst")
  {
    sigma = (795.2); // Pbp > atau
  }


  Double_t L_int = 10;
  Double_t N_events = nentries;
  Double_t scale1 = ((sigma)/(N_events));



  if(proc == 1 && collider =="central" ){
    cout << "================================================================="<<endl;
    cout << "Raw events"<<endl;
    cout << "================================================================="<<endl;
    cout << "N_events(Total Xs)= " << h_nevents_cuts->GetBinContent(1) <<endl;
    cout << "N_events(pT(a,tau1) > 2 GeV) = " << h_nevents_cuts->GetBinContent(2) <<endl;
    cout << "N_events(CMS Selection (-2.5 < eta < 2.5)) = " << h_nevents_cuts->GetBinContent(3) <<endl;
    cout << "N_events(CMS Selection (-2.5 < eta < 2.5)-CT-PPS) = " << h_nevents_cuts->GetBinContent(4) <<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section after cuts (Diffractive selection)"<<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section(Total Xs) = " << (h_nevents_cuts->GetBinContent(1)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(pT(tau1e2) > 2 GeV) =  " << (h_nevents_cuts->GetBinContent(2)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(CMS(-2.5 < eta < 2.5)) = " << (h_nevents_cuts->GetBinContent(3)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(CMS(-2.5 < eta < 2.5)-CT-PPS) = " << (h_nevents_cuts->GetBinContent(4)/chain->GetEntries())*sigma <<endl;
  }

  if(proc == 1 && collider =="frontal" ){
    cout << "================================================================="<<endl;
    cout << "Raw events"<<endl;
    cout << "================================================================="<<endl;
    cout << "N_events(Total Xs)= " << h_nevents_cuts->GetBinContent(1) <<endl;
    cout << "N_events(pT(a,tau1) > 2 GeV) = " << h_nevents_cuts->GetBinContent(2) <<endl;
    cout << "N_events(CMS Selection (2.0 < eta < 4.5)) = " << h_nevents_cuts->GetBinContent(3) <<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section after cuts (Diffractive selection)"<<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section(Total Xs) = " << (h_nevents_cuts->GetBinContent(1)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(pT(tau1e2) > 2 GeV) =  " << (h_nevents_cuts->GetBinContent(2)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(CMS(2.0 < eta < 4.5)) = " << (h_nevents_cuts->GetBinContent(3)/chain->GetEntries())*sigma <<endl;
  }


  if(proc == 2 && collider =="central" ){
    cout << "================================================================="<<endl;
    cout << "Raw events"<<endl;
    cout << "================================================================="<<endl;
    cout << "N_events(Total Xs)= " << h_nevents_cuts->GetBinContent(1) <<endl;
    cout << "N_events(pT(a,tau1) > 2 GeV) = " << h_nevents_cuts->GetBinContent(2) <<endl;
    cout << "N_events(CMS Selection (-2.5 < eta < 2.5)) = " << h_nevents_cuts->GetBinContent(3) <<endl;
    cout << "N_events(CMS Selection (-2.5 < eta < 2.5)-CT-PPS) = " << h_nevents_cuts->GetBinContent(4) <<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section after cuts (Diffractive selection)"<<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section(Total Xs) = " << (h_nevents_cuts->GetBinContent(1)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(pT(tau1e2) > 2 GeV) =  " << (h_nevents_cuts->GetBinContent(2)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(CMS(-2.5 < eta < 2.5)) = " << (h_nevents_cuts->GetBinContent(3)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(CMS(-2.5 < eta < 2.5)-CT-PPS) = " << (h_nevents_cuts->GetBinContent(4)/chain->GetEntries())*sigma <<endl;

  }

  if(proc == 2 && collider =="frontal" ){
    cout << "================================================================="<<endl;
    cout << "Raw events"<<endl;
    cout << "================================================================="<<endl;
    cout << "N_events(Total Xs)= " << h_nevents_cuts->GetBinContent(1) <<endl;
    cout << "N_events(pT(a,tau1) > 2 GeV) = " << h_nevents_cuts->GetBinContent(2) <<endl;
    cout << "N_events(CMS Selection (2.0 < eta < 4.5)) = " << h_nevents_cuts->GetBinContent(3) <<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section after cuts (Diffractive selection)"<<endl;
    cout << "================================================================="<<endl;
    cout << "Cross Section(Total Xs) = " << (h_nevents_cuts->GetBinContent(1)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(pT(tau1e2) > 2 GeV) =  " << (h_nevents_cuts->GetBinContent(2)/chain->GetEntries())*sigma <<endl;
    cout << "Cross Section(CMS(2.0 < eta < 4.5)) = " << (h_nevents_cuts->GetBinContent(3)/chain->GetEntries())*sigma <<endl;
  }




  h_ntau1->Scale(scale1);
  h_tau1_px->Scale(scale1);
  h_tau1_py->Scale(scale1);
  h_tau1_pz->Scale(scale1);
  h_tau1_pt->Scale(scale1);
  h_tau1_eta->Scale(scale1);
  h_tau1_phi->Scale(scale1);
  h_tau1_energy->Scale(scale1);
  h_tau1_charge->Scale(scale1);
  h_tau1_eta_central->Scale(scale1);
  h_tau1_eta_frontal->Scale(scale1);


  h_na->Scale(scale1);
  h_a_px->Scale(scale1);
  h_a_py->Scale(scale1);
  h_a_pz->Scale(scale1);
  h_a_pt->Scale(scale1);
  h_a_eta->Scale(scale1);
  h_a_phi->Scale(scale1);
  h_a_energy->Scale(scale1);
  h_a_charge->Scale(scale1);

  h_a_eta_central_efixed->Scale(scale1);
  h_tau_eta_central_afixed->Scale(scale1);
  h_a_eta_frontal_efixed->Scale(scale1);
  h_tau_eta_frontal_afixed->Scale(scale1);

  h_pttaua_eta_central->Scale(scale1);
  h_pttaua_eta_frontal->Scale(scale1);

  h_atau_eta_frontal->Scale(scale1);
  h_atau_eta_central->Scale(scale1);

  h_aco_taua->Scale(scale1);
  h_aco_taua_central->Scale(scale1);
  h_aco_taua_frontal->Scale(scale1);

  h_atau_m_inv->Scale(scale1);
  h_atau_m_inv_central->Scale(scale1);
  h_atau_m_inv_frontal->Scale(scale1);

  h_tau1_pt_central->Scale(scale1);
  h_tau1_pt_frontal->Scale(scale1);

  //========================================================================
  //Output File
  TFile* output = 0;

  if(proc == 1 && collider == "central" && pdf=="luxlep")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_atau_atau_13TeV_central_luxlep.root","RECREATE");
  }
  if(proc == 1 && collider == "frontal" && pdf=="luxlep")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_atau_atau_13TeV_frontal_luxlep.root","RECREATE");
  }
  if(proc == 1 && collider == "central" && pdf=="nn23nlo")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_atau_atau_13TeV_central_nn23nlo.root","RECREATE");
  }
  if(proc == 1 && collider == "frontal" && pdf=="nn23nlo")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_atau_atau_13TeV_frontal_nn23nlo.root","RECREATE");
  }
  if(proc == 1 && collider == "central" && pdf=="mrst")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_atau_atau_13TeV_central_mrst.root","RECREATE");
  }
  if(proc == 1 && collider == "frontal" && pdf=="mrst")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_atau_atau_13TeV_frontal_mrst.root","RECREATE");
  }




  if(proc == 2 && collider == "central" && pdf=="luxlep")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_Pbp_atau_atau_8TeV_central_luxlep.root","RECREATE");
  }
  if(proc == 2 && collider == "frontal" && pdf=="luxlep")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_Pbp_atau_atau_8TeV_frontal_luxlep.root","RECREATE");
  }
  if(proc == 2 && collider == "central" && pdf=="nn23nlo")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_Pbp_atau_atau_8TeV_central_nn23nlo.root","RECREATE");
  }
  if(proc == 2 && collider == "frontal" && pdf=="nn23nlo")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_Pbp_atau_atau_8TeV_frontal_nn23nlo.root","RECREATE");
  }
  if(proc == 2 && collider == "central" && pdf=="mrst")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_Pbp_atau_atau_8TeV_central_mrst.root","RECREATE");
  }
  if(proc == 2 && collider == "frontal" && pdf=="mrst")
  {
    output = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_Pbp_atau_atau_8TeV_frontal_mrst.root","RECREATE");
  }

  //========================================================================



  h_ntau1->Write();
  h_tau1_pt->Write();
  h_tau1_px->Write();
  h_tau1_py->Write();
  h_tau1_pz->Write();
  h_tau1_energy->Write();
  h_tau1_eta->Write();
  h_tau1_phi->Write();

  h_tau1_eta_central->Write();
  h_tau1_eta_frontal->Write();

  h_na->Write();
  h_a_pt->Write();
  h_a_px->Write();
  h_a_py->Write();
  h_a_pz->Write();
  h_a_energy->Write();
  h_a_eta->Write();
  h_a_phi->Write();

  h_a_eta_central_efixed->Write();
  h_tau_eta_central_afixed->Write();
  h_a_eta_frontal_efixed->Write();
  h_tau_eta_frontal_afixed->Write();

  h_pttaua_eta_central->Write();
  h_pttaua_eta_frontal->Write();

  h_atau_eta_frontal->Write();
  h_atau_eta_central->Write();

  h_aco_taua->Write();
  h_aco_taua_central->Write();
  h_aco_taua_frontal->Write();

  h_atau_m_inv->Write();
  h_atau_m_inv_central->Write();
  h_atau_m_inv_frontal->Write();

  h_tau1_pt_central->Write();
  h_tau1_pt_frontal->Write();

  h_pte_pta_2D_correlation_central->Write();
  hme_ma_2D_correlation_central->Write();

  output->Close();

  /*
  FILE* PbPb_ee_analysis_analysis;
  PbPb_ee_analysis_analysis = fopen ("out_PbPb_ee_analysis_analysis.txt","a");
  fprintf(PbPb_ee_analysis_analysis,"%s\n\n","================================================================================================================================");
  fprintf(PbPb_ee_analysis_analysis,"%s","PbPb_ee_analysis -> PbPb > dieletron analysis\n");
  fprintf(PbPb_ee_analysis_analysis,"%s\n","--------------------------------------------------------------------------------------------------------------------------------");
  fprintf(PbPb_ee_analysis_analysis,"%s\t\t %f\n","Cross Section(Total Xs) = : ",(h_nevents_cuts->GetBinContent(1)/chain->GetEntries())*sigma);
  fprintf(PbPb_ee_analysis_analysis,"%s\t\t %f\n","Cross Section(m(e,a2) > 5 GeV) = : ",(h_nevents_cuts->GetBinContent(2)/chain->GetEntries())*sigma);
  fprintf(PbPb_ee_analysis_analysis,"%s\t\t %f\n","Cross Section(pT(ea2) > 3 GeV) = : ",(h_nevents_cuts->GetBinContent(3)/chain->GetEntries())*sigma);
  fprintf(PbPb_ee_analysis_analysis,"%s\t\t %f\n","Cross Section(CMS)) = : ",(h_nevents_cuts->GetBinContent(4)/chain->GetEntries())*sigma);
  fprintf(PbPb_ee_analysis_analysis,"%s\n\n","================================================================================================================================");
  fprintf(PbPb_ee_analysis_analysis,"%s\n","");
  fclose(PbPb_ee_analysis_analysis);
  //--------------------------------------------------------------------------
  */

}
