#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
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
#include "TLatex.h"


void hist_ae_pp(string collider ="",string variable ="")
{

  //==============================================
  //Estilo do histograma
  //=================================================
  //     gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x lhcbStyle.C");

  TFile *_f0=0;
  TFile *_f1=0;
  TFile *_f2=0;

  TH1F *histo1=0;
  TH1F *histo2=0;
  TH1F *histo3=0;

  if(collider == "central" && variable =="eta"){
    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_nn23nlo.root");

    histo1 = (TH1F*)_f0->Get("e1_eta_central");
    histo2 = (TH1F*)_f1->Get("e1_eta_central");
    histo3 = (TH1F*)_f2->Get("e1_eta_central");
  }

  if(collider == "frontal"&& variable =="eta"){

    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_nn23nlo.root");

    histo1 = (TH1F*)_f0->Get("e1_eta_frontal");
    histo2 = (TH1F*)_f1->Get("e1_eta_frontal");
    histo3 = (TH1F*)_f2->Get("e1_eta_frontal");
  }


  if(collider == "central" && variable =="pt"){
    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_nn23nlo.root");

    histo1 = (TH1F*)_f0->Get("ptea_eta_central");
    histo2 = (TH1F*)_f1->Get("ptea_eta_central");
    histo3 = (TH1F*)_f2->Get("ptea_eta_central");
  }

  if(collider == "frontal" && variable =="pt"){
    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_nn23nlo.root");

    histo1 = (TH1F*)_f0->Get("ptea_eta_frontal");
    histo2 = (TH1F*)_f1->Get("ptea_eta_frontal");
    histo3 = (TH1F*)_f2->Get("ptea_eta_frontal");
  }

  if(collider == "central" && variable =="pt_e"){
    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_nn23nlo.root");

    histo1 = (TH1F*)_f0->Get("e1_pt_central");
    histo2 = (TH1F*)_f1->Get("e1_pt_central");
    histo3 = (TH1F*)_f2->Get("e1_pt_central");
  }

  if(collider == "frontal" && variable =="pt_e"){
    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_nn23nlo.root");

    histo1 = (TH1F*)_f0->Get("e1_pt_frontal");
    histo2 = (TH1F*)_f1->Get("e1_pt_frontal");
    histo3 = (TH1F*)_f2->Get("e1_pt_frontal");
  }

  //=================================================================================================
  //Inv Mass
  //=================================================================================================

  if(collider == "central" && variable =="ae_minv"){
    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_central_nn23nlo.root");


    histo1 = (TH1F*)_f0->Get("ae_minv_central");
    histo2 = (TH1F*)_f1->Get("ae_minv_central");
    histo3 = (TH1F*)_f2->Get("ae_minv_central");
  }

  if(collider == "frontal" && variable =="ae_minv"){
    TFile *_f0 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_luxlep.root");
    TFile *_f1 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_mrst.root");
    TFile *_f2 = new TFile("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/analysis_pp_ae_ae_13TeV_frontal_nn23nlo.root");

    histo1 = (TH1F*)_f0->Get("ae_minv_frontal");
    histo2 = (TH1F*)_f1->Get("ae_minv_frontal");
    histo3 = (TH1F*)_f2->Get("ae_minv_frontal");
  }





  //Mundando o título

  histo1->SetLabelSize(0.050,"xy");

  TCanvas *canvas1 = new TCanvas("plot1","plot1",800,600); canvas1->Range(0,0,25,18);
  canvas1->SetRightMargin(0.1917293);
  canvas1->Range(1.637102,-2.886239,5.302827,1.247769);
  canvas1->SetLeftMargin(0.15); canvas1->SetRightMargin(0.2); canvas1->SetTopMargin(0.1); canvas1->SetBottomMargin(0.15);
  canvas1->SetTicky(0);

  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(2);

  histo1->SetLineColor(1);
  histo1->SetLineWidth(4);
  histo1->SetLineStyle(1);

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(4);
  histo2->SetLineStyle(9);

  histo3->SetLineColor(kBlue);
  histo3->SetLineWidth(4);
  histo3->SetLineStyle(2);

  if(variable=="pt_e"){
    histo1->Rebin(4);
    histo2->Rebin(4);
    histo3->Rebin(4);
  }
  if(variable=="eta"){
    histo1->Rebin(10);
    histo2->Rebin(10);
    histo3->Rebin(10);
  }

  if(variable=="ae_minv"){
    histo1->Rebin(40);
    histo2->Rebin(40);
    histo3->Rebin(40);
  }

  histo1->Draw("HIST");
  histo2->Draw("HIST SAME");
  histo3->Draw("HIST SAME");

  if(collider =="central" && variable =="eta"){
    histo1->GetXaxis()->SetRangeUser(-3,3);
    histo1->GetYaxis()->SetRangeUser(0.0001,0.2);
  }
  if(collider =="frontal" && variable =="eta"){
    histo1->GetXaxis()->SetRangeUser(1.5,5);
    histo1->GetYaxis()->SetRangeUser(0.0001,0.2);
  }
  if(collider =="frontal"  && variable =="pt"){
    histo1->GetXaxis()->SetRangeUser(0,15);
    histo1->GetYaxis()->SetRangeUser(0.001,0.1);
  }

  if( collider =="central" && variable =="pt"){
    histo1->GetXaxis()->SetRangeUser(0,15);
    histo1->GetYaxis()->SetRangeUser(0.001,1);
  }

  if(collider =="central" && variable =="pt_e"){
    histo1->GetXaxis()->SetRangeUser(1,15);
    histo1->GetYaxis()->SetRangeUser(0.001,1);
  }
  if(collider =="frontal"  && variable =="pt_e"){
    histo1->GetXaxis()->SetRangeUser(1,15);
    histo1->GetYaxis()->SetRangeUser(0.001,1);
  }

  if( collider =="central" && variable =="ae_minv"){
    histo1->GetXaxis()->SetRangeUser(0,150);
    histo1->GetYaxis()->SetRangeUser(0.0001,1);
  }
  if( collider =="frontal" && variable =="ae_minv"){
    histo1->GetXaxis()->SetRangeUser(0,150);
    histo1->GetYaxis()->SetRangeUser(0.0001,1);
  }

  if(variable =="eta"){
    histo1->GetYaxis()->SetTitle("d#sigma/d#eta_{e}[pb]");
    histo1->GetXaxis()->SetTitle("#eta_{e}");
  }
  if(variable =="pt"){
    histo1->GetYaxis()->SetTitle("d#sigma/dp_T(e#gamma)[pb/GeV]");
    histo1->GetXaxis()->SetTitle("pT(e#gamma)[GeV]");
  }
  if(variable =="pt_e"){
    histo1->GetYaxis()->SetTitle("d#sigma/dp_{T}^{e}[pb/GeV]");
    histo1->GetXaxis()->SetTitle("p_{T}^{e}[GeV]");
  }
  if(variable =="ae_minv"){
    histo1->GetYaxis()->SetTitle("d#sigma/dm_{e#gamma}[pb/GeV]");
    histo1->GetXaxis()->SetTitle("m_{e#gamma}[GeV]");
  }


  histo1->GetYaxis()->SetTitleSize(0.07);
  histo1->GetXaxis()->SetTitleSize(0.07);
  histo1->GetYaxis()->SetNdivisions(510);
  histo1->GetXaxis()->SetNdivisions(510);

 if(variable=="pt_e"){
  TLegend *leg1 = new TLegend(0.45,0.6,0.85,0.8,NULL,"brNDC"); //x1,y1,x2,y2   (x1,x2 comp) (y1,y2 alt)
  leg1->SetFillColor(0); leg1->SetFillStyle(0); leg1->SetBorderSize(0);
  leg1->SetFillColor(0); leg1->SetFillStyle(0); leg1->SetBorderSize(0);

  leg1->AddEntry(histo1,"LUXLep ","l");
  leg1->AddEntry(histo2,"APFEL (MRST)","l");
  leg1->AddEntry(histo3,"APFEL (NN23NLO)","l");
  leg1->SetFillColor(0); // white background
  leg1->SetBorderSize(0); // get rid of the box
  leg1->SetTextSize(0.035); // set text size
  leg1->Draw();
}

  TLatex* beam = new TLatex(0.17,0.92,"pp collision - #sqrt{s} = 13 TeV");
  beam->SetNDC();
  beam->SetTextSize(0.05);
  beam->Draw("same");


    TLatex* proc = new TLatex(0.5,0.85,"#gamma e^{+} #rightarrow #gamma e^{+}");
    proc->SetNDC();
    proc->SetTextSize(0.045);
//    proc->Draw("same");

  if(collider =="central"){
    TLatex* central = new TLatex(0.5,0.80,"Central selection");
    central->SetNDC();
    central->SetTextSize(0.045);
    central->Draw("same");
  }

  if(collider =="frontal"){
    TLatex* frontal = new TLatex(0.5,0.80,"Forward selection");
    frontal->SetNDC();
    frontal->SetTextSize(0.045);
    frontal->Draw("same");
  }


  //CMcms->Draw("same");
  canvas1->SetLogy();
  canvas1->Update();

  if(collider =="central" && variable=="eta"){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/eta_e_central_pp.pdf");
  }
  if(collider =="frontal" && variable=="eta" ){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/eta_e_frontal_pp.pdf");
  }
  if(collider =="central" && variable=="pt"){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/pt_ea_central_pp.pdf");
  }
  if(collider =="frontal" && variable=="pt"){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/pt_ea_frontal_pp.pdf");
  }
  if(collider =="central" && variable=="pt_e"){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/pt_e_central_pp.pdf");
  }
  if(collider =="frontal" && variable=="pt_e"){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/pt_e_frontal_pp.pdf");
  }
  if(collider =="central" && variable=="ae_minv"){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/minv_ea_central_pp.pdf");
  }
  if(collider =="frontal" && variable=="ae_minv"){
    canvas1->Print("/home/danbia/workspace/build_hepmc2/examples/ae_pdf/plots/minv_ea_frontal_pp.pdf");
  }


}
