#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

const int ncut = 13;
const int max_ptbin = 100;
const char* ptbin_filename = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/draw_Rgamma/ptbin.txt";
int npt = 0;
double pt_low[max_ptbin] = {0};
double pt_high[max_ptbin] = {0};
double pt_mid[max_ptbin] = {0};
double pt_diff[max_ptbin] = {0};

void read_ptbin(){
  std::ifstream fin_ptbin(ptbin_filename);
  if(!fin_ptbin){
    std::cout << "Cannot open " << ptbin_filename << std::endl;
    exit(1);
  }

  fin_ptbin >> npt;
  if(npt <= 0 || npt > max_ptbin){
    std::cout << "Bad ptbin number in " << ptbin_filename << std::endl;
    exit(1);
  }

  for(int ipt = 0; ipt < npt; ++ipt){
    fin_ptbin >> pt_low[ipt] >> pt_high[ipt];
    if(!fin_ptbin || pt_high[ipt] <= pt_low[ipt]){
      std::cout << "Bad ptbin line " << ipt + 1 << " in " << ptbin_filename << std::endl;
      exit(1);
    }
    pt_mid[ipt] = 0.5 * (pt_low[ipt] + pt_high[ipt]);
    pt_diff[ipt] = pt_high[ipt] - pt_low[ipt];
  }

  fin_ptbin.close();
}

//const char* cut_name[ncut] = {"cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9", "cut10", "cut11", "cut12"};
//const char* tex_cut[ncut] = {"cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9", "cut10", "cut11", "cut12"};

const char* cut_name[ncut] = {
  "nominal",
  "conv_solution",
  "eid0",
  "eid1",
  "eid3",
  "pte02",
  "pte04",
  "ecore03",
  "ecore05",
  "chi2_4",
  "chi2_5",
  "escale099",
  "escale101"
};
const char* tex_cut[ncut] = {
  "nominal",
  "conv_solution",
  "eid0",
  "eid1",
  "eid3",
  "pte02",
  "pte04",
  "ecore03",
  "ecore05",
  "chi2_4",
  "chi2_5",
  "escale099",
  "escale101"
};


void plot(TH1D* h1_ee, TH1D* h1_pi0, int icut, int ipt, double pt1, double pt2){

  TCanvas* c1 = new TCanvas(Form("c_ee_%s_pt%d", cut_name[icut], ipt + 1), "", 900, 700);
  c1->cd();

  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.05);
  gPad->SetTicks(1, 1);

  h1_ee->SetLineColor(4);
  h1_ee->SetMarkerColor(4);
  h1_ee->SetMarkerStyle(2);
  h1_ee->SetMarkerSize(1.0);
  h1_ee->SetLineWidth(1);

  h1_ee->GetXaxis()->SetRangeUser(0.0, 0.35);
  h1_ee->GetXaxis()->SetTitle("m_{ee} [GeV]");
  h1_ee->GetYaxis()->SetTitle("N_{ee}");
  h1_ee->GetXaxis()->SetTitleSize(0.05);
  h1_ee->GetYaxis()->SetTitleSize(0.05);
  h1_ee->GetXaxis()->SetLabelSize(0.045);
  h1_ee->GetYaxis()->SetLabelSize(0.045);
  h1_ee->GetYaxis()->SetTitleOffset(1.0);
  h1_ee->SetTitle("");

  h1_ee->Draw("E0");

  TLegend* leg1 = new TLegend(0.73, 0.56, 0.93, 0.66);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.04);
  leg1->AddEntry(h1_ee, "ee pairs", "lep");
  leg1->Draw();

  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.045);
  latex1.DrawLatex(0.64, 0.83, tex_cut[icut]);
  latex1.DrawLatex(0.64, 0.78, "Pure MC");
  latex1.DrawLatex(0.64, 0.73, Form("%.1f < p_{T}^{ee} < %.1f GeV", pt1, pt2));

  c1->SaveAs(Form("./plot/ee_%s_pt%d.pdf", cut_name[icut], ipt + 1));

  TCanvas* c2 = new TCanvas(Form("c_pi0_%s_pt%d", cut_name[icut], ipt + 1), "", 900, 700);
  c2->cd();

  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.05);
  gPad->SetTicks(1, 1);

  h1_pi0->SetLineColor(4);
  h1_pi0->SetMarkerColor(4);
  h1_pi0->SetMarkerStyle(2);
  h1_pi0->SetMarkerSize(1.0);
  h1_pi0->SetLineWidth(1);

  h1_pi0->GetXaxis()->SetRangeUser(0.0, 0.35);
  h1_pi0->GetXaxis()->SetTitle("m_{ee#gamma} [GeV]");
  h1_pi0->GetYaxis()->SetTitle("N_{ee#gamma}");
  h1_pi0->GetXaxis()->SetTitleSize(0.05);
  h1_pi0->GetYaxis()->SetTitleSize(0.05);
  h1_pi0->GetXaxis()->SetLabelSize(0.045);
  h1_pi0->GetYaxis()->SetLabelSize(0.045);
  h1_pi0->GetYaxis()->SetTitleOffset(1.0);
  h1_pi0->SetTitle("");

  h1_pi0->Draw("E0");

  TLegend* leg2 = new TLegend(0.73, 0.56, 0.93, 0.66);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(h1_pi0, "ee#gamma pairs", "lep");
  leg2->Draw();

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.045);
  latex2.DrawLatex(0.64, 0.83, "Default cuts");
  latex2.DrawLatex(0.64, 0.78, "Pure MC");
  latex2.DrawLatex(0.64, 0.73, Form("%.1f < p_{T}^{ee} < %.1f GeV", pt1, pt2));

  c2->SaveAs(Form("./plot/pi0_%s_pt%d.pdf", cut_name[icut], ipt + 1));

  delete c1;
  delete c2;
  delete leg1;
  delete leg2;
}

void make_1d(){

  read_ptbin();

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);
  gSystem->mkdir("plot", kTRUE);

  TFile* fin = new TFile("2d_mass_pt.root", "READ");
  TFile* fout = new TFile("1d.root", "RECREATE");

  for (int icut = 0; icut < ncut; ++icut) {

    TH2F* h2_ee = (TH2F*) fin->Get(Form("eemass_%s_cent0_sect0", cut_name[icut]));
    TH2F* h2_pi0 = (TH2F*) fin->Get(Form("eegmass_%s_cent0_sect0", cut_name[icut]));
    //TH2F* h2_ee = (TH2F*) fin->Get(Form("eemass_unw_%s_cent0_sect0", cut_name[icut]));
    //TH2F* h2_pi0 = (TH2F*) fin->Get(Form("eegmass_unw_%s_cent0_sect0", cut_name[icut]));

    for (int ipt = 0; ipt < npt; ++ipt) {

      int biny1_ee = h2_ee->GetYaxis()->FindBin(pt_low[ipt] + 1e-6);
      int biny2_ee = h2_ee->GetYaxis()->FindBin(pt_high[ipt] - 1e-6);

      int biny1_pi0 = h2_pi0->GetYaxis()->FindBin(pt_low[ipt] + 1e-6);
      int biny2_pi0 = h2_pi0->GetYaxis()->FindBin(pt_high[ipt] - 1e-6);

      int nbinsx_ee = h2_ee->GetNbinsX();
      double xlow_ee = h2_ee->GetXaxis()->GetXmin();
      double xup_ee = h2_ee->GetXaxis()->GetXmax();
      TH1D* h1_ee = new TH1D(Form("h1_ee_%s_pt%d", cut_name[icut], ipt + 1), "", nbinsx_ee, xlow_ee, xup_ee);
      h1_ee->SetDirectory(0);
      h1_ee->Sumw2();

      for (int binx = 1; binx <= nbinsx_ee; ++binx) {
        double content = 0.0;
        double err2 = 0.0;
        for (int biny = biny1_ee; biny <= biny2_ee; ++biny) {
          double c = h2_ee->GetBinContent(binx, biny);
          double e = h2_ee->GetBinError(binx, biny);
          content += c;
          err2 += e * e;
        }
        h1_ee->SetBinContent(binx, content);
        h1_ee->SetBinError(binx, TMath::Sqrt(err2));
      }

      int nbinsx_pi0 = h2_pi0->GetNbinsX();
      double xlow_pi0 = h2_pi0->GetXaxis()->GetXmin();
      double xup_pi0 = h2_pi0->GetXaxis()->GetXmax();
      TH1D* h1_pi0 = new TH1D(Form("h1_pi0_%s_pt%d", cut_name[icut], ipt + 1), "", nbinsx_pi0, xlow_pi0, xup_pi0);
      h1_pi0->SetDirectory(0);
      h1_pi0->Sumw2();

      for (int binx = 1; binx <= nbinsx_pi0; ++binx) {
        double content = 0.0;
        double err2 = 0.0;
        for (int biny = biny1_pi0; biny <= biny2_pi0; ++biny) {
          double c = h2_pi0->GetBinContent(binx, biny);
          double e = h2_pi0->GetBinError(binx, biny);
          content += c;
          err2 += e * e;
        }
        h1_pi0->SetBinContent(binx, content);
        h1_pi0->SetBinError(binx, TMath::Sqrt(err2));
      }

      fout->cd();
      h1_ee->Write();
      h1_pi0->Write();

      plot(h1_ee, h1_pi0, icut, ipt, pt_low[ipt], pt_high[ipt]);

      delete h1_ee;
      delete h1_pi0;
    }
  }

  fout->Close();
  fin->Close();
}
