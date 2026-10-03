#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include <fstream>
#include <iomanip>
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

const char* input_cut_name[ncut] = {
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

const double x_eeg_min = 0.09;
const double x_eeg_max = 0.18;

void ratio()
{
  read_ptbin();
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  TFile* fin = new TFile("1d.root", "READ");
  TFile* fout = new TFile("graph.root", "RECREATE");

  int color[ncut] = {1, 2, 4, 6, 8, 9, 28, 46, 38, 30, 41, 42, 49};

  double x_ee_min = 0.01;
  double x_ee_max = 0.15;

  TGraphErrors* gr_all[ncut];

  std::ofstream txt("ratio_values.txt");
  txt << std::setprecision(10);
  txt << "cut ipt pt ratio ratio_err" << std::endl;

  for (int icut = 0; icut < ncut; ++icut) {
    double x[max_ptbin];
    double ex[max_ptbin];
    double y[max_ptbin];
    double ey[max_ptbin];

    for (int ipt = 0; ipt < npt; ++ipt) {
      double pt_center = 0.5 * (pt_low[ipt] + pt_high[ipt]);
      //double pt_err = 0.5 * (pt_high[ipt] - pt_low[ipt]);
      double pt_err = 0.0;

      x[ipt] = pt_center;
      ex[ipt] = pt_err;
      y[ipt] = 0.0;
      ey[ipt] = 0.0;

      TH1D* h1_ee = (TH1D*) fin->Get(Form("h1_ee_%s_pt%d", input_cut_name[icut], ipt + 1));
      TH1D* h1_pi0 = (TH1D*) fin->Get(Form("h1_pi0_%s_pt%d", input_cut_name[icut], ipt + 1));

      if(!h1_ee || !h1_pi0){
        std::cout << "Missing histogram for cut " << input_cut_name[icut]
                  << ", pt bin " << ipt + 1 << std::endl;
        continue;
      }

      double err_ee = 0.0;
      double err_eeg = 0.0;

      int bin1_ee = h1_ee->GetXaxis()->FindBin(x_ee_min + 1e-6);
      int bin2_ee = h1_ee->GetXaxis()->FindBin(x_ee_max - 1e-6);

      int bin1_eeg = h1_pi0->GetXaxis()->FindBin(x_eeg_min + 1e-6);
      int bin2_eeg = h1_pi0->GetXaxis()->FindBin(x_eeg_max - 1e-6);

      double yield_ee = h1_ee->IntegralAndError(bin1_ee, bin2_ee, err_ee);
      double yield_eeg = h1_pi0->IntegralAndError(bin1_eeg, bin2_eeg, err_eeg);

      if (yield_ee > 0.0 && yield_eeg > 0.0) {
        double ratio_val = yield_eeg / yield_ee;
        double ratio_err = ratio_val * TMath::Sqrt(
          (err_eeg / yield_eeg) * (err_eeg / yield_eeg) +
          (err_ee / yield_ee) * (err_ee / yield_ee)
        );
        y[ipt] = ratio_val;
        ey[ipt] = ratio_err;
      }
    }

    gr_all[icut] = new TGraphErrors(npt, x, y, ex, ey);
    gr_all[icut]->SetName(Form("gr_ratio_%s", cut_name[icut]));
    gr_all[icut]->SetTitle("");

    gr_all[icut]->SetLineColor(color[icut]);
    gr_all[icut]->SetMarkerColor(color[icut]);
    gr_all[icut]->SetMarkerStyle(20);
    gr_all[icut]->SetMarkerSize(1.7);
    gr_all[icut]->SetLineWidth(1);

    gr_all[icut]->GetXaxis()->SetTitle("p_{T}^{ee} [GeV]");
    gr_all[icut]->GetYaxis()->SetTitle("#LT#epsilon_{f}#GT");
    gr_all[icut]->GetXaxis()->SetTitleSize(0.05);
    gr_all[icut]->GetYaxis()->SetTitleSize(0.05);
    gr_all[icut]->GetXaxis()->SetLabelSize(0.045);
    gr_all[icut]->GetYaxis()->SetLabelSize(0.045);
    gr_all[icut]->GetYaxis()->SetTitleOffset(0.95);

    for (int ipt = 0; ipt < npt; ++ipt) {
      txt << cut_name[icut] << " "
          << ipt + 1 << " "
          << x[ipt] << " "
          << y[ipt] << " "
          << ey[ipt] << std::endl;
    }
    txt << std::endl;

    fout->cd();
    gr_all[icut]->Write();

    TCanvas* c = new TCanvas(Form("c_ratio_%s", cut_name[icut]), "", 900, 700);
    c->cd();

    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);
    gPad->SetTicks(1, 1);

    gr_all[icut]->GetXaxis()->SetLimits(0.0, 6.0);
    gr_all[icut]->Draw("AP E1");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.25, 0.86, "ee#gamma / ee yield ratio, Pure MC");
    latex.DrawLatex(0.25, 0.81, "#LT#epsilon_{f}#GT:  Run12 Cu + Au @ 200 GeV");
    latex.DrawLatex(0.25, 0.76, tex_cut[icut]);

    c->SaveAs(Form("ratio_%s.pdf", cut_name[icut]));

    delete c;
  }

  TCanvas* c_all = new TCanvas("c_ratio_allcuts", "", 900, 700);
  c_all->cd();

  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.05);
  gPad->SetTicks(1, 1);

  gr_all[0]->GetXaxis()->SetTitle("p_{T}^{ee} [GeV]");
  gr_all[0]->GetYaxis()->SetTitle("#LT#epsilon_{f}#GT");
  gr_all[0]->GetXaxis()->SetLimits(0.0, 6.0);
  gr_all[0]->SetMinimum(0.0);
  gr_all[0]->SetMaximum(0.5);
  gr_all[0]->Draw("AP E1");

  for (int icut = 1; icut < ncut; ++icut) {
    gr_all[icut]->Draw("P E1 SAME");
  }

  TLegend* leg = new TLegend(0.18, 0.55, 0.80, 0.78);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  leg->SetNColumns(3);
  leg->SetMargin(0.18);
  leg->SetColumnSeparation(0.05);
  leg->SetEntrySeparation(0.10);
  for (int icut = 0; icut < ncut; ++icut) {
    leg->AddEntry(gr_all[icut], tex_cut[icut], "lep");
  }
  leg->Draw();

  TLatex latex_all;
  latex_all.SetNDC();
  latex_all.SetTextSize(0.04);
  latex_all.DrawLatex(0.25, 0.86, "ee#gamma / ee yield ratio, Pure MC");
  latex_all.DrawLatex(0.25, 0.81, "#LT#epsilon_{f}#GT:  Run12 Cu + Au @ 200 GeV");

  c_all->SaveAs("allcuts.pdf");

  fout->cd();
  c_all->Write();

  delete leg;
  delete c_all;

  txt.close();
  fout->Close();
  fin->Close();
}
