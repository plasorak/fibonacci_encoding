//____________________________________________________________________ 
//  
// $Id$ 
// Author: Pierre Lasorak <plasorak@Pierres-MacBook-Pro.local>
// Update: 2020-01-21 11:38:31+0000
// Copyright: 2020 (C) Pierre Lasorak
//
//
#ifndef __CINT__
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#endif

int Analyse()
{
  int max_noise = 20;
  TH2D*     compression_factor_th2  = new TH2D    ("comp_fact_th2" , ";Exponential noise normalisation;Compression factor", max_noise, 1, max_noise+1, 20, 0, 5);
  TProfile* compression_factor_prof = new TProfile("comp_fact_prof", ";Exponential noise normalisation;Compression factor", max_noise, 1, max_noise+1, "S");
  TCanvas c;
  
  for (int i=3; i<max_noise; ++i) {
    TFile f(Form("~/Desktop/empty_gen_g4_noise%i_fibonacci_encoding_hist.root", i), "READ");
    TH1D* histo = (TH1D*)f.Get("daq/Compression_factor");
    for (int ibin=1; ibin<=histo->GetNbinsX(); ++ibin) {
      double b_center  = histo->GetBinCenter(ibin);
      double b_content = histo->GetBinContent(ibin);
      if (b_content>0) {
        compression_factor_th2 ->Fill(i, b_center, b_content);
        compression_factor_prof->Fill(i, b_center, b_content);
      }
    }
    f.Close();
  }
  
  compression_factor_th2 ->SetStats(0);
  compression_factor_prof->SetStats(0);
  compression_factor_prof->SetLineColor(kBlack);
  compression_factor_prof->SetLineWidth(2);
  c.Print("compression_noise.pdf[");
  compression_factor_th2->Draw("COLZ");
  c.Print("compression_noise.pdf");
  compression_factor_prof->Draw("");
  c.Print("compression_noise.pdf");
  compression_factor_th2 ->Draw("COLZ");
  compression_factor_prof->Draw("SAME");
  c.Print("compression_noise.pdf");
  c.Print("compression_noise.pdf]");
  return 0;
}

#ifndef __CINT__
int main(int argc, char** argv)
{
  return Analyse();
}
#endif

//____________________________________________________________________ 
//  
// EOF
//
