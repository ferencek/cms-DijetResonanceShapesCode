#include "LineShape_pdf_gg.h"

void GetResonanceShape_gg()
{

    double resmass;

    TH1D * h_shape[36];
    TH1D * h_pdf[36];
    TH1D * h_cdf[36];

    TFile *output = new TFile("Resonance_Shapes_gg.root","RECREATE");
    double m0[36]= {500., 600.,  700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000.,2100.,2200.,2300.,2400.,2500.,2600.,2700.,2800.,2900.,3000.,3100.,3200.,3300.,3400.,3500.,3600.,3700.,3800.,3900.,4000.};

    for(int j=0; j<36; j++){
       resmass = m0[j];
       cout<<resmass<<endl;
       char *histname = new char[100];
       sprintf(histname,"h_gg_%d",resmass);

       h_shape[j] = new TH1D(histname,"gg Resonance Shape",nMassBins,massBoundaries);

       LineShape_pdf(resmass);

       for(int i=0; i<h_shape[j]->GetNbinsX(); i++){
               double mss    =  h_shape[j]->GetBinCenter(i+1);
               double prob   =  FastQstarBinnedProb(mss) ;
               h_shape[j]-> SetBinContent(i+1,prob );
       }

       h_shape[j]-> SetXTitle("Dijet Mass (GeV)");
       h_shape[j]-> SetYTitle("Probability");

       output->cd();
       h_shape[j]-> Write();

       // make the histograms
       h_pdf[j]=new TH1D(TString(histname)+"_pdf", "gg Resonance Shape", 14000, 0, 14000);
       h_cdf[j]=new TH1D(TString(histname)+"_cdf", "gg Resonance Shape CDF", 2800, 0, 14000);

       for(int i=1; i<=h_shape[j]->GetNbinsX(); i++){

         int bin_min = h_pdf[j]->GetXaxis()->FindBin(h_shape[j]->GetXaxis()->GetBinLowEdge(i)+0.5);
         int bin_max = h_pdf[j]->GetXaxis()->FindBin(h_shape[j]->GetXaxis()->GetBinUpEdge(i)-0.5);
         double bin_content = h_shape[j]->GetBinContent(i)/double(bin_max-bin_min+1);
         for(int b=bin_min; b<=bin_max; b++){
            h_pdf[j]->SetBinContent(b, bin_content);
         }
       }

       for(int i=1; i<=h_cdf[j]->GetNbinsX(); i++){

         int bin_min = h_pdf[j]->GetXaxis()->FindBin(h_cdf[j]->GetXaxis()->GetBinLowEdge(i)+0.5);
         int bin_max = h_pdf[j]->GetXaxis()->FindBin(h_cdf[j]->GetXaxis()->GetBinUpEdge(i)-0.5);

         double curr = 0;
         for(int b=bin_min; b<=bin_max; b++){
            curr+=h_pdf[j]->GetBinContent(b);
         }

         double prev=h_cdf[j]->GetBinContent(i-1);

         h_cdf[j]->SetBinContent(i, prev+curr);
       }

       output->cd();
       h_cdf[j]->Write();
    }

    output->Close();
 }

