#include "LineShape_bb_pdf.h"

void GetResonanceShape_bb()
{

    double resmass;

    int N_SHAPE_BINS = 2800;
    double SHAPE_BINS_MIN = 0;
    double SHAPE_BINS_STEPS = 5;

    TH1D * h_shape[36];
    TH1D * hist[36];
    TH1D * hist_cdf[36];
    TH1D * hist_new[36];
    TH1D * hist_cdf_new[36];

    TFile *output = new TFile("Resonance_Shapes_bb_test.root","RECREATE");
    double m0[36]= {500., 600.,  700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000.,2100.,2200.,2300.,2400.,2500.,2600.,2700.,2800.,2900.,3000.,3100.,3200.,3300.,3400.,3500.,3600.,3700.,3800.,3900.,4000.};

    for(int j=0; j<36; j++){
        resmass = m0[j];
        cout<<resmass<<endl;
        char *histname = new char[100];
        sprintf(histname,"h_bb_%d",resmass);

        h_shape[j] = new TH1D(histname,"b#bar{b} Resonance Shape",nMassBins,massBoundaries);

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
       hist[j]=new TH1D(TString(histname)+"_shape", "b#bar{b} Resonance Shape", N_SHAPE_BINS, SHAPE_BINS_MIN, SHAPE_BINS_MIN+N_SHAPE_BINS*SHAPE_BINS_STEPS);
       hist_cdf[j]=new TH1D(TString(histname)+"_cdf", "b#bar{b} Resonance Shape CDF", N_SHAPE_BINS, SHAPE_BINS_MIN, SHAPE_BINS_MIN+N_SHAPE_BINS*SHAPE_BINS_STEPS);

       for(int i=1; i<=N_SHAPE_BINS; i++) {
         double mss = hist[j]->GetBinCenter(i);
         double prob = FastQstarBinnedProb(mss);
         hist[j]->SetBinContent(i, prob);
       }
       // re-normalize the histogram
       double integral=hist[j]->GetSumOfWeights();
       hist[j]->Scale(1/integral);

       // compute the cdf
       for(int i=1; i<=N_SHAPE_BINS; i++) {
         double prev=hist_cdf[j]->GetBinContent(i-1);
         double curr=hist[j]->GetBinContent(i);
         hist_cdf[j]->SetBinContent(i, prev+curr);
       }

       output->cd();
       hist_cdf[j]->Write();

       // make the histograms
       hist_new[j]=new TH1D(TString(histname)+"_shape_new", "b#bar{b} Resonance Shape", 14000, 0, 14000);
       hist_cdf_new[j]=new TH1D(TString(histname)+"_cdf_new", "b#bar{b} Resonance Shape CDF", N_SHAPE_BINS, SHAPE_BINS_MIN, SHAPE_BINS_MIN+N_SHAPE_BINS*SHAPE_BINS_STEPS);
       
       for(int i=1; i<=h_shape[j]->GetNbinsX(); i++){

         int bin_min = hist_new[j]->GetXaxis()->FindBin(h_shape[j]->GetXaxis()->GetBinLowEdge(i)+0.5);
         int bin_max = hist_new[j]->GetXaxis()->FindBin(h_shape[j]->GetXaxis()->GetBinUpEdge(i)-0.5);
         double bin_content = h_shape[j]->GetBinContent(i)/double(bin_max-bin_min+1);
         for(int b=bin_min; b<=bin_max; b++){
            hist_new[j]->SetBinContent(b, bin_content);
         }
       }

       //output->cd();
       //hist_new[j]->Write();

       for(int i=1; i<=hist_cdf_new[j]->GetNbinsX(); i++){

         int bin_min = hist_new[j]->GetXaxis()->FindBin(hist_cdf_new[j]->GetXaxis()->GetBinLowEdge(i)+0.5);
         int bin_max = hist_new[j]->GetXaxis()->FindBin(hist_cdf_new[j]->GetXaxis()->GetBinUpEdge(i)-0.5);
         
         double curr = 0;
         for(int b=bin_min; b<=bin_max; b++){
            curr+=hist_new[j]->GetBinContent(b);
         }

         double prev=hist_cdf_new[j]->GetBinContent(i-1);
         
         hist_cdf_new[j]->SetBinContent(i, prev+curr);
       }

       output->cd();
       hist_cdf_new[j]->Write();
    }
    
    output->Close();
 }

