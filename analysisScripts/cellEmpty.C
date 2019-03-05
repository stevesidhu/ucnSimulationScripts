#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include "TGraphErrors.h"
#include "TH1.h"

/*
This script calculates collection efficiency from cell emptying
1.  Make sure the script is inputting the proper .root files from cell emptying simulations
2.  Adjust scale to total simulated / number of energy bins (ucn in each bin)
3.  Confirm solidend cut for detetectors in draw function

*/


void cellEmpty(){
    std::vector<std::string> fileNames ={"bottomCell_emptying.root","topCell_emptying.root"};

    double scale = 20000; //was 20000  total ucn in each energy bin simulated. 
    TCut detectors = "solidend > 500";
    //Load file and trees.
    for(std::string file : fileNames){
        TFile f1(file.c_str());
        TTree* neutronend =(TTree*) f1.Get("neutronend");
        TTree* neutronsnapshot =(TTree*) f1.Get("neutronsnapshot");
        
        //Simulation file name for plot titles
        std::string fileNameString = file.substr(0, file.size()-5);

        TFile outfile((fileNameString + "_DetEff.root").c_str(),"RECREATE");

        TCanvas *t = new TCanvas("d", "d", 800, 600);
        neutronend->Draw("Hstart:tend>>emptyEff(250,0,250,50,0,250e-9)",detectors && "stopID>0");
	    TH2 *emptyEff = (TH2*) gDirectory->Get("emptyEff");
	    emptyEff->Scale(1./scale);
        emptyEff->Draw("LEGO");
       // empty->FitSlicesX(&emptyFit,0,-1,0,"",&arr);
       // TH1* tendVSh = static_cast<TH1*>(arr[1]);
       // t->SetLogy(1);
       // tendVSh->SetTitle(Form("Emptying time vs Energy for %s",fileNameString.c_str()));
       // tendVSh->GetXaxis()->SetTitle("Energy (eV)");
       // tendVSh->GetYaxis()->SetTitle("Emptying Time (s)");
       // tendVSh->Draw();
        t->Print((fileNameString + "_tendVsHstart.eps").c_str());


        delete t;
///////////////////////////////// 1-d hist        
        TCanvas *d = new TCanvas("d", "d", 800, 600);
        neutronend->Draw("Hstart>>detEff(50,0,250e-9)",detectors && "stopID>0");
	    TH1 *detEff = (TH1*) gDirectory->Get("detEff");
	    detEff->Scale(1./scale);
       // empty->FitSlicesX(&emptyFit,0,-1,0,"",&arr);
       // TH1* tendVSh = static_cast<TH1*>(arr[1]);
       // t->SetLogy(1);
        detEff->SetTitle(Form("Collection efficiency of detector vs Energy for %s",fileNameString.c_str()));
        detEff->GetXaxis()->SetTitle("Energy (eV)");
        detEff->GetYaxis()->SetTitle("Collection efficiency");
        detEff->Draw();
        d->Print((fileNameString + "_det_eff.C").c_str());

        delete d;
      
/////////////// Collection rate per second for detector analysis
//
        double normalize;
        TCanvas *b = new TCanvas("b", "b", 800, 600);
        neutronend->Draw("tend>>countRate(200,0,200)",detectors && "stopID>0");
	    TH1 *countRate = (TH1*) gDirectory->Get("countRate");
        countRate->SetTitle(Form("DetectorCountRate %s",fileNameString.c_str()));
        countRate->GetXaxis()->SetTitle("Time (s)");
        countRate->GetYaxis()->SetTitle("Counts");
        normalize = 1.45E6 / countRate->GetEntries();
        countRate->Scale(normalize);
        countRate->Draw("HIST");
        b->Print((fileNameString + "_CountRate.C").c_str());
        b->Print((fileNameString + "_CountRate.eps").c_str());

        delete b;
               
        outfile.Write();
    }
    std::cout<< "finished "<< std::endl;
}

