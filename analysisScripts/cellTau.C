#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include "TGraphErrors.h"
#include "TH1.h"

void cellTau(){
    int binsize = 5;
    std::vector<std::string> fileNames ={"topCell.root","bottomCell.root"};

    //Load file and trees.
    for(std::string file : fileNames){
        TFile f1(file.c_str());
        TTree* neutronend =(TTree*) f1.Get("neutronend");
        TTree* neutronsnapshot =(TTree*) f1.Get("neutronsnapshot");
        
        //Simulation file name for plot titles
        std::string fileNameString = file.substr(0, file.size()-5);

        TFile outfile((fileNameString + "_hist.root").c_str(),"RECREATE");
        TObjArray arr;
        TF1 tauFit("tauFit","[0]*exp(-x/[1])");
        tauFit.SetParameters(100,50);
 
        TF1 dualFit("dualFit","[0]*exp(-x/[1]) + [2]*exp(-x/[3])");
        dualFit.SetParameters(10000,100,10000,30);

        TCanvas *d = new TCanvas("d", "d", 800, 600);
        TH1D* fallout = new TH1D("fallout","Histogram Title", 55, -25E-9, 250E-9);
        neutronend->Draw("Hend:tend >> fallout" ,"stopID >0");
        fallout->Draw();
		d->Print((fileNameString + "_totalLoss.eps").c_str());
		delete d;

        

        TCanvas *t = new TCanvas("d", "d", 800, 600);
        neutronend->Draw("Hend:tend>>lifetime(200,0,200,50,0,250e-9)","stopID>0 || stopID ==-4","LEGO");
	    TH2 *lifetime = (TH2*) gDirectory->Get("lifetime");
        lifetime->FitSlicesX(&tauFit,0,-1,0,"",&arr);
        TH1* tauVSh = static_cast<TH1*>(arr[1]);
        TH1* tauVSh_const = static_cast<TH1*>(arr[0]);
        t->SetLogy(1);
        tauVSh->SetTitle(Form("Storage lifetime vs Energy for %s",fileNameString.c_str()));
        tauVSh->GetXaxis()->SetTitle("Energy (eV)");
        tauVSh->GetYaxis()->SetTitle("Storage lifetime (s)");
        tauVSh->Draw();
        TLine *line = new TLine(0,100,250E-9,100);
        line->SetLineColor(kRed);
        line->Draw();
        TLine *line2 = new TLine(0,80,250E-9,80);
        line2->SetLineColor(kBlue);
        line2->Draw();
        TLine *line3 = new TLine(0,120,250E-9,120);
        line3->SetLineColor(kGreen);
        line3->Draw();
        t->Print((fileNameString + "_TauVsHend.C").c_str());
        tauVSh_const->Draw();
        t->Print((fileNameString + "_TauConstVsHend.eps").c_str());


        delete t;
        
        
        
        
        TCanvas *e = new TCanvas("d", "d", 800, 600);
        neutronend->Draw("tend>>time(200,0,200)","stopID>0 || stopID==-4");
        TH1* time = (TH1*) gDirectory->Get("time");
        time -> Fit(&dualFit,"","",20,199);
        double slope = dualFit.GetParameter(1);
        std::cout << 1/slope << std::endl;
		e->Print((fileNameString + "_time.eps").c_str());
        outfile.Write();
    }
    std::cout<< "finished "<< std::endl;
}
