#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include <TH3.h>
#include <algorithm>
#include "SimResult.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include <numeric>
#include "TROOT.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include <fstream>
#include "Rtypes.h"
#include "TLegend.h"
#include <chrono> 
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

    int counter = 0;

    double beamHeating = 10; //W
    double UCN_production = 2.62E6; //UCN/W
    double cellVolume = 16286 *2; // singular cell volume in cm^3

//Experimental constants used in simulation and analysis
//    double emptyingTime = 55;  //cell emptying time !DEFAULT VALUE!
//    double t_edm = 132;  // Ramsey storage time in edm cell: NOT CELL STORAGE LIFETIME!   !DEFAULT VALUE!
    double degaussingTime = 1800; //seconds = 30 mins.
    double stableField = 16.14; //hr stable field period per day
    double fillsPerCycle = 8;  //
    double polarityTime = 150; //polarity change time. seconds
    double hbar = 6.71E-16; //eV*s
    double polarization = 0.5; // (when only simulating high-field seekers)
    double EField = 12000; // V/cm
//    double t_Ramsey = 124; //s
//    double tau_walls = 120;
    double tau_Xe = 6566.92;
    double spinTrans = 0.5;
//    double transToDet = 0.9; //transport efficiency to detector (collection efficiency)
    double detEff = 0.9; //Detector efficiency
    double alpha = 0.95; //Initial polarization (alpha_0)
    double T2 = 500; //Transverse relaxation time /
    double t_wait = 2; //Time before Pi/2 pulse
    double t_pulse = 2; //duration of Pi/2 pulse
    double T1 = 1000; //longitudinal relaxation time
    double detSpinTrans = 0.95; //Polarization loss during emptying
    double Panalyzer = 0.9; //analyzing power of analyzer
//	double cellUpperZBound = 0.2; //Highest z coordinate within the UPPER EDM cell, relative to the origin, in m
//	double cellLowerZBound = -0.2; //Lowest z coordinate within the LOWER EDM cell, relative to the origin, in m
	double cellCenter = 0; //The z coordinate of the center of the two edm cells !DEFAULT VALUE!
	Double_t cellLowerXBound = 5.3; //Minimum x coordinate corresponding to the EDM cell, relative to the origin, in m. This is the distance from the origin to the start of the cell
	double valveOpenTime = 100;
	double activeTime = 200; //Source active time (simulation)
    
    
	//Tcuts for cells for simulations with EDM Cell fill geometry
	TCut topCellCut = TCut(Form("solidend>274 &&solidend<276 && zend>%f", cellCenter));
	TCut bottomCellCut = TCut(Form("solidend>274 &&solidend<276 && zend<%f", cellCenter));
	TCut cellsCut = (topCellCut || bottomCellCut);	
    TCut sourceTime = TCut(Form("tend<%.0f", valveOpenTime)); // time cut before valve open time
    TCut valveOpenTimeCut = TCut(Form("tend==%.0f", valveOpenTime)); // time cut for snapshot at valve open time


//This structure contains all the .root files and important parameters for a study ex. KinkHeight.
struct Study{
    std::string parameter_name;
    std::vector<double> parameters;
    std::vector<std::string> filenames; 
    std::vector<SimResult> results;
   // std::vector< std::pair<double,double>> cellCuts;  //top cell max height, bottom cell min height
    std::vector<double> cellCenter;
};

//Type definition to allow the SimResult struct to pass on an individual member of a vector.  This is used for making parameter plots with makeGraph.
typedef double (SimResult::*ResultMember);



//Function declariations
/////////////////////////
// Plots different members of a SimResult for all simulations in a study.
TGraphErrors* makeGraph(Study &study,std::string xName,std::string yName ,ResultMember mem_xs, ResultMember mem_ys, ResultMember mem_xerr, ResultMember mem_yerr);

//Plots energy binned historgrams for Cell filling, source storage lifetime, and filling efficiency for different operational modes in the filling simulations
void transportHistograms(std::string FileName, SimResult &results);

//Returns number of days to reach
double days(int mode, SimResult &result,std::string FileName,TH1F* toptauhist, TH1F* bottomtauhist,TH2* topemptyinghist,TH2* bottomemptyinghist);


//SimulationResults: Function that analyzes each individual .root file in a study.  
//
//After analyzing transport and Days to reach, it optimizes days to reach for t_optimalEDM, t_optimalEmptying,
void analyzeSim(SimResult &results, std::string FileName, int mode, double cellCenter = 0);

//Total number in cell
std::vector<double> cellTotal(int mode, SimResult &results);

//Numerical error calculator
//double numError(int mode, SimResult tempResult, std::string FileName);

//Stolen function: used to better fit the ideal filling time.  This fixes double peaks and takes global maximum.
double LargerProbability(int value, int ref);
std::array<double, 3> GetMaxBinRange(const TH1* hist, const double confidence_interval = 0.6827);

//Draw different histograms for cell finish and emptying
void DrawCellHistogram(SimResult &result, std::string fileName);




//Main function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void testSimulationAnalysisV4()
{

    ///
    //Set up timer to measure run time
    auto start = std::chrono::high_resolution_clock::now();


    gROOT->SetBatch(1);
    int mode = 1;   //Options are 0:batch, 1: steady beam, 2, steady state.

//////////////////    
// Create a study    
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Create an instance of different studies, simulation files grouped together with varying parameters.
    //Tester study
    Study test_study1 = {"Tester",{45},{"kinkAng2V45deg.root"},{},{0}};

    //Kink Angle
     Study KinkAngle_study{"KinkAngle",{45,50,55,60,65,70,75,80,85,90},{"kinkAng45deg.root","kinkAng50deg.root","kinkAng55deg.root","kinkAng60deg.root","kinkAng65deg.root","kinkAng70deg.root","kinkAng75deg.root","kinkAng80deg.root","kinkAng85deg.root","kinkAng90deg.root"},{},{0,0,0,0,0,0,0,0,0,0}};

     Study KinkAngle2V_study{"KinkAngle2V",{45,50,55,60,65,70,75,80,85,90},{"kinkAng2V45deg.root","kinkAng2V50deg.root","kinkAng2V55deg.root","kinkAng2V60deg.root","kinkAng2V65deg.root","kinkAng2V70deg.root","kinkAng2V75deg.root","kinkAng2V80deg.root","kinkAng2V85deg.root","kinkAng2V90deg.root"},{},{0,0,0,0,0,0,0,0,0,0}};

///////////////////////////////////////////////////////////////////////////////////////////////////////////    
    
   //Tester function 
   std::vector<Study> allStudies = {test_study1};    

   //Specific studies
  // std::vector<Study> allStudies = {KinkAngle_study,KinkAngle2V_study};
  // std::vector<Study> allStudies = {KinkAngle2V_study};


    //Run analysis functions on all simulations in each study and plot graphs
    for(Study study : allStudies)
    {
        SimResult tempResult;
        for(uint i = 0; i < study.filenames.size(); ++i) 
        {
            tempResult.parameter = study.parameters[i];
            analyzeSim(tempResult, study.filenames[i], mode, study.cellCenter[i]);
            study.results.push_back(tempResult);
        }


   //Plots of whole study vs parameters.
        //Days vs P
       // makeGraph(study,"Parameter","DaysToReach", &SimResult::parameter, &SimResult::daysToReach, NULL, &SimResult::errDaysToReach);
        makeGraph(study,"Parameter","DaysToReach", &SimResult::parameter, &SimResult::daysToReach, NULL, NULL);

        //Tau vs Days
//        makeGraph(study,"Tau_source","DaysToReach", &SimResult::sourceStorageLifetime, &SimResult::daysToReach, &SimResult::sourceStorageLifetimeError, &SimResult::errDaysToReach);

        //Tau vs P
//        makeGraph(study,"Parameter","Tau_source", &SimResult::parameter, &SimResult::sourceStorageLifetime, NULL, &SimResult::sourceStorageLifetimeError);

        //Days vs Extraction
       // makeGraph(study,"ExtractionEfficiency","DaysToReach", &SimResult::efficiency, &SimResult::daysToReach, &SimResult::efficiencyError, &SimResult::errDaysToReach);
        makeGraph(study,"ExtractionEfficiency","DaysToReach", &SimResult::efficiency, &SimResult::daysToReach, &SimResult::efficiencyError, NULL);

        //Total Number in Cell
        makeGraph(study,"Parameter","UCNinCell", &SimResult::parameter, &SimResult::numInCellReal, NULL, &SimResult::uNumInCellReal);
        
        //System storage lifetime 
//        makeGraph(study,"Parameter","SystemTau", &SimResult::parameter, &SimResult::systemTau, NULL, &SimResult::uSystemTau);
        
        //Transport efficiency
//        makeGraph(study,"Parameter","Transport_efficiency", &SimResult::parameter, &SimResult::transportEff, NULL, &SimResult::uTransportEff);

        //Avg Tau cells 
        makeGraph(study,"Parameter","AvgTauCells", &SimResult::parameter, &SimResult::avgTauCells, NULL, &SimResult::uavgTauCells);

        //Weighted Collection Efficiency 
        makeGraph(study,"Parameter","WeightedDetEff", &SimResult::parameter, &SimResult::WDetEff, NULL, &SimResult::uWDetEff);
       
        //Optimal t_edm 
        makeGraph(study,"Parameter","Optimal_Tedm", &SimResult::parameter, &SimResult::optimalTedm, NULL, NULL);
        makeGraph(study,"Parameter","Optimal_FillTime", &SimResult::parameter, &SimResult::fillTime, NULL, NULL);
        makeGraph(study,"Parameter","Optimal_SourcePumpingTime", &SimResult::parameter, &SimResult::sourcePumpingTime, NULL, NULL);

       //Optimal emtpying time 
        makeGraph(study,"Parameter","Optimal_EmptyingTime", &SimResult::parameter, &SimResult::emptyingTime, NULL, NULL);

       //Total collected 
        makeGraph(study,"Parameter","TotalCollectedUCN", &SimResult::parameter, &SimResult::totalCollected , NULL, NULL);

/*
    double scaleTop =0.;
    double scaleBottom =0.;
/////////////////////////////// Combined Filling plots on one canvas ////////////////////
    EColor colors[] = {kRed,kBlue,kBlack,kViolet,kGreen,kPink,kOrange,kTeal,kMagenta,kGray}; 
    gStyle->SetOptFit(1111);
//////Cell Filling    
    TCanvas *can = new TCanvas("can" , "can", 400, 900);
    can->Divide(2,7);

    auto *topCellEHist = new TH1F("topCellEHist","topCellMeanEnergy",50,0,250e-9);
    gStyle->SetOptStat(0);

    for(uint i = 0; i < study.filenames.size(); ++i){
     can->cd(2*i+1);
        scaleTop = 1E6 / study.results[i].topCellESpec[0].GetEntries();
       // study.results[i].topCellESpec[0].Scale(1.0/study.results[i].efficiency);
        study.results[i].topCellESpec[0].Scale(scaleTop);
        study.results[i].topCellESpec[0].SetLineColor(colors[i]);
        study.results[i].topCellESpec[0].GetYaxis()->SetRangeUser(0, 3E8);
        study.results[i].topCellESpec[0].GetXaxis()->SetTitle("Energy (neV)");
        study.results[i].topCellESpec[0].GetYaxis()->SetTitle("UCN counts in cell");
        study.results[i].topCellESpec[0].Draw("HIST same");
      
      can->cd(2*i + 2);   
        scaleBottom = 1E6 / study.results[i].bottomCellESpec[0].GetEntries();
        study.results[i].bottomCellESpec[0].Scale(scaleBottom);
        study.results[i].bottomCellESpec[0].SetLineColor(colors[i]);
        study.results[i].bottomCellESpec[0].GetYaxis()->SetRangeUser(0, 3E8);
        study.results[i].bottomCellESpec[0].GetXaxis()->SetTitle("Energy (neV)");
        study.results[i].bottomCellESpec[0].GetYaxis()->SetTitle("UCN counts in cell");
        study.results[i].bottomCellESpec[0].Draw("HIST same");
    }
    can->Print(Form("%s_CellFillingEnergy.png",study.parameter_name.c_str()));
    delete can;
      
///////////////////////////////////////////////////////////////////// Finish (Ramsey surviors) on one plot /////////////
    TCanvas *can2 = new TCanvas("can2" , "can2", 400, 900);
    can2->Divide(2,7);

    auto *topCellEHistFinish = new TH1F("topCellEHistFinish","topCellMeanEnergyFinish",50,0,250e-9);
    gStyle->SetOptStat(0);

    //auto leg2 = new TLegend(0.7, 0.7, 0.89, 0.89);
    for(uint i = 0; i < study.filenames.size(); ++i){
        can2->cd(2*i +1);
//        scale = 1E7 / study.results[i].topCellESpecFinish.GetEntries();
        study.results[i].topCellESpecFinish.Scale(scaleTop);
        study.results[i].topCellESpecFinish.SetLineColor(colors[i]);
        study.results[i].topCellESpecFinish.GetYaxis()->SetRangeUser(0, 3E8);
        study.results[i].topCellESpecFinish.GetXaxis()->SetTitle("Energy (neV)");
        study.results[i].topCellESpecFinish.GetYaxis()->SetTitle("UCN surviving");
        study.results[i].topCellESpecFinish.Draw("HIST same");  
        // leg2->AddEntry(&study.results[i].topCellESpecFinish,Form("%s %i",study.parameter_name.c_str(), (int)study.parameters[i]));
    
        can2->cd(2*i +2);
  //      scale = 1E7 / study.results[i].bottomCellESpecFinish.GetEntries();
        study.results[i].bottomCellESpecFinish.Scale(scaleBottom);
        study.results[i].bottomCellESpecFinish.SetLineColor(colors[i]);
        study.results[i].bottomCellESpecFinish.GetYaxis()->SetRangeUser(0, 3E8);
        study.results[i].bottomCellESpecFinish.GetXaxis()->SetTitle("Energy (neV)");
        study.results[i].bottomCellESpecFinish.GetYaxis()->SetTitle("UCN surviving");
        study.results[i].bottomCellESpecFinish.Draw("HIST same"); 
    }
    can2->Print(Form("%s_CellFinishEnergy.eps",study.parameter_name.c_str()));
    delete can2;

////////////////////////////////////////////////////////////////////// Detected histograms on one plot ////////////////
//Top Cell  
    TCanvas *can3 = new TCanvas("can3" , "can3", 400, 900);
    can3->Divide(2,7);
    auto *topCellEHistDetected = new TH1F("topCellEHistDetected","topCellMeanEnergyDetected",50,0,250e-9);
    gStyle->SetOptStat(0);
    for(uint i = 0; i < study.filenames.size(); ++i){
        can3->cd(2*i+1);
    //    scale = 1E7 / study.results[i].topCellESpecFinish.GetEntries();
        study.results[i].topCellESpecDetected.Scale(scaleTop);
        study.results[i].topCellESpecDetected.SetLineColor(colors[i]);
        study.results[i].topCellESpecDetected.GetYaxis()->SetRangeUser(0, 2E8);
        study.results[i].topCellESpecDetected.GetXaxis()->SetTitle("Energy (neV)");
        study.results[i].topCellESpecDetected.GetYaxis()->SetTitle("UCN collected");
        study.results[i].topCellESpecDetected.Draw("HIST same"); 

        can3->cd(2*i +2);
      //  scale = 1E7 / study.results[i].bottomCellESpecFinish.GetEntries();
        study.results[i].bottomCellESpecDetected.Scale(scaleBottom);
        study.results[i].bottomCellESpecDetected.SetLineColor(colors[i]);
        study.results[i].bottomCellESpecDetected.GetYaxis()->SetRangeUser(0, 2E8);
        study.results[i].bottomCellESpecDetected.GetXaxis()->SetTitle("Energy (neV)");
        study.results[i].bottomCellESpecDetected.GetYaxis()->SetTitle("UCN collected");
        study.results[i].bottomCellESpecDetected.Draw("HIST same"); 
    }
    can3->Print(Form("%s_CellDetectedEnergy.eps",study.parameter_name.c_str()));
    delete can3;
*/
///////////////////////

        std::ofstream resultsFile;
        resultsFile.open(Form("%s_results.txt",study.parameter_name.c_str()));
        resultsFile <<"file, parameter,totalSimulated,totalProducedReal, fillTime, fillTimeError, efficiency, efficiencyError, sourceStorageLifetime, sourceStorageLifetimeError, systemTau, uSystemTau, totalUCNinSim, uTotalUCNinSim, numInCellSim, uNumInCellSim, numInCellReal, uNumInCellReal, UCNdensityInCells, transportEff, uTransportEff, avgTauCells, uavgTauCells, wTopCellDetEff, uWTopCellDetEff, wBottomCellDetEff, uWBottomCellDetEff, weightedDetEff, uWeightedDetEff, topCellRamseySurvivorsSim, bottomCellRamseySurvivorsSim, topCellFilledReal, topMeanEFill, uTopMeanEFill, bottomCellFilledReal, bottomMeanEFill, uBottomMeanEFill, topCellFinishedReal, topMeanEFinish, uTopMeanEFinish, bottomCellFinishedReal, bottomMeanEFinish, uBottomMeanEFinish, topCellDetectedReal, topMeanEDetected, uTopMeanEDetected, bottomCellDetectedReal, bottomMeanEDetected, uBottomMeanEDetected, totalDetected, perFilledFromProduced, perDetectedFromProduced, perDetectedFromFilled, perTopCellDetectedFromProduced, perBottomCellDetectedFromProduced, optimalT_edm, optimalEmptyingTime, daysToReach, errDaysToReach \n";
        
        for(uint i = 0; i < study.filenames.size(); ++i) {
            resultsFile <<study.filenames[i] << ", "<< study.results[i].parameter <<", "<<study.results[i].totalSimulated<<", " << study.results[i].totalProducedReal <<", "<< study.results[i].fillTime << ", " <<study.results[i].fillTimeError << ", " <<study.results[i].efficiency << ", " << study.results[i].efficiencyError << ", " << study.results[i].sourceStorageLifetime << ", " << study.results[i].sourceStorageLifetimeError << ", " << study.results[i].systemTau << ", " << study.results[i].uSystemTau << ", " << study.results[i].totalUCN << ", " << study.results[i].uTotalUCN << ", "<<study.results[i].numInCellSim << ", "<<study.results[i].uNumInCellSim<<", " << study.results[i].numInCellReal << ", " << study.results[i].uNumInCellReal << ", " << study.results[i].numInCellReal/cellVolume << ", "<< study.results[i].transportEff << ", " << study.results[i].uTransportEff << ", " << study.results[i].avgTauCells << ", " << study.results[i].uavgTauCells << ", "<< study.results[i].wTopDetEff << ", " <<study.results[i].uWTopDetEff << ", "<< study.results[i].wBottomDetEff << ", " <<study.results[i].uWBottomDetEff << ", "  <<study.results[i].WDetEff<< ", " << study.results[i].uWDetEff<< ", " <<study.results[i].topCellRamsey << ", " <<study.results[i].bottomCellRamsey << ", "  <<study.results[i].topCellTotalFilled << ", "  <<study.results[i].topMeanFill<< ", " << study.results[i].uTopMeanFill << ", " <<study.results[i].bottomCellTotalFilled << ", " <<study.results[i].bottomMeanFill << ", "<<study.results[i].uBottomMeanFill << ", " <<study.results[i].topCellTotalFinished << ", "<<study.results[i].topMeanFinish << ", "<<study.results[i].uTopMeanFinish << ", "<<study.results[i].bottomCellTotalFinished << ", "<<study.results[i].bottomMeanFinish << ", "<<study.results[i].uBottomMeanFinish << ", "<<study.results[i].topCellTotalDetected << ", "<<study.results[i].topMeanDetected << ", "<<study.results[i].uTopMeanDetected << ", " <<study.results[i].bottomCellTotalDetected << ", "<<study.results[i].bottomMeanDetected << ", "<<study.results[i].uBottomMeanDetected << ", "<< study.results[i].topCellTotalDetected + study.results[i].bottomCellTotalDetected <<", "<<study.results[i].numInCellReal/ study.results[i].totalProducedReal << ", "<< (study.results[i].topCellTotalDetected + study.results[i].bottomCellTotalDetected)/study.results[i].totalProducedReal <<", "<< (study.results[i].topCellTotalDetected + study.results[i].bottomCellTotalDetected)/study.results[i].numInCellReal<<", "<< study.results[i].topCellTotalDetected/study.results[i].totalProducedReal <<", "<< study.results[i].bottomCellTotalDetected/study.results[i].totalProducedReal<<", " << study.results[i].optimalTedm << ", " << study.results[i].emptyingTime << ", "<< study.results[i].daysToReach << ", " << study.results[i].errDaysToReach << "\n";
        }
        resultsFile.close();      
    }

    ///
    //Stop timer
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout<< "Time is " << duration.count() <<" s" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////
//Fucntions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void transportHistograms(std::string fileName, SimResult &results)
{
	//Open file//
	TFile f1(fileName.c_str());
        TTree *neutronend =(TTree*) f1.Get("neutronend");
        TTree *neutronsnapshot =(TTree*) f1.Get("neutronsnapshot");
    results.totalSimulated = neutronend->GetEntries();    
    results.productionRateSim = results.totalSimulated/200.;

//////////////////////////////
	//New set up: make a 3-D histogram box plot of UCN filling the cells with time on one axis (1s binning) and Hend on other axis (5 nev binning).
//////////////////////////////
	//Get optimal filling time from minimizing days to reach function.
	//This gives energy spectrum and tend for top and bottom cells
    
    //Used for assymetric feeders or when cell center is not located at z = 0
    //Histogram draw TString for top, bottom, both cells with energy shifts
	double zeroPoint = cellCenter * 1.025E-7;  // used to shift histograms if cell center is not z=0 in simulations
    std::cout<< "zeropoint energy is : "<<zeroPoint << std::endl;
	TString drawBoth = TString::Format("Hend-%g: tend: tstart>>bothCell(200,0,200,200,0,200,50,0,250e-9)",zeroPoint);
	TString drawAll = TString::Format("Hend-%g: tend: tstart>>All(200,0,200,200,0,200,50,0,250e-9)",zeroPoint);
	TString drawTop = TString::Format("Hend-%g: tend: tstart>>topCell(200,0,200,200,0,200,50,0,250e-9)",zeroPoint);
	TString drawBottom = TString::Format("Hend-%g: tend: tstart>>bottomCell(200,0,200,200,0,200,50,0,250e-9)",zeroPoint);
	
	//Both cells
	TCanvas *f = new TCanvas("f", "f", 800, 600);
    neutronsnapshot->Draw(drawBoth,"xend>5.3" && cellsCut); // xend>5.3 is relic from previous simulations without EDM cell fill volume
    TH3D *bothCell = (TH3D*) gDirectory->Get("bothCell");
    results.bothCellFill = *bothCell;
    delete f;

    //Top cell
	TCanvas *f3 = new TCanvas("f", "f", 800, 600);
    neutronsnapshot->Draw(drawTop,"xend>5.3" && topCellCut);
    TH3D* topCell = (TH3D*) gDirectory->Get("topCell");
    results.topCellFill = *topCell;
    delete f3;
    
    //Bottom cell
	TCanvas *f2 = new TCanvas("f", "f", 800, 600);
    neutronsnapshot->Draw(drawBottom,"xend>5.3" && bottomCellCut);
    TH3D *bottomCell = (TH3D*) gDirectory->Get("bottomCell");
    results.bottomCellFill = *bottomCell;
    delete f2;
///////////////////////////////    
///////////////////////////////
//Plots for energy dependent source storage lifetime
///////////////////////////////
/*	TCanvas *c = new TCanvas("c1", "c1", 800, 600);
	neutronsnapshot->Draw("tend-tstart :Hend : tstart>>sourceLifetime(200,0,200,50,0,250E-9,200,0,200)");
	TH3D *sourceLifetime = (TH3D*) gDirectory->Get("sourceLifetime");

    results.sourceLifetime = *sourceLifetime;
    delete c;
*/
///////////////////////////////
///////////////////////////////
//Scatter plot vs Time
///////////////////////////////
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    neutronsnapshot->Draw(drawAll); // xend>5.3 is relic from previous simulations without EDM cell fill volume
	TH3D *scatter = (TH3D*) gDirectory->Get("All");
	scatter->SetDirectory(0);
    results.scatter = scatter;
    delete c1;

}



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//Days to reach function.  
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
double days(int mode, SimResult &result,std::string FileName,TH1F* toptauhist, TH1F* bottomtauhist,TH2* topemptyinghist,TH2* bottomemptyinghist)
{
    counter++;
	
	//start by calculating fill time from histogram.  Then optimize for fill time
	//look at time projectection of fill time and run fill time maximization
	//from that pull energy and plug into source storage lifetime histogram (maybe)
	//set up source pumping time as a parameter
	//look at efficiency by energy
	//plug all this into NiTop, NiBottom
	
    TH1F *topDetEff = (TH1F*)topemptyinghist->ProjectionY("_detEff", 0,topemptyinghist->GetXaxis()->FindBin(result.emptyingTime));
    TH1F *bottomDetEff = (TH1F*)bottomemptyinghist->ProjectionY("_detEff",0,bottomemptyinghist->GetXaxis()->FindBin(result.emptyingTime));
	
	int fillBin = result.topCellFill.GetXaxis()->FindBin(result.fillTime);
	int startBin = result.topCellFill.GetXaxis()->FindBin(result.fillTime - result.sourcePumpingTime);
	int endBin = result.topCellFill.GetYaxis()->FindBin(result.fillTime);
	std::printf("FillTime: %i, startBin: %i, emptyingTime: %f, t_edm: %f \n",fillBin,startBin,result.emptyingTime,result.optimalTedm);

	//TCut fillTimeCut = TCut(Form("tend == %.0f", result.fillTime)); // TCut for end of fillTime
	//TCut modeCuts = TCut(Form("tstart>%f && tstart < %f", 100 - result.sourcePumpingTime, result.fillTime));

    result.topCellFill.GetXaxis()->SetRange(startBin,fillBin); //tstart axis.  need to check this.  This would mean we start producing UCN at this source pumping time.
    result.topCellFill.GetYaxis()->SetRange(endBin,endBin);  //tend axis
    result.topCellESpec =(TH1F*)result.topCellFill.Project3D("z");
    result.topCellESpec->Scale(1E6/result.totalSimulated);

	fillBin = result.bottomCellFill.GetXaxis()->FindBin(result.fillTime);
	startBin = result.bottomCellFill.GetXaxis()->FindBin(result.fillTime - result.sourcePumpingTime);
	endBin = result.bottomCellFill.GetYaxis()->FindBin(result.fillTime);
    result.bottomCellFill.GetXaxis()->SetRange(startBin,fillBin);
    result.bottomCellFill.GetYaxis()->SetRange(endBin,endBin);
    result.bottomCellESpec =(TH1F*)result.bottomCellFill.Project3D("z");
    result.bottomCellESpec->Scale(1E6/result.totalSimulated);

    std::cout<< "Total UCN that reach the cell at the end of filling time: " << result.topCellESpec->GetEntries() + result.bottomCellESpec->GetEntries() <<std::endl;


    //Create a results.scatter plot with same binning to calculate total ucn created in pumping time.
//	fillBin = result.scatter->GetXaxis()->FindBin(result.fillTime);
//	startBin = result.scatter->GetXaxis()->FindBin(result.fillTime - result.sourcePumpingTime);
//	endBin = result.scatter->GetYaxis()->FindBin(result.fillTime);
//    result.scatter->GetXaxis()->SetRange(startBin,fillBin);
//    result.scatter->GetYaxis()->SetRange(endBin,endBin);
//    TH1F *totalCount =(TH1F*)result.scatter->Project3D("z");
    //result.bottomCellESpec->Scale(1E6/result.totalSimulated);




	//Step 4 - Number of UCN in Cell 
//	result.numInCellSim = neutronsnapshot->GetEntries(cellsCut && modeCuts && fillTimeCut); //use same for batch and steady BEAM


	result.numInCellSim = result.topCellESpec->GetEntries() + result.bottomCellESpec->GetEntries(); //use same for batch and steady BEAM
	//result.numInCellSim = neutronsnapshot->GetEntries(cellsCut && fillTimeCut); //use same for batch and steady BEAM	
	result.uNumInCellSim = sqrt(result.numInCellSim);

    std::cout<< "Total UCN that reach the cell at the end of filling time (neutronend): " << result.numInCellSim <<std::endl;

	//Step 5 - Total number of UCN//
	//Total UCN created in the source that could make it in the same time
	result.totalUCN = result.productionRateSim * result.sourcePumpingTime; 
//	result.uTotalUCN =sqrt(pow(productionRate*uSourcePumpingTime,2)+ pow(sqrt(produced),2));  //includes error in production rate
	result.uTotalUCN =sqrt(result.totalUCN);  //includes error in production rate

/*    //Transport efficiency
//	if(mode == 0 || mode == 1) {
    //double numbTransported = neutronsnapshot->GetEntries(modeCuts && valveOpenTimeCut);
    double numTransported = neutronsnapshot->GetEntries(valveOpenTimeCut);
    double uNumTransported = sqrt(numTransported);
	result.transportEff = result.numInCellSim/numTransported; // Number in cell at the end of ideal filling time divided by ucn produced within last 2.4* tau_storage of the source.  roughly 90 percent of them would survive up scatter
	result.uTransportEff = abs(result.numInCellSim/numTransported) * sqrt(pow(result.uNumInCellSim/result.numInCellSim,2) + pow(uNumTransported/numTransported,2));  //confirmed by SS: 2019-02-28
*/	
	//Step 6 - Calculate Efficiencies//
	//total cell UCN divided by the total created in the source in that mode.
    result.efficiency = result.numInCellSim/result.totalUCN;
    result.efficiencyError = result.efficiency * sqrt(pow(result.uNumInCellSim/result.numInCellSim,2) + pow(result.uTotalUCN/result.totalUCN,2));  //confirmed by SS: 2019-02-28		
    std::cout << "Efficiency: " << result.efficiency << std::endl;

/*
    //Step 7 - System Storage Lifetime//
	TCut systemStorageCut = TCut(Form("tstart < %f && tstart > %f",valveOpenTime+5, valveOpenTime));	
	neutronsnapshot->Draw("tend-tstart>>systemStorage", systemStorageCut); //(Form("tend-tstart>>systemStorage(%f, %f, %f)", activeTime - valveOpenTime + 5, valveOpenTime+5, activeTime), systemStorageCut);
        TF1 systemStorFit("fit","expo", 0, activeTime-valveOpenTime-5); //,valveOpenTime+5, activeTime);
        systemStorage->Fit(&systemStorFit,"R");
        Double_t slope2 = systemStorFit.GetParameter(1);
        Double_t slope2Error = systemStorFit.GetParError(1);
        result.systemTau = -1/slope2;
        result.uSystemTau = slope2Error/(slope2*slope2);
/////////////////////////////////////////////////////////////////////////////////////////////////
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////
//Estimates UCN surviving the Ramsey cycle for top and bottom cell individually, Parameter is used for asymmetric feeders to refer to different emptying.root files for each cell.  The parameter should be in he name of the root file 
//UCN remaining: look at N_0 that fill the cell and estimate how many survive after Ramsey cycle, emptying, collection as a function of energy (5 nev binning)
    // Top and bottom spectrum come from transport() function.  This is the energy spectrum histogram of EDM cell filling broken into 5nev bins.  These contain simulation numbers of filling.
    // TopTauFile and bottomTauFile contains a histogram of storage lifetime of ucn in edm cell broken up by ucn energy (broken into 5nev bins).
    // Top and bottom emptying file contains histogram of collection efficiency vs engery spectrum of ucn (broken into 5nev bins).

    //Initial ucn in both cells  NOT USED.  INTEGRATE INSTEAD
    double N0 = result.topCellESpec->GetEntries() + result.bottomCellESpec->GetEntries();


    double nFCellTop = 0.; // Top cell Ramsey survivors in simulation
    double nFCellBottom = 0.;
    double survivors = 0.; // Total ucn that get collected by detectors in simulation
    double dsurvivors = 0.;
    double wTopDet = 0.;  //Total number of UCN collected from top cell
    double wBottomDet =0.;

    //We now have 6 1-D histograms(top and bottom filling, collection efficiency, ramsey survivors), broken up into 5nev energy bins, that we can now multiply bin by bin. 
    //Filling histogram gives us simulated number of UCN in each cell.  Collection efficiency and ramsey survivor histograms give us probabilities.

    // looping over each energy bin (5 nev bins from 0 to 250 nev, i think)
    for (int i = 1; i < toptauhist->GetNbinsX() - 1; ++i)
    {
        double Nitop = result.topCellESpec->GetBinContent(i);  //Initial number of UCN, at end of filling, in that energy bin
        double dNitop = result.topCellESpec->GetBinError(i);
        double tauitop = toptauhist->GetBinContent(i); //tau storage for that energy bin
        double dtauitop = toptauhist->GetBinError(i);
        double detEffTop = topDetEff->GetBinContent(i);
        double uDetEffTop = topDetEff->GetBinError(i);  //Emptying and collection efficiency of that energy bin
        double Nibottom = result.bottomCellESpec->GetBinContent(i);
        double dNibottom = result.bottomCellESpec->GetBinError(i);
        double tauibottom = bottomtauhist->GetBinContent(i);
        double dtauibottom = bottomtauhist->GetBinError(i);
        double detEffBottom = bottomDetEff->GetBinContent(i);
        double uDetEffBottom = bottomDetEff->GetBinError(i);

// Just used to give weighted mean averages to be displayed and quoted.  Not used in UCNremaining calculation.
        //Number of UCN of that energy bin collected
        wTopDet += detEffTop * Nitop*exp(-result.optimalTedm/tauitop); 
        wBottomDet+= detEffBottom * Nibottom*exp(-result.optimalTedm/tauibottom);   
        //Number that survive Ramsey cycle (right before when the emptying valve opens)
        nFCellTop += Nitop*exp(-result.optimalTedm/tauitop);
        nFCellBottom += Nibottom*exp(-result.optimalTedm/tauibottom);

//        //Most important
        //Total number of collected ucn from each energy bin in simulated number of ucn.
        survivors += detEffTop * Nitop*exp(-result.optimalTedm/tauitop) + detEffBottom * Nibottom*exp(-result.optimalTedm/tauibottom);   
        if (dtauitop != 0)
            dsurvivors += pow(uDetEffTop*Nitop*exp(-result.optimalTedm/tauitop), 2) + pow(dNitop*detEffTop*exp(-result.optimalTedm/tauitop), 2) + pow(Nitop*detEffTop*exp(-result.optimalTedm/tauitop)/tauitop/tauitop*dtauitop, 2);
        if (dtauibottom != 0)
            dsurvivors += pow(uDetEffBottom* Nibottom *exp(-result.optimalTedm/tauibottom), 2) + pow(dNibottom *detEffBottom *exp(-result.optimalTedm/tauibottom), 2) + pow(Nibottom*detEffBottom*exp(-result.optimalTedm/tauibottom)/tauibottom/tauibottom*dtauibottom, 2);   
    }  //end for loop looping over each bin.
    

    //taking the square root of all of the uncertainties
    dsurvivors = sqrt(dsurvivors);

    //Used to calculed weight def eff
    double NF = nFCellTop + nFCellBottom;	// Number that finish the Ramsey cycle
    result.avgTauCells = result.optimalTedm / log(N0/NF);
    result.uavgTauCells = result.optimalTedm * (sqrt(1./N0 + 1./NF) * NF/N0) ; 
    
    //Used for display purposes only
    double topDetMean = wTopDet / nFCellTop;
    double bottomDetMean = wBottomDet / nFCellBottom;
    double utopDetMean =  sqrt(1.0/nFCellTop);
    double ubottomDetMean = sqrt(1.0/nFCellBottom);

    //Used for display purposes only
    result.wTopDetEff = (topDetMean /pow(utopDetMean,2))/(1./pow(utopDetMean,2));
    result.uWTopDetEff = sqrt(1/( 1./pow(ubottomDetMean,2)));
    result.wBottomDetEff = ( bottomDetMean /pow(ubottomDetMean,2))/(1./pow(ubottomDetMean,2));
    result.uWBottomDetEff = sqrt(1/(1./pow(ubottomDetMean,2)));
    result.WDetEff = (topDetMean /pow(utopDetMean,2) + bottomDetMean /pow(ubottomDetMean,2))/(1./pow(utopDetMean,2) + 1./pow(ubottomDetMean,2));
    result.uWDetEff = sqrt(1/(1./pow(ubottomDetMean,2) + 1./pow(ubottomDetMean,2)));


    result.survivalprob  = survivors/N0;  // survival probability
    std::cout << "Survivors: " << survivors << ", survivalProb: " << result.survivalprob << std::endl;
    result.dSurvivalprob = dsurvivors/N0;
    result.topCellRamsey = nFCellTop;
    result.bottomCellRamsey = nFCellBottom;
    
/*    
    //Determining the time we irriadiate for
    double sourcePumpingTime; // function: t_irradiation in RP spreadsheet. 
    if (mode ==0)
        sourcePumpingTime = result.sourceStorageLifetime * 2.4;
    else if (mode ==1)
        sourcePumpingTime = result.sourceStorageLifetime * 2.4 + result.fillTime;
    else
        sourcePumpingTime = 100;
       // sourcePumpingTime = result.fillingTime;
*/
    double cyclesToReach;
    double cyclesPerDay;
    double sensitivityPerFill;
   

    //cyclesPerDay  (unit of time is in minutes)
    if (mode ==1)
        cyclesPerDay = (stableField * 60)/((fillsPerCycle * (result.emptyingTime + result.optimalTedm + result.sourcePumpingTime) + 2*polarityTime + degaussingTime/10)/60);
        
        //beam on at same time as emptying and only switching E field once, degaussing only on off hours and once per day.
        //cyclesPerDay = (stableField * 60)/((fillsPerCycle * (std::max(result.emptyingTime, result.sourcePumpingTime) +result.optimalTedm) + polarityTime)/60) -1 ;
    else
        cyclesPerDay = (stableField * 60)/((fillsPerCycle * (result.fillTime + result.emptyingTime + result.optimalTedm + result.sourcePumpingTime) + 2*polarityTime + degaussingTime/10)/60);
   
    std::cout<< "Cycles per day: "<<cyclesPerDay <<std::endl;
    //cyclesToReach
    //sensitivityPerFill
    double N_0;  //Initial number of UCN
    double N_det;  //Number of UCN that get detected
    double alpha_after;  //Visibility after Ramsey sequence

    //Number of actual UCN in cell at start of Ramsey cycle.  Total production * time to produce * time to fill * filling efficiency
    N_0 = spinTrans * beamHeating * UCN_production * result.sourcePumpingTime * result.efficiency;

    //Number of actual UCN collected
    N_det = N_0 * result.survivalprob * exp(-result.optimalTedm/tau_Xe) * detEff;  //detEff = efficiency of the detector, global constant

    alpha_after = alpha * exp(-((result.optimalTedm- t_wait - 2*t_pulse) / T2) - (t_wait + 2* t_pulse)/T1) * detSpinTrans * Panalyzer;

    sensitivityPerFill = hbar / (2 * (result.optimalTedm - t_wait - 2*t_pulse) * EField * sqrt(N_det) * alpha_after);
    std::cout<<"N0: " << N_0 << ", N_det: " << N_det << ", Sensitivity per fill: " <<sensitivityPerFill <<std::endl;

    cyclesToReach = pow(sensitivityPerFill/sqrt(fillsPerCycle)/(1E-27),2);

    double daysToReach =  cyclesToReach / cyclesPerDay;

    std::cout<< "Temp days to reach is " << daysToReach<< "\n"<<std::endl;
    return daysToReach;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyzeSim(SimResult &result, std::string FileName, int mode, double cellCenter)
{
         std::cout << FileName <<" , " << mode << std::endl;  
         std::string fileNameString = FileName.substr(0, FileName.size()-5); //Strips .root from file name.
	     std::vector<double> cellNumber;
         
         //function call for transport
	     transportHistograms(FileName, result);

         //initial values:
         result.fillTime = 120.;
         result.sourcePumpingTime = 80.;
         result.emptyingTime = 60.;
         result.optimalTedm = 130.;
	     std::printf("Initial Values are FillTime: %f, sourcePumpingTime: %f, emptyingTime: %f, t_edm: %f \n",result.fillTime,result.sourcePumpingTime,result.emptyingTime,result.optimalTedm);


         /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
           //// topCell = old cell shape with dPS walls.  kinkHtTop (v1) is new cell shape with NiP cell walls.  kinkHtTop2V is new cell shape with dPS walls.
           //  TFile toptaufile("topCell_hist.root");
           //  TFile bottomtaufile("bottomCell_hist.root");
          //   TFile toptaufile("kinkHtTopTau_hist.root");
          //   TFile bottomtaufile("kinkHtBottomTau_hist.root");
             TFile toptaufile("kinkHtTopTau2V_hist.root");
             TFile bottomtaufile("kinkHtBottomTau2V_hist.root");
             TH1F *toptauhist = (TH1F*)toptaufile.Get("lifetime_1");
             TH1F *bottomtauhist =(TH1F*)bottomtaufile.Get("lifetime_1");
             toptauhist->SetDirectory(0);
             bottomtauhist->SetDirectory(0);
         //    TFile topemptyingfile("topCell_emptying_DetEff.root");
         //    TFile bottomemptyingfile("bottomCell_emptying_DetEff.root");
         //    TFile topemptyingfile("kinkHtTopEmptying_DetEff.root");
         //    TFile bottomemptyingfile("kinkHtBottomEmptying_DetEff.root");
             TFile topemptyingfile("kinkHtTopEmptying2V_DetEff.root");
             TFile bottomemptyingfile("kinkHtBottomEmptying2V_DetEff.root");
         //    TFile topemptyingfile(Form("asymFeedEmpt2VTop%icm_DetEff.root",(int)result.parameter));
         //    TFile bottomemptyingfile(Form("asymFeedEmpt2VBottom%icm_DetEff.root",(int)result.parameter));
            
             TH2 *topemptyinghist = (TH2*)topemptyingfile.Get("emptyEff");
             TH2 *bottomemptyinghist = (TH2*)bottomemptyingfile.Get("emptyEff");
             bottomemptyinghist->SetDirectory(0);
             topemptyinghist->SetDirectory(0);
         
             //Collection efficiency limiting the x-axis by optimal emptying time. Not including ucn detected after optimal emptying time.
         
             toptaufile.Close();
             bottomtaufile.Close();
             topemptyingfile.Close();
             bottomemptyingfile.Close();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////         
    TH1F *topDetEff = (TH1F*)topemptyinghist->ProjectionY("_detEff", 0,topemptyinghist->GetXaxis()->FindBin(result.emptyingTime));
    TH1F *bottomDetEff = (TH1F*)bottomemptyinghist->ProjectionY("_detEff",0,bottomemptyinghist->GetXaxis()->FindBin(result.emptyingTime));

    for(int i = 101; i<150; i++)
    {
    	int fillBin = result.topCellFill.GetXaxis()->FindBin(i);
    	int startBin = result.topCellFill.GetXaxis()->FindBin(10);
    	int endBin = result.topCellFill.GetYaxis()->FindBin(i);

        TH1F *fillingTop, *fillingBottom;

        TCanvas top("top","top");
        result.topCellFill.GetXaxis()->SetRange(startBin,fillBin); //tstart axis.  need to check this.  This would mean we start producing UCN at this source pumping time.
        result.topCellFill.GetYaxis()->SetRange(endBin,endBin);  //tend axis
        fillingTop =(TH1F*)result.topCellFill.Project3D("z");
        fillingTop->Scale(1E6/result.totalSimulated);
        fillingTop->Draw("HIST");
        top.Print(Form("FillingTop_%i.png",i));

        TCanvas bottom("bottom","bottom");
        result.bottomCellFill.GetXaxis()->SetRange(startBin,fillBin); //tstart axis.  need to check this.  This would mean we start producing UCN at this source pumping time.
        result.bottomCellFill.GetYaxis()->SetRange(endBin,endBin);  //tend axis
        fillingBottom =(TH1F*)result.bottomCellFill.Project3D("z");
        fillingBottom->Scale(1E6/result.totalSimulated);
        fillingBottom->Draw("HIST");
        bottom.Print(Form("FillingBottom_%i.png",i));
     }
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////         
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////         
    
        //Starting multivariate minimization of days() function
         ROOT::Math::Functor f([mode,&result,FileName,toptauhist,bottomtauhist,topemptyinghist,bottomemptyinghist](const double *x){
                 result.fillTime = x[0];
                 result.sourcePumpingTime = x[1];
                 result.emptyingTime = x[2];
                 result.optimalTedm = x[3];
                 double d = days(mode, result, FileName,toptauhist,bottomtauhist,topemptyinghist,bottomemptyinghist);
                 return d;
                 },4);

         ROOT::Math::GSLSimAnMinimizer min;
         min.SetTolerance(5);
         min.SetMaxFunctionCalls(1);
         min.SetMaxIterations(1);
         min.SetFunction(f);
         min.SetVariable(0, "fillTime", 120, 2.0);
         min.SetVariable(1, "sourcePumpingTime", 80, 2.0);
         min.SetVariable(2, "emptyingTime", 60, 2.0);
         min.SetVariable(3, "t_edm", 130, 2.0);
         min.SetVariableLimits(0,115,145);
         min.SetVariableLimits(1,60,100);
         min.SetVariableLimits(2,30,75);
         min.SetVariableLimits(3,80,150);

         std::cout<< "minimizing" << std::endl;
         min.Minimize();

         std::cout << "Found minimum: f(" << min.X()[0] << ", " << min.X()[1] <<", " << min.X()[2] << ", " << min.X()[3] << ") = " << min.MinValue() <<  std::endl;
         std::cout << "Found minimum error:" << min.Errors() <<  std::endl;

         std::cout << "Should be f(" << result.optimalTedm << ", " << result.emptyingTime << ", " << result.survivalprob << ")" << std::endl;
         std::cout << "counter value: " << counter <<std::endl;

         result.fillTime = min.X()[0];
         result.sourcePumpingTime = min.X()[1];
         result.optimalTedm = min.X()[3];
         result.emptyingTime = min.X()[2];
       //  result.uEmptyingTime = min.X()[];
        // result.uOptimalTedm = min.X()[];
        
         //After finding the optimal emptyingTime and t_edm, run the Days to reach function
         double daysToReach = days(mode, result, FileName,toptauhist,bottomtauhist,topemptyinghist,bottomemptyinghist);
         result.daysToReach = daysToReach;

         //Numerical approximation in the error in days to reach
      //   double errDaysToReach = numError(mode, result, FileName);
      //   result.errDaysToReach = errDaysToReach;

         std::cout << "Days to reach: " <<  result.daysToReach << " +/- "<< result.errDaysToReach <<"\n\n";
 
         //Estimating real UCN in actual source   
         cellNumber = cellTotal(mode, result);
    	 result.numInCellReal = cellNumber[0];
    	 result.uNumInCellReal = cellNumber[1];
/*     
         result.fillTime = 120.;
         result.sourcePumpingTime = 100.;
         result.emptyingTime = 60.;
         result.optimalTedm = 130.;
         result.daysToReach = days(mode, result, FileName,toptauhist,bottomtauhist,topemptyinghist,bottomemptyinghist);
         std::cout << "Days to reach: " <<  result.daysToReach << " +/- "<< result.errDaysToReach <<"\n\n";
         //Drawing filling, surviving, emptying energy spectra for each cell for each simulation
      //   DrawCellHistogram(result, fileNameString);
*/

}

/*
std::vector<double> SA_prepareInitialState(int seed)
{
    return {rand() % 30 + 1894, rand() % 30 + 1231, rand() % 20 + 1893, rand() % 1294};
}

double SA_deltaEnergy(int mode, SimResult &result, std::string FileName,TH1F* toptauhist, TH1F* bottomtauhist,TH1F* topDetEff,TH1F* bottomDetEff)
{

}
*/


// Total number of UCN in cell approximation
/////////////////////////////////////////////////////////////////
std::vector<double> cellTotal(int mode, SimResult &results)
{
    double UCNproductionPers = UCN_production * beamHeating;
    double productionTimeSteadyState = 100;

    double totalSourceUCN;
    std::vector<double> cellUCN;

    if(mode ==0)
        totalSourceUCN = results.sourceStorageLifetime * 2.4 * UCNproductionPers;
    else if(mode ==1)
        totalSourceUCN = (results.sourceStorageLifetime * 2.4 + results.fillTime) * UCNproductionPers;
    else if(mode ==2)
        totalSourceUCN = UCNproductionPers * productionTimeSteadyState;
    else{
        std::cout <<"You have entered incorrect mode" << std::endl;
        return cellUCN;
    }
   
    results.totalProducedReal = totalSourceUCN;
    cellUCN.push_back(totalSourceUCN * results.efficiency * polarization);
    cellUCN.push_back(cellUCN[0] * results.efficiencyError / results.efficiency);

    return cellUCN;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors* makeGraph(Study &study,std::string xName, std::string yName, ResultMember mem_xs, ResultMember mem_ys, ResultMember mem_xerr, ResultMember mem_yerr)
{
    int simNumber =0;

    std::vector<double> xs;
    std::vector<double> xerrs;
    std::vector<double> ys;
    std::vector<double> yerrs;
    for (auto it = study.results.begin(); it != study.results.end(); ++it) {
        xs.push_back((*it) .* mem_xs);
        xerrs.push_back(mem_xerr == NULL ? 0.0 : (*it) .* mem_xerr);
        ys.push_back((*it) .* mem_ys);
        yerrs.push_back(mem_yerr == NULL ? 0.0 : (*it) .* mem_yerr);
        simNumber++;
    }
    gStyle->SetOptFit(1111);
    TCanvas *can = new TCanvas("can" , "can", 1200, 900);
    can->SetWindowSize(1200, 900);
       TGraphErrors *gr = new TGraphErrors (simNumber, &(xs[0]), &(ys[0]), &(xerrs[0]) , &(yerrs[0]));
       gr -> SetMarkerStyle(20);
       gr -> SetTitle(Form("%s: %s vs %s", study.parameter_name.c_str(), yName.c_str(), xName.c_str())); 
       gr -> GetXaxis()-> SetTitle(xName.c_str());
       gr -> GetYaxis()-> SetTitle(yName.c_str());
       gr -> GetXaxis() -> SetTitleSize(0.05);
       gr -> GetXaxis() -> SetTitleOffset(0.75);
       gr -> GetYaxis() -> SetTitleSize(0.05);
       gr -> GetYaxis() -> SetTitleOffset(0.9);
       gr -> Draw("Ap");
    can->Print(Form("%s_%s_vs_%s.png",study.parameter_name.c_str(), yName.c_str() ,xName.c_str()));
    delete can;
    return gr;

}

/*
//Numerical error in days to reach function.
double numError(int mode, SimResult tempResult, std::string FileName)
{
    SimResult result = tempResult;
    const int n = 5; //number of parameters    

    std::string names[n] = {"efficiency","fillTime","Tau_source","emptyingTime","RamseyTime"};
    double error = 0;
    double errorPar[n] ={result.efficiencyError, result.fillTimeError, result.sourceStorageLifetimeError,result.uEmptyingTime,result.uOptimalTedm}; 
   // double x = days(mode, result, FileName, toptauhist, bottomtauhist, topDetEff, bottomDetEff); 
    double x = days(mode, result, FileName); 
    
    TF1 f("f",[mode,&result,FileName](double *x, double *p){
                result.efficiency = p[0];
                result.fillTime = p[1];
                result.sourceStorageLifetime = p[2];
                result.emptyingTime = p[3];
                result.optimalTedm = p[4];
               // return days(mode, result, FileName, toptauhist, bottomtauhist, topDetEff, bottomDetEff);
                return days(mode, result, FileName);
            },0,10000,n);
   
    //Set parameters so GradiantPar has a step size for optimization.  Uses 0.01*parameterError for step size
    f.SetParameters(result.efficiency, result.fillTime, result.sourceStorageLifetime, result.emptyingTime, result.optimalTedm);

    for (int i = 0; i < n; ++i)
    {  
        f.SetParError(i, errorPar[i]);
        error = error + abs(errorPar[i]*f.GradientPar(i, &x));
        std::cout << "Error in "<< names[i] <<" is: " << errorPar[i] << " times gradient " << f.GradientPar(i, &x) << " = "<< error <<std::endl;
    }
    error = error + abs(result.dSurvivalprob*x/result.survivalprob);
    std::cout << "Error in survivalprob is: " << result.dSurvivalprob << " times gradient " << x/result.survivalprob << " = "<< error <<std::endl;
    return error;
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// return probability that poisson process that produced result "value" has a higher rate than the poisson process that produced result "ref"
// see "conditional test" from http://www.ucs.louisiana.edu/~kxk4695/JSPI-04.pdf
double LargerProbability(int value, int ref){
  double prob = 0.;
  for (int i = value; i <= value + ref; ++i){
    prob += ROOT::Math::binomial_pdf(i, 0.5, value + ref);
  }
  return 1. - prob;
}

// return range of bins that have cumulative probabilites within the confidence interval
// first returned value is lower edge of leftmost bin
// second returned value is center of maximum bin
// third returned value is upper edge of rightmost bin
std::array<double, 3> GetMaxBinRange(const TH1* hist, const double confidence_interval){
  std::vector<double> probs;
  for (int i = 1; i < hist->GetNbinsX(); ++i){ // loop over bins of hist (excluding underflow and overflow bins)
    double prod = 1.;
    for (int j = 1; j < hist->GetNbinsX(); ++j){ // for each bin calculate probability that it is produced by poisson process with a rate larger than every other bin
      if (j != i)
        prod = prod*LargerProbability(hist->GetBinContent(i), hist->GetBinContent(j)); // accumulate probabilities
    }
    probs.push_back(prod); // add probabilites to list
  }
  double integral = std::accumulate(probs.begin(), probs.end(), 0.); // calculate normalization (sum of all probabilities)
  std::vector<int> maxbins;
  while (std::accumulate(probs.begin(), probs.end(), 0.) > (1. - confidence_interval)*integral){
    // remove bin with max probability and add its bin number to maxbin until sum of all probabilities falls below (1 - confidence interval)*normalization
    auto maxbin = std::max_element(probs.begin(), probs.end());
    maxbins.push_back(std::distance(probs.begin(), maxbin) + 1);
    probs.erase(maxbin);
  }

  // first element in maxbins contains overall maximum bin, min and max element give range of bins
  return {hist->GetBinLowEdge(*std::min_element(maxbins.begin(), maxbins.end())),
          hist->GetBinCenter(maxbins.front()),
          hist->GetBinLowEdge(*std::max_element(maxbins.begin(), maxbins.end()) + 1)};
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////
void DrawCellHistogram(SimResult &result, std::string fileName)
{

    double scale;
//    TFile toptaufile("topCell_hist.root");
//    TFile bottomtaufile("bottomCell_hist.root");
//    TFile toptaufile("kinkHtTopTau_hist.root");
//    TFile bottomtaufile("kinkHtBottomTau_hist.root");
    TFile toptaufile("kinkHtTopTau2V_hist.root");
    TFile bottomtaufile("kinkHtBottomTau2V_hist.root");
    TH1F *toptauhist = (TH1F*)toptaufile.Get("lifetime_1");
    TH1F *bottomtauhist = (TH1F*)bottomtaufile.Get("lifetime_1");

//    TFile topemptyingfile("topCell_emptying_DetEff.root");
//    TFile bottomemptyingfile("bottomCell_emptying_DetEff.root");
//    TFile topemptyingfile("kinkHtTopEmptying_DetEff.root");
//    TFile bottomemptyingfile("kinkHtBottomEmptying_DetEff.root");
    TFile topemptyingfile("kinkHtTopEmptying2V_DetEff.root");
    TFile bottomemptyingfile("kinkHtBottomEmptying2V_DetEff.root");

   //// used for emptying files when a study has multiple emptying geometries
//    TFile topemptyingfile(Form("asymFeedEmpt2VTop%icm_DetEff.root",(int)result.parameter));      
//    TFile bottomemptyingfile(Form("asymFeedEmpt2VBottom%icm_DetEff.root",(int)result.parameter));
    TH2 *topemptyinghist = (TH2*)topemptyingfile.Get("emptyEff");
    TH2 *bottomemptyinghist = (TH2*)bottomemptyingfile.Get("emptyEff");

    TH1D *topDetEff = (TH1D*)topemptyingfile.Get("detEff");
    TH1D *bottomDetEff = (TH1D*)bottomemptyingfile.Get("detEff");
///// Top Cell histograms 

    TCanvas Can("Can","TopFill",3000,2000);
    ////////
    scale = UCN_production*beamHeating * polarization * activeTime /result.totalSimulated; // ratio of production rates actual:simulated,  activeTime = simulation active time for source
    ////////
    result.topCellESpec->Scale(scale);
    result.topCellESpec->Draw("HIST");
    result.topCellESpec->SetTitle(Form("%s Top cell filling", fileName.c_str())); 
    result.topCellESpec->GetXaxis()->SetTitle("Energy (neV)");
    result.topCellESpec->GetYaxis()->SetTitle("UCN counts");
    result.topMeanFill = result.topCellESpec->GetMean();
    result.uTopMeanFill = result.topCellESpec->GetStdDev();
    result.topCellTotalFilled = result.topCellESpec->Integral();
    Can.Print(Form("%s_TopCellFilling.eps",fileName.c_str()));

    TH1F topFinish;
    TCanvas *can11 = new TCanvas("can11","TopFinish",3000,2000);
    topFinish = (*toptauhist) * (*result.topCellESpec);
    scale = ( result.topCellESpec->GetEntries()/ result.topCellRamsey) / topFinish.GetEntries();  //Scaling surviving UCN to Filled, then normalizing plot 
    topFinish.Scale(scale);
    topFinish.Draw("HIST same");
       topFinish.SetTitle(Form("%s: top cell end of Ramsey", fileName.c_str())); 
       topFinish.GetXaxis()->SetTitle("Energy (neV)");
//      // topFinish.GetYaxis().SetTitle("");
    can11->Update();
    can11->Print(Form("%s_TopCellFinish.eps",fileName.c_str()));
    topFinish.Copy(result.topCellESpecFinish);
    result.topMeanFinish = topFinish.GetMean();
    result.uTopMeanFinish = topFinish.GetStdDev();
    result.topCellTotalFinished = topFinish.Integral();
    delete can11;

    TH1D h1D;
    topFinish.Copy(h1D);
    TH1D topFinishP = h1D;

    TCanvas *can12 = new TCanvas("can12","TopDetected",3000,2000);
    topFinishP.Multiply(topDetEff);
    topFinishP.Draw("Hist");
      topFinishP.SetTitle(Form("%s: collected from top cell", fileName.c_str())); 
      topFinishP.GetXaxis()-> SetTitle("Energy (neV)");
//      // topFinishP.GetYaxis()-> SetTitle("");
    can12->Update();
    can12->Print(Form("%s_TopCellDetect.eps",fileName.c_str()));
    topFinishP.Copy(result.topCellESpecDetected);
    result.topMeanDetected = topFinishP.GetMean();
    result.uTopMeanDetected = topFinishP.GetMean();
    result.topCellTotalDetected = topFinishP.Integral();
    delete can12;

///// Bottom fill histograms

    TCanvas Can1("Can1","BottomFill",3000,2000);
    scale = UCN_production*beamHeating * 0.5 * 200/result.totalSimulated;
    result.bottomCellESpec->Scale(scale);
    result.bottomCellESpec->Draw("HIST");
    result.bottomCellESpec->SetTitle(Form("%s Bottom cell filling", fileName.c_str())); 
    result.bottomCellESpec->GetXaxis()->SetTitle("Energy (neV)");
    result.bottomCellESpec->GetYaxis()->SetTitle("UCN counts");
    result.bottomMeanFill = result.bottomCellESpec->GetMean();
    result.uBottomMeanFill = result.bottomCellESpec->GetStdDev();
    result.bottomCellTotalFilled = result.bottomCellESpec->Integral();
    Can1.Print(Form("%s_BottomCellFilling.eps",fileName.c_str()));
   
    
    TCanvas *can21 = new TCanvas("can21","BottomFinish",3000,2000);
    TH1F bottomFinish;
    bottomFinish = (*bottomtauhist) * (*result.bottomCellESpec);
    scale = ( result.bottomCellESpec->GetEntries()/ result.bottomCellRamsey) / bottomFinish.GetEntries(); //Scaling surviving UCN to Filled, then normalizing plot
    bottomFinish.Scale(scale);
    bottomFinish.Draw("HIST same");
       bottomFinish.SetTitle(Form("%s: bottom cell end of Ramsey", fileName.c_str())); 
       bottomFinish.GetXaxis()->SetTitle("Energy (neV)");
//      // bottomFinishGetYaxis().SetTitle("");
    can21->Update();
    can21->Print(Form("%s_BottomCellFinish.eps",fileName.c_str()));
    bottomFinish.Copy(result.bottomCellESpecFinish);
    result.bottomMeanFinish = bottomFinish.GetMean();
    result.uBottomMeanFinish = bottomFinish.GetStdDev();
    result.bottomCellTotalFinished = bottomFinish.Integral();
    delete can21;

    TH1D h2D;
    bottomFinish.Copy(h2D);
    TH1D bottomFinishP = h2D;
    
    TCanvas *can22 = new TCanvas("can12","BottomDetected",3000,2000);
    bottomFinishP.Multiply(bottomDetEff);
    bottomFinishP.Draw("Hist");
       bottomFinishP.SetTitle(Form("%s: collected from bottom cell", fileName.c_str())); 
       bottomFinishP.GetXaxis()-> SetTitle("Energy (neV)");
      // bottomFinishP.GetYaxis()-> SetTitle("");
    can22->Update();
    can22->Print(Form("%s_BottomCellDetect.eps",fileName.c_str()));
    bottomFinishP.Copy(result.bottomCellESpecDetected);
    result.bottomMeanDetected =  bottomFinishP.GetMean();
    result.uBottomMeanDetected = bottomFinishP.GetStdDev();
    result.bottomCellTotalDetected = bottomFinishP.Integral();
    delete can22;

    result.totalCollected = result.bottomCellTotalDetected + result.topCellTotalDetected;
    toptaufile.Close();
    bottomtaufile.Close();
    topemptyingfile.Close();
    bottomemptyingfile.Close();
}



