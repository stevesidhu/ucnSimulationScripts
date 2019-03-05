#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include "TGraphErrors.h"
#include "TH1.h"
#include <algorithm>
#include "SimResult.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include <numeric>
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include <fstream>
#include "Rtypes.h"
#include "TLegend.h"

    double beamHeating = 10; //W
    double UCN_production = 2.62E6; //UCN/W
    double cellVolume = 16286 *2; // singular cell volume in cm^3

//Experimental constants used in simulation and analysis
    double emptyingTime = 55;  //cell emptying time !DEFAULT VALUE!
    double t_edm = 132;  // Ramsey storage time in edm cell: NOT CELL STORAGE LIFETIME!   !DEFAULT VALUE!
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
	double activeTime = 200; //Source active time, usually same as simulation time
    
    

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

//Calculates extraction efficiency, storage lifetime, filling time for different operational modes and creates filling energy spectrum histograms
void transport(std::string FileName, int modeParameter, SimResult &results, double cellCenter);

//Returns number of days to reach
double days(int mode, SimResult &result);

//SimulationResults: Function that analyzes each individual .root file in a study.  
//
//After analyzing transport and Days to reach, it optimizes days to reach for t_optimalEDM, t_optimalEmptying,
void analyzeSim(SimResult &results, std::string fileName, int mode, double cellCenter = 0);

//Total number in cell
std::vector<double> cellTotal(int mode, SimResult &results);

//Numerical error calculator
double numError(int mode, SimResult tempResult);

//Stolen function: used to better fit the ideal filling time.  This fixes double peaks and takes global maximum.
double LargerProbability(int value, int ref);
std::array<double, 3> GetMaxBinRange(const TH1* hist, const double confidence_interval = 0.6827);

//Estimates UCN surviving the Ramsey cycle for top and bottom cell individually
std::vector<double> UCNremainingAfterCycle(TH1* topSpectrum, TH1 *bottomSpectrum, double parameter);

//Draw different histograms for cell finish and emptying
void DrawCellHistogram(SimResult &result, std::string fileName);




//Main function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void testSimulationAnalysis()
{
    
    gROOT->SetBatch(1);
    int mode = 1;   //Options are 0:batch, 1: steady beam, 2, steady state.

//////////////////    
// Create a study    
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Create an instance of different studies, simulation files grouped together with varying parameters.
    //Tester study
    Study test_study1 = {"Tester",{45},{"asymFunKinkHt45cm.root"},{},{0}};

    Study asymFeeder_study{"AsymmetricFeeders2v",{0.25,4,8,12},{"asymFeed2V0cm.root","asymFeed2V4cm.root","asymFeed2V8cm.root","asymFeed2V12cm.root"},{},{0.18,0.18,0.18,0.18}};


///////////////////////////////////////////////////////////////////////////////////////////////////////////    
    
   //Tester function 
   // std::vector<Study> allStudies = {test_study1};    

   //Specific studies
   std::vector<Study> allStudies = {asymFeeder_study};


    //Run analysis functions on all simulations in each study and plot graphs
    for(Study study : allStudies){
        SimResult tempResult;
        for(int i = 0; i < study.filenames.size(); ++i) {
            tempResult.parameter = study.parameters[i];
                analyzeSim(tempResult, study.filenames[i], mode, study.cellCenter[i]);
            }    
            tempResult.parameter = study.parameters[i];
            study.results.push_back(tempResult);
            }

   //Plots of whole study vs parameters.
        //Days vs P
        makeGraph(study,"Parameter","DaysToReach", &SimResult::parameter, &SimResult::daysToReach, NULL, &SimResult::errDaysToReach);

        //Tau vs Days
        makeGraph(study,"Tau_source","DaysToReach", &SimResult::sourceStorageLifetime, &SimResult::daysToReach, &SimResult::sourceStorageLifetimeError, &SimResult::errDaysToReach);

        //Tau vs P
        makeGraph(study,"Parameter","Tau_source", &SimResult::parameter, &SimResult::sourceStorageLifetime, NULL, &SimResult::sourceStorageLifetimeError);

        //Days vs Extraction
        makeGraph(study,"ExtractionEfficiency","DaysToReach", &SimResult::efficiency, &SimResult::daysToReach, &SimResult::efficiencyError, &SimResult::errDaysToReach);

        //Total Number in Cell
        makeGraph(study,"Parameter","UCNinCell", &SimResult::parameter, &SimResult::numInCellReal, NULL, &SimResult::uNumInCellReal);
        
        //System storage lifetime 
        makeGraph(study,"Parameter","SystemTau", &SimResult::parameter, &SimResult::systemTau, NULL, &SimResult::uSystemTau);
        
        //Transport efficiency
        makeGraph(study,"Parameter","Transport_efficiency", &SimResult::parameter, &SimResult::transportEff, NULL, &SimResult::uTransportEff);
       //Avg Tau cells 
        makeGraph(study,"Parameter","AvgTauCells", &SimResult::parameter, &SimResult::avgTauCells, NULL, &SimResult::uavgTauCells);

       //Avg Tau cells 
        makeGraph(study,"Parameter","WeightedDetEff", &SimResult::parameter, &SimResult::WDetEff, NULL, &SimResult::uWDetEff);
       
        //Avg Optimal t_edm 
        makeGraph(study,"Parameter","Optimal_Tedm", &SimResult::parameter, &SimResult::optimalTedm, NULL, NULL);

       //Avg Tau cells 
        makeGraph(study,"Parameter","Optimal_EmptyingTime", &SimResult::parameter, &SimResult::optimalEmptyingTime, NULL, NULL);

       //Total collected 
        makeGraph(study,"Parameter","TotalCollectedUCN", &SimResult::parameter, &SimResult::totalCollected , NULL, NULL);


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

    for(int i = 0; i < study.filenames.size(); ++i){
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
    for(int i = 0; i < study.filenames.size(); ++i){
        can->cd(2*i +1);
//        scale = 1E7 / study.results[i].topCellESpecFinish.GetEntries();
        study.results[i].topCellESpecFinish.Scale(scaleTop);
        study.results[i].topCellESpecFinish.SetLineColor(colors[i]);
        study.results[i].topCellESpecFinish.GetYaxis()->SetRangeUser(0, 3E8);
        study.results[i].topCellESpecFinish.GetXaxis()->SetTitle("Energy (neV)");
        study.results[i].topCellESpecFinish.GetYaxis()->SetTitle("UCN surviving");
        study.results[i].topCellESpecFinish.Draw("HIST same");  
        // leg2->AddEntry(&study.results[i].topCellESpecFinish,Form("%s %i",study.parameter_name.c_str(), (int)study.parameters[i]));
    
        can->cd(2*i +2);
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
    for(int i = 0; i < study.filenames.size(); ++i){
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

///////////////////////

        std::ofstream resultsFile;
        resultsFile.open(Form("%s_results.txt",study.parameter_name.c_str()));
        resultsFile <<"file, parameter,totalSimulated,totalProducedReal, fillTime, fillTimeError, efficiency, efficiencyError, sourceStorageLifetime, sourceStorageLifetimeError, systemTau, uSystemTau, totalUCNinSim, uTotalUCNinSim, numInCellSim, uNumInCellSim, numInCellReal, uNumInCellReal, UCNdensityInCells, transportEff, uTransportEff, avgTauCells, uavgTauCells, weightedDetEff, uWeightedDetEff, topCellRamseySurvivorsSim, bottomCellRamseySurvivorsSim, topCellFilledReal, topMeanEFill, uTopMeanEFill, bottomCellFilledReal, bottomMeanEFill, uBottomMeanEFill, topCellFinishedReal, topMeanEFinish, uTopMeanEFinish, bottomCellFinishedReal, bottomMeanEFinish, uBottomMeanEFinish, topCellDetectedReal, topMeanEDetected, uTopMeanEDetected, bottomCellDetectedReal, bottomMeanEDetected, uBottomMeanEDetected, totalDetected, perFilledFromProduced, perDetectedFromProduced, perDetectedFromFilled, perTopCellDetectedFromProduced, perBottomCellDetectedFromProduced, optimalT_edm, optimalEmptyingTime, daysToReach, errDaysToReach \n";
        
        for(int i = 0; i < study.filenames.size(); ++i) {
            resultsFile <<study.filenames[i] << ", "<< study.results[i].parameter <<", "<<study.results[i].totalSimulated<<", " << study.results[i].totalProducedReal <<", "<< study.results[i].fillTime << ", " <<study.results[i].fillTimeError << ", " <<study.results[i].efficiency << ", " << study.results[i].efficiencyError << ", " << study.results[i].sourceStorageLifetime << ", " << study.results[i].sourceStorageLifetimeError << ", " << study.results[i].systemTau << ", " << study.results[i].uSystemTau << ", " << study.results[i].totalUCN << ", " << study.results[i].uTotalUCN << ", "<<study.results[i].numInCellSim << ", "<<study.results[i].uNumInCellSim<<", " << study.results[i].numInCellReal << ", " << study.results[i].uNumInCellReal << ", " << study.results[i].numInCellReal/cellVolume << ", "<< study.results[i].transportEff << ", " << study.results[i].uTransportEff << ", " << study.results[i].avgTauCells << ", " << study.results[i].uavgTauCells << ", "<<study.results[i].WDetEff<< ", " << study.results[i].uWDetEff<< ", " <<study.results[i].topCellRamsey << ", " <<study.results[i].bottomCellRamsey << ", "  <<study.results[i].topCellTotalFilled << ", "  <<study.results[i].topMeanFill<< ", " << study.results[i].uTopMeanFill << ", " <<study.results[i].bottomCellTotalFilled << ", " <<study.results[i].bottomMeanFill << ", "<<study.results[i].uBottomMeanFill << ", " <<study.results[i].topCellTotalFinished << ", "<<study.results[i].topMeanFinish << ", "<<study.results[i].uTopMeanFinish << ", "<<study.results[i].bottomCellTotalFinished << ", "<<study.results[i].bottomMeanFinish << ", "<<study.results[i].uBottomMeanFinish << ", "<<study.results[i].topCellTotalDetected << ", "<<study.results[i].topMeanDetected << ", "<<study.results[i].uTopMeanDetected << ", " <<study.results[i].bottomCellTotalDetected << ", "<<study.results[i].bottomMeanDetected << ", "<<study.results[i].uBottomMeanDetected << ", "<< study.results[i].topCellTotalDetected + study.results[i].bottomCellTotalDetected <<", "<<study.results[i].numInCellReal/ study.results[i].totalProducedReal << ", "<< (study.results[i].topCellTotalDetected + study.results[i].bottomCellTotalDetected)/study.results[i].totalProducedReal <<", "<< (study.results[i].topCellTotalDetected + study.results[i].bottomCellTotalDetected)/study.results[i].numInCellReal<<", "<< study.results[i].topCellTotalDetected/study.results[i].totalProducedReal <<", "<< study.results[i].bottomCellTotalDetected/study.results[i].totalProducedReal<<", " << study.results[i].optimalTedm << ", " << study.results[i].optimalEmptyingTime << ", "<< study.results[i].daysToReach << ", " << study.results[i].errDaysToReach << "\n";
        }
        resultsFile.close();      
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/////
//Fucntions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void transport(std::string fileName, int modeParameter,SimResult &results, double cellCenter)
{

    int mode = modeParameter; //Operational mode.
	
	//Tcuts for cells for simulations with EDM Cell fill geometry
	TCut topCell = TCut(Form("solidend>274 &&solidend<276 && zend>%f", cellCenter));
	TCut bottomCell = TCut(Form("solidend>274 &&solidend<276 && zend<%f", cellCenter));
	TCut cell = (topCell || bottomCell);	
        

	//Open file//
	TFile f1(fileName.c_str());
        TTree *neutronend =(TTree*) f1.Get("neutronend");
        TTree *neutronsnapshot =(TTree*) f1.Get("neutronsnapshot");

    results.totalSimulated = neutronend->GetEntries();    
    std::string fileNameString = fileName.substr(0, fileName.size()-5);
	
	//Names used in plot names//
	std::string Name[3];
        Name[0]= "Batch";
        Name[1]= "Steady Beam";
        Name[2]= "Steady State";

	double fillTime;
	double fillTimeError;
	double sourceStorageLifetime;	
	double sourceStorageLifetimeError;
	double totalUCN;
	double uTotalUCN;
	double numInCell;
	double uNumInCell;
	double efficiency;
	double efficiencyError;
	double transportEff;
	double uTransportEff;
	
	//Creating blank instance of various TCuts
	TCut fillCondition;
	string fillEvolutionFileName;
	TCut tendCut;
	TCut valveOpenTimeCut = TCut(Form("tend==%.0f", valveOpenTime));
        string filledCellFileName;
	string SourceLifetimeFileName;
	TCut modeCuts;


	//Step 1 - Filling Time//
	//initializing the fillCondition TCut based on operational mode.
	//This only considers UCN created in the source relavent to that operational mode
	if(mode == 0) {
		fillCondition = TCut(Form("tstart < %f", valveOpenTime));   
	}
	else if(mode == 1) {
		fillCondition = TCut(Form("tstart < %f", valveOpenTime)); //no condition
	}
	else if(mode == 2) {
		fillCondition = TCut(Form("tstart > %f", valveOpenTime));
	}	

	//Plot snapshot of UCN in cell
	TCanvas *f = new TCanvas("f", "f", 800, 600);
	// creates a tend histogram of ~ 100 bins, t = 100 to t = 200 
	neutronsnapshot->Draw(Form("tend >> filling(%f, %f, %f)", activeTime-valveOpenTime, valveOpenTime, activeTime), cell && fillCondition); //use form for valve open time
	TH1S *filling = (TH1S*) gDirectory->Get("filling");
    filling->Scale(1E6/results.totalSimulated);
	filling->Draw();
   
    //fit the ideal filling time
    auto range = GetMaxBinRange(filling, 0.9545); // get range of bins that contain maximum with 95% probability
    fillTime = (range[0] + range[2])/2.0;
    fillTimeError = (fabs(fillTime - range[0]) + fabs(fillTime - range[2]))/2.0;
 
    //plot Filling time based on mode.
    filling->SetTitle("EDM Cell Filling");
        filling->GetXaxis()->SetTitle("Time [s]");
        filling->GetYaxis()->SetTitle("Number of UCN in EDM cells");
        filling->GetYaxis()->SetTitleOffset(1.3);
        filling->SetStats(0);
	filling->Draw();
	
        if(mode == 0) {
		f->Print((fileNameString + "_fillEvolution_batch.pdf").c_str());
        }
        else if (mode == 1) {
		f->Print((fileNameString + "_fillEvolution_steadyBeam.pdf").c_str());
        }
        else if (mode == 2) {
		f->Print((fileNameString + "_fillEvolution_steadyState.pdf").c_str());
        }
	delete f;
	if (mode == 2) {
		fillTime = activeTime - 1;
	}
	// initialize TCut for the end of filling time ie. Closing EDM cell valves
	tendCut = TCut(Form("tend == %.0f", fillTime)); // TCut for end of fillTime

    
    
    //Step 2 - Cell Energy Spectrum//		
	double zeroPoint = cellCenter * 1.025E-7;  // used to shift histograms if cell center is not z=0 in simulations
    std::cout<< "zeropoint energy is : "<<zeroPoint << std::endl;

    //Histogram draw TString for top, bottom, both cells with energy shifts
	TString drawBoth = TString::Format("Hend - %g >>cellESpectrum(50,0,250E-9)",zeroPoint);
	TString drawTop = TString::Format("Hend - %g >>topcellESpectrum(50,0,250E-9)",zeroPoint);
	TString drawBottom = TString::Format("Hend -%g >>bottomcellESpectrum(50,0,250E-9)",zeroPoint); // shift bottom spectrum to the right a bit. 

///////// Drawing energy spectrum in both cells combined.
	TCanvas *d = new TCanvas("d", "d", 800, 600);
        neutronsnapshot->Draw(drawBoth,"Estart<233.5e-9" && cell && tendCut);  // Draw histogram based on ideal filling time.
	    TH1 *cellESpectrum = (TH1*) gDirectory->Get("cellESpectrum");
	    cellESpectrum->SetTitle("Energy distribution of UCN at the end of filling");
	    cellESpectrum->Scale(1E6 / results.totalSimulated);
	    cellESpectrum->GetXaxis()->SetTitle("UCN energy (eV)");
	    cellESpectrum->GetYaxis()->SetTitle("Number of UCN");
	    cellESpectrum->GetYaxis()->SetTitleOffset(1.3);
	        cellESpectrum->ResetStats();
	    //cut off ends of the spectrum for a cleaner histogram
	    double quantiles[2];
	    double cutOff[2] = {0.05, 0.95};
	    cellESpectrum->GetQuantiles(2, quantiles, cutOff);
        cellESpectrum->Draw("HIST");
        double lowerBound = quantiles[0];
        double upperBound = quantiles[1];

     ///Drawing espec   
	if(mode == 0) {
        	d->Print((fileNameString + "_EDCellFilled_batch.eps").c_str());
        }
	else if(mode == 1) {
		d->Print((fileNameString + "_EDCellFilled_steadyBeam.eps").c_str());
	}
	if(mode == 2) {
		d->Print((fileNameString + "_EDCellFilled_steadyState.eps").c_str());
        }
	delete d;

    //Drawing filling energy spectrum for top and bottom cells seperately.  This are used in the rest of analysis
    double scale = 1E6 / neutronend->GetEntries();
	d = new TCanvas("d", "d", 800, 600);
	neutronsnapshot->Draw(drawTop,"Estart<233.5e-9"&& topCell && tendCut);
	TH1F *topcellESpectrum = (TH1F*) gDirectory->Get("topcellESpectrum");
	topcellESpectrum->SetDirectory(0);
	neutronsnapshot->Draw(drawBottom,"Estart<233.5e-9"&& bottomCell && tendCut);
	TH1F *bottomcellESpectrum = (TH1F*) gDirectory->Get("bottomcellESpectrum");
	bottomcellESpectrum->SetDirectory(0);
	topcellESpectrum->Scale(scale);
	bottomcellESpectrum->Scale(scale);
	topcellESpectrum->Draw("HIST");
	bottomcellESpectrum->Draw("HIST");
    results.topCellESpec = topcellESpectrum;
    results.bottomCellESpec = bottomcellESpectrum;
		d->Print((fileNameString + "_testFill.eps").c_str());
    delete d;



	//Step 3 - Source Storage Lifetime//
	//Plot decay rate of UCN in side the source before valve opens.  
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	neutronsnapshot->Draw(Form("tend-tstart>>lifetime(%.0f,5,%.0f)",valveOpenTime-10, valveOpenTime-5), Form("tstart < 5 && Hend > %e && Hend < %e", lowerBound, upperBound));
	TH1D *lifetime = (TH1D*) gDirectory->Get("lifetime");
	lifetime->Scale(1E6 / results.totalSimulated);
	TF1 g2("fit","expo",5,valveOpenTime - 5);  //exponential fit to calculate tau
	lifetime->Fit(&g2,"R");
	        
	Double_t slope = g2.GetParameter(1);
	Double_t slopeError = g2.GetParError(1);
	sourceStorageLifetime = -1/slope;
	sourceStorageLifetimeError = slopeError/(slope*slope);  //confirmed by SS and SM
	       
	lifetime->SetTitle("Storage lifetime of source");
	lifetime->GetXaxis()->SetTitle("UCN Lifetime [s]");
	lifetime->GetYaxis()->SetTitle("Number of UCN");
	lifetime->GetYaxis()->SetTitleOffset(1.3);
	c1->SetLogy(1);
	lifetime->Draw();
	 
	//Plotting source lifetime 
	if(mode == 0) {
		c1->Print((fileNameString + "_sourceLifetime_batch.pdf").c_str());	
	}
	else if (mode == 1) {
		c1->Print((fileNameString + "_sourceLifetime_steadyBeam.pdf").c_str()); 
	}
	else if (mode == 2) {
		c1->Print((fileNameString + "_sourceLifetime_steadyState.pdf").c_str()); 	
	}		
	delete c1;

    // Used for as a TCut.  This gives the start time of UCN creation for steady beam and steady state modes
    // All the UCN that have the ability to make it to the cell in a mode are created in this start time.
    // e^2.4 = .9
	double ninetyPer = valveOpenTime - 2.4*sourceStorageLifetime; // 100 - optimal irradiation time
	
	//mode cuts for source production start time
	if(mode == 0) {
		 modeCuts = TCut(Form("tstart<%.0f && tstart > %f", valveOpenTime, ninetyPer)); //batchMode;
        }
        else if(mode == 1) {
            modeCuts = TCut(Form("tstart>%f && tstart < %f", ninetyPer, fillTime));
        }
        else if(mode == 2) {
                modeCuts = TCut(Form("tstart>%.0f", valveOpenTime)); //steadyStateMode;
        }
	

	//Step 4 - Number of UCN in Cell 
	//Number in cell at the end of ideal filling time
	TCanvas *u = new TCanvas("u", "u", 800, 600);
	neutronsnapshot->Draw("zend:yend:xend>>inCell", cell && modeCuts && tendCut); //use same for batch and steady BEAM
	TH1 *inCell = (TH1*) gDirectory->Get("inCell");
	numInCell = inCell->GetEntries();
	uNumInCell = sqrt(numInCell);
	delete u;


	//Step 5 - Total number of UCN//
	//Total UCN created in the source that could make it in the same time
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	neutronend->Draw("zend:yend:xend >> total", modeCuts);
	TH1 *total = (TH1*) gDirectory->Get("total");
	totalUCN = total->GetEntries();
	uTotalUCN = sqrt(totalUCN);
	delete c2;
	
	//if batch, take snapshot at tend == valveOpenTime, plot with modeCut and tend cut, getEntries, get transport efficiency (numInCell/getEntries)
	
	if(mode == 0 || mode == 1) {
		TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
		neutronsnapshot->Draw("zend:yend:xend>>transport", modeCuts && valveOpenTimeCut,"stopID != -2");
	    TH1 *transport = (TH1*) gDirectory->Get("transport");
		transport->SetTitle("UCNs transported");
		double numTransported = transport->GetEntries();
		double uNumTransported = sqrt(numTransported);
		delete c3;
		transportEff = numInCell/numTransported; // Number in cell at the end of ideal filling time divided by ucn produced within last 2.4* tau_storage of the source.  roughly 90 percent of them would survive up scatter
		uTransportEff = sqrt((uNumInCell/numTransported)*(uNumInCell/numTransported) + (numInCell*uNumTransported/numTransported/numTransported)*(numInCell*uNumTransported/numTransported/numTransported));
		}
    else{
    transportEff =0;
    uTransportEff=0;
    }

	//Step 6 - Calculate Efficiencies//
	//total cell UCN divided by the total created in the source in that mode.
    efficiency = numInCell/totalUCN;
	efficiencyError = sqrt( (uNumInCell/totalUCN)*(uNumInCell/totalUCN) + (numInCell*uTotalUCN/totalUCN/totalUCN)*(numInCell*uTotalUCN/totalUCN/totalUCN) );		
	

	//Step 7 - System Storage Lifetime//
	TCut systemStorageCut = TCut(Form("tstart < %f && tstart > %f",valveOpenTime+5, valveOpenTime));	
	TCanvas *s = new TCanvas("s", "s", 800, 600);
	neutronsnapshot->Draw("tend-tstart>>systemStorage", systemStorageCut); //(Form("tend-tstart>>systemStorage(%f, %f, %f)", activeTime - valveOpenTime + 5, valveOpenTime+5, activeTime), systemStorageCut);
	TH1D *systemStorage =(TH1D*) gDirectory->Get("systemStorage");
        TF1 systemStorFit("fit","expo", 0, activeTime-valveOpenTime-5); //,valveOpenTime+5, activeTime);
        systemStorage->Fit(&systemStorFit,"R");

        Double_t slope2 = systemStorFit.GetParameter(1);
        Double_t slope2Error = systemStorFit.GetParError(1);
        double systemStorageLifetime = -1/slope2;
        double systemStorageLifetimeError = slope2Error/(slope2*slope2);
	delete s;


//Need to change outputs to a return function:
    results.efficiency = efficiency; //number, not %
    results.efficiencyError = efficiencyError;
    results.fillTime = fillTime -100; //proper fill time
    results.sourceStorageLifetime = sourceStorageLifetime;
    results.sourceStorageLifetimeError = sourceStorageLifetimeError;
    results.totalUCN = totalUCN;
    results.uTotalUCN = uTotalUCN;
    results.numInCellSim = numInCell;
    results.uNumInCellSim = uNumInCell;
    results.systemTau = systemStorageLifetime;
    results.uSystemTau = systemStorageLifetimeError;
    results.transportEff = transportEff;
    results.uTransportEff = uTransportEff;
    results.fillTimeError = fillTimeError;
}




//UCN remaining: look at N_0 that fill the cell and estimate how many survive after Ramsey cycle, emptying, collection as a function of energy (5 nev binning)
std::vector<double> UCNremainingAfterCycle(TH1F *topspectrum, TH1F *bottomspectrum, double parameter)
{
    TFile toptaufile("topCell_hist.root");
    TFile bottomtaufile("bottomCell_hist.root");
    TH1F *toptauhist = (TH1F*)toptaufile.Get("lifetime_1");
    TH1F *bottomtauhist = (TH1F*)bottomtaufile.Get("lifetime_1");

    TFile topemptyingfile(Form("asymFeedEmpt2VTop%icm_DetEff.root",(int)parameter));
    TFile bottomemptyingfile(Form("asymFeedEmpt2VBottom%icm_DetEff.root",(int)parameter));
    TH2 *topemptyinghist = (TH2*)topemptyingfile.Get("emptyEff");
    TH2 *bottomemptyinghist = (TH2*)bottomemptyingfile.Get("emptyEff");
    TH1 *topDetEff = topemptyinghist->ProjectionY("_detEff", 0,topemptyinghist->GetXaxis()->FindBin(emptyingTime));
    TH1 *bottomDetEff = bottomemptyinghist->ProjectionY("_detEff",0,bottomemptyinghist->GetXaxis()->FindBin(emptyingTime));


    //Initial ucn in both cells
    double n0Cell = topspectrum->GetEntries() + bottomspectrum->GetEntries();

    double nFCellTop = 0.;
    double nFCellBottom = 0.;
    double survivors = 0.;
    double dsurvivors = 0.;
    double wTopDet = 0.;
    double wBottomDet =0.;

    // looping over each energy bin
    for (int i = 1; i < toptauhist->GetNbinsX() - 1; ++i){
        double Nitop = topspectrum->GetBinContent(i);  //Initial number of UCN at end of filling, in that energy bin
        double dNitop = topspectrum->GetBinError(i);
        double tauitop = toptauhist->GetBinContent(i); //tau storage for that energy bin
        double dtauitop = toptauhist->GetBinError(i);
        double detEffTop = topDetEff->GetBinContent(i);
        double uDetEffTop = topDetEff->GetBinError(i);  //Emptying and collection efficiency of that energy bin
        double Nibottom = bottomspectrum->GetBinContent(i);
        double dNibottom = bottomspectrum->GetBinError(i);
        double tauibottom = bottomtauhist->GetBinContent(i);
        double dtauibottom = bottomtauhist->GetBinError(i);
        double detEffBottom = bottomDetEff->GetBinContent(i);
        double uDetEffBotom = bottomDetEff->GetBinError(i);

        //Number of UCN of that energy bin detected
        wTopDet += detEffTop * Nitop*exp(-t_edm/tauitop); 
        wBottomDet+= detEffBottom * Nibottom*exp(-t_edm/tauibottom);
        
        //Number that survive Ramsey cycle (right before when the emptying valve opens)
        nFCellTop += Nitop*exp(-t_edm/tauitop);
        nFCellBottom += Nibottom*exp(-t_edm/tauibottom);
        
        //Most important
        //Total survivors for each energy bin 
        survivors += detEffTop * Nitop*exp(-t_edm/tauitop) + detEffBottom * Nibottom*exp(-t_edm/tauibottom);
        if (dtauitop != 0)
            dsurvivors += pow(udetEffTop*Nitop*exp(-t_edm/tauitop), 2) + pow(dNitop*detEffTop*exp(-t_edm/tauitop), 2) + pow(Nitop*detEffTop*exp(-t_edm/tauitop)/tauitop/tauitop*dtauitop, 2);
        if (dtauibottom != 0)
            dsurvivors += pow(udetEffBottom* Nibottom *exp(-t_edm/tauibottom), 2) + pow(dNibottom *detEffBottom *exp(-t_edm/tauibottom), 2) + pow(Nibottom*detEffBottom*exp(-t_edm/tauibottom)/tauibottom/tauibottom*dtauibottom, 2);   
    }
   
    //Used to calculed weight def eff
    double NF = nFCellTop + nFCellBottom;
    double N0 = topspectrum->Integral() + bottomspectrum->Integral();
    double avgTauCell = 130.0 / log(N0/NF);
    double uavgTauCell = 130.0 * (sqrt(1./N0 + 1./NF) * NF/N0) ; 
    
    double topDetMean = wTopDet / nFCellTop;
    double bottomDetMean = wBottomDet / nFCellBottom;
    double utopDetMean =  sqrt(1.0/nFCellTop);
    double ubottomDetMean = sqrt(1.0/nFCellBottom);

    double wMeanDetEff = (topDetMean /pow(utopDetMean,2) + bottomDetMean /pow(ubottomDetMean,2))/(1./pow(utopDetMean,2) + 1./pow(ubottomDetMean,2));
    double uwMeanDetEff = sqrt(1/(1./pow(ubottomDetMean,2) + 1./pow(ubottomDetMean,2)));
    ///

    toptaufile.Close();
    bottomtaufile.Close();
    topemptyingfile.Close();
    bottomemptyingfile.Close();
    
    return {survivors/N0, sqrt(dsurvivors)/N0, avgTauCell, uavgTauCell,wMeanDetEff,uwMeanDetEff,nFCellTop,nFCellBottom};
}


//////////////////////////////////////////////////////////////////////////////////


double days(int mode, SimResult &result)
{
    double daysToReach;

    double sourcePumpingTime; // function: t_irradiation in RP spreadsheet. 
    if (mode ==0)
        sourcePumpingTime = result.sourceStorageLifetime * 2.4;
    else if (mode ==1)
        sourcePumpingTime = result.sourceStorageLifetime * 2.4 + result.fillTime;
    else
        sourcePumpingTime = 100;  //need to change this later to fillingTime of mode =2
    
    double cyclesToReach;
    double cyclesPerDay;
    double sensitivityPerFill;
   
    std::vector<double> survprob = UCNremainingAfterCycle(result.topCellESpec, result.bottomCellESpec, result.parameter);
    result.survivalprob  = survprob[0];
    result.dSurvivalprob = survprob[1];
    result.avgTauCells = survprob[2];
    result.uavgTauCells = survprob[3];
    result.WDetEff = survprob[4];
    result.uWDetEff = survprob[5];
    result.topCellRamsey = survprob[6];
    result.bottomCellRamsey = survprob[7];

    //cyclesPerDay
    if (mode ==1)
        cyclesPerDay = (stableField * 60)/((fillsPerCycle * (emptyingTime + t_edm +sourcePumpingTime) + 2*polarityTime + degaussingTime/10)/60);
        
        //beam on at same time as emptying and only switching E field once, degaussing only on off hours and once per day.
        //cyclesPerDay = (stableField * 60)/((fillsPerCycle * (std::max(emptyingTime, sourcePumpingTime) +t_edm) + polarityTime)/60) -1 ;
    else
        cyclesPerDay = (stableField * 60)/((fillsPerCycle * (result.fillTime + emptyingTime + t_edm +sourcePumpingTime) + 2*polarityTime + degaussingTime/10)/60);
   
    
    //cyclesToReach
    
    //sensitivityPerFill
    double N_0;  //Initial number of UCN
    double N_det;  //Number of UCN that get detected
    double alpha_after;  //Visibility after Ramsey sequence

    N_0 = spinTrans * beamHeating * UCN_production * sourcePumpingTime * result.efficiency;

    N_det = N_0 * result.survivalprob * exp(-t_edm/tau_Xe) * detEff;

    alpha_after = alpha * exp(-((t_edm- t_wait - 2*t_pulse) / T2) - (t_wait + 2* t_pulse)/T1) * detSpinTrans * Panalyzer;

    sensitivityPerFill = hbar / (2 * (t_edm - t_wait - 2*t_pulse) * EField * sqrt(N_det) * alpha_after);

    cyclesToReach = 100 * (sensitivityPerFill/sqrt(fillsPerCycle)/(1E-26)) * (sensitivityPerFill/sqrt(fillsPerCycle)/(1E-26));

//    std::cout<< "cycles per day "<< cyclesPerDay << std::endl;
    daysToReach =  cyclesToReach / cyclesPerDay;
    
    return daysToReach;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyzeSim(SimResult &result, std::string fileName, int mode, double cellCenter)
{
         std::cout << fileName <<" , " << mode << std::endl;  
         std::string fileNameString = fileName.substr(0, fileName.size()-5); //Strips .root from file name.
	     std::vector<double> cellNumber;
         
         //function call for transport
	     transport(fileName, mode, result,cellCenter);

        //2019-02-27:  SS and WS just realized that we're optimizing t_EDM and t_emptying and then plugging those values into days function
        /////////////

         TF1 optimizeddays_tEDM("optdays_tEDM",[mode, &result](double *x, double *p){
                 double old_t_edm = t_edm;
                 t_edm = x[0];
                 double d = days(mode, result);
                 t_edm = old_t_edm;
                 return d;
            }, 20, 500, 0);

	     TCanvas *u = new TCanvas("u", "u", 800, 600);
         optimizeddays_tEDM.Draw();
         u->Print((fileNameString + "_optimized_day_tEDM.eps").c_str());
         t_edm = optimizeddays_tEDM.GetMinimumX(20,500);
         std::cout << "Optimal t_edm: " <<  t_edm << "s\n";
         delete u;
         
         result.optimalTedm = t_edm;
       
         TF1 optimizeddays_empty("optdays_empty",[mode, &result](double *x, double *p){
                 double old_emptyingTime = t_edm;
                 emptyingTime = x[0];
                 double d = days(mode, result);
                 emptyingTime = old_emptyingTime;
                 return d;
            }, 10, 200, 0);

	     TCanvas *v = new TCanvas("v", "v", 800, 600);
         optimizeddays_empty.Draw();
         v->Print((fileNameString + "_optimized_day_empty.eps").c_str());
         emptyingTime = optimizeddays_empty.GetMinimumX(10,200);
         std::cout << "Optimal emptyingTime: " <<  emptyingTime << "s\n";
         delete v;
         result.optimalEmptyingTime = emptyingTime;
         
         double daysToReach = days(mode, result);
         result.daysToReach = daysToReach;

         std::vector<double> pDays = {result.efficiency, result.fillTime, result.sourceStorageLifetime,result.survivalprob};
         std::vector<double> errPDays = {result.efficiencyError, result.fillTimeError, result.sourceStorageLifetimeError, result.dSurvivalprob};
         double errDaysToReach = numError(mode, result);
         result.errDaysToReach = errDaysToReach;
         cellNumber = cellTotal(mode, result);
    	 result.numInCellReal = cellNumber[0];
    	 result.uNumInCellReal = cellNumber[1];
     
         DrawCellHistogram(result, fileNameString);
}




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

double numError(int mode, SimResult tempResult)
{
    SimResult result = tempResult;
    
    double error = 0;
    double errorPar[3] ={result.efficiencyError, result.fillTimeError, result.sourceStorageLifetimeError}; 
    double x = days(mode, result); 
    
    TF1 f("f",[mode,&result](double *x, double *p){
                result.efficiency = p[0];
                result.fillTime = p[1];
                result.sourceStorageLifetime = p[2];
                return days(mode, result);
            },0,10000,3);
    f.SetParameters(result.efficiency, result.fillTime, result.sourceStorageLifetime);
    f.SetParError(0,result.efficiencyError);
    f.SetParError(1,result.fillTimeError);
    f.SetParError(2,result.sourceStorageLifetimeError);

    for (int i = 0; i < 3; ++i)
    { 
        error = error + pow(errorPar[i]*f.GradientPar(i, &x), 2);
    }
    return error = sqrt(error);
    
}

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


void DrawCellHistogram(SimResult &result, std::string fileName)
{

    double scale;
    TFile toptaufile("topCell_hist.root");
    TFile bottomtaufile("bottomCell_hist.root");
    TH1F *toptauhist = (TH1F*)toptaufile.Get("lifetime_1");
    TH1F *bottomtauhist = (TH1F*)bottomtaufile.Get("lifetime_1");

    std::cout<< Form("asymFeedEmptBottom%icm_DetEff.root",(int)result.parameter) << std::endl;

    TFile topemptyingfile(Form("asymFeedEmpt2VTop%icm_DetEff.root",(int)result.parameter));
    TFile bottomemptyingfile(Form("asymFeedEmpt2VBottom%icm_DetEff.root",(int)result.parameter));
    TH2 *topemptyinghist = (TH2*)topemptyingfile.Get("emptyEff");
    TH2 *bottomemptyinghist = (TH2*)bottomemptyingfile.Get("emptyEff");

    TH1D *topDetEff = (TH1D*)topemptyingfile.Get("detEff");
    TH1D *bottomDetEff = (TH1D*)bottomemptyingfile.Get("detEff");
///// Top Cell histograms 

    TCanvas Can("Can","TopFill",3000,2000);
    scale = UCN_production*beamHeating * 0.5 * 200/result.totalSimulated; // ratio of production rates actual:simulated
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



