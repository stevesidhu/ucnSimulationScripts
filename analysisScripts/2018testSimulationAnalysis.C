#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include "TGraphErrors.h"
#include "TH1.h"
#include <algorithm>


//Experimental constants used in simulation and analysis
    double emptyingTime = 30;  //cell emptying time
    double t_edm = 132;  // Ramsey storage time in edm cell: NOT CELL STORAGE LIFETIME!
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
    double beamHeating = 10; //W
    double UCN_production = 2.62E6; //UCN/W
//    double transToDet = 0.9; //transport efficiency to detector (collection efficiency)
    double detEff = 0.9; //Detector efficiency
    double alpha = 0.95; //Initial polarization (alpha_0)
    double T2 = 500; //Transverse relaxation time /
    double t_wait = 2; //Time before Pi/2 pulse
    double t_pulse = 2; //duration of Pi/2 pulse
    double T1 = 1000; //longitudinal relaxation time
    double detSpinTrans = 0.95; //Polarization loss during emptying
    double Panalyzer = 0.9; //analyzing power of analyzer

	double cellUpperZBound = 0.2; //Highest z coordinate within the UPPER EDM cell, relative to the origin, in m
	double cellLowerZBound = -0.2; //Lowest z coordinate within the LOWER EDM cell, relative to the origin, in m
	Double_t cellLowerXBound = 5.3; //Minimum x coordinate corresponding to the EDM cell, relative to the origin, in m. This is the distance from the origin to the start of the cell
	double valveOpenTime = 100;
	double activeTime = 200; //Source active time, usually same as simulation time
    
    

struct SimResult{
    double parameter;
    double fillTime;
    double fillTimeError;
    double efficiency;
    double efficiencyError;
    double sourceStorageLifetime;
    double sourceStorageLifetimeError;
    double systemTau;
    double uSystemTau;
    double totalUCN;
    double uTotalUCN;
    double numInCell;
    double uNumInCell;
    double transportEff;
    double uTransportEff;
    double survivalprob;
    double dSurvivalprob;
    double daysToReach;
    double errDaysToReach;
    double avgTauCells;
    double uavgTauCells;
    TH1 *topCellESpec;
    TH1 *bottomCellESpec;
    double WDetEff;
    double uWDetEff;
    double optimalTedm;
    double optimalEmptyingTime;
};

struct Study{
    std::string parameter_name;
    std::vector<double> parameters;
    std::vector<std::string> filenames; 
    std::vector<SimResult> results;
};

typedef double (SimResult::*ResultMember);



//Function declariations
/////////////////////////
// Plots a plot.
void makeGraph(Study &study,std::string xName,std::string yName ,ResultMember mem_xs, ResultMember mem_ys, ResultMember mem_xerr, ResultMember mem_yerr);

//Calculates extraction efficiency, storage lifetime, filling time for different operational modes
void transport(std::string FileName, int modeParameter, SimResult &results);

//Returns number of days to reach
double days(int mode, SimResult &result);

//SimulationResults
void analyzeSim(SimResult &results, std::string fileName, int mode);

//Total number in cell
std::vector<double> cellTotal(int mode, SimResult results);

//Numerical error calculator
double numError(int mode, SimResult tempResult);

//
double LargerProbability(int value, int ref);
//
std::array<double, 3> GetMaxBinRange(const TH1* hist, const double confidence_interval = 0.6827);

std::vector<double> UCNremainingAfterCycle(TH1* topSpectrum, TH1 *bottomSpectrum);



//Main function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void testSimulationAnalysis()
{
    
    gROOT->SetBatch(1);
    int mode = 1;   
    
    TFile outfile("AllResults.root","RECREATE");
    TTree studyTree("Study","A study");


// Create a study    
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Create an instance of different studies, simulation files grouped together with varying parameters.
    //Tester study
    Study test_study = {"Tester",{1},{"kinkHeight35.root"},{}};

    //V3 cryostat; standard (0.05 mm thick foil; solidend ==3 )
    Study extraction_study = {"Extraction_height",{15,7.5,7.5,0.5},{"standard.root","double_kink.root","middle_extraction.root","low_extraction.root"},{}};
    
    //V3 cryostat; standard (0.05 mm thick foil; solidend == 3)
    Study tube_study ={"LHe_guide_diameter",{10,12.5,15,18},{"10cmTube.root","12.5cmTube.root","standard.root","18cmTube.root"},{}} ;
   
    //V3 cryostat; standard (0.05 mm thick foil; solidend == 3)
    Study HEX_study ={"HEX_diameter",{12.5,15,18,20,12.5},{"HEX12cm_75.root","standard.root","HEX18cm.root","HEX20.root","HEX_12cm.root"},{}};
    
    //V4 cryostat; cedarV4standard (0.1 mm thick foil; solidend ==370) 
    Study v4funnel_study = {"FunnelShape",{12,12,15,15,18,18,45},{"funnel_bottom_12deg.root","funnel_top_12deg.root","funnel_bottom_15deg.root","funnel_top_15deg.root","funnel_bottom_18deg.root","cedarV4standard.root","funnel_bottom_45deg.root"},{}};
   
    //V4 cryostat; cedarV4standard (0.1 mm thick foil; solidend == 370)
    Study kinkDiameter_study = {"KinkDiameter",{8.5,10,10,8.5,10,10,15,10},{"funnel_bottom_18deg.root","kinkDiameter10.root","kinkDiameterOnward_10cm.root","cedarV4standard.root","afterKink100_85_100.root","kink100_100_85_85.root","guideDiameter15.root","shortBore_8_onward.root"},{}}; 
    
    //V4 cryostat; (0.1 mm thick foil; solidend == 370)
    Study kinkHeight_study{"oldKinkHeight",{50,75,35,100,60,90,82.5,45,65,70,80},{"kinkHeight50.root", "kinkHeight75.root","kinkHeight35.root","kinkHeight100.root","kinkHeight60.root","kinkHeight90.root","kinkHeight82.root","kinkHeight45.root","kinkHeight65.root","kinkHeight70.root","kinkHeight80.root"},{}};
    
    //kink height study: 10cm overfill and 50cm long HEX
    Study newkinkHeight_study{"KinkHeight",{35,45,55,65,70,75,80,85,95},{"newkinkHeight35.root", "newkinkHeight45.root","newkinkHeight55.root","newkinkHeight65.root","newkinkHeight70.root","newkinkHeight75.root","newkinkHeight80.root","newkinkHeight85.root","newkinkHeight95.root"},{}};
    
    //no bfield but foil vs kink height: 10cm overfill and 50cm long HEX
    Study noBkinkHeight_study{"NoBKinkHeight",{35,45,55,65,70,75,80,85,95},{"noBkinkHeight35.root", "noBkinkHeight45.root","noBkinkHeight55.root","noBkinkHeight65.root","noBkinkHeight70.root","noBkinkHeight75.root","noBkinkHeight80.root","noBkinkHeight85.root","noBkinkHeight95.root"},{}};
   
    // no bfield or foil vs kink height:10cm overfill and 50cm long HEX
    Study noBFkinkHeight_study{"NoBFKinkHeight",{35,45,55,65,70,75,80,85,95},{"noBFkinkHeight35.root", "noBFkinkHeight45.root","noBFkinkHeight55.root","noBFkinkHeight65.root","noBFkinkHeight70.root","noBFkinkHeight75.root","noBFkinkHeight80.root","noBFkinkHeight85.root","noBFkinkHeight95.root"},{}};

    //double bfield coild diamater vs kink height: 10cm overfill and 50cm long HEX
    Study DoubleBKinkHeight_study{"2BKinkHeight",{35,45,55,65,70,75,85,95},{"2BkinkHeight35.root", "2BkinkHeight45.root","2BkinkHeight55.root","2BkinkHeight65.root","2BkinkHeight70.root","2BkinkHeight75.root","2BkinkHeight85.root","2BkinkHeight95.root"},{}};
    
    //double bfield coil vs bore diameter: 10cm overfill and 50cm long HEX
    Study DoubleBBoreDiameter_study{"2BBoreDiameter",{60,70,75,85,100,125,150},{"2BboreDiameter60.root","2BboreDiameter70.root","2BboreDiameter75.root","2BboreDiameter85.root","2BboreDiameter100.root","2BboreDiameter125.root","2BboreDiameter150.root"},{}};
    
    
    //V3 cryostat; standard_hipriorFoil (0.1 mm thick foil; solidend == 91)
    Study afterKink_study{"AfterKink",{10,15,18,8.5},{"afterKink10.root","afterKink15.root","afterKink18.root","standard_hipriorFoil.root"},{}};

    //V4 cryostat; 1 config.in (0.1 mm thick foil; solidend ==370) 
    Study foilLD_study{"Foil_LD",{20,5,0.5},{"kinkDiameterOnward_10cm.root","kinkDiameterOnward_10_foilLD5.root","kinkDiameterOnward_10_foilLD0.root"},{}};
   

    //V4 cryostat; (0.1 mm thick foil)
    Study boreDiameter_study = {"boreDiameter",{8.5,10,8.5,7},{"highRes_kinkHeight35.root","bore10_scaledB.root","bore85.root","bore70.root"},{}}; 

    //V5 Cryostat;
    Study hardBoreDiameter_study = {"hardBoreDiameter",{7,8.5,10,12.5},{"70HardBore.root","85HardBore.root","100HardBore.root","125HardBore.root"},{}};

    //V5 Cryostat;
    Study coil100HardBoreDiameter_study = {"100coilHardBoreDiameter",{7,8.5,10,12.5},{"70HardBore_100coil.root","85HardBore_100coil.root","100HardBore_100coil.root","125HardBore_100coil.root"},{}};
    
    //V5 cryostat
    Study hardVertical_study = {"HardVerticalStudy",{212,347},{"VerticalNiPBore.root","VerticalHardBore.root"},{}};

    //v5 cryostat
    Study brokenHeight_study = {"BrokenHeightStudy",{51,75},{"bottle51to75.root","bottle75to75.root"},{}};
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////    
    
   //Tester function 
    std::vector<Study> allStudies = {test_study};    

   //All studies
   // std::vector<Study> allStudies= {extraction_study, tube_study, HEX_study, v4funnel_study, kinkDiameter_study, kinkHeight_study, afterKink_study, boreDiameter_study,hardBoreDiameter_study};
   
   //Specific studies
   // std::vector<Study> allStudies = {newkinkHeight_study};
   // std::vector<Study> allStudies = {DoubleBBoreDiameter_study,DoubleBKinkHeight_study};


    //Run analysis functions on all simulations in each study and plot graphs
    for(Study study : allStudies){
        for(int i = 0; i < study.filenames.size(); ++i) {
            TTree studyTree(study.parameter_name.c_str(),"RECREATE");
            SimResult tempResult;
            analyzeSim(tempResult, study.filenames[i], mode);
            tempResult.parameter = study.parameters[i];
          //  study.results.push_back(tempResult);
           // studyTree.Branch("DaysToReach",&SimResult::daysToReach,"DaysToReach/D");
            }

        //Days vs P
         makeGraph(study,"Parameter","DaysToReach", &SimResult::parameter, &SimResult::daysToReach, NULL, &SimResult::errDaysToReach);
             
        //Tau vs Days
        makeGraph(study,"Tau_source","DaysToReach", &SimResult::sourceStorageLifetime, &SimResult::daysToReach, &SimResult::sourceStorageLifetimeError, &SimResult::errDaysToReach);

        //Tau vs P
        makeGraph(study,"Parameter","Tau_source", &SimResult::parameter, &SimResult::sourceStorageLifetime, NULL, &SimResult::sourceStorageLifetimeError);

        //Days vs Extraction
        makeGraph(study,"ExtractionEfficiency","DaysToReach", &SimResult::efficiency, &SimResult::daysToReach, &SimResult::efficiencyError, &SimResult::errDaysToReach);

        //Total Number in Cell
        makeGraph(study,"Parameter","UCNinCell", &SimResult::parameter, &SimResult::numInCell, NULL, &SimResult::uNumInCell);
        
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
       
        
        
        std::ofstream resultsFile;
        resultsFile.open(Form("%s_results.txt",study.parameter_name.c_str()));
        resultsFile <<"file, parameter, fillTime, fillTimeError, efficiency, efficiencyError, sourceStorageLifetime, sourceStorageLifetimeError, systemTau, uSystemTau, totalUCN, uTotalUCN, numInCell, uNumInCell, transportEff, uTransportEff, avgTauCells, uavgTauCells, weightedDetEff, uWeightedDetEff, optimalT_edm, optimalEmptyingTime, daysToReach, errDaysToReach \n";
        for(int i = 0; i < study.filenames.size(); ++i){
            resultsFile <<study.filenames[i] << ", "<< study.results[i].parameter <<", " << study.results[i].fillTime << ", " <<study.results[i].fillTimeError << ", " <<study.results[i].efficiency << ", " << study.results[i].efficiencyError << ", " << study.results[i].sourceStorageLifetime << ", " << study.results[i].sourceStorageLifetimeError << ", " << study.results[i].systemTau << ", " << study.results[i].uSystemTau << ", " << study.results[i].totalUCN << ", " << study.results[i].uTotalUCN << ", " << study.results[i].numInCell << ", " << study.results[i].uNumInCell << ", " << study.results[i].transportEff << ", " << study.results[i].uTransportEff << ", " << study.results[i].avgTauCells << ", " << study.results[i].uavgTauCells << ", "<<study.results[i].WDetEff<< ", " << study.results[i].uWDetEff<< ", " <<study.results[i].optimalTedm << ", " << study.results[i].optimalEmptyingTime << ", "<< study.results[i].daysToReach << ", " << study.results[i].errDaysToReach << "\n";
        }

        resultsFile.close();     
    }
    outfile.Write();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/////
//Fucntions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void transport(std::string fileName, int modeParameter,SimResult &results)
{

    int mode = modeParameter;
	double neutronsSimulated = 1000000;	
	
	
	TCut topCell = TCut(Form("zend<%f&&zend>0&&xend>%f", cellUpperZBound, cellLowerXBound));
	TCut bottomCell = TCut(Form("zend<0&&zend>%f&&xend>%f", cellLowerZBound, cellLowerXBound));	
	TCut cell = (topCell || bottomCell);	

	//Open file//
	TFile f1(fileName.c_str());
        TTree *neutronend =(TTree*) f1.Get("neutronend");
        TTree *neutronsnapshot =(TTree*) f1.Get("neutronsnapshot");
	
    std::string fileNameString = fileName.substr(0, fileName.size()-5);
	
	//Outputs//
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
	
	TCut fillCondition;
	string fillEvolutionFileName;
	TCut tendCut;
	TCut valveOpenTimeCut = TCut(Form("tend==%.0f", valveOpenTime));
        string filledCellFileName;
	string SourceLifetimeFileName;
	TCut modeCuts;

	//Step 1 - Filling Time//	
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
	neutronsnapshot->Draw(Form("tend >> filling(%f, %f, %f)", activeTime-valveOpenTime, valveOpenTime, activeTime), cell && fillCondition); //use form for valve open time
	TH1S *filling = (TH1S*) gDirectory->Get("filling");
    
//	TF1 fillTimeFit("fit", "pol3", 105, 140);  
//    filling->Fit(&fillTimeFit, "R");
	filling->Draw();
    
    auto range = GetMaxBinRange(filling, 0.9545); // get range of bins that contain maximum with 95% probability
    
    fillTime = (range[0] + range[2])/2.0;
    fillTimeError = (fabs(fillTime - range[0]) + fabs(fillTime - range[2]))/2.0;
  
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
	tendCut = TCut(Form("tend == %.0f", fillTime));
	//Step 2 - Cell Energy Spectrum//
	
	TCanvas *d = new TCanvas("d", "d", 800, 600);
        neutronsnapshot->Draw("Hend>>cellESpectrum","Hstart<233.5e-9 && Hend < 250E-9" &&  cell && tendCut);
	    TH1 *cellESpectrum = (TH1*) gDirectory->Get("cellESpectrum");
	    cellESpectrum->SetTitle("Energy distribution of UCN at the end of filling");
	    cellESpectrum->GetXaxis()->SetTitle("UCN energy (eV)");
	    cellESpectrum->GetYaxis()->SetTitle("Number of UCN");
	    cellESpectrum->GetYaxis()->SetTitleOffset(1.3);
	    double quantiles[2];
	    double cutOff[2] = {0.05, 0.95};
	    cellESpectrum->GetQuantiles(2, quantiles, cutOff);
        cellESpectrum->Draw();
        double lowerBound = quantiles[0];
        double upperBound = quantiles[1];
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

	d = new TCanvas("d", "d", 800, 600);
	neutronsnapshot->Draw("Hend>>topcellESpectrum(50,0,250e-9)","Hend<250e-9" && topCell && tendCut);
	TH1 *topcellESpectrum = (TH1*) gDirectory->Get("topcellESpectrum");
	topcellESpectrum->SetDirectory(0);
	neutronsnapshot->Draw("Hend>>bottomcellESpectrum(50,0,250e-9)","Hend<250e-9" && bottomCell && tendCut);
	TH1 *bottomcellESpectrum = (TH1*) gDirectory->Get("bottomcellESpectrum");
	bottomcellESpectrum->SetDirectory(0);
    delete d;
    results.topCellESpec = topcellESpectrum;
    results.bottomCellESpec = bottomcellESpectrum;
	
	//Step 3 - Source Storage Lifetime//
		
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	neutronsnapshot->Draw(Form("tend-tstart>>lifetime(%.0f,5,%.0f)",valveOpenTime-10, valveOpenTime-5), Form("tstart < 5 && Hend > %e && Hend < %e", lowerBound, upperBound));
	TH1D *lifetime = (TH1D*) gDirectory->Get("lifetime");
	TF1 g2("fit","expo",5,valveOpenTime - 5);
	lifetime->Fit(&g2,"R");
	        
	Double_t slope = g2.GetParameter(1);
	Double_t slopeError = g2.GetParError(1);
	sourceStorageLifetime = -1/slope;
	sourceStorageLifetimeError = slopeError/(slope*slope);
	       
	lifetime->SetTitle("Storage lifetime of source");
	lifetime->GetXaxis()->SetTitle("UCN Lifetime [s]");
	lifetime->GetYaxis()->SetTitle("Number of UCN");
	lifetime->GetYaxis()->SetTitleOffset(1.3);
	c1->SetLogy(1);
	lifetime->Draw();
	       
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

	double ninetyPer = valveOpenTime - 2.4*sourceStorageLifetime; // 100 - optimal irradiation time
	
	if(mode == 0) {
		 modeCuts = TCut(Form("tstart<%.0f && tstart > %f", valveOpenTime, ninetyPer)); //batchMode;
        }
        else if(mode == 1) {
            modeCuts = TCut(Form("tstart>%f && tstart < %f", ninetyPer, fillTime));
        }
        else if(mode == 2) {
                modeCuts = TCut(Form("tstart>%.0f", valveOpenTime)); //steadyStateMode;
        }
	
	//Step 4 - Number of UCN in Cell//  Number in cell at the end of ideal filling time
	TCanvas *u = new TCanvas("u", "u", 800, 600);
	neutronsnapshot->Draw("zend:yend:xend>>inCell", cell && modeCuts && tendCut); //use same for batch and steady BEAM
	TH1 *inCell = (TH1*) gDirectory->Get("inCell");
	numInCell = inCell->GetEntries();
	uNumInCell = sqrt(numInCell);
	delete u;


	//Step 5 - Total number of UCN//
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
		transportEff = numInCell/numTransported;	// Number in cell at the end of ideal filling time divided by ucn produced within last 2.4* tau_storage of the source.  roughly 90 percent of them would survive up scatter
		uTransportEff = sqrt((uNumInCell/numTransported)*(uNumInCell/numTransported) + (numInCell*uNumTransported/numTransported/numTransported)*(numInCell*uNumTransported/numTransported/numTransported));
		}
	
	//Step 6 - Calculate Efficiencies//
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
    results.numInCell = numInCell;
    results.uNumInCell = uNumInCell;
    results.systemTau = systemStorageLifetime;
    results.uSystemTau = systemStorageLifetimeError;
    results.transportEff = transportEff;
    results.uTransportEff = uTransportEff;
    results.fillTimeError = fillTimeError;
}




//UCN remaining: look at N_0 that fill the cell and estimate how many survive after Ramsey cycle, emptying, collection as a function of energy (5 nev binning)
std::vector<double> UCNremainingAfterCycle(TH1 *topspectrum, TH1 *bottomspectrum){

    TFile toptaufile("topCell_hist.root");
    TFile bottomtaufile("bottomCell_hist.root");
    TH1 *toptauhist = (TH1*)toptaufile.Get("lifetime_1");
    TH1 *bottomtauhist = (TH1*)bottomtaufile.Get("lifetime_1");

    TFile topemptyingfile("topCell_emptying_DetEff.root");
    TFile bottomemptyingfile("bottomCell_emptying_DetEff.root");
    TH2 *topemptyinghist = (TH2*)topemptyingfile.Get("emptyEff");
    TH2 *bottomemptyinghist = (TH2*)bottomemptyingfile.Get("emptyEff");
    TH1 *topDetEff = topemptyinghist->ProjectionY("_detEff", 0,topemptyinghist->GetXaxis()->FindBin(emptyingTime));
    TH1 *bottomDetEff = bottomemptyinghist->ProjectionY("_detEff",0,bottomemptyinghist->GetXaxis()->FindBin(emptyingTime));
//    TH1 *topDetEff = (TH1*)topemptyingfile.Get("detEff");
//    TH1 *bottomDetEff = (TH1*)bottomemptyingfile.Get("detEff");

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
            dsurvivors += pow(dNitop*exp(-t_edm/tauitop), 2) + pow(Nitop*exp(-t_edm/tauitop)/tauitop/tauitop*dtauitop, 2);
        if (dtauibottom != 0)
            dsurvivors += pow(dNibottom*exp(-t_edm/tauibottom), 2) + pow(Nibottom*exp(-t_edm/tauibottom)/tauibottom/tauibottom*dtauibottom, 2);   
    }
    
    double NF = nFCellTop + nFCellBottom;
    double N0 = topspectrum->Integral() + bottomspectrum->Integral();
    double avgTauCell = 130.0 / log(N0/NF);
   // double uavgTauCell = 130.0 * sqrt(log(N0/NF)*(1/N0 + 1/NF));  //CHECK THIS!!!!!!!!!!!!!!!!!!!
    double uavgTauCell = 130.0 * (sqrt(1./N0 + 1./NF) * NF/N0) ;  //CHECK THIS!!!!!!!!!!!!!!!!!!!
    
    double topDetMean = wTopDet / nFCellTop;
    double bottomDetMean = wBottomDet / nFCellBottom;
    double utopDetMean =  sqrt(1.0/nFCellTop);
    double ubottomDetMean = sqrt(1.0/nFCellBottom);

    double wMeanDetEff = (topDetMean /pow(utopDetMean,2) + bottomDetMean /pow(ubottomDetMean,2))/(1./pow(utopDetMean,2) + 1./pow(ubottomDetMean,2));

    double uwMeanDetEff = sqrt(1/(1./pow(ubottomDetMean,2) + 1./pow(ubottomDetMean,2)));
    
    return {survivors/N0, sqrt(dsurvivors)/N0, avgTauCell, uavgTauCell,wMeanDetEff,uwMeanDetEff};
}


//////////////////////////////////////////////////////////////////////////////////


double days(int mode, SimResult &result)
{
    double daysToReach;

    double sourcePumpingTime; // function: t_irradiation 
    if (mode ==0)
        sourcePumpingTime = result.sourceStorageLifetime * 2.4;
    else if (mode ==1)
        sourcePumpingTime = result.sourceStorageLifetime * 2.4 + result.fillTime;
    else if (mode ==2)
        sourcePumpingTime = 100;  //need to change this later to fillingTime of mode =2
    
    double cyclesToReach;
    double cyclesPerDay;
    double sensitivityPerFill;
   
    std::vector<double> survprob = UCNremainingAfterCycle(result.topCellESpec, result.bottomCellESpec);
    result.survivalprob  = survprob[0];
    result.dSurvivalprob = survprob[1];
    result.avgTauCells = survprob[2];
    result.uavgTauCells = survprob[3];
    result.WDetEff = survprob[4];
    result.uWDetEff = survprob[5];


    //cyclesPerDay
    if (mode ==1)
        cyclesPerDay = (stableField * 60)/((fillsPerCycle * (emptyingTime + t_edm +sourcePumpingTime) + 2*polarityTime + degaussingTime/10)/60);
    else
        cyclesPerDay = (stableField * 60)/((fillsPerCycle * (result.fillTime + emptyingTime + t_edm +sourcePumpingTime) + 2*polarityTime + degaussingTime/10)/60);
   
    
    //cyclesToReach
    
    //sensitivityPerFill
    double N_0;  //Initial number of UCN
    double N_det;  //Number of UCN that get detected
    double alpha_after;  //Visibility after Ramsey sequence

    N_0 = spinTrans * beamHeating * UCN_production * sourcePumpingTime * result.efficiency;

//    N_det = N_0 * exp(-t_edm/tau_walls -t_edm/tau_Xe) * transToDet * detEff;
    // N_det = N_0 * result.survivalprob * exp(-t_edm/tau_Xe) * transToDet * detEff;
      N_det = N_0 * result.survivalprob * exp(-t_edm/tau_Xe) * detEff;

    alpha_after = alpha * exp(-((t_edm- t_wait - 2*t_pulse) / T2) - (t_wait + 2* t_pulse)/T1) * detSpinTrans * Panalyzer;

    sensitivityPerFill = hbar / (2 * (t_edm - t_wait - 2*t_pulse) * EField * sqrt(N_det) * alpha_after);

    cyclesToReach = 100 * (sensitivityPerFill/sqrt(fillsPerCycle)/(1E-26)) * (sensitivityPerFill/sqrt(fillsPerCycle)/(1E-26));

    daysToReach =  cyclesToReach / cyclesPerDay;
    
    return daysToReach;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyzeSim(SimResult &result, std::string fileName, int mode)
{
         std::cout << fileName <<" , " << mode << std::endl;
         std::string fileNameString = fileName.substr(0, fileName.size()-5);
	     std::vector<double> cellNumber;
	     transport(fileName, mode, result);

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
    	 result.numInCell = cellNumber[0];
    	 result.uNumInCell = cellNumber[1];
}




// Total number of UCN in cell approximation
/////////////////////////////////////////////////////////////////
std::vector<double> cellTotal(int mode, SimResult results)
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
    
    cellUCN.push_back(totalSourceUCN * results.efficiency * polarization);
    cellUCN.push_back(cellUCN[0] * results.efficiencyError / results.efficiency);

    return cellUCN;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void makeGraph(Study &study,std::string xName, std::string yName, ResultMember mem_xs, ResultMember mem_ys, ResultMember mem_xerr, ResultMember mem_yerr)
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
   //    gr -> Fit("pol0","F");
   //
       gr -> SetTitle(Form("%s: %s vs %s", study.parameter_name.c_str(), yName.c_str(), xName.c_str())); 
       gr -> GetXaxis()-> SetTitle(xName.c_str());
       gr -> GetYaxis()-> SetTitle(yName.c_str());
   //    gr -> GetYaxis()-> SetRangeUser(100, 1000000);
   //    gr -> GetXaxis()-> SetRangeUser(0, 50);
       gr -> GetXaxis() -> SetTitleSize(0.05);
       gr -> GetXaxis() -> SetTitleOffset(0.75);
       gr -> GetYaxis() -> SetTitleSize(0.05);
       gr -> GetYaxis() -> SetTitleOffset(0.9);
       gr -> Draw("Ap");
    can->Print(Form("%s_%s_vs_%s.eps",study.parameter_name.c_str(), yName.c_str() ,xName.c_str()));
    delete gr;
    delete can;

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
std::array<double, 3> GetMaxBinRange(const TH1* hist, const double confidence_interval = 0.6827){
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
