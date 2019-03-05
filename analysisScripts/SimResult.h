
#include "TH1.h"

struct SimResult{
    TH1F *topCellESpec;
    TH1F *bottomCellESpec;
    TH1F topCellESpecFinish;
    TH1F bottomCellESpecFinish;
    TH1D topCellESpecDetected;
    TH1D bottomCellESpecDetected;
    double parameter;
    double totalSimulated;
    double totalProducedReal;
    double fillTime;
    double fillTimeError;
    double efficiency;
    double efficiencyError;
    double sourceStorageLifetime;
    double sourceStorageLifetimeError;
    double systemTau;
    double uSystemTau;
    double totalUCN;  //That reach EDM cells in the simulation
    double uTotalUCN;
    double numInCellSim;
    double uNumInCellSim;
    double numInCellReal; //Estimated number uf UCN in cells in real experiment based on estimate of source production
    double uNumInCellReal;
    double transportEff;
    double uTransportEff;
    double survivalprob;
    double dSurvivalprob;
    double daysToReach;
    double errDaysToReach;
    double avgTauCells;
    double uavgTauCells;
    double wTopDetEff;
    double uWTopDetEff;
    double wBottomDetEff;
    double uWBottomDetEff;
    double WDetEff;
    double uWDetEff;
    double optimalTedm;
    double emptyingTime;
    double uOptimalTedm;
    double uEmptyingTime;

    double topMeanFill;
    double uTopMeanFill;
    double bottomMeanFill;
    double uBottomMeanFill;
    double topMeanFinish;
    double uTopMeanFinish;
    double bottomMeanFinish;
    double uBottomMeanFinish;
    double topMeanDetected;
    double uTopMeanDetected;
    double bottomMeanDetected;
    double uBottomMeanDetected;
    //Estimated survivors after Ramsey cyclefrom survivors equation from calculation
    double topCellRamsey;
    double bottomCellRamsey; //
    // integrated totals from histograms
    double topCellTotalFilled;
    double bottomCellTotalFilled;
    double topCellTotalFinished;
    double bottomCellTotalFinished;
    double topCellTotalDetected;
    double bottomCellTotalDetected;
    double totalCollected;
};

