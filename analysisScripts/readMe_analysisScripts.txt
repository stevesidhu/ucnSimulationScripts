2018testSimulationAnalysis.C
    Oldest version of analysis script


All here need SimResult.h
2019testSimulationAnalysisV1.C
    Used for old simulations with cell cuts = zend cuts

2019testSimulationAnalysisV2.C
    Contains dual functionality to use EDMcell vacuum volume in the cells

testSimulationAnalysisV3.C
    Fixed many bugs, fixed plotting.  Updated error in days to reach.  Added more comments throught.  Still need to add more comments for Days function

Needs SimResult1.h
testSimulationAnalysisV4.C
    Elinimated error analysis for days to reach since it was wrong.  Added simulated annealing minimization instead.  Added 3D histograms so how sourcePumpingTime goes with time as well

daystoreach.py
    New python script to replace C scripts.





cellTau.C 
    Script used for determining the storage lifetime per energy bin in the EDM cells

cellEmpty.C 
    Script used for determining collection efficiency as a function of energy for EDM cells (to the detectors)

