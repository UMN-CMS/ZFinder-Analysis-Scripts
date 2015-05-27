GETTING EVERYTHING WORKING.

Assuming you downloaded everything properly the first step is to change the path inside rootlogon.C to the new location of RooUnfold-1.1.1. Next you can cmsenv in a release that you want to try. I have been using CMSSW_7_1_1 and  haven't been having any issues though that may change. Next you start root. Currently almost all the variables are hard coded in so you need to change them inside the code and recompile

MakeFinalPlots.C
Assuming you want to make final plots you can change what is graphed b changing the string type to either elec, muon or combined.  
