//going to graph gen Y to reco Y
// Standard Library
#include <stdexcept>
#include <iostream>
#include <sstream>


#include <TFile.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLeaf.h>

#include "TH1.h"
#include <TTree.h>
#include "TSystem.h"

#include <algorithm>
#include <set>
#include <vector>

using namespace std; // yeah I am evil, get over it

TTree* GetTTree(const std::string TFILE, const std::string TTREE) {
    // Open the TFiles
    TFile* mc_tfile = new TFile(TFILE.c_str(), "READ");
    if (!mc_tfile || mc_tfile->IsZombie()) {
        const std::string ERR = "Could not open the file " + TFILE;
        throw std::runtime_error(ERR);
    }

    // Load the tree
    TTree* tree = new TTree();
    mc_tfile->GetObject(TTREE.c_str(), tree);
    if (!tree || tree->IsZombie()) {
        const std::string ERR = "Could not open the TTree " + TTREE;
        throw std::runtime_error(ERR);
    }

    return tree;
}

void GenVSRecoY() {
    const string FileName = "";
    const string TreeName = "Combined Single Reco";
    const string RecoBranch = "reco";
    const string TruthBranch = "truth";
    const string RapidityLeaf = "z_y";
    TTree* Tree = GetTTree(FileName, TreeName);
    TBranch* RBranch = Tree->GetBranch(RecoBranch);//reco branch
    TBranch* TBranch = Tree->GetBranch(TruthBranch);//truth branch dumbass
    TLeaf* RecoRap = RBranch->GetLeaf();
    TLeaf* TruthRap = TBranch->GetLeaf();
    
    TH2* GenVSReco_Rap= new TH2F("Z_Y","Z_Y",100,-2,2,100,-2,2); 
    for (int i = 0; i < Tree->GetEntries(); i++) {
        Tree->GetEntry(i);
        GenVSReco_Rap->Fill(TruthRap->GetValue(),RecoRap->GetValue());
    }
    
    TFile *theFile = new TFile("GenvsRecoY.root", "RECREATE");
    GenVSReco_Rap->Write();
    delete GenVSReco_Rap;
    
}