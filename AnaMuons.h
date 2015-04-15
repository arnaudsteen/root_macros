//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 11 15:17:56 2014 by ROOT version 5.34/05
// from TTree tree/Shower variables
// found on file: TDHCAL_715747.root
//////////////////////////////////////////////////////////

#ifndef AnaMuons_h
#define AnaMuons_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <map>

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnaMuons {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  std::map<int,float> mulmap;
  std::map<int,float> mulsquaremap;
  std::map<int,int> mulcount;
  std::map<int,float> effmap;
  std::map<int,float> effsquaremap;
  std::map<int,int> effcount;
  std::map<int,float> eff2map;
  std::map<int,float> eff2squaremap;
  std::map<int,float> eff3map;
  std::map<int,float> eff3squaremap;

  TH1D *hTheta;
  TH2D *hMul;
  TH2D *hEff;
  TH2D *hMulLayer;
  TH2D *hEffLayer;
  TH2D *hEff2Layer;
  TH2D *hEff3Layer;
  TF1* fit_func;

  // Declaration of leaf types
  ULong64_t       eventTime;
  ULong64_t       spillEventTime;
  Float_t         transversRatio;
  Int_t           Nhit1;
  Int_t           Nhit2;
  Int_t           Nhit3;
  Int_t           Nlayer;
  Int_t           Ebeam;
  Float_t         Efficiency1[48];
  Float_t         Efficiency2[48];
  Float_t         Efficiency3[48];
  Float_t         Multiplicity[48];
  Float_t         CorrectedMultiplicity[48];
  vector<int>     *ClusterSize;
  Float_t         Chi2[48];
  Int_t           evtNum;
  Float_t         effGlobal;
  Float_t         mulGlobal;
  Float_t         chi2Global;
  Int_t           trackEnd;
  Float_t         trackParams[4];

  // List of branches
  TBranch        *b_eventTime;   //!
  TBranch        *b_spillEventTime;   //!
  TBranch        *b_transversRatio;   //!
  TBranch        *b_Nhit1;   //!
  TBranch        *b_Nhit2;   //!
  TBranch        *b_Nhit3;   //!
  TBranch        *b_Nlayer;   //!
  TBranch        *b_Ebeam;   //!
  TBranch        *b_Efficiency1;   //!
  TBranch        *b_Efficiency2;   //!
  TBranch        *b_Efficiency3;   //!
  TBranch        *b_Multiplicity;   //!
  TBranch        *b_CorrectedMultiplicity;   //!
  TBranch        *b_ClusterSize;   //!
  TBranch        *b_Chi2;   //!
  TBranch        *b_evtNum;   //!
  TBranch        *b_effGlobal;   //!
  TBranch        *b_mulGlobal;   //!
  TBranch        *b_chi2Global;   //!
  TBranch        *b_trackEnd;   //!
  TBranch        *b_trackParams;   //!

  AnaMuons(const char* input);
  virtual ~AnaMuons();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  TFile* file;
};

#endif

#ifdef AnaMuons_cxx
AnaMuons::AnaMuons(const char* input) : fChain(0) 
{
  file = new TFile(input);
  TTree *tree = (TTree*)file->Get("tree");
  Init(tree);
  mulmap.clear();
  mulsquaremap.clear();
  effmap.clear();
  effsquaremap.clear();
  effcount.clear();
  mulcount.clear();
}

AnaMuons::~AnaMuons()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t AnaMuons::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t AnaMuons::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void AnaMuons::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  ClusterSize = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("eventTime", &eventTime, &b_eventTime);
  fChain->SetBranchAddress("spillEventTime", &spillEventTime, &b_spillEventTime);
  fChain->SetBranchAddress("transversRatio", &transversRatio, &b_transversRatio);
  fChain->SetBranchAddress("Nhit1", &Nhit1, &b_Nhit1);
  fChain->SetBranchAddress("Nhit2", &Nhit2, &b_Nhit2);
  fChain->SetBranchAddress("Nhit3", &Nhit3, &b_Nhit3);
  fChain->SetBranchAddress("Nlayer", &Nlayer, &b_Nlayer);
  fChain->SetBranchAddress("Ebeam", &Ebeam, &b_Ebeam);
  fChain->SetBranchAddress("Efficiency1", Efficiency1, &b_Efficiency1);
  fChain->SetBranchAddress("Efficiency2", Efficiency2, &b_Efficiency2);
  fChain->SetBranchAddress("Efficiency3", Efficiency3, &b_Efficiency3);
  fChain->SetBranchAddress("Multiplicity", Multiplicity, &b_Multiplicity);
  fChain->SetBranchAddress("CorrectedMultiplicity", CorrectedMultiplicity, &b_CorrectedMultiplicity);
  fChain->SetBranchAddress("ClusterSize", &ClusterSize, &b_ClusterSize);
  fChain->SetBranchAddress("Chi2", Chi2, &b_Chi2);
  fChain->SetBranchAddress("evtNum", &evtNum, &b_evtNum);
  fChain->SetBranchAddress("effGlobal", &effGlobal, &b_effGlobal);
  fChain->SetBranchAddress("mulGlobal", &mulGlobal, &b_mulGlobal);
  fChain->SetBranchAddress("chi2Global", &chi2Global, &b_chi2Global);
  fChain->SetBranchAddress("trackEnd", &trackEnd, &b_trackEnd);
  fChain->SetBranchAddress("trackParams", trackParams, &b_trackParams);
  Notify();
}

Bool_t AnaMuons::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void AnaMuons::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t AnaMuons::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void CaliceStyle()
{
  /*CALICE style for figure: use in a ROOT macro like this:*/
  //gROOT->ProcessLine(".L ~/RootStuff/CaliceStyle.C");
  //CaliceStyle();

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetStatFont(42);

  gStyle->SetFrameFillColor(kWhite);
  //  gStyle->SetFrameLineWidth(2);
  gStyle->SetCanvasColor(kWhite);  
  gStyle->SetOptStat(0); /*don't show statistics box*/
  gStyle->SetTitleSize(0.05, "xyz"); 
  gStyle->SetLegendBorderSize();

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gROOT->ForceStyle();
}

#endif // #ifdef AnaMuons_cxx
