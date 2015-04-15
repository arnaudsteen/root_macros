//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 15 15:08:31 2014 by ROOT version 5.34/05
// from TTree tree/Shower variables
// found on file: /home/steen/resultRootFile/tb_data/DHCAL_715747.root
//////////////////////////////////////////////////////////

#ifndef NclusterPion_h
#define NclusterPion_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TText.h>
#include <TLegend.h>

// Header file for the classes stored in the TTree if any.
#include <iostream>
#include <vector>
#include <vector>
#include <string>
#include <limits>

// Fixed size dimensions of array or collections stored in the TTree if any.

class NclusterPion {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  Float_t cutBary[4];
  std::vector<int> ncluster;
  std::vector<float> time;
  
  struct radialProfile{
    int ring[96];
  };
  struct longitudinalProfile{
    int layer[48];
  };
  std::vector<radialProfile> radialVec;
  std::vector<longitudinalProfile> longitudinalVec;
  // Declaration of leaf types
  ULong64_t       eventTime;
  ULong64_t       spillEventTime;
  Int_t           NShowers;
  Int_t           Nhit;
  Int_t           Nhit1;
  Int_t           Nhit2;
  Int_t           Nhit3;
  Int_t           Nhough1;
  Int_t           Nhough2;
  Int_t           Nhough3;
  Int_t           Nlayer;
  Int_t           NInteractinglayer;
  Int_t           Begin;
  Int_t           End;
  Float_t         Radius;
  Int_t           LongiProfile[48];
  Int_t           LongiProfileBis[48];
  Int_t           RadialProfile[96];
  Int_t           ClusterRadialProfile[96];
  Int_t           RadialProfileBis[96];
  Float_t         MaxRadius;
  Int_t           Hole;
  Float_t         CoG[4];
  Float_t         CentralRatio;
  Float_t         F3D;
  vector<int>     *Density;
  vector<double>  *TrackLength;
  Int_t           TrackMultiplicity;
  vector<double>  *EfficiencyPerLayer;
  vector<double>  *MultiplicityPerLayer;
  Float_t         MeanClusterSize;
  Int_t           Nclusters;
  vector<int>     *TrackClusterSize;
  vector<int>     *TrackClusterNumber;
  Int_t           ClusterMip;
  Int_t           ClusterEM;
  Int_t           ClusterIsolated;
  Int_t           Neutral;
  Int_t           Single;
  Float_t         TransverseRatio;
  Float_t         IncidentParticleCosTheta;
  Float_t         ReconstructedCosTheta;

  // List of branches
  TBranch        *b_eventTime;   //!
  TBranch        *b_spillEventTime;   //!
  TBranch        *b_NShowers;   //!
  TBranch        *b_Nhit;   //!
  TBranch        *b_Nhit1;   //!
  TBranch        *b_Nhit2;   //!
  TBranch        *b_Nhit3;   //!
  TBranch        *b_Nhough1;   //!
  TBranch        *b_Nhough2;   //!
  TBranch        *b_Nhough3;   //!
  TBranch        *b_Nlayer;   //!
  TBranch        *b_NInteractinglayer;   //!
  TBranch        *b_Begin;   //!
  TBranch        *b_End;   //!
  TBranch        *b_Radius;   //!
  TBranch        *b_LongiProfile;   //!
  TBranch        *b_LongiProfileBis;   //!
  TBranch        *b_RadialProfile;   //!
  TBranch        *b_ClusterRadialProfile;   //!
  TBranch        *b_RadialProfileBis;   //!
  TBranch        *b_MaxRadius;   //!
  TBranch        *b_Hole;   //!
  TBranch        *b_CoG;   //!
  TBranch        *b_CentralRatio;   //!
  TBranch        *b_F3D;   //!
  TBranch        *b_Density;   //!
  TBranch        *b_TrackLength;   //!
  TBranch        *b_TrackMultiplicity;   //!
  TBranch        *b_EfficiencyPerLayer;   //!
  TBranch        *b_MultiplicityPerLayer;   //!
  TBranch        *b_MeanClusterSize;   //!
  TBranch        *b_Nclusters;   //!
  TBranch        *b_TrackClusterSize;   //!
  TBranch        *b_TrackClusterNumber;   //!
  TBranch        *b_ClusterMip;   //!
  TBranch        *b_ClusterEM;   //!
  TBranch        *b_ClusterIsolated;   //!
  TBranch        *b_Neutral;   //!
  TBranch        *b_Single;   //!
  TBranch        *b_TransverseRatio;   //!
  TBranch        *b_IncidentParticleCosTheta;   //!
  TBranch        *b_ReconstructedCosTheta;   //!

  NclusterPion(const char* input);
  virtual ~NclusterPion();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual void ShowerBarycenterCut(std::string runPeriod);
  TFile *file;
};

#endif

#ifdef NclusterPion_cxx
NclusterPion::NclusterPion(const char*input) : fChain(0) 
{
  file = new TFile(input);
  TTree *tree = (TTree*)file->Get("tree");
  Init(tree);
  cutBary[0]=-300;
  cutBary[1]=300;
  cutBary[2]=-300;
  cutBary[3]=300;
  ncluster.clear();
  time.clear();
  radialVec.clear();
  longitudinalVec.clear();
  //const char* badFile="/home/steen/Back-Up/TB_Data/SPS_2012/RootFiles/DHCAL_716305.root";
  //if(strcmp (input,badFile) == 0)
  //  maxentry=356870;
  //else maxentry=std::numeric_limits<int>::max();
  
}

NclusterPion::~NclusterPion()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t NclusterPion::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t NclusterPion::LoadTree(Long64_t entry)
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

void NclusterPion::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  Density = 0;
  TrackLength = 0;
  EfficiencyPerLayer = 0;
  MultiplicityPerLayer = 0;
  TrackClusterSize = 0;
  TrackClusterNumber = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("eventTime", &eventTime, &b_eventTime);
  fChain->SetBranchAddress("spillEventTime", &spillEventTime, &b_spillEventTime);
  fChain->SetBranchAddress("NShowers", &NShowers, &b_NShowers);
  fChain->SetBranchAddress("Nhit", &Nhit, &b_Nhit);
  fChain->SetBranchAddress("Nhit1", &Nhit1, &b_Nhit1);
  fChain->SetBranchAddress("Nhit2", &Nhit2, &b_Nhit2);
  fChain->SetBranchAddress("Nhit3", &Nhit3, &b_Nhit3);
  fChain->SetBranchAddress("Nhough1", &Nhough1, &b_Nhough1);
  fChain->SetBranchAddress("Nhough2", &Nhough2, &b_Nhough2);
  fChain->SetBranchAddress("Nhough3", &Nhough3, &b_Nhough3);
  fChain->SetBranchAddress("Nlayer", &Nlayer, &b_Nlayer);
  fChain->SetBranchAddress("NInteractinglayer", &NInteractinglayer, &b_NInteractinglayer);
  fChain->SetBranchAddress("Begin", &Begin, &b_Begin);
  fChain->SetBranchAddress("End", &End, &b_End);
  fChain->SetBranchAddress("Radius", &Radius, &b_Radius);
  fChain->SetBranchAddress("LongiProfile", LongiProfile, &b_LongiProfile);
  fChain->SetBranchAddress("LongiProfileBis", LongiProfileBis, &b_LongiProfileBis);
  fChain->SetBranchAddress("RadialProfile", RadialProfile, &b_RadialProfile);
  fChain->SetBranchAddress("ClusterRadialProfile", ClusterRadialProfile, &b_ClusterRadialProfile);
  fChain->SetBranchAddress("RadialProfileBis", RadialProfileBis, &b_RadialProfileBis);
  fChain->SetBranchAddress("MaxRadius", &MaxRadius, &b_MaxRadius);
  fChain->SetBranchAddress("Hole", &Hole, &b_Hole);
  fChain->SetBranchAddress("CoG", CoG, &b_CoG);
  fChain->SetBranchAddress("CentralRatio", &CentralRatio, &b_CentralRatio);
  fChain->SetBranchAddress("F3D", &F3D, &b_F3D);
  fChain->SetBranchAddress("Density", &Density, &b_Density);
  fChain->SetBranchAddress("TrackLength", &TrackLength, &b_TrackLength);
  fChain->SetBranchAddress("TrackMultiplicity", &TrackMultiplicity, &b_TrackMultiplicity);
  fChain->SetBranchAddress("EfficiencyPerLayer", &EfficiencyPerLayer, &b_EfficiencyPerLayer);
  fChain->SetBranchAddress("MultiplicityPerLayer", &MultiplicityPerLayer, &b_MultiplicityPerLayer);
  fChain->SetBranchAddress("MeanClusterSize", &MeanClusterSize, &b_MeanClusterSize);
  fChain->SetBranchAddress("Nclusters", &Nclusters, &b_Nclusters);
  fChain->SetBranchAddress("TrackClusterSize", &TrackClusterSize, &b_TrackClusterSize);
  fChain->SetBranchAddress("TrackClusterNumber", &TrackClusterNumber, &b_TrackClusterNumber);
  fChain->SetBranchAddress("ClusterMip", &ClusterMip, &b_ClusterMip);
  fChain->SetBranchAddress("ClusterEM", &ClusterEM, &b_ClusterEM);
  fChain->SetBranchAddress("ClusterIsolated", &ClusterIsolated, &b_ClusterIsolated);
  fChain->SetBranchAddress("Neutral", &Neutral, &b_Neutral);
  fChain->SetBranchAddress("Single", &Single, &b_Single);
  fChain->SetBranchAddress("TransverseRatio", &TransverseRatio, &b_TransverseRatio);
  fChain->SetBranchAddress("IncidentParticleCosTheta", &IncidentParticleCosTheta, &b_IncidentParticleCosTheta);
  fChain->SetBranchAddress("ReconstructedCosTheta", &ReconstructedCosTheta, &b_ReconstructedCosTheta);
  Notify();
}

Bool_t NclusterPion::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void NclusterPion::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t NclusterPion::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

int Xmin(int E)
{
  if(E==5) return 0;
  else if(E==10) return 0;
  else if(E==15) return 0;
  else if(E==20) return 0;
  else if(E==25) return 0;
  else if(E==30) return 0;
  else if(E==40) return 0;
  else if(E==50) return 0;
  else if(E==60) return 0;
  else if(E==70) return 0;
  else if(E==80) return 0;
  else
    std::cout << "Energy = " << E << " GeV IS NOT A CORRECT ENERGY" << std::endl;
  return 0;
}

int Xmax(int E)
{
  if(E==5) return 100;
  else if(E==10) return 150;
  else if(E==15) return 180;
  else if(E==20) return 200;
  else if(E==25) return 240;
  else if(E==30) return 270;
  else if(E==40) return 330;
  else if(E==50) return 380;
  else if(E==60) return 400;
  else if(E==70) return 440;
  else if(E==80) return 450;
  else
    std::cout << "Energy = " << E << " GeV IS NOT A CORRECT ENERGY" << std::endl;
  return 0;
}

int Nbin(int E)
{
  if(E==5) return 20;
  else if(E==10) return 30;
  else if(E==15) return 36;
  else if(E==20) return 40;
  else if(E==25) return 48;
  else if(E==30) return 54;
  else if(E==40) return 55;
  else if(E==50) return 63;
  else if(E==60) return 66;
  else if(E==70) return 73;
  else if(E==80) return 75;
  else
    std::cout << "Energy = " << E << " GeV IS NOT A CORRECT ENERGY" << std::endl;
  return 0;
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
  gStyle->SetFrameLineWidth(1);
  gStyle->SetCanvasColor(kWhite);  
  gStyle->SetOptStat(0); /*don't show statistics box*/
  gStyle->SetTitleSize(0.05, "xyz"); 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetOptTitle(0);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.05);

  gROOT->ForceStyle();
}
#endif // #ifdef NclusterPion_cxx
