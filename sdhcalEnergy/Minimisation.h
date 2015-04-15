//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 27 10:01:29 2013 by ROOT version 5.34/01
// from TTree tree/Shower variables
// found on file: TMVA/ROOTFiles/single_pi-_20GeV.root
//////////////////////////////////////////////////////////

#ifndef Minimisation_h
#define Minimisation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <TFitter.h>
#include <cmath>
#include <iostream>
#include <TMinuit.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

void Minimize(int nEpoch);
double myFunc(double x1, double y1, double z1, 
	      double x2, double y2, double z2, 
	      double x3, double y3, double z3,
	      double cste);

std::vector<int> n0;
std::vector<int> n1;
std::vector<int> n2;
std::vector<int> n3;
std::vector<int> ho1;
std::vector<int> ho2;
std::vector<int> ho3;
std::vector<int> energy;
double BestW[10];
void minuitFunction(int& nDim, double* gout,double& result, 
		    double par[],int flg);

class Minimisation {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  float emean;
  float eres;
  float err_emean;
  float err_eres;
  Float_t cutFD;
  Float_t cutFD_bis;
  std::vector<int> nhit;
  std::vector<int> nhit1;
  std::vector<int> nhit2;
  std::vector<int> nhit3;
  std::vector<int> hough1;
  std::vector<int> hough2;
  std::vector<int> hough3;
  std::vector<int> time;
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
  Int_t           DeadHitNumber;
  Int_t           Nlayer;
  Int_t           NInteractinglayer;
  Int_t           Begin;
  Int_t           End;
  Float_t         Radius;
  Int_t           LongiProfile[48];
  Int_t           LongiProfileBis[48];
  Int_t           RadialProfile[96];
  Float_t         MaxRadius;
  Int_t           Lasthit;
  Int_t           Hole;
  Float_t         CoG[3];
  Float_t         FirstLayerRatio;
  Float_t         CentralRatio;
  Float_t         F3D;
  vector<double>  *TrackLength;
  Int_t           TrackMultiplicity;
  vector<double>  *EfficiencyPerLayer;
  vector<double>  *MultiplicityPerLayer;
  Float_t         MeanClusterSize;
  Int_t           Nclusters;
  Int_t           NclusterLayer[48];
  vector<int>     *TrackClusterSize;
  vector<int>     *TrackClusterNumber;
  Int_t           ClusterMip;
  Int_t           ClusterEM;
  Int_t           Edge;
  Int_t           Neutral;
  Int_t           Single;

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
  TBranch        *b_DeadHitNumber;   //!
  TBranch        *b_Nlayer;   //!
  TBranch        *b_NInteractinglayer;   //!
  TBranch        *b_Begin;   //!
  TBranch        *b_End;   //!
  TBranch        *b_Radius;   //!
  TBranch        *b_LongiProfile;   //!
  TBranch        *b_LongiProfileBis;   //!
  TBranch        *b_RadialProfile;   //!
  TBranch        *b_MaxRadius;   //!
  TBranch        *b_Lasthit;   //!
  TBranch        *b_Hole;   //!
  TBranch        *b_CoG;   //!
  TBranch        *b_FirstLayerRatio;   //!
  TBranch        *b_CentralRatio;   //!
  TBranch        *b_F3D;   //!
  TBranch        *b_TrackLength;   //!
  TBranch        *b_TrackMultiplicity;   //!
  TBranch        *b_EfficiencyPerLayer;   //!
  TBranch        *b_MultiplicityPerLayer;   //!
  TBranch        *b_MeanClusterSize;   //!
  TBranch        *b_Nclusters;   //!
  TBranch        *b_NclusterLayer;   //!
  TBranch        *b_TrackClusterSize;   //!
  TBranch        *b_TrackClusterNumber;   //!
  TBranch        *b_ClusterMip;   //!
  TBranch        *b_ClusterEM;   //!
  TBranch        *b_Edge;   //!
  TBranch        *b_Neutral;   //!
  TBranch        *b_Single;   //!

  Minimisation(const char* input);
  virtual ~Minimisation();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual void     Erec(int ebeam);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
 private:
  TFile *file;
};

#endif

#ifdef Minimisation_cxx
Minimisation::Minimisation(const char* input) : fChain(0)
{
  file = new TFile(input);
  TTree *tree = (TTree*)file->Get("tree");

  Init(tree);
  cutFD=0;
  cutFD_bis=0;
  nhit.clear();
  nhit1.clear();
  nhit2.clear();
  nhit3.clear();
  hough1.clear();
  hough2.clear();
  hough3.clear();
  time.clear();
}

Minimisation::~Minimisation()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t Minimisation::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t Minimisation::LoadTree(Long64_t entry)
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

void Minimisation::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
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
  fChain->SetBranchAddress("DeadHitNumber", &DeadHitNumber, &b_DeadHitNumber);
  fChain->SetBranchAddress("Nlayer", &Nlayer, &b_Nlayer);
  fChain->SetBranchAddress("NInteractinglayer", &NInteractinglayer, &b_NInteractinglayer);
  fChain->SetBranchAddress("Begin", &Begin, &b_Begin);
  fChain->SetBranchAddress("End", &End, &b_End);
  fChain->SetBranchAddress("Radius", &Radius, &b_Radius);
  fChain->SetBranchAddress("LongiProfile", LongiProfile, &b_LongiProfile);
  fChain->SetBranchAddress("LongiProfileBis", LongiProfileBis, &b_LongiProfileBis);
  fChain->SetBranchAddress("RadialProfile", RadialProfile, &b_RadialProfile);
  fChain->SetBranchAddress("MaxRadius", &MaxRadius, &b_MaxRadius);
  fChain->SetBranchAddress("Lasthit", &Lasthit, &b_Lasthit);
  fChain->SetBranchAddress("Hole", &Hole, &b_Hole);
  fChain->SetBranchAddress("CoG", CoG, &b_CoG);
  fChain->SetBranchAddress("FirstLayerRatio", &FirstLayerRatio, &b_FirstLayerRatio);
  fChain->SetBranchAddress("CentralRatio", &CentralRatio, &b_CentralRatio);
  fChain->SetBranchAddress("F3D", &F3D, &b_F3D);
  fChain->SetBranchAddress("TrackLength", &TrackLength, &b_TrackLength);
  fChain->SetBranchAddress("TrackMultiplicity", &TrackMultiplicity, &b_TrackMultiplicity);
  fChain->SetBranchAddress("EfficiencyPerLayer", &EfficiencyPerLayer, &b_EfficiencyPerLayer);
  fChain->SetBranchAddress("MultiplicityPerLayer", &MultiplicityPerLayer, &b_MultiplicityPerLayer);
  fChain->SetBranchAddress("MeanClusterSize", &MeanClusterSize, &b_MeanClusterSize);
  fChain->SetBranchAddress("Nclusters", &Nclusters, &b_Nclusters);
  fChain->SetBranchAddress("NclusterLayer", NclusterLayer, &b_NclusterLayer);
  fChain->SetBranchAddress("TrackClusterSize", &TrackClusterSize, &b_TrackClusterSize);
  fChain->SetBranchAddress("TrackClusterNumber", &TrackClusterNumber, &b_TrackClusterNumber);
  fChain->SetBranchAddress("ClusterMip", &ClusterMip, &b_ClusterMip);
  fChain->SetBranchAddress("ClusterEM", &ClusterEM, &b_ClusterEM);
  fChain->SetBranchAddress("Edge", &Edge, &b_Edge);
  fChain->SetBranchAddress("Neutral", &Neutral, &b_Neutral);
  fChain->SetBranchAddress("Single", &Single, &b_Single);
  Notify();
}

Bool_t Minimisation::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void Minimisation::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t Minimisation::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef Minimisation_cxx
