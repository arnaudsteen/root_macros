#define LinearSimEreco_cxx
#include "LinearSimEreco.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "iostream"
#include "TF1.h"
#include <TText.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include "cmath"
#include "fstream"
#include "algorithm"
#include "numeric"
#include "string"

void LinearSimEreco::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if( 
       (float)Nhit/Nlayer>3 && 
       Begin>=0 && 
       (float)NInteractinglayer/Nlayer>0.01
	){
      nhit.push_back(Nhit);
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
    }
  }
}

void writeHisto(TH1D* h,std::string fileName,std::string option)
{
  TFile *file=new TFile(fileName.c_str(),option.c_str());
  file->cd();
  h->Write();
  file->Write();
  file->Close();
}

double myLinearFunc(double w0, double w1, double w2)
{
  unsigned int size = n0.size();
  float test = 0;
  float energy_reco;
  for(unsigned int i=0; i<size; i++){
    energy_reco = w0*n1[i] + w1*n2[i] + w2*n3[i];
		   
    //    test+=fabs((energy[i]-energy_reco)/energy[i]);
    test = test + ( (energy[i]-energy_reco)*
		    (energy[i]-energy_reco)/energy[i] );
  }
  return test/size;
}

void minuitLinearFunction(int& nDim, double* gout,double& result, 
			  double *par,int flg){
  result = myLinearFunc(par[0], par[1], par[2]);
}

void LinearMinimize(int nEpoch)
{
  double bestW[3];
  bestW[0] = 0.1;//0289365; 
  bestW[1] = 0.1;//112674; 
  bestW[2] = 0.1;//320099; 

  double tryW[3];
  for(int i=0; i<3; i++) tryW[i]=bestW[i];
  
  double min=0;
  min = myLinearFunc(bestW[0], bestW[1], bestW[2]);

  std::cout << "first min: " << min << std::endl;

  for(int i=0; i<nEpoch; i++){
    TFitter* minimizer = new TFitter(10);
    double p1 = -1; minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);   
    minimizer->SetFCN(minuitLinearFunction);
    minimizer->SetParameter(0,"W0",tryW[0],1,0,0);
    minimizer->SetParameter(1,"W1",tryW[1],1,0,0);
    minimizer->SetParameter(2,"W2",tryW[2],1,0,0);

    //minimizer->ExecuteCommand("SIMPLEX",0,0);
    //minimizer->ExecuteCommand("MIGRAD",0,0);
    //
    minimizer->ExecuteCommand("MINOS",0,0);
    
    double newW[3];
    newW[0] = minimizer->GetParameter(0); 
    newW[1] = minimizer->GetParameter(1); 
    newW[2] = minimizer->GetParameter(2);
    
    double newmin = myLinearFunc(newW[0],newW[1],newW[2]);
    std::cout << "iteration: " << i << ", new min: " << newmin << std::endl;
    if(newmin<min){
      min = newmin;
      for(int j=0; j<3; j++){
	bestW[j]=newW[j];
      }
    }
    for(int j=0; j<3; j++){
      tryW[j]=newW[j];
    }
  }
  std::cout << min << std::endl;
  for(int j=0; j<3; j++){
    char out[200];
    sprintf(out,"%s%d%s","bestW[",j,"] = ");
    std::cout << out << bestW[j] << "; ";
    std::cout << "" << std::endl;
    BestW[j]=bestW[j];
  }
}

TH1D* LinearSimEreco::FindLinearErec(int ebeam)
{
  char hName[100];
  sprintf(hName,"%s%d%s","Erec_",ebeam,"GeV");
  TH1D* h=new TH1D(hName,hName,300,0,150);
  unsigned int nevent=nhit.size();
  float ERECO;
  for(unsigned int j=0; j<nevent; j++){
    ERECO=nhit1.at(j)*BestW[0] + nhit2.at(j)*BestW[1] + nhit3.at(j)*BestW[2];
    h->Fill(ERECO);
  }
  if( ebeam==5 )
    writeHisto(h,std::string("./histo/LinearErec.root"),"RECREATE");
  else 
    writeHisto(h,std::string("./histo/LinearErec.root"),"UPDATE");

  TF1* func=new TF1("fGaus","gaus",0,150);
  h->Fit(func,"NQ");
  emean=func->GetParameter(1);
  eres=func->GetParameter(2);
  err_emean=func->GetParError(1);
  err_eres=sqrt( pow(func->GetParError(2)/func->GetParameter(1),2) + pow(func->GetParError(1)*func->GetParameter(2)/func->GetParameter(1)*func->GetParameter(1),2) );
  return h;
}
void Process()
{
  n0.clear();
  n1.clear();
  n2.clear();
  n3.clear();
  energy.clear();
  int Ebeam[]={5,10,15,20,25,30,40,50,60,70,80};
  char input[200];
  for(unsigned int i=0; i<sizeof(Ebeam)/sizeof(int); i++){
    sprintf(input,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/ftfp_bert/ild_digit/single_pi-_",Ebeam[i],"GeV.root");
    std::cout << input << std::endl;
    LinearSimEreco *minim = new LinearSimEreco(input);
    minim->Loop();
    unsigned int nevent=minim->nhit.size();
    for(unsigned int j=0; j<nevent; j++){
      n0.push_back(minim->nhit.at(j));
      n1.push_back(minim->nhit1.at(j));
      n2.push_back(minim->nhit2.at(j));
      n3.push_back(minim->nhit3.at(j));
      energy.push_back(Ebeam[i]);
    }
  }
  std::cout << n0.size() << std::endl;
  LinearMinimize(20);

  int EbeamR[]={5,10,15,20,25,30,40,50,60,70,80};
  for(unsigned int i=0; i<sizeof(EbeamR)/sizeof(int); i++){
    sprintf(input,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/qgsp_bert/ild_digit/single_pi-_",EbeamR[i],"GeV.root");
    LinearSimEreco *ereco = new LinearSimEreco(input);
    ereco->Loop();
    TH1D* h=ereco->FindLinearErec(EbeamR[i]);
    delete h;
  }
}
