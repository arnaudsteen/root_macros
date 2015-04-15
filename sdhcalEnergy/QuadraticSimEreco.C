#define QuadraticSimEreco_cxx
#include "QuadraticSimEreco.h"
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

void QuadraticSimEreco::Loop(bool withMax)
{
  if (fChain == 0) return;
  int ncount=0;
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if( (float)Nhit/Nlayer>3 && Begin>=0 && (float)NInteractinglayer/Nlayer>0.01 ){
      nhit.push_back(Nhit);
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      ncount++;
    }
    if(withMax&&ncount>2000)break;
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

double myQuadraticFunc(double w0, double w1, double w2, 
		    double w3, double w4, double w5,
		    double w6, double w7, double w8)
{
  unsigned int size = n0.size();
  float test = 0;
  float energy_reco;
  for(unsigned int i=0; i<size; i++){
    energy_reco = ((w0 + w1*n0[i] + w2*n0[i]*n0[i])*n1[i]+
		   (w3 + w4*n0[i] + w5*n0[i]*n0[i])*n2[i]+
		   (w6 + w7*n0[i] + w8*n0[i]*n0[i])*n3[i]);
		   
    test+=fabs((energy[i]-energy_reco)/energy[i]);
    //test = test + ( (energy[i]-energy_reco)*
    //		    (energy[i]-energy_reco)/energy[i]/energy[i] );
  }
  return test/size;
}

void minuitQuadraticFunction(int& nDim, double* gout,double& result, 
			  double *par,int flg){
  result = myQuadraticFunc(par[0], par[1], par[2], 
			   par[3], par[4], par[5], 
			   par[6], par[7], par[8]);
}

//void QuadraticMinimize(int nEpoch)
//{
//  double bestW[9];
//  bestW[0] = 0.0395046; bestW[1] = 3.84302e-05; bestW[2] =-2.30225e-08; 
//  bestW[3] = 0.0781585; bestW[4] =-5.96698e-05; bestW[5] =-7.04039e-09; 
//  bestW[6] = 0.1152960; bestW[7] = 1.59277e-05; bestW[8] = 2.39387e-09;  
//  double tryW[9];
//  for(int i=0; i<9; i++) tryW[i]=bestW[i];
//  
//  double min=0;
//  min = myFunc(bestW[0], bestW[1], bestW[2], 
//	       bestW[3], bestW[4], bestW[5], 
//	       bestW[6], bestW[7], bestW[8]);
//
//  std::cout << "first min: " << min << std::endl;
//
//  for(int i=0; i<nEpoch; i++){
//    TFitter* minimizer = new TFitter(10);
//    double p1 = -1; minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);   
//    minimizer->SetFCN(minuitFunction);
//    minimizer->SetParameter(0,"W0",tryW[0],.1,0,0);
//    minimizer->SetParameter(1,"W1",tryW[1],.1,0,0);
//    minimizer->SetParameter(2,"W2",tryW[2],.1,0,0);
//    minimizer->SetParameter(3,"W3",tryW[3],.1,0,0);
//    minimizer->SetParameter(4,"W4",tryW[4],.1,0,0);
//    minimizer->SetParameter(5,"W5",tryW[5],.1,0,0);
//    minimizer->SetParameter(6,"W6",tryW[6],.1,0,0);  
//    minimizer->SetParameter(7,"W7",tryW[7],.1,0,0);
//    minimizer->SetParameter(8,"W8",tryW[8],.0000000001,0,0);
//
//    
//    minimizer->ExecuteCommand("SIMPLEX",0,0);
//    minimizer->ExecuteCommand("MIGRAD",0,0);
//    //
//    //    minimizer->ExecuteCommand("MINOS",0,0);
//    
//    double newW[10];
//    newW[0] = minimizer->GetParameter(0); 
//    newW[1] = minimizer->GetParameter(1); 
//    newW[2] = minimizer->GetParameter(2);
//    newW[3] = minimizer->GetParameter(3); 
//    newW[4] = minimizer->GetParameter(4); 
//    newW[5] = minimizer->GetParameter(5);
//    newW[6] = minimizer->GetParameter(6);
//    newW[7] = minimizer->GetParameter(7);
//    newW[8] = minimizer->GetParameter(8);
//    
//    double newmin = myFunc(newW[0],newW[1],newW[2],
//			   newW[3],newW[4],newW[5],
//			   newW[6],newW[7],newW[8]);
//    std::cout << "iteration: " << i << ", new min: " << newmin << std::endl;
//    if(newmin<min){
//      min = newmin;
//      for(int j=0; j<9; j++){
//	bestW[j]=newW[j];
//      }
//    }
//    for(int j=0; j<9; j++){
//      tryW[j]=newW[j];
//    }
//  }
//  std::cout << min << std::endl;
//  for(int j=0; j<9; j++){
//    char out[200];
//    sprintf(out,"%s%d%s","bestW[",j,"] = ");
//    if(j==3 || j==6 ) std::cout << "" << std::endl;;
//    std::cout << out << bestW[j] << "; ";
//    if(j==8) std::cout << "" << std::endl;
//    BestW[j]=bestW[j];
//  }
//}

void QuadraticMinimize(int nEpoch)
{
  double bestW[9];
  bestW[0] = 0.0433517; 
  bestW[1] = 2.32656e-05; 
  bestW[2] = 3.09656e-08; 
  bestW[3] = 0.0722234; 
  bestW[4] = 1.81082e-05; 
  bestW[5] = 9.14195e-09; 
  bestW[6] = 0.224; 
  bestW[7] = 0.000187388; 
  bestW[8] = 7.71522e-08; 
  //bestW[0] = 0.0432348; 
  //bestW[1] = 3.6678e-05; 
  //bestW[2] = -3.90726e-08; 
  //bestW[3] = 0.0693766; 
  //bestW[4] = 2.18316e-05; 
  //bestW[5] = -1.32885e-08; 
  //bestW[6] = 0.175724; 
  //bestW[7] = 0.000151068; 
  //bestW[8] = 1.61142e-07; 
  
  //bestW[0] = 0.0395046; bestW[1] = 3.84302e-05; bestW[2] = 2.30225e-08; 
  //bestW[3] = 0.0781585; bestW[4] = 2.96698e-05; bestW[5] = 2.04039e-09; 
  //bestW[6] = 0.1152960; bestW[7] = 1.59277e-05; bestW[8] = 2.39387e-09;  

  double tryW[9];
  for(int i=0; i<9; i++) tryW[i]=bestW[i];
  
  double min=0;
  min = myQuadraticFunc(bestW[0], bestW[1], bestW[2], 
			bestW[3], bestW[4], bestW[5], 
			bestW[6], bestW[7], bestW[8]);

  std::cout << "first min: " << min << std::endl;

  for(int i=0; i<nEpoch; i++){
    TFitter* minimizer = new TFitter(10);
    double p1 = -1; minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);   
    minimizer->SetFCN(minuitQuadraticFunction);
    minimizer->SetParameter(0,"W0",tryW[0],.1,0,0);
    minimizer->SetParameter(1,"W1",tryW[1],.1,0,0);
    minimizer->SetParameter(2,"W2",tryW[2],.1,0,0);
    minimizer->SetParameter(3,"W3",tryW[3],.1,0,0);
    minimizer->SetParameter(4,"W4",tryW[4],.1,0,0);
    minimizer->SetParameter(5,"W5",tryW[5],.1,0,0);
    minimizer->SetParameter(6,"W6",tryW[6],.1,0,0);  
    minimizer->SetParameter(7,"W7",tryW[7],.1,0,0);
    minimizer->SetParameter(8,"W8",tryW[8],.0000001,0,0);

    //minimizer->ExecuteCommand("SIMPLEX",0,0);
    minimizer->ExecuteCommand("MIGRAD",0,0);
    //
    //    minimizer->ExecuteCommand("MINOS",0,0);
    
    double newW[9];
    newW[0] = minimizer->GetParameter(0); 
    newW[1] = minimizer->GetParameter(1); 
    newW[2] = minimizer->GetParameter(2);
    newW[3] = minimizer->GetParameter(3); 
    newW[4] = minimizer->GetParameter(4); 
    newW[5] = minimizer->GetParameter(5);
    newW[6] = minimizer->GetParameter(6);
    newW[7] = minimizer->GetParameter(7);
    newW[8] = minimizer->GetParameter(8);
    
    double newmin = myQuadraticFunc(newW[0], newW[1], newW[2], 
				    newW[3], newW[4], newW[5], 
				    newW[6], newW[7], newW[8]);

    std::cout << "iteration: " << i << ", new min: " << newmin << std::endl;
    if(newmin<min){
      min = newmin;
      for(int j=0; j<9; j++){
	bestW[j]=newW[j];
      }
    }
    for(int j=0; j<9; j++){
      tryW[j]=newW[j];
    }
  }
  std::cout << min << std::endl;
  for(int j=0; j<9; j++){
    char out[200];
    sprintf(out,"%s%d%s","bestW[",j,"] = ");
    std::cout << out << bestW[j] << "; ";
    std::cout << "" << std::endl;
    BestW[j]=bestW[j];
  }
}

TH1D* QuadraticSimEreco::FindQuadraticErec(int ebeam)
{
  char hName[100];
  sprintf(hName,"%s%d%s","Erec_",ebeam,"GeV");
  TH1D* h=new TH1D(hName,hName,300,0,150);
  unsigned int nevent=nhit.size();
  float ERECO;
  for(unsigned int j=0; j<nevent; j++){
    ERECO=nhit1.at(j)*(BestW[0] + BestW[1]*nhit.at(j) + BestW[2]*nhit.at(j)*nhit.at(j))+
      nhit2.at(j)*(BestW[3] + BestW[4]*nhit.at(j) + BestW[5]*nhit.at(j)*nhit.at(j))+
      nhit3.at(j)*(BestW[6] + BestW[7]*nhit.at(j) + BestW[8]*nhit.at(j)*nhit.at(j));
    h->Fill(ERECO);
  }
  if( ebeam==5 )
    writeHisto(h,std::string("./histo/QuadraticErec.root"),"RECREATE");
  else 
    writeHisto(h,std::string("./histo/QuadraticErec.root"),"UPDATE");

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
    QuadraticSimEreco *minim = new QuadraticSimEreco(input);
    minim->Loop(0);
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
  QuadraticMinimize(14);

  int EbeamR[]={5,10,15,20,25,30,40,50,60,70,80};
  for(unsigned int i=0; i<sizeof(EbeamR)/sizeof(int); i++){
    sprintf(input,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/qgsp_bert/ild_digit/single_pi-_",EbeamR[i],"GeV.root");
    std::cout << input << std::endl;
    QuadraticSimEreco *ereco = new QuadraticSimEreco(input);
    ereco->Loop(0);
    TH1D* h=ereco->FindQuadraticErec(EbeamR[i]);
    delete h;
  }
}
