#define Minimisation_cxx
#include "Minimisation.h"
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

void Minimisation::Loop()
{
  if (fChain == 0) return;
  int nEvent=0;
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(float(Nhit)/Nlayer>3&&
       Begin>=0&&
       float(NInteractinglayer)/Nlayer>0.2//&&
       //Begin<=5
       ){
    nhit.push_back(Nhit);
    nhit1.push_back(Nhit1);
    nhit2.push_back(Nhit2);
    nhit3.push_back(Nhit3);
    hough1.push_back(0);
    hough2.push_back(0);
    hough3.push_back(0);
    time.push_back(spillEventTime);
    nEvent++;
    }
    // if(nEvent>=1000) return;
  }
}

void Minimisation::Erec(int ebeam)
{
  if (fChain == 0) return;
  TFile *resfile;
  if(ebeam==5)
    resfile=new TFile("resultQua.root","RECREATE");
  else
    resfile=new TFile("resultQua.root","UPDATE");
  char output[200];
  sprintf(output,"%s%d%s","h_",ebeam,"GeV");
  TH1D *h=new TH1D(output,"",100,0,(int)2*ebeam);

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(float(Nhit)/Nlayer>3&&
       Begin>=0&&
       float(NInteractinglayer)/Nlayer>0.2&&
       Begin<=10){
    //LongiProfileBis[47]<=10&&LongiProfileBis[46]<=10){
      h->Fill((BestW[0] + BestW[1]*Nhit + BestW[2]*Nhit*Nhit)*(Nhit1-Nhough1)+
	      (BestW[3] + BestW[4]*Nhit + BestW[5]*Nhit*Nhit)*(Nhit2-Nhough2)+
	      (BestW[6] + BestW[7]*Nhit + BestW[8]*Nhit*Nhit)*(Nhit3-Nhough3)+
	      BestW[9]*(Nhough1+Nhough2+Nhough3));
    }
  }
  h->Write();
  TCanvas *cc=new TCanvas();
  cc->cd();
  h->Draw();
  //  cc->WaitPrimitive();
  resfile->Write();
  TF1 *f=new TF1("f","gaus",0,120);
  h->Draw();
  h->Fit(f,"","",h->GetMean()-2*h->GetRMS(),h->GetMean()+2*h->GetRMS());
  std::cout << "Ebeam = " << ebeam << "\t Resolution = " << f->GetParameter(2)/f->GetParameter(1) << std::endl;
  emean=f->GetParameter(1);
  eres=f->GetParameter(2);
  err_emean=f->GetParError(1);
  err_eres=f->GetParError(2);
}

double myFunc(double w0, double w1, double w2, 
	      double w3, double w4, double w5,
	      double w6, double w7, double w8,
	      double w9)
{
//std::cout << w0 << "\t" << w1 << "\t" << w2 << "\n"
//	    << w3 << "\t" << w4 << "\t" << w5 << "\n"
//	    << w6 << "\t" << w7 << "\t" << w8 << std::endl;

  unsigned int size = n0.size();
  float test = 0;
  float energy_reco;
  for(unsigned int i=0; i<size; i++){
    energy_reco = ((w0 + w1*n0[i] + w2*n0[i]*n0[i])*(n1[i]-ho1[i])+
		   (w3 + w4*n0[i] + w5*n0[i]*n0[i])*(n2[i]-ho2[i])+
		   (w6 + w7*n0[i] + w8*n0[i]*n0[i])*(n3[i]-ho3[i])+
		   w9*(ho1[i]+ho2[i]+ho3[i]));
		   
    test = test + (
		   (energy[i]-energy_reco)*
		   (energy[i]-energy_reco)/energy[i]
		   );
    //if(energy_reco<0) std::cout << "Bug: Negative Reconstructed Energy " << std::endl;
    //    std::cout << energy[i] << " " << energy_reco << std::endl;
  }
  return test/(size-1);
}

void minuitFunction(int& nDim, double* gout,double& result, 
		    double *par,int flg){
  result = myFunc(par[0], par[1], par[2], 
		  par[3], par[4], par[5], 
		  par[6], par[7], par[8],
		  par[9]);
}


void Minimize(int nEpoch)
{
  double bestW[10];
  bestW[0] = 0.0395046; bestW[1] = 3.84302e-05; bestW[2] =-2.30225e-08; 
  bestW[3] = 0.0781585; bestW[4] =-5.96698e-05; bestW[5] =-7.04039e-09; 
  bestW[6] = 0.1152960; bestW[7] = 1.59277e-05; bestW[8] = 2.39387e-09;  
  bestW[9] = 0.04499;
  double tryW[10];
  for(int i=0; i<10; i++) tryW[i]=bestW[i];
  
  double min=0;
  min = myFunc(bestW[0], bestW[1], bestW[2], 
	       bestW[3], bestW[4], bestW[5], 
	       bestW[6], bestW[7], bestW[8],
	       bestW[9]);

  std::cout << "first min: " << min << std::endl;

  for(int i=0; i<nEpoch; i++){
    TFitter* minimizer = new TFitter(10);
    double p1 = -1; minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);   
    minimizer->SetFCN(minuitFunction);
    minimizer->SetParameter(0,"W0",tryW[0],.1,0,0);
    minimizer->SetParameter(1,"W1",tryW[1],.1,0,0);
    minimizer->SetParameter(2,"W2",tryW[2],.1,0,0);
    minimizer->SetParameter(3,"W3",tryW[3],.1,0,0);
    minimizer->SetParameter(4,"W4",tryW[4],.1,0,0);
    minimizer->SetParameter(5,"W5",tryW[5],.1,0,0);
    minimizer->SetParameter(6,"W6",tryW[6],.1,0,0);  
    minimizer->SetParameter(7,"W7",tryW[7],.1,0,0);
    minimizer->SetParameter(8,"W8",tryW[8],.1,0,0);
    minimizer->SetParameter(9,"W9",tryW[9],.1,0,0);

    
    minimizer->ExecuteCommand("SIMPLEX",0,0);
    minimizer->ExecuteCommand("MIGRAD",0,0);
    //
    //    minimizer->ExecuteCommand("MINOS",0,0);
    
    double newW[10];
    newW[0] = minimizer->GetParameter(0); 
    newW[1] = minimizer->GetParameter(1); 
    newW[2] = minimizer->GetParameter(2);
    newW[3] = minimizer->GetParameter(3); 
    newW[4] = minimizer->GetParameter(4); 
    newW[5] = minimizer->GetParameter(5);
    newW[6] = minimizer->GetParameter(6);
    newW[7] = minimizer->GetParameter(7);
    newW[8] = minimizer->GetParameter(8);
    newW[9] = minimizer->GetParameter(9);
    
    double newmin = myFunc(newW[0],newW[1],newW[2],
			   newW[3],newW[4],newW[5],
			   newW[6],newW[7],newW[8],
			   newW[9]);
    std::cout << "iteration: " << i << ", new min: " << newmin << std::endl;
    if(newmin<min){
      min = newmin;
      for(int j=0; j<10; j++){
	bestW[j]=newW[j];
      }
    }
    for(int j=0; j<10; j++){
      tryW[j]=newW[j];
    }
  }
  std::cout << min << std::endl;
  for(int j=0; j<10; j++){
    char out[200];
    sprintf(out,"%s%d%s","bestW[",j,"] = ");
    if(j==3 || j==6 || j==9) std::cout << "" << std::endl;;
    std::cout << out << bestW[j] << "; ";
    if(j==9) std::cout << "" << std::endl;
    BestW[j]=bestW[j];
  }
}
void Process()
{
  n0.clear();
  n1.clear();
  n2.clear();
  n3.clear();
  energy.clear();
  int Ebeam[]={5,10,15,20,25,30,40,50,60,70,80};
  //int Ebeam[]={5,10,15,20,25,30,40};
  //  energy={Ebeam,Ebeam+sizeof(Ebeam)/sizeof(int)};
  char input[200];
  for(unsigned int i=0; i<sizeof(Ebeam)/sizeof(int); i++){
    sprintf(input,"%s%d%s","/home/steen/v02-06_analysis/sim_data/single_pi-_",Ebeam[i],"GeV.root");
    //    sprintf(input,"%s%d%s","/home/steen/v02-06_analysis/sim_data/single_pi-_",Ebeam[i],"GeV_ideal.root");
    printf(input);
    std::cout << " " << std::endl;
    Minimisation *m = new Minimisation(input);
    m->Loop();
    unsigned int nevent=m->nhit.size();
    for(unsigned int j=0; j<nevent; j++){
      n0.push_back(m->nhit[j]);
      n1.push_back(m->nhit1[j]);
      n2.push_back(m->nhit2[j]);
      n3.push_back(m->nhit3[j]);
      ho1.push_back(m->hough1[j]);
      ho2.push_back(m->hough2[j]);
      ho3.push_back(m->hough3[j]);
      energy.push_back(Ebeam[i]);
    }
  }
  std::cout << n0.size() << std::endl;
  Minimize(1);

  TFile *file=new TFile("result.root","recreate");
  //int Energy[]={5,10,15,20,25,30,40,50,60,70,80};
  ////int size=sizeof(Energy)/sizeof(int);
  TGraph *erec=new TGraph();
  erec->SetMarkerStyle(20);
  erec->SetMarkerSize(0.8);
  TGraphErrors *eres=new TGraphErrors();
  eres->SetMarkerStyle(20);
  eres->SetMarkerSize(0.8);
  for(unsigned int i=0; i<sizeof(Ebeam)/sizeof(int); i++){
    sprintf(input,"%s%d%s","/home/steen/v02-06_analysis/sim_data/single_pi-_",Ebeam[i],"GeV.root");
    //sprintf(input,"%s%d%s","/home/steen/v02-06_analysis/sim_data/single_pi-_",Ebeam[i],"GeV_ideal.root");
    printf(input);
    std::cout << " " << std::endl;
    Minimisation *n = new Minimisation(input);
    n->Erec(Ebeam[i]);
    erec->SetPoint(i,Ebeam[i],n->emean);
    eres->SetPoint(i,Ebeam[i],n->eres/n->emean);
    eres->SetPointError(i,0,sqrt((n->err_eres/n->emean)*(n->err_eres/n->emean) + (n->err_emean*n->eres/(n->emean*n->emean))*(n->err_emean*n->eres/(n->emean*n->emean))) );
  }
  TH2D *h1=new TH2D("h1","",100,0,100,100,0,100);
  TH2D *h2=new TH2D("h2","",100,0,100,100,0,0.50);
  TCanvas *can=new TCanvas();
  can->SetWindowSize(1000,600);
  can->Divide(2,1);
  can->cd(1);
  TF1 *lin1=new TF1("lin1","[0]*x",0,100);
  h1->Draw();
  erec->Draw("p");
  erec->Fit(lin1,"","",0,75);
  
  
  can->cd(2);
  //TF1 *func = new TF1("func",myFunc,0,50,3);
  //func->SetLineColor(1);
  //func->SetLineWidth(2);
  //func->SetParameter(0,0.6);
  //func->SetParameter(1,0.06);
  //func->SetParameter(2,0.);
  h2->Draw();
  eres->Draw("p");
  //eres->Fit(func);

  
  
  //for(unsigned int i=0; i<sizeof(Ebeam)/sizeof(int); i++){
  //  sprintf(input,"%s%d%s","/home/arnaud/Work/TB_Data/SPS_2012/RootFiles/DHCAL_",Run(Ebeam[i]),".root");
  //  printf(input);
  //  std::cout << " " << std::endl;
  //  Minimisation *n = new Minimisation(input);
  //  n->Erec(Ebeam[i]);
  //}
  
  //file->Close();
}
