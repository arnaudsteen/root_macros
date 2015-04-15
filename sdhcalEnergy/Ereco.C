#define Ereco_cxx
#include "Ereco.h"
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

int Run(double E,std::string period)
{
  if(period==std::string("AugSep2012")){
    if(E==80) return 715756;
    else if(E==70) return 715754;
    else if(E==60) return 715753;
    else if(E==50) return 715751;
    else if(E==40) return 715748;
    else if(E==30) return 715747;
    else if(E==25) return 715703;
    else if(E==20) return 715675;
    else if(E==15) return 715699;
    else if(E==10) return 715693;
    else if(E==5) return 715694;
    else {
      std::cout << "no run number at beam energy = " << E << " GeV for the period " << period << std::endl;
      return 0;
    }
  }
  else if(period==std::string("Nov2012")){
    if(E==80) return 716319;
    else if(E==70) return 716290;
    else if(E==60) return 716298;
    else if(E==50) return 716305;
    else if(E==40) return 716307;
    else if(E==30) return 716308;
    else if(E==20) return 716315;
    else if(E==10) return 716321;
    else {
      std::cout << "no run number at beam energy = " << E << " GeV for the period " << period << std::endl;
      return 0;
    }
  }
  std::cout << period << " is not a good period" << std::endl;
  return 0;
}

void Ereco::ShowerBarycenterCut(std::string runPeriod)
{
  if(runPeriod==std::string("AugSep2012")){
    cutBary[0]=490;
    cutBary[1]=690;
    cutBary[2]=400;
    cutBary[3]=600;
  }
  if(runPeriod==std::string("Nov2012")){
    cutBary[0]=490;
    cutBary[1]=690;
    cutBary[2]=400;
    cutBary[3]=600;
  }
  if(runPeriod==std::string("Dec2014")){
    cutBary[0]=490;
    cutBary[1]=690;
    cutBary[2]=400;
    cutBary[3]=600;
  }
}


void Ereco::Loop()
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
       (float)NInteractinglayer/Nlayer>0.2&&
       Single==1 &&
       Neutral==0 &&
       Hole==0 &&
       (float)spillEventTime*200/pow(10.,9)<15 &&
       CoG[0]>cutBary[0] &&
       CoG[0]<cutBary[1] &&
       CoG[2]>cutBary[2] &&
       CoG[2]<cutBary[3] &&
       ReconstructedCosTheta>.9&&
       TrackMultiplicity>0
       //TrackClusterNumber->at(0)>5&&
       //TransverseRatio>0.01&&
       //(float)(Nhough1+Nhough2+Nhough3)/Nhit<0.5) 
	){
      nhit.push_back(Nhit);
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      time.push_back(spillEventTime*200/pow(10.,9));
    }
  }
}

void Ereco::Loop(int NEVENT)
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
    if( 
       (float)Nhit/Nlayer>3 && 
       Begin>=0 && 
       (float)NInteractinglayer/Nlayer>0.2&&
       Single==1 &&
       Neutral==0 &&
       Hole==0 &&
       (float)spillEventTime*200/pow(10.,9)<15 &&
       CoG[0]>cutBary[0] &&
       CoG[0]<cutBary[1] &&
       CoG[2]>cutBary[2] &&
       CoG[2]<cutBary[3] &&
       ReconstructedCosTheta>.9&&
       TrackMultiplicity>0
       //TrackClusterNumber->at(0)>5&&
       //TransverseRatio>0.01&&
       //(float)(Nhough1+Nhough2+Nhough3)/Nhit<0.5) 
	){
      nhit.push_back(Nhit);
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      time.push_back(spillEventTime*200/pow(10.,9));
      ncount++;
    }
    if(ncount>=NEVENT) break;
  }
}

//void Ereco::Erec(int ebeam)
//{
//  if (fChain == 0) return;
//  TFile *resfile;
//  if(ebeam==5)
//    resfile=new TFile("resultQua.root","RECREATE");
//  else
//    resfile=new TFile("resultQua.root","UPDATE");
//  char output[200];
//  sprintf(output,"%s%d%s","h_",ebeam,"GeV");
//  TH1D *h=new TH1D(output,"",100,0,(int)2*ebeam);
//
//  Long64_t nentries = fChain->GetEntriesFast();
//  
//  Long64_t nbytes = 0, nb = 0;
//  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//    Long64_t ientry = LoadTree(jentry);
//    if (ientry < 0) break;
//    nb = fChain->GetEntry(jentry);   nbytes += nb;
//    // if (Cut(ientry) < 0) continue;
//    if( 
//       (((float)Nhit/Nlayer>3 && Begin>=0 && (float)NInteractinglayer/Nlayer>0.2)||(Nlayer<40))&&
//       Single==1 &&
//       Neutral==0 &&
//       Hole==0 &&
//       (float)spillEventTime*200/pow(10.,9)<15 &&
//       CoG[0]>cutBary[0] &&
//       CoG[0]<cutBary[1] &&
//       CoG[2]>cutBary[2] &&
//       CoG[2]<cutBary[3] &&
//       ReconstructedCosTheta>.9&&
//       TrackMultiplicity>0&&
//       TrackClusterNumber->at(0)>5&&
//       TransverseRatio>0.01&&
//       (float)(Nhough1+Nhough2+Nhough3)/Nhit<0.5) {
//      h->Fill((BestW[0] + BestW[1]*Nhit + BestW[2]*Nhit*Nhit)*Nhit1+
//	      (BestW[3] + BestW[4]*Nhit + BestW[5]*Nhit*Nhit)*Nhit2+
//	      (BestW[6] + BestW[7]*Nhit + BestW[8]*Nhit*Nhit)*Nhit3);
//    }
//  }
//  h->Write();
//  TCanvas *cc=new TCanvas();
//  cc->cd();
//  h->Draw();
//  resfile->Write();
//  TF1 *f=new TF1("f","gaus",0,120);
//  h->Draw();
//  h->Fit(f,"","",h->GetMean()-2*h->GetRMS(),h->GetMean()+2*h->GetRMS());
//  std::cout << "Ebeam = " << ebeam << "\t Resolution = " << f->GetParameter(2)/f->GetParameter(1) << std::endl;
//  emean=f->GetParameter(1);
//  eres=f->GetParameter(2);
//  err_emean=f->GetParError(1);
//  err_eres=f->GetParError(2);
//}

double myFunc(double w0, double w1, double w2, 
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
		   
    test = test + (
		   (energy[i]-energy_reco)*
		   (energy[i]-energy_reco)/energy[i]
		   );
  }
  return test/(size-1);
}

void minuitFunction(int& nDim, double* gout,double& result, 
		    double *par,int flg){
  result = myFunc(par[0], par[1], par[2], 
		  par[3], par[4], par[5], 
		  par[6], par[7], par[8]);
}


void Minimize(int nEpoch)
{
  double bestW[9];
  bestW[0] = 0.0395046; bestW[1] = 3.84302e-05; bestW[2] =-2.30225e-08; 
  bestW[3] = 0.0781585; bestW[4] =-5.96698e-05; bestW[5] =-7.04039e-09; 
  bestW[6] = 0.1152960; bestW[7] = 1.59277e-05; bestW[8] = 2.39387e-09;  
  double tryW[9];
  for(int i=0; i<9; i++) tryW[i]=bestW[i];
  
  double min=0;
  min = myFunc(bestW[0], bestW[1], bestW[2], 
	       bestW[3], bestW[4], bestW[5], 
	       bestW[6], bestW[7], bestW[8]);

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
    minimizer->SetParameter(8,"W8",tryW[8],.0000000001,0,0);

    
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
    
    double newmin = myFunc(newW[0],newW[1],newW[2],
			   newW[3],newW[4],newW[5],
			   newW[6],newW[7],newW[8]);
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
    if(j==3 || j==6 ) std::cout << "" << std::endl;;
    std::cout << out << bestW[j] << "; ";
    if(j==8) std::cout << "" << std::endl;
    BestW[j]=bestW[j];
  }
}

std::vector<float> readCalib(int run)
{
  std::vector<float> calib;
  char inputCalib[200];
  sprintf(inputCalib,"%s%d%s","/home/steen/timeCalib/Nhit/calib_",run,".txt");
  fstream in;
  in.open(inputCalib);
  if(!in.is_open()){ std::cout << inputCalib << "\t NO SUCH FILE OR DIRECTORY" << std::endl; throw;}
  float coeff[3];
  float coeffError[3];
  while(1){
    if(!in.good()) break;
    in >> coeff[0] >> coeffError[0] >> coeff[1] >> coeffError[1] >> coeff[2] >> coeffError[2] ;
    calib.push_back(coeff[0]);
    calib.push_back(coeff[1]);
    calib.push_back(coeff[2]);
    //calib.push_back(0);
    //calib.push_back(0);
    //calib.push_back(0);
  }
  return calib;
}

TH1D* Ereco::FindErec(int ebeam)
{
  TH1D* h=new TH1D("h","",300,0,150);
  std::vector<float> calibrations=readCalib(Run(ebeam,std::string("AugSep2012")));
  float calib1[3];
  float calib2[3];
  float calib3[3];
  for(unsigned int j=0; j<3; j++){
    calib1[j]=calibrations.at(j);
    calib2[j]=calibrations.at(j+3);
    calib3[j]=calibrations.at(j+6);
  }
  unsigned int nevent=nhit.size();
  float nhitCalib;
  float nhit1Calib;
  float nhit2Calib;
  float nhit3Calib;
  float coeff1=calib1[1]+calib2[1]+calib3[1];
  float coeff2=calib1[2]+calib2[2]+calib3[2];
  float ERECO;
  for(unsigned int j=0; j<nevent; j++){
    nhitCalib=( nhit1.at(j)+nhit2.at(j)+nhit3.at(j) - coeff1*time.at(j) - coeff2*time.at(j)*time.at(j) );
    nhit1Calib=nhit1.at(j)-calib1[1]*time.at(j)-calib1[2]*time.at(j)*time.at(j);
    nhit2Calib=nhit2.at(j)-calib2[1]*time.at(j)-calib2[2]*time.at(j)*time.at(j);
    nhit3Calib=nhit3.at(j)-calib3[1]*time.at(j)-calib3[2]*time.at(j)*time.at(j);
    ERECO=nhit1Calib*(BestW[0] + BestW[1]*nhitCalib + BestW[2]*nhitCalib*nhitCalib)+
      nhit2Calib*(BestW[3] + BestW[4]*nhitCalib + BestW[5]*nhitCalib*nhitCalib)+
      nhit3Calib*(BestW[6] + BestW[7]*nhitCalib + BestW[8]*nhitCalib*nhitCalib);
    h->Fill(ERECO);
  }
  TCanvas* cc=new TCanvas();
  TF1* func=new TF1("fGaus","gaus",0,150);
  h->Fit(func,"NQ");
  h->Fit(func,"","",func->GetParameter(1)-1.5*func->GetParameter(2),func->GetParameter(1)+1.5*func->GetParameter(2));
  emean=func->GetParameter(1);
  eres=func->GetParameter(2);
  err_emean=func->GetParError(1);
  err_eres=sqrt( pow(func->GetParError(2)/func->GetParameter(1),2) + pow(func->GetParError(1)*func->GetParameter(2)/func->GetParameter(1)*func->GetParameter(1),2) );
  //h->Draw();
  //  cc->WaitPrimitive();
  return h;
}

void Process()
{
  int nmax=10000;
  n0.clear();
  n1.clear();
  n2.clear();
  n3.clear();
  energy.clear();
  //int Ebeam[]={5,10,15,20,25,30,40,50,60,70,80};
  int Ebeam[]={10,30,50,80};
  char input[200];
  for(unsigned int i=0; i<sizeof(Ebeam)/sizeof(int); i++){
    sprintf(input,"%s%d%s","/home/steen/resultRootFile/tb_data/DHCAL_",Run(Ebeam[i],std::string("AugSep2012")),".root");
    std::cout << input << std::endl;
    Ereco *minim = new Ereco(input);
    minim->ShowerBarycenterCut(std::string("AugSep2012"));
    minim->Loop(nmax);
    std::vector<float> calibrations=readCalib(Run(Ebeam[i],std::string("AugSep2012")));
    float calib1[3];
    float calib2[3];
    float calib3[3];
    for(unsigned int j=0; j<3; j++){
      calib1[j]=calibrations.at(j);
      calib2[j]=calibrations.at(j+3);
      calib3[j]=calibrations.at(j+6);
    }
    unsigned int nevent=minim->nhit.size();
    float nhitCalib;
    float coeff1=calib1[1]+calib2[1]+calib3[1];
    float coeff2=calib1[2]+calib2[2]+calib3[2];
    for(unsigned int j=0; j<nevent; j++){
      nhitCalib=( minim->nhit1.at(j)+minim->nhit2.at(j)+minim->nhit3.at(j) - 
		  coeff1*minim->time.at(j) -
		  coeff2*minim->time.at(j)*minim->time.at(j) );
      n0.push_back(nhitCalib);
      n1.push_back(minim->nhit1[j]-calib1[1]*minim->time.at(j)-calib1[2]*minim->time.at(j)*minim->time.at(j));
      n2.push_back(minim->nhit2[j]-calib2[1]*minim->time.at(j)-calib2[2]*minim->time.at(j)*minim->time.at(j));
      n3.push_back(minim->nhit3[j]-calib3[1]*minim->time.at(j)-calib3[2]*minim->time.at(j)*minim->time.at(j));
      energy.push_back(Ebeam[i]);
    }
    //    std::cout << std::accumulate(n0.begin(),n0.end(),0.)/n0.size() << std::endl;
    //    n0.clear();
  }
  std::cout << n0.size() << std::endl;
  Minimize(1);

  //TFile *file=new TFile("result.root","recreate");
  int EbeamR[]={5,10,15,20,25,30,40,50,60,70,80};
  TGraphErrors *erec=new TGraphErrors();
  erec->SetMarkerStyle(20);
  erec->SetMarkerSize(0.8);
  TGraphErrors *eres=new TGraphErrors();
  eres->SetMarkerStyle(20);
  eres->SetMarkerSize(0.8);
  for(unsigned int i=0; i<sizeof(EbeamR)/sizeof(int); i++){
    sprintf(input,"%s%d%s","/home/steen/resultRootFile/tb_data/DHCAL_",Run(EbeamR[i],std::string("AugSep2012")),".root");
    std::cout << input << std::endl;
    Ereco *ereco = new Ereco(input);
    ereco->ShowerBarycenterCut(std::string("AugSep2012"));
    ereco->Loop();
    TH1D* h=ereco->FindErec(EbeamR[i]);
    //h->Draw();
    //cc->WaitPrimitive();
    //ereco->Erec(EbeamR[i]);
    erec->SetPoint(i,EbeamR[i],ereco->emean);
    erec->SetPointError(i,0,ereco->err_emean);
    eres->SetPoint(i,EbeamR[i],ereco->eres/ereco->emean);
    eres->SetPointError(i,0,sqrt(ereco->err_eres/ereco->emean*ereco->err_eres/ereco->emean+ereco->err_emean*ereco->eres/ereco->emean*ereco->err_emean*ereco->eres/ereco->emean));
    //    std::cout << EbeamR[i] << " " << ereco->eres << std::endl;
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
  erec->Fit(lin1);
  
  can->cd(2);
  h2->Draw();
  //  eres->Draw("p");
  TF1 *resoF=new TF1("resoF","sqrt([0]*[0]/x +[1]*[1] +[2]*[2]/x**2)",0,100);
  resoF->SetParameters(0.6,0.05,0.01);
  eres->Draw("p");
  eres->Fit(resoF,"","",5,80);
  
  
}
