#include <iostream>
#include <TROOT.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <cmath>
#include <TMath.h>
#include <TRandom.h>

void DrawGausCosTheta(/*int nevent=10000,float mean=0.0,*/float sigma=0.1)
{
  int nevent=10000;
  float mean=0.0;
  TF1* gaus=new TF1("func","gaus",-1,1);
  gaus->SetParameters(1,mean,sigma);
  TF1* gaus2=new TF1("func2","gaus",-1,1);
  gaus2->SetParameters(1,mean,sigma*10);
  TH1D* hCosT=new TH1D("hcosT","",100,.9,1);
  hCosT->SetLineColor(kBlack);
  hCosT->SetLineWidth(2);
  for(int i=0; i<nevent;i++){
    double px=gaus->GetRandom();//+0.1*gaus2->GetRandom();
    double py=gaus->GetRandom();//+0.1*gaus2->GetRandom();
    double pz=1.;
    TVector3 temp(px,py,pz);
    TVector3 v(temp.x()/temp.Mag(),temp.y()/temp.Mag(),temp.z()/temp.Mag());
    hCosT->Fill(v.CosTheta());
  }
  TCanvas *cc=new TCanvas();
  cc->SetWindowSize(600,600);
  hCosT->Draw();
}

void DrawSolidAngleCosTheta(int nevent=10000,float solidAngleX0=2.0,float solidAngleRad=0.05)
{
  double R0 = std::sqrt(solidAngleRad*solidAngleRad/4+solidAngleRad*solidAngleRad/4);  
  double rndm1, rndm2;  
  double px, py, pz, projx, projy; 
  double MinTheta, MaxTheta, MinPhi, MaxPhi; 
  double Phi;
    
  MinTheta = 0.; 
  MaxTheta = std::atan(R0/solidAngleX0);   
  MinPhi = 0.; 
  MaxPhi = 2*TMath::Pi(); 		
  double sintheta, sinphi, costheta, cosphi, tantheta; 
  TRandom *rand=new TRandom();
  TH1D* hCosT=new TH1D("hcosT","",100,.8,1);
  hCosT->SetLineColor(kBlack);
  hCosT->SetLineWidth(2);
  for(int i=0; i<nevent; i++){
    do{
      rndm1 = rand->Uniform();  
      costheta = std::cos(MinTheta) - rndm1 * (std::cos(MinTheta) - std::cos(MaxTheta));
      sintheta = std::sqrt(1. - costheta*costheta);  
      tantheta = sintheta/costheta;  
      rndm2 = rand->Uniform();    
      Phi = MinPhi + (MaxPhi - MinPhi) * rndm2;  
      sinphi = std::sin(Phi); 
      cosphi = std::cos(Phi); 
      px = sintheta * cosphi;  
      py = sintheta * sinphi;   
      pz = costheta; 
      projx = solidAngleX0*tantheta*cosphi;  
      projy = solidAngleX0*tantheta*sinphi;   
    }while(sqrt(projx*projx)>500||sqrt(projy*projy)>500);
    TVector3 v(px,py,pz);
    hCosT->Fill(v.CosTheta());
  }
  TCanvas *cc=new TCanvas();
  cc->SetWindowSize(600,600);
  hCosT->Draw();
}
