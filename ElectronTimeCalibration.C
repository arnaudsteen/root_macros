#define ElectronTimeCalibration_cxx
#include "ElectronTimeCalibration.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstring>
#include <stdio.h>
#include <cmath>
#include <fstream>

int Run(double E,std::string period)
{
  if(period==std::string("AugSep2012")){
    if(E==50) return 715716;
    else if(E==40) return 715713;
    else if(E==30) return 715714;
    else if(E==20) return 715715;
    else if(E==10) return 715725;
    else {
      std::cout << "no run number at beam energy = " << E << " GeV for the period " << period << std::endl;
      return 0;
    }
  }
  std::cout << period << " is not a good period" << std::endl;
  return 0;
}

void ElectronTimeCalibration::ShowerBarycenterCut(std::string runPeriod)
{
  if(runPeriod==std::string("AugSep2012")){
    cutBary[0]=450;
    cutBary[1]=700;
    cutBary[2]=350;
    cutBary[3]=650;
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

void ElectronTimeCalibration::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //TH2D* htemp=new TH2D("htemp","",100,0,10,1500,0,1500);
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if( (float)Nhit/Nlayer>3 &&
	Begin>=0 &&
	Begin<5 &&
	Nlayer<30 &&
	TrackMultiplicity==0 &&
	Single==1 &&
	Neutral==0 &&
	(float)spillEventTime*200/pow(10.,9)<15 &&
	(float)spillEventTime*200/pow(10.,9)>0.3 &&
	CoG[0]>cutBary[0] &&
	CoG[0]<cutBary[1] &&
	CoG[2]>cutBary[2] &&
	CoG[2]<cutBary[3] &&
	ReconstructedCosTheta>.9){
      //htemp->Fill((float)spillEventTime*200/pow(10.,9),Nhit1+Nhit2+Nhit3);
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      nclusters.push_back(Nclusters);
      time.push_back( (float)spillEventTime*200/pow(10.,9) );
      radialProfile radPro;
      for(unsigned int i=0; i<96; i++)
      	radPro.ring[i]=RadialProfile[i];
      radialVec.push_back(radPro);
      longitudinalProfile longiPro;
      for(unsigned int i=0; i<48; i++)
      	longiPro.layer[i]=LongiProfile[i];
      longitudinalVec.push_back(longiPro);
    }      
  }
  //TCanvas *cc=new TCanvas();
  //htemp->Draw("colz");
  //cc->WaitPrimitive();
  //htemp->ProfileX()->Fit("pol2");
  //cc->WaitPrimitive();
  //delete htemp;
  //delete cc;
}

void NhitCalibration(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  ElectronTimeCalibration* timeCalib=new ElectronTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  TH2D* hN1=new TH2D("hN1","",100,0,10,10000,0,2000);
  TH2D* hN2=new TH2D("hN2","",100,0,10,5000,0,1000);
  TH2D* hN3=new TH2D("hN3","",100,0,10,2500,0,500);
  for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
    hN1->Fill(timeCalib->time.at(i),timeCalib->nhit1.at(i));
    hN2->Fill(timeCalib->time.at(i),timeCalib->nhit2.at(i));
    hN3->Fill(timeCalib->time.at(i),timeCalib->nhit3.at(i));
  }
  //TCanvas *cc=new TCanvas();
  TProfile* profile1=new TProfile();
  TProfile* profile2=new TProfile();
  TProfile* profile3=new TProfile();
  profile1=hN1->ProfileX();
  profile2=hN2->ProfileX();
  profile3=hN3->ProfileX();

  TF1 *f1=new TF1("f1","pol4",0,10);
  TF1 *f2=new TF1("f2","pol4",0,10);
  TF1 *f3=new TF1("f3","pol4",0,10);
  profile1->Fit(f1,"");
  //cc->WaitPrimitive();
  profile2->Fit(f2,"");
  //cc->WaitPrimitive();
  profile3->Fit(f3,"");
  //cc->WaitPrimitive();
  
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/Nhit/calib_",Run(E,runPeriod),".txt");
  std::cout << "outputFile = " << outputFile << std::endl;
  fstream out;
  out.open(outputFile,ios::out);
  out << f1->GetParameter(0) << " " << f1->GetParError(0) << " "
      << f1->GetParameter(1) << " " << f1->GetParError(1) << " "
      << f1->GetParameter(2) << " " << f1->GetParError(2) << " " 
      << f1->GetParameter(3) << " " << f1->GetParError(3) << " " 
      << f1->GetParameter(4) << " " << f1->GetParError(4) << std::endl;
  out << f2->GetParameter(0) << " " << f2->GetParError(0) << " "
      << f2->GetParameter(1) << " " << f2->GetParError(1) << " "
      << f2->GetParameter(2) << " " << f2->GetParError(2) << " " 
      << f2->GetParameter(3) << " " << f2->GetParError(3) << " " 
      << f2->GetParameter(4) << " " << f2->GetParError(4) << std::endl;
  out << f3->GetParameter(0) << " " << f3->GetParError(0) << " "
      << f3->GetParameter(1) << " " << f3->GetParError(1) << " "
      << f3->GetParameter(2) << " " << f3->GetParError(2) << " " 
      << f3->GetParameter(3) << " " << f3->GetParError(3) << " " 
      << f3->GetParameter(4) << " " << f3->GetParError(4) << std::endl;
  out.close();
  out.open("testChi2",ios::out|ios::app);
  out << E << "\t" << f1->GetChisquare()/f1->GetNDF() << "\t" << f2->GetChisquare()/f2->GetNDF() << "\t" << f3->GetChisquare()/f3->GetNDF() << std::endl;
  out.close();
}

void RadialProfile(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  ElectronTimeCalibration* timeCalib=new ElectronTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  std::cout<< timeCalib->radialVec.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/RadProfile/calib_",Run(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  for(unsigned int r=0; r<48; r++){
    TH2D* hRadPro=new TH2D("hName","",50,0,10,4000,0,200);
    func=new TF1("func","pol2",0,10);
    for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
      hRadPro->Fill(timeCalib->time.at(i),timeCalib->radialVec.at(i).ring[r]);
    }
    hRadPro->ProfileX()->Fit(func,"N");
    out << r 
	<< " " << func->GetParameter(0) << " " << func->GetParError(0) 
	<< " " << func->GetParameter(1) << " " << func->GetParError(1)
	<< " " << func->GetParameter(2) << " " << func->GetParError(2)
	<< " " << func->GetChisquare()/func->GetNDF() << std::endl;
    delete hRadPro;
  }
  out.close();
  delete func;
  delete timeCalib;
}

void NclustersCalibration(int E, std::string runPeriod)
{
  //TCanvas *can=new TCanvas();
  //can->SetWindowSize(600,600);
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  ElectronTimeCalibration* timeCalib=new ElectronTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  TH2D* hN1=new TH2D("hN1","",100,0,10,100,0,100);
  for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
    hN1->Fill(timeCalib->time.at(i),timeCalib->nclusters.at(i));
  }
  TProfile* profile1=new TProfile();
  profile1=hN1->ProfileX();
  profile1->Draw();
  TF1 *f1=new TF1("f1","pol1",0,10);
  profile1->Fit(f1,"");
  //can->Update();
  //can->WaitPrimitive();
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/Nclusters/calib_",Run(E,runPeriod),".txt");
  std::cout << "outputFile = " << outputFile << std::endl;
  fstream out;
  out.open(outputFile,ios::out);
  out << f1->GetParameter(0) << " " << f1->GetParError(0) << " "
      << f1->GetParameter(1) << " " << f1->GetParError(1) << std::endl;
  out.close();
}

void LongiProfile(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  ElectronTimeCalibration* timeCalib=new ElectronTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  std::cout<< timeCalib->longitudinalVec.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/LongiProfile/calib_",Run(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  TCanvas *cc=new TCanvas();
  for(unsigned int k=0; k<48; k++){
    TH2D* hRadPro=new TH2D("hName","",100,0,10,4000,0,100);
    func=new TF1("func","pol2",0,10);
    for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
      hRadPro->Fill(timeCalib->time.at(i),timeCalib->longitudinalVec.at(i).layer[k]);
    }
    hRadPro->ProfileX()->Fit(func,"NQ");
    //cc->WaitPrimitive();
    out << k 
	<< " " << func->GetParameter(0) << " " << func->GetParError(0) 
	<< " " << func->GetParameter(1) << " " << func->GetParError(1)
	<< " " << func->GetParameter(2) << " " << func->GetParError(2)
	<< " " << func->GetChisquare()/func->GetNDF() << std::endl;
    delete hRadPro;
  }
  out.close();
  delete func;
  delete timeCalib;
}

void ProcessElectron()
{
  int energy[]={10,20,30,40,50};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    //LongiProfile(energy[i],std::string("AugSep2012"));
    NhitCalibration(energy[i],std::string("AugSep2012"));
    //NclustersCalibration(energy[i],std::string("AugSep2012"));
    //RadialProfile(energy[i], std::string("AugSep2012"));
  }
}
