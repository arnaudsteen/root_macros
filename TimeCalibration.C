#define TimeCalibration_cxx
#include "TimeCalibration.h"
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

int RunElectron(double E,std::string period)
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

void TimeCalibration::FractalDimCut(int energy)
{
  if(energy==5)       cutFD=0.20;
  else if(energy==10) cutFD=0.18;
  else if(energy==15) cutFD=0.18;
  else if(energy==20) cutFD=0.17;
  else if(energy==25) cutFD=0.17;
  else if(energy==30) cutFD=0.16;
  else if(energy==40) cutFD=0.15;
  else if(energy==50) cutFD=0.13;
  else if(energy==60) cutFD=0.13;
  else if(energy==70) cutFD=0.13;
  else if(energy==80) cutFD=0.12;
}

void TimeCalibration::ShowerBarycenterCut(std::string runPeriod)
{
  if(runPeriod==std::string("Simulation")){
    cutBary[0]=40;
    cutBary[1]=55;
    cutBary[2]=40;
    cutBary[3]=55;
  }
  if(runPeriod==std::string("AugSep2012")){
    cutBary[0]=45;
    cutBary[1]=65;
    cutBary[2]=35;
    cutBary[3]=60;
  }
  if(runPeriod==std::string("Nov2012")){
    cutBary[0]=40;
    cutBary[1]=60;
    cutBary[2]=45;
    cutBary[3]=65;    
  }
#ifndef NO_DEBUG
  std::cout << cutBary[0] << " " << cutBary[1] << " " << cutBary[2] << " " << cutBary[3] << " " << std::endl;
#endif
}

void TimeCalibration::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(jentry>=maxentry) return;
    if( (float)Nhit/Nlayer>3 &&
	Begin>=0 &&
	(float)NInteractinglayer/Nlayer>0.2 &&
	Single==1 &&
	Neutral==0 &&
	Hole==0 &&
	(float)spillEventTime*200/pow(10.,9)<15 &&
	//F3D/log(Nhit)*CentralRatio<cutFD &&
	CoG[0]>cutBary[0] &&
	CoG[0]<cutBary[1] &&
	CoG[2]>cutBary[2] &&
	CoG[2]<cutBary[3] ){
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      time.push_back( (float)spillEventTime*200/pow(10.,9) );
      //radialProfile radPro;
      //for(unsigned int i=0; i<96; i++)
      //	radPro.ring[i]=RadialProfile[i];
      //radialVec.push_back(radPro);
      //longitudinalProfile longiPro;
      //for(unsigned int i=0; i<48; i++)
      //	longiPro.layer[i]=RadialProfile[i];
      //longitudinalVec.push_back(longiPro);
    }      
  }
}

void TimeCalibration::LoopElectron()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if( (float)Nhit/Nlayer>3 &&
	Begin>=0 &&
	Begin<=5 &&
	Nlayer<=30 &&
	(float)NInteractinglayer/Nlayer>0.2 &&
	Single==1 &&
	Neutral==0 &&
	Hole==0 &&
	(float)spillEventTime*200/pow(10.,9)<15 &&
	F3D/log(Nhit)*CentralRatio>cutFD &&
	CoG[0]>cutBary[0] &&
	CoG[0]<cutBary[1] &&
	CoG[2]>cutBary[2] &&
	CoG[2]<cutBary[3] ){
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
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
}

void CalibrationPion(int E, std::string runPeriod)
{
  char inputFile[200];
  //  sprintf(inputFile,"%s%d%s","/home/arnaud/Work/TB_Data/SPS_2012/RootFiles/DHCAL_",Run(E,runPeriod),".root");
  sprintf(inputFile,"%s%d%s","/home/steen/Back-Up/TB_Data/SPS_2012/RootFiles/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  TimeCalibration* timeCalib=new TimeCalibration(inputFile);
  timeCalib->FractalDimCut(E);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  TH2D* hN1=new TH2D("hN1","",100,0,10,750,0,1500);
  TH2D* hN2=new TH2D("hN2","",100,0,10,200,0,400);
  TH2D* hN3=new TH2D("hN3","",100,0,10,100,0,200);
  for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
    hN1->Fill(timeCalib->time.at(i),timeCalib->nhit1.at(i));
    hN2->Fill(timeCalib->time.at(i),timeCalib->nhit2.at(i));
    hN3->Fill(timeCalib->time.at(i),timeCalib->nhit3.at(i));
  }
  TProfile* profile1=new TProfile();
  TProfile* profile2=new TProfile();
  TProfile* profile3=new TProfile();
  profile1=hN1->ProfileX();
  profile2=hN2->ProfileX();
  profile3=hN3->ProfileX();

  TF1 *f1=new TF1("f1","pol2",0,10);
  TF1 *f2=new TF1("f2","pol2",0,10);
  TF1 *f3=new TF1("f3","pol2",0,10);
  profile1->Fit(f1,"NQ");
  profile2->Fit(f2,"NQ");
  profile3->Fit(f3,"NQ");
  
  char outputFile[200];
  //  sprintf(outputFile,"%s%d%s","/home/arnaud/Work/TB_Data/SPS_2012/Calibration/calib_",Run(E,runPeriod),"_error.txt");
  sprintf(outputFile,"%s%d%s","/home/steen/Back-Up/TB_Data/SPS_2012/Calibration/calib_",Run(E,runPeriod),"_error.txt");
  //  printf(outputFile);
  fstream out;
  out.open(outputFile,ios::out);
  out << f1->GetParameter(0) << " " << f1->GetParError(0) << " "
      << f1->GetParameter(1) << " " << f1->GetParError(1) << " "
      << f1->GetParameter(2) << " " << f1->GetParError(2) << std::endl;
  out << f2->GetParameter(0) << " " << f2->GetParError(0) << " "
      << f2->GetParameter(1) << " " << f2->GetParError(1) << " "
      << f2->GetParameter(2) << " " << f2->GetParError(2) << std::endl;
  out << f3->GetParameter(0) << " " << f3->GetParError(0) << " "
      << f3->GetParameter(1) << " " << f3->GetParError(1) << " "
      << f3->GetParameter(2) << " " << f3->GetParError(2) << std::endl;
  out.close();
  //out << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << std::endl;
  //out << f3->GetParameter(0) << " " << f3->GetParameter(1) << " " << f3->GetParameter(2) << std::endl;
}

void RadialProfileCalibrationPion(int E, std::string runPeriod)
{
  char inputFile[200];
  //sprintf(inputFile,"%s%d%s","/home/arnaud/Work/TB_Data/SPS_2012/RootFiles/DHCAL_",Run(E,runPeriod),".root");
  sprintf(inputFile,"%s%d%s","/home/steen/Back-Up/TB_Data/SPS_2012/RootFiles/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  TimeCalibration* timeCalib=new TimeCalibration(inputFile);
  timeCalib->FractalDimCut(E);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  std::cout<< timeCalib->radialVec.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  //  sprintf(outputFile,"%s%d%s","/home/arnaud/Work/TB_Data/SPS_2012/Calibration/radialCalib_",Run(E,runPeriod),".txt");
  sprintf(outputFile,"%s%d%s","/home/steen/Back-Up/TB_Data/SPS_2012/Calibration/radialCalib_",Run(E,runPeriod),"_error.txt");
  fstream out;
  out.open(outputFile,ios::out);
  //  TCanvas *cc=new TCanvas();
  for(unsigned int r=0; r<48; r++){
    TH2D* hRadPro=new TH2D("hName","",50,0,10,4000,0,200);
    func=new TF1("func","pol2",0,10);
    for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
      hRadPro->Fill(timeCalib->time.at(i),timeCalib->radialVec.at(i).ring[r]);
    }
    //   std::cout << r << std::endl;
    hRadPro->ProfileX()->Fit(func,"NQ");
    //cc->WaitPrimitive();
    out << r 
	<< " " << func->GetParameter(0) << " " << func->GetParError(0) 
	<< " " << func->GetParameter(1) << " " << func->GetParError(1)
	<< " " << func->GetParameter(2) << " " << func->GetParError(2)
	<< " " << func->GetChisquare()/func->GetNDF() << std::endl;
    //    out << r << " " << func->GetParameter(0) << " " << func->GetParameter(1) << " " << func->GetParameter(2) << " " << func->GetChisquare()/func->GetNDF() << std::endl;
    delete hRadPro;
  }
  delete func;
  delete timeCalib;
}

void ProcessPion()
{
  int energy[]={5,10,15,20,25,30,40,50,60,70,80};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    //CalibrationPion(energy[i], std::string("AugSep2012"));
    RadialProfileCalibrationPion(energy[i], std::string("AugSep2012"));
  }
  int energyNov[]={10,20,30,40,50,60,70,80};
  unsigned int eSizeNov=sizeof(energyNov)/sizeof(int);
  for(unsigned int i=0; i<eSizeNov; i++){
    //CalibrationPion(energyNov[i], std::string("Nov2012"));
    RadialProfileCalibrationPion(energyNov[i], std::string("Nov2012"));
  }
}

void CalibrationElectron(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/home/steen/v02-06_analysis/tb_data/RootFiles/DHCAL_",RunElectron(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  TimeCalibration* timeCalib=new TimeCalibration(inputFile);
  timeCalib->FractalDimCut(E);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->LoopElectron();
  
  TH2D* hN1=new TH2D("hN1","",100,0,10,300,0,600);
  TH2D* hN2=new TH2D("hN2","",100,0,10,90,0,180);
  TH2D* hN3=new TH2D("hN3","",100,0,10,60,0,120);
  for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
    hN1->Fill(timeCalib->time.at(i),timeCalib->nhit1.at(i));
    hN2->Fill(timeCalib->time.at(i),timeCalib->nhit2.at(i));
    hN3->Fill(timeCalib->time.at(i),timeCalib->nhit3.at(i));
  }
  TProfile* profile1=new TProfile();
  TProfile* profile2=new TProfile();
  TProfile* profile3=new TProfile();
  profile1=hN1->ProfileX();
  profile2=hN2->ProfileX();
  profile3=hN3->ProfileX();

  TF1 *f1=new TF1("f1","pol3",0,10);
  TF1 *f2=new TF1("f2","pol3",0,10);
  TF1 *f3=new TF1("f3","pol3",0,10);
  profile1->Fit(f1,"NQ");
  profile2->Fit(f2,"NQ");
  profile3->Fit(f3,"NQ");
  
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/v02-06_analysis/tb_data/Calibration/calib_",RunElectron(E,runPeriod),"_error.txt");
  fstream out;
  out.open(outputFile,ios::out);
  out << f1->GetParameter(0) << " " << f1->GetParError(0) << " "
      << f1->GetParameter(1) << " " << f1->GetParError(1) << " "
      << f1->GetParameter(2) << " " << f1->GetParError(2) << " "
      << f1->GetParameter(3) << " " << f1->GetParError(3) << std::endl;
  out << f2->GetParameter(0) << " " << f2->GetParError(0) << " "
      << f2->GetParameter(1) << " " << f2->GetParError(1) << " "
      << f2->GetParameter(2) << " " << f2->GetParError(2) << " "
      << f2->GetParameter(3) << " " << f2->GetParError(3) << std::endl;
  out << f3->GetParameter(0) << " " << f3->GetParError(0) << " "
      << f3->GetParameter(1) << " " << f3->GetParError(1) << " "
      << f3->GetParameter(2) << " " << f3->GetParError(2) << " "
      << f3->GetParameter(3) << " " << f3->GetParError(3) << std::endl;
  out.close();
  //out << f1->GetParameter(0) << " " << f1->GetParameter(1) << " " << f1->GetParameter(2) << " " << f1->GetParameter(3) << std::endl;
  //out << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << " " << f2->GetParameter(3) << std::endl;
  //out << f3->GetParameter(0) << " " << f3->GetParameter(1) << " " << f3->GetParameter(2) << " " << f3->GetParameter(3) << std::endl;
}

void RadialProfileCalibrationElectron(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","./RootFiles/DHCAL_",RunElectron(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  TimeCalibration* timeCalib=new TimeCalibration(inputFile);
  timeCalib->FractalDimCut(E);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->LoopElectron();

  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","./TXTFile/radialCalib_",RunElectron(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  //TCanvas *cc=new TCanvas;
  for(unsigned int r=0; r<48; r++){
    func=new TF1("func","pol2",0,10);
    TH2D* hRadPro=new TH2D("hName","",100,0,10,2000,0,200);
    for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
      hRadPro->Fill(timeCalib->time.at(i),timeCalib->radialVec.at(i).ring[r]);
    }
    //TCanvas *cc=new TCanvas();
    hRadPro->ProfileX()->Fit(func,"");
    //cc->WaitPrimitive();
    out << r 
	<< " " << func->GetParameter(0) << " " << func->GetParError(0) 
	<< " " << func->GetParameter(1) << " " << func->GetParError(1)
	<< " " << func->GetParameter(2) << " " << func->GetParError(2)
	<< " " << func->GetParameter(3) << " " << func->GetParError(3) 
	<< " " << func->GetChisquare()/func->GetNDF() << std::endl;
    delete hRadPro;
  }
  delete func;
  delete timeCalib;
}

void LongitudinalProfileCalibrationElectron(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/home/steen/Back-Up/TB_Data/SPS_2012/RootFiles/DHCAL_",RunElectron(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  TimeCalibration* timeCalib=new TimeCalibration(inputFile);
  timeCalib->FractalDimCut(E);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->LoopElectron();

  TCanvas *cc=new TCanvas();
  cc->cd();
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/Back-Up/TB_Data/SPS_2012/Calibration/longitudinalCalib_",RunElectron(E,runPeriod),".txt");
  std::cout << "outputFile = " << outputFile << std::endl;
  fstream out;
  out.open(outputFile,ios::out);
  for(unsigned int r=0; r<48; r++){
    TH2D* hLongPro=new TH2D("hName","",100,0,10,10000,-100,100);
    func=new TF1("func","pol3",0,10);
    for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
      hLongPro->Fill(timeCalib->time.at(i),timeCalib->longitudinalVec.at(i).layer[r]);
    }
    hLongPro->Draw("colz");
    //cc->WaitPrimitive();
    hLongPro->ProfileX()->Fit(func,"NQ");
    //cc->WaitPrimitive();
    out << r << " " << func->GetParameter(0) << " " << func->GetParameter(1) << " " << func->GetParameter(2) << " " << func->GetParameter(3) << std::endl;
    delete hLongPro;
  }
  delete func;
}

void ProcessElectron()
{
 int energy[]={10,20,30,40,50};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    //CalibrationElectron(energy[i], std::string("AugSep2012"));
    RadialProfileCalibrationElectron(energy[i], std::string("AugSep2012"));
    //LongitudinalProfileCalibrationElectron(energy[i], std::string("AugSep2012"));
  }
}

