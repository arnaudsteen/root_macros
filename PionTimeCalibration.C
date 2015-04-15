#define PionTimeCalibration_cxx
#include "PionTimeCalibration.h"
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
    else if(E==70) return 715493;
    else if(E==60) return 715531;
    else if(E==50) return 715751;
    else if(E==40) return 715651;
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

void PionTimeCalibration::ShowerBarycenterCut(std::string runPeriod)
{
  if(runPeriod==std::string("AugSep2012")){
    cutBary[0]=400;
    cutBary[1]=800;
    cutBary[2]=250;
    cutBary[3]=750;
  }
  if(runPeriod==std::string("Nov2012")){
    cutBary[0]=400;
    cutBary[1]=750;
    cutBary[2]=350;
    cutBary[3]=700;
  }
  if(runPeriod==std::string("Dec2014")){
    cutBary[0]=490;
    cutBary[1]=690;
    cutBary[2]=400;
    cutBary[3]=600;
  }
}

void PionTimeCalibration::Loop()
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
	((float)NInteractinglayer/Nlayer>0.1||Nlayer<35) &&
	Single==1 &&
	Neutral==0 &&
	Hole==0 &&
	(float)spillEventTime*200/pow(10.,9)<15 &&
	CoG[0]>cutBary[0] &&
	CoG[0]<cutBary[1] &&
	CoG[2]>cutBary[2] &&
	CoG[2]<cutBary[3] &&
	ReconstructedCosTheta>.9&&
	(Nlayer>=30||Begin>=5||TrackMultiplicity>0)){
      //htemp->Fill((float)spillEventTime*200/pow(10.,9),Nhit1+Nhit2+Nhit3);
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      nhit2by2.push_back(Nhit2by2);
      nhit3by3.push_back(Nhit3by3);
      nhit4by4.push_back(Nhit4by4);
      nhit5by5.push_back(Nhit5by5);
      ncluster.push_back(Nclusters);
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
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
  //timeCalib->FractalDimCut(E);
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
  profile1->Fit(f1,"");
  profile2->Fit(f2,"");
  profile3->Fit(f3,"");  
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/Nhit/calib_",Run(E,runPeriod),".txt");
  std::cout << "outputFile = " << outputFile << std::endl;
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
}

void RadialProfile(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
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
    hRadPro->ProfileX()->Fit(func,"NQ");
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
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();  
  TH2D* hN1=new TH2D("hN1","",100,0,10,300,0,300);
  for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
    hN1->Fill(timeCalib->time.at(i),timeCalib->ncluster.at(i));
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
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  std::cout<< timeCalib->longitudinalVec.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/LongiProfile/calib_",Run(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  for(unsigned int k=0; k<48; k++){
    TH2D* hRadPro=new TH2D("hName","",100,0,10,4000,0,200);
    func=new TF1("func","pol2",0,10);
    for(unsigned int i=0; i<timeCalib->nhit1.size(); i++){
      hRadPro->Fill(timeCalib->time.at(i),timeCalib->longitudinalVec.at(i).layer[k]);
    }
    hRadPro->ProfileX()->Fit(func,"NQ");
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

void Nhit2by2Calibration(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  std::cout<< timeCalib->nhit2by2.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/Nhit2by2/calib_",Run(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  TH2D* hRadPro=new TH2D("hName","",100,0,10,1000,0,1000);
  func=new TF1("func","pol1",0,10);
  for(unsigned int i=0; i<timeCalib->nhit2by2.size(); i++){
    hRadPro->Fill(timeCalib->time.at(i),timeCalib->nhit2by2.at(i));
  }
  //TCanvas *cc=new TCanvas();
  //hRadPro->Draw("colz");
  //cc->Update();
  //cc->WaitPrimitive();
  hRadPro->ProfileX()->Fit(func,"NQ");
  //cc->Update();
  //cc->WaitPrimitive();
  out << func->GetParameter(0) << " " << func->GetParError(0) 
      << " " << func->GetParameter(1) << " " << func->GetParError(1) 
      << " " << func->GetChisquare()/func->GetNDF() << std::endl;
  delete hRadPro;
  out.close();
  delete func;
  delete timeCalib;
  //delete cc;
}

void Nhit3by3Calibration(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  std::cout<< timeCalib->nhit3by3.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/Nhit3by3/calib_",Run(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  TH2D* hRadPro=new TH2D("hName","",100,0,10,800,0,800);
  func=new TF1("func","pol1",0,10);
  for(unsigned int i=0; i<timeCalib->nhit3by3.size(); i++){
    hRadPro->Fill(timeCalib->time.at(i),timeCalib->nhit3by3.at(i));
  }
  //TCanvas *cc=new TCanvas();
  //hRadPro->Draw("colz");
  //cc->Update();
  //cc->WaitPrimitive();
  hRadPro->ProfileX()->Fit(func,"NQ");
  //cc->Update();
  //cc->WaitPrimitive();
  out << func->GetParameter(0) << " " << func->GetParError(0) 
      << " " << func->GetParameter(1) << " " << func->GetParError(1) 
      << " " << func->GetChisquare()/func->GetNDF() << std::endl;
  delete hRadPro;
  out.close();
  delete func;
  delete timeCalib;
  //delete cc;
}

void Nhit4by4Calibration(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  std::cout<< timeCalib->nhit4by4.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/Nhit4by4/calib_",Run(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  TH2D* hRadPro=new TH2D("hName","",100,0,10,700,0,700);
  func=new TF1("func","pol1",0,10);
  for(unsigned int i=0; i<timeCalib->nhit4by4.size(); i++){
    hRadPro->Fill(timeCalib->time.at(i),timeCalib->nhit4by4.at(i));
  }
  //TCanvas *cc=new TCanvas();
  //hRadPro->Draw("colz");
  //cc->Update();
  //cc->WaitPrimitive();
  hRadPro->ProfileX()->Fit(func,"NQ");
  //cc->Update();
  //cc->WaitPrimitive();
  out << func->GetParameter(0) << " " << func->GetParError(0) 
      << " " << func->GetParameter(1) << " " << func->GetParError(1) 
      << " " << func->GetChisquare()/func->GetNDF() << std::endl;
  delete hRadPro;
  out.close();
  delete func;
  delete timeCalib;
  //delete cc;
}

void Nhit5by5Calibration(int E, std::string runPeriod)
{
  char inputFile[200];
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  PionTimeCalibration* timeCalib=new PionTimeCalibration(inputFile);
  timeCalib->ShowerBarycenterCut(runPeriod);
  timeCalib->Loop();
  
  std::cout<< timeCalib->nhit5by5.size()<<std::endl;
  TF1 *func;
  char outputFile[200];
  sprintf(outputFile,"%s%d%s","/home/steen/timeCalib/Nhit5by5/calib_",Run(E,runPeriod),".txt");
  fstream out;
  out.open(outputFile,ios::out);
  TH2D* hRadPro=new TH2D("hName","",100,0,10,600,0,600);
  func=new TF1("func","pol1",0,10);
  for(unsigned int i=0; i<timeCalib->nhit5by5.size(); i++){
    hRadPro->Fill(timeCalib->time.at(i),timeCalib->nhit5by5.at(i));
  }
  //TCanvas *cc=new TCanvas();
  //hRadPro->Draw("colz");
  //cc->Update();
  //cc->WaitPrimitive();
  hRadPro->ProfileX()->Fit(func,"NQ");
  //cc->Update();
  //cc->WaitPrimitive();
  out << func->GetParameter(0) << " " << func->GetParError(0) 
      << " " << func->GetParameter(1) << " " << func->GetParError(1) 
      << " " << func->GetChisquare()/func->GetNDF() << std::endl;
  delete hRadPro;
  out.close();
  delete func;
  delete timeCalib;
  //delete cc;
}

void ProcessPion()
{
  int energy[]={5,10,15,20,25,30,40,50,60,70,80};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    Nhit2by2Calibration(energy[i],std::string("AugSep2012"));
    Nhit3by3Calibration(energy[i],std::string("AugSep2012"));
    Nhit4by4Calibration(energy[i],std::string("AugSep2012"));
    Nhit5by5Calibration(energy[i],std::string("AugSep2012"));
    //LongiProfile(energy[i],std::string("AugSep2012"));
    //NclustersCalibration(energy[i],std::string("AugSep2012"));
    //NhitCalibration(energy[i],std::string("AugSep2012"));
    //RadialProfile(energy[i], std::string("AugSep2012"));
  }
  int energyNov[]={10,20,30,40,50,60,70,80};
  unsigned int eSizeNov=sizeof(energyNov)/sizeof(int);
  for(unsigned int i=0; i<eSizeNov; i++){
    //LongiProfile(energyNov[i],std::string("Nov2012"));
    //NclustersCalibration(energyNov[i],std::string("Nov2012"));
    //NhitCalibration(energyNov[i], std::string("Nov2012"));
    //RadialProfile(energyNov[i], std::string("Nov2012"));
  }
}
