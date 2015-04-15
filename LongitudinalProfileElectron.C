#define LongitudinalProfileElectron_cxx
#include "LongitudinalProfileElectron.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstring>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <TLegend.h>

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

void LongitudinalProfileElectron::ShowerBarycenterCut(std::string runPeriod)
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

void LongitudinalProfileElectron::FractalDimCut(int energy)
{
  if(energy==10) cutFD=0.18;
  else if(energy==20) cutFD=0.17;
  else if(energy==30) cutFD=0.16;
  else if(energy==40) cutFD=0.15;
  else if(energy==50) cutFD=0.13;
  else cutFD=0;
}

void LongitudinalProfileElectron::Loop()
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
	Begin<5 &&
	Nlayer<30 &&
	TrackMultiplicity==0 &&
	Single==1 &&
	Neutral==0 &&
	Hole==0 &&
	(float)spillEventTime*200/pow(10.,9)<15 &&
	CoG[0]>cutBary[0] &&
	CoG[0]<cutBary[1] &&
	CoG[2]>cutBary[2] &&
	CoG[2]<cutBary[3] &&
	ReconstructedCosTheta>.9){
      time.push_back( (float)spillEventTime*200/pow(10.,9) );
      longiProfile longiPro;
      for(unsigned int i=0; i<48; i++)
      	longiPro.layer[i]=LongiProfile[i];
      longiVec.push_back(longiPro);
    }      
  }
}

std::vector<LongitudinalCalibration> readCalib(int run)
{
  std::vector<LongitudinalCalibration> calibVec;
  char inputCalib[200];
  sprintf(inputCalib,"%s%d%s","/home/steen/timeCalib/LongiProfile/calib_",run,".txt");
  fstream in;
  in.open(inputCalib);
  if(!in.is_open()){ std::cout << inputCalib << "\t NO SUCH FILE OR DIRECTORY" << std::endl; throw;}
  double coeff[3];
  double coeffError[3];
  double chi2;
  int r;
  while(1){
    if(!in.good()) break;
    in >> r >> coeff[0] >> coeffError[0] >> coeff[1] >> coeffError[1] >> coeff[2] >> coeffError[2] >> chi2;
    LongitudinalCalibration longiCalib;
    longiCalib.layer=r;
    longiCalib.coeff0=coeff[0];
    longiCalib.coeff1=coeff[1];
    longiCalib.coeff2=coeff[2];
    longiCalib.coefferror0=coeffError[0];
    longiCalib.coefferror1=coeffError[1];
    longiCalib.coefferror2=coeffError[2];
    longiCalib.chi2=chi2;
    calibVec.push_back(longiCalib);
    /*std::cout << longiCalib.radius << " " 
	      << longiCalib.coeff0 << " " 
	      << longiCalib.coeff1 << " " 
	      << longiCalib.coeff2 << std::endl;*/
  }
  return calibVec;
}

void SimulationProfileStyleDefinition(TProfile* prof,std::string model)
{
  if( model==std::string("ftfp_bert_hp") ){
    prof->SetLineWidth(1);
    prof->SetFillColor(kRed-3);
    prof->SetLineColor(kRed-3);
    prof->SetFillStyle(1001);
  }
  else if( model==std::string("qgsp_bert_hp") ){
    prof->SetLineWidth(1);
    prof->SetFillColor(kBlue-6);
    prof->SetLineColor(kBlue-6);
    prof->SetFillStyle(1001);
  }
  else 
    std::cout << "WRONG CHOSEN MODEL OPTION" << std::endl;
}

void LongitudinalProfile(int E, std::string runPeriod,std::string model)
{
  char inputFile[200];
  char text[100];
  char plotName[100];
  char fDataName[200];
  char fSimuName[200];
  float simMean=0;
  float simRMS=0;
  float simCount=0;
  float simMeanError=0;
  float simRMSError=0;
  float dataMean=0;
  float dataRMS=0;
  float dataMeanError=0;
  float dataRMSError=0;
  float dataCount=0;

  CaliceStyle();
  int run=Run(E,runPeriod);
  std::vector<LongitudinalCalibration> calibVec=readCalib(run);
  sprintf(inputFile,"%s%d%s","/home/steen/resultRootFile/tb_data/DHCAL_",run,".root");
  std::cout << inputFile << std::endl;
  LongitudinalProfileElectron* data=new LongitudinalProfileElectron(inputFile);
  data->ShowerBarycenterCut(runPeriod);
  data->FractalDimCut(E);
  data->Loop();  
  TH2D* hData = new TH2D("hData","",40,0,40,1000,0,1);
  std::vector<float> dataNorm;
  std::cout << data->time.size() << " data " << std::endl;
  float correction=0;
  for(unsigned int evt=0; evt<data->time.size(); evt++){
    dataNorm.push_back(0);
    for(unsigned int r=0; r<calibVec.size(); r++){
      correction=calibVec.at(r).coeff1*data->time.at(evt)+
  	calibVec.at(r).coeff2*data->time.at(evt)*data->time.at(evt) ;
      if(correction>0)correction=0;
      dataNorm.back()+=(data->longiVec.at(evt).layer[r]-correction);
    }
    for(unsigned int r=0; r<calibVec.size(); r++){
      correction = calibVec.at(r).coeff1*data->time.at(evt)+
	calibVec.at(r).coeff2*data->time.at(evt)*data->time.at(evt) ;
      if(correction>0)correction=0;
      hData->Fill(r,(data->longiVec.at(evt).layer[r]-correction)/dataNorm.at(evt));
      dataMean+=r*(data->longiVec.at(evt).layer[r]-correction)/dataNorm.at(evt);
      dataRMS+=r*r*(data->longiVec.at(evt).layer[r]-correction)/dataNorm.at(evt);
      dataCount+=(data->longiVec.at(evt).layer[r]-correction)/dataNorm.at(evt);
    }
  }
  TProfile* profileData=new TProfile();
  profileData=hData->ProfileX();
  profileData->SetMarkerStyle(20);
  profileData->SetMarkerSize(1);
  profileData->SetLineWidth(2);
  dataMean=dataMean/dataCount;
  dataRMS=sqrt(dataRMS/dataCount-dataMean*dataMean);  
  dataMeanError=dataRMS/sqrt(dataCount);
  dataRMSError=sqrt(2/dataCount)*dataRMS;
  std::cout << "dataMean = " << dataMean << "\t dataRMS = " << dataRMS << std::endl; 

  sprintf(inputFile,"%s%s%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/",model.c_str(),"/single_e-_",E,"GeV.root");
  std::cout << inputFile << std::endl;
  LongitudinalProfileElectron* simu=new LongitudinalProfileElectron(inputFile);
  simu->FractalDimCut(E);
  simu->Loop();  
  TH2D* hSimu = new TH2D("hSimu","",40,0,40,1000,0,1);
  std::vector<float> simuNorm;
  for(unsigned int evt=0; evt<simu->longiVec.size(); evt++){
    simuNorm.push_back(0);
    for(unsigned int r=0; r<calibVec.size(); r++){
      simuNorm.back()+=simu->longiVec.at(evt).layer[r];
    }
    for(unsigned int r=0; r<49; r++){
      hSimu->Fill(r,simu->longiVec.at(evt).layer[r]/simuNorm.at(evt));
      simMean+=r*simu->longiVec.at(evt).layer[r]/simuNorm.at(evt);
      simRMS+=r*r*simu->longiVec.at(evt).layer[r]/simuNorm.at(evt);
      simCount+=simu->longiVec.at(evt).layer[r]/simuNorm.at(evt);
    }
  }
  TProfile* profileSimu=new TProfile();
  profileSimu=hSimu->ProfileX();  
  SimulationProfileStyleDefinition(profileSimu,model);
  profileSimu->GetXaxis()->SetTitle("distance from shower axis [cm]");
  profileSimu->GetYaxis()->SetTitle("N_{hit}");
  simMean=simMean/simCount;
  simRMS=sqrt(simRMS/simCount-simMean*simMean);
  simMeanError=simRMS/sqrt(simCount);
  simRMSError=sqrt(2/simCount)*simRMS;
  std::cout << "simMean = " << simMean << "\t simRMS = " << simRMS << std::endl; 

  TCanvas *can=new TCanvas();
  can->SetWindowSize(600,600);
  if( profileSimu->GetBinContent(profileSimu->GetMaximumBin()) > profileData->GetBinContent(profileData->GetMaximumBin()) )
    profileSimu->Draw("axis");
  else 
    profileData->Draw("axis");
  profileSimu->Draw("histsame");
  profileData->Draw("same");
  if( profileSimu->GetBinContent(profileSimu->GetMaximumBin()) > profileData->GetBinContent(profileData->GetMaximumBin()) )
    profileSimu->Draw("axissame");
  else 
    profileData->Draw("axissame");
  TText *tex=new TText();
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.46,0.96,"CALICE Fe-SDHCAL Preliminary");
  sprintf(text,"%s%d%s","E = ",E," GeV");
  tex->DrawTextNDC(0.11,0.96,text);
  can->Update();
  
  TLegend *leg=new TLegend(0.52,0.73,0.91,0.9);
  leg->SetFillStyle(0);
  std::transform(model.begin(),model.end(),model.begin(),::toupper);
  sprintf(text,"%s%s%s","MC (",model.c_str(),")");
  leg->AddEntry(profileSimu,text,"f");
  leg->AddEntry(profileData,"Data","p");
  leg->Draw();
  can->Update();

  std::transform(model.begin(),model.end(),model.begin(),::tolower);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProf_e-_",E,"GeV_",runPeriod.c_str(),".C");
  can->SaveAs(plotName);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProf_e-_",E,"GeV_",runPeriod.c_str(),".pdf");
  can->SaveAs(plotName);
  
  //can->SetLogy();
  //can->Update();
  //sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProfLog_e-_",E,"GeV_",runPeriod.c_str(),".C");
  //can->SaveAs(plotName);
  //sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProfLog_e-_",E,"GeV_",runPeriod.c_str(),".pdf");
  //can->SaveAs(plotName);
  
  fstream out;
  sprintf(fDataName,"%s%s%s","./txtFile/longiProf_e-_",runPeriod.c_str(),".txt");
  if( model==std::string("ftfp_bert_hp") && E==10&&runPeriod==std::string("AugSep2012") )
    out.open(fDataName,ios::out);
  else if(model==std::string("ftfp_bert_hp"))
    out.open(fDataName,ios::out|ios::app);
  if(out.is_open()){
    out << E << " " << dataMean << " " << dataMeanError << " " << dataRMS << " " << dataRMSError << std::endl;
    out.close();
  }
  
  sprintf(fSimuName,"%s%s%s","./txtFile/longiProf_e-_",model.c_str(),".txt");
  if( E==10 && runPeriod=="AugSep2012")
    out.open(fSimuName,ios::out);
  else if(runPeriod=="AugSep2012")
    out.open(fSimuName,ios::out|ios::app);
  if(out.is_open()){
    out << E << " " << simMean << " " << simMeanError << " " << simRMS << " " << simRMSError << std::endl;
    out.close();
  }
}

void ProcessElectron()
{
  int energy[]={10,20,30,40,50};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    LongitudinalProfile(energy[i],std::string("AugSep2012"),std::string("ftfp_bert_hp"));
    LongitudinalProfile(energy[i],std::string("AugSep2012"),std::string("qgsp_bert_hp"));
  }
}
