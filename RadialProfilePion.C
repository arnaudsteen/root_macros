#define RadialProfilePion_cxx
#include "RadialProfilePion.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstring>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <numeric>
#include <TLegend.h>

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

void RadialProfilePion::ShowerBarycenterCut(std::string runPeriod)
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

void RadialProfilePion::Loop()
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
	//(float)NInteractinglayer/Nlayer>0.2 &&
	TransverseRatio>0.01&&
	Single==1 &&
	Neutral==0 &&
	Hole==0 &&
	(float)spillEventTime*200/pow(10.,9)<15 &&
	CoG[0]>cutBary[0] &&
	CoG[0]<cutBary[1] &&
	CoG[2]>cutBary[2] &&
	CoG[2]<cutBary[3] &&
	ReconstructedCosTheta>.9&&
	TrackMultiplicity>0){
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
//       //       TrackClusterNumber->at(0)>5&&
//       TransverseRatio>0.01&&
//       (float)(Nhough1+Nhough2+Nhough3)/Nhit<0.5) {
      time.push_back( (float)spillEventTime*200/pow(10.,9) );
      radialProfile radPro;
      for(unsigned int i=0; i<96; i++)
      	radPro.ring[i]=RadialProfileBis[i];
      radialVec.push_back(radPro);
    }      
  }
}

std::vector<RadialCalibration> readCalib(int run)
{
  std::vector<RadialCalibration> calibVec;
  char inputCalib[200];
  sprintf(inputCalib,"%s%d%s","/home/steen/timeCalib/RadProfile/calib_",run,".txt");
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
    RadialCalibration radialCalib;
    radialCalib.radius=r;
    radialCalib.coeff0=coeff[0];
    radialCalib.coeff1=coeff[1];
    radialCalib.coeff2=coeff[2];
    radialCalib.coefferror0=coeffError[0];
    radialCalib.coefferror1=coeffError[1];
    radialCalib.coefferror2=coeffError[2];
    radialCalib.chi2=chi2;
    calibVec.push_back(radialCalib);
    /*std::cout << radialCalib.radius << " " 
      << radialCalib.coeff0 << " " 
      << radialCalib.coeff1 << " " 
      << radialCalib.coeff2 << std::endl;*/
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
  else if( model==std::string("ftfp_bert") ){
    prof->SetFillColor(kRed+1);
    prof->SetFillStyle(1001);
    prof->SetLineColor(kRed+1);
    prof->SetLineWidth(1);
  }
  else if( model==std::string("qgsp_bert") ){
    prof->SetFillColor(kBlue+2);
    prof->SetFillStyle(1001);
    prof->SetLineColor(kBlue+2);
    prof->SetLineWidth(1);
  }
  else 
    std::cout << "WRONG CHOSEN MODEL OPTION" << std::endl;
}

void RadialProfile(int E, std::string runPeriod,std::string model)
{
  char inputFile[200];
  char text[100];
  char plotName[100];
  char fDataName[200];
  char fSimuName[200];
  float simMean=0;
  float simRMS=0;
  float simMeanError=0;
  float simRMSError=0;
  float simCount=0;
  float dataMean=0;
  float dataRMS=0;
  float dataMeanError=0;
  float dataRMSError=0;
  float dataCount=0;
  float correction;

  CaliceStyle();
  int run=Run(E,runPeriod);
  std::vector<RadialCalibration> calibVec=readCalib(run);
  //sprintf(inputFile,"%s%d%s","/home/steen/resultRootFile/tb_data/DHCAL_",run,".root");
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",run,".root");
  std::cout << inputFile << std::endl;
  RadialProfilePion* data=new RadialProfilePion(inputFile);
  data->ShowerBarycenterCut(runPeriod);
  data->Loop();  
  TH2D* hData = new TH2D("hData","",50,0,50,800,0,1);
  std::vector<float> dataNorm;
  for(unsigned int evt=0; evt<data->time.size(); evt++){
    dataNorm.push_back(0);
    for(unsigned int r=0; r<calibVec.size(); r++){
      correction=calibVec.at(r).coeff1*data->time.at(evt)+
  	calibVec.at(r).coeff2*data->time.at(evt)*data->time.at(evt) ;
      if(correction>0)correction=0;
      dataNorm.back()+=(data->radialVec.at(evt).ring[r]-correction);
    }
  }
  for(unsigned int evt=0; evt<data->time.size(); evt++){
    for(unsigned int r=0; r<calibVec.size(); r++){
      correction = calibVec.at(r).coeff1*data->time.at(evt)+
  	calibVec.at(r).coeff2*data->time.at(evt)*data->time.at(evt) ;
      if(correction>0)
	correction=0;
      hData->Fill(r,(data->radialVec.at(evt).ring[r]-correction)/dataNorm.at(evt));
      //if(r>20)continue;
      dataMean+=r*(data->radialVec.at(evt).ring[r]-correction)/dataNorm.at(evt);
      dataRMS+=r*r*(data->radialVec.at(evt).ring[r]-correction)/dataNorm.at(evt);
      dataCount+=(data->radialVec.at(evt).ring[r]-correction)/dataNorm.at(evt);
    }
  }
  TProfile* profileData=new TProfile();
  profileData=hData->ProfileX();
  profileData->SetMarkerStyle(20);
  profileData->SetMarkerSize(1);
  profileData->SetLineWidth(2);
  profileData->GetXaxis()->SetTitle("distance from shower axis [cm]");
  profileData->GetYaxis()->SetTitle("# hits (normalised to unity)");
  profileData->GetYaxis()->SetTitleOffset(1.2);
  dataMean=dataMean/dataCount;
  dataRMS=sqrt(dataRMS/dataCount-dataMean*dataMean);
  dataMeanError=dataRMS/sqrt(dataCount);
  dataRMSError=sqrt(2/dataCount)*dataRMS;
  
  //sprintf(inputFile,"%s%s%s%d%s","/home/steen/resultRootFile/sim_data/",model.c_str(),"/single_pi-_",E,"GeV.root");
  //sprintf(inputFile,"%s%d%s","/home/steen/sdhcal_analysis/single_pi-_",E,"GeV.root");
  sprintf(inputFile,"%s%s%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/",model.c_str(),"/timecut/single_pi-_",E,"GeV.root");
  std::cout << inputFile << std::endl;
  RadialProfilePion* simu=new RadialProfilePion(inputFile);
  simu->Loop();  
  TH2D* hSimu = new TH2D("hSimu","",50,0,50,800,0,1);
  std::vector<float> simuNorm;
  for(unsigned int evt=0; evt<simu->time.size(); evt++){
    simuNorm.push_back(0);
    for(unsigned int r=0; r<calibVec.size(); r++){
      simuNorm.back()+=simu->radialVec.at(evt).ring[r];
    }
  }
  for(unsigned int evt=0; evt<simu->radialVec.size(); evt++){
    for(unsigned int r=0; r<49; r++){
      hSimu->Fill(r,simu->radialVec.at(evt).ring[r]/simuNorm.at(evt));
      simMean+=r*simu->radialVec.at(evt).ring[r]/simuNorm.at(evt);
      simRMS+=r*r*simu->radialVec.at(evt).ring[r]/simuNorm.at(evt);
      simCount+=simu->radialVec.at(evt).ring[r]/simuNorm.at(evt);
    }
  }
  TProfile* profileSimu=new TProfile();
  profileSimu=hSimu->ProfileX();  
  SimulationProfileStyleDefinition(profileSimu,model);
  profileSimu->GetXaxis()->SetTitle("distance from shower axis [cm]");
  profileSimu->GetYaxis()->SetTitle("# hits (normalised to unity)");
  profileSimu->GetYaxis()->SetTitleOffset(1.2);
  simMean=simMean/simCount;
  simRMS=sqrt(simRMS/simCount-simMean*simMean);
  simMeanError=simRMS/sqrt(simCount);
  simRMSError=sqrt(2/simCount)*simRMS;

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
  tex->DrawTextNDC(0.45,0.96,"CALICE-Fe SDHCAL Preliminary");
  sprintf(text,"%s%d%s","E = ",E," GeV");
  tex->DrawTextNDC(0.11,0.96,text);
  //can->Update();
  
  TLegend *leg=new TLegend(0.52,0.73,0.91,0.9);
  leg->SetFillStyle(0);
  std::transform(model.begin(),model.end(),model.begin(),::toupper);
  sprintf(text,"%s%s%s","MC (",model.c_str(),")");
  leg->AddEntry(profileSimu,text,"f");
  leg->AddEntry(profileData,"Data","p");
  leg->Draw();
  //can->Update();
  
  std::transform(model.begin(),model.end(),model.begin(),::tolower);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProf_pi-_",E,"GeV_",runPeriod.c_str(),".C");
  can->SaveAs(plotName);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProf_pi-_",E,"GeV_",runPeriod.c_str(),".pdf");
  can->SaveAs(plotName);
  
  can->SetLogy();
  can->Update();
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProfLog_pi-_",E,"GeV_",runPeriod.c_str(),".C");
  can->SaveAs(plotName);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/radProfLog_pi-_",E,"GeV_",runPeriod.c_str(),".pdf");
  can->SaveAs(plotName);
  
  //TH1D *hsim=new TH1D("hsim","",50,0,50);
  //TH1D *hdata=new TH1D("hdata","",50,0,50);
  //for(unsigned int r=0; r<=49; r++){
  //  hsim->Fill(r,profileSimu->GetBinContent(profileSimu->GetBin(r+1)));
  //  hdata->Fill(r,profileData->GetBinContent(profileData->GetBin(r+1)));
  //}
  //dataMean=hdata->GetMean();
  //dataMeanError=hdata->GetMeanError()/sqrt(data->radialVec.size());
  //dataRMS=hdata->GetRMS();
  //dataRMSError=hdata->GetRMSError()/sqrt(data->radialVec.size());
  //simMean=hsim->GetMean();
  //simMeanError=hsim->GetMeanError()/sqrt(simu->radialVec.size());
  //simRMS=hsim->GetRMS();
  //simRMSError=hsim->GetRMSError()/sqrt(simu->radialVec.size());
  
  fstream out;
  sprintf(fDataName,"%s%s%s","./txtFile/radProf_pi-_",runPeriod.c_str(),".txt");
  if( model==std::string("ftfp_bert_hp") && ((E==5&&runPeriod==std::string("AugSep2012")) || (E==10&&runPeriod==std::string("Nov2012"))) )
    out.open(fDataName,ios::out);
  else if(model==std::string("ftfp_bert_hp"))
    out.open(fDataName,ios::out|ios::app);
  if(out.is_open()){
    out << E << " " << dataMean << " " << dataMeanError << " " << dataRMS << " " << dataRMSError << std::endl;
    out.close();
  }
  
  sprintf(fSimuName,"%s%s%s","./txtFile/radProf_pi-_",model.c_str(),".txt");
  if( E==5 && runPeriod=="AugSep2012")
    out.open(fSimuName,ios::out);
  else if(runPeriod=="AugSep2012")
    out.open(fSimuName,ios::out|ios::app);
  if(out.is_open()){
    out << E << " " << simMean << " " << simMeanError << " " << simRMS << " " << simRMSError << std::endl;
    out.close();
  }

   //std::cout << "simu" << std::endl;
   //std::cout << simMean << " #pm " << simMeanError << std::endl;
   //std::cout << simRMS << " #pm " << simRMSError << std::endl;
  
  
//   std::cout << "data" << std::endl;
//   std::cout << dataMean << " #pm " << dataMeanError << std::endl;
//   std::cout << dataRMS << " #pm " << dataRMSError << std::endl;
  
  
// std::cout << "mean relative deviation = " << (simMean-dataMean)/dataMean << "\t"
// 	  << "mean error on deviation = " << 1/dataMean*sqrt(simMeanError*simMeanError+dataMeanError*dataMeanError*simMean*simMean/dataMean/dataMean)
// 	  << std::endl;
 
//  std::cout << "rms relative deviation = " << (simRMS-dataRMS)/dataRMS << "\t"
// 	   << "rms error on deviation = " << 1/dataRMS*sqrt(simRMSError*simRMSError+dataRMSError*dataRMSError*simRMS*simRMS/dataRMS/dataRMS)
// 	   << std::endl;

}

void ProcessPion()
{
  int energy[]={5,10,15,20,25,30,40,50,60,70,80};
  //int energy[]={10,20,60,70};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    RadialProfile(energy[i],std::string("AugSep2012"),std::string("ftfp_bert_hp"));
    RadialProfile(energy[i],std::string("AugSep2012"),std::string("qgsp_bert_hp"));
    //RadialProfile(energy[i],std::string("AugSep2012"),std::string("ftfp_bert"));
    //RadialProfile(energy[i],std::string("AugSep2012"),std::string("qgsp_bert"));
  }
  int energyNov[]={10,20,30,40,50,60,70,80};
  unsigned int eSizeNov=sizeof(energyNov)/sizeof(int);
  for(unsigned int i=0; i<eSizeNov; i++){
    //RadialProfile(energyNov[i],std::string("Nov2012"),std::string("ftfp_bert_hp"));
    //RadialProfile(energyNov[i],std::string("Nov2012"),std::string("qgsp_bert_hp"));
    //RadialProfile(energyNov[i],std::string("Nov2012"),std::string("ftfp_bert"));
    //RadialProfile(energyNov[i],std::string("Nov2012"),std::string("qgsp_bert"));
  }
}
