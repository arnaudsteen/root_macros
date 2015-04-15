
#define NhitPion_cxx
#include "NhitPion.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstring>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>

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

void NhitPion::ShowerBarycenterCut(std::string runPeriod)
{
  if(runPeriod==std::string("AugSep2012")){
    cutBary[0]=400;
    cutBary[1]=800;
    cutBary[2]=300;
    cutBary[3]=700;
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

void NhitPion::Loop()
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
	TransverseRatio>0.01&&
	(float)NInteractinglayer/Nlayer>0.01&&
	Single==1 &&
	Neutral==0 &&
	(float)spillEventTime*200/pow(10.,9)<15 &&
	CoG[0]>cutBary[0] &&
	CoG[0]<cutBary[1] &&
	CoG[2]>cutBary[2] &&
	CoG[2]<cutBary[3] &&
	ReconstructedCosTheta>.9&&
	(Nlayer>=30||Begin>=5||TrackMultiplicity>0)){
//    if( 
//       (float)Nhit/Nlayer>3&&
//       (float)NInteractinglayer/Nlayer>.2&&
//       Begin>=0&&
//       //(((float)Nhit/Nlayer>3 && Begin>=0 && (float)NInteractinglayer/Nlayer>0.2)||(Nlayer<40))&&
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
//       //TrackClusterNumber->at(0)>5&&
//       TransverseRatio>0.02
//       //(float)(Nhough1+Nhough2+Nhough3)/Nhit<0.5) 
//	){
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      time.push_back( (float)spillEventTime*200/pow(10.,9) );
    }      
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
  }
  return calib;
}

void writeHisto(TH1D* h,std::string fileName,std::string option)
{
  TFile *file=new TFile(fileName.c_str(),option.c_str());
  file->cd();
  h->Write();
  file->Write();
  file->Close();
}

void SimulationHistoStyleDefinition(TH1D* h,std::string model)
{
  if( model==std::string("ftfp_bert_hp") ){
    h->SetFillColor(kRed-3);
    h->SetFillStyle(1001);
    h->SetLineColor(kRed-3);
    h->SetLineWidth(1);
  }
  else if( model==std::string("qgsp_bert_hp") ){
    h->SetFillColor(kBlue-6);
    h->SetFillStyle(1001);
    h->SetLineColor(kBlue-6);
    h->SetLineWidth(1);
  }
  else if( model==std::string("ftfp_bert") ){
    h->SetFillColor(kRed+1);
    h->SetFillStyle(1001);
    h->SetLineColor(kRed+1);
    h->SetLineWidth(1);
  }
  else if( model==std::string("qgsp_bert") ){
    h->SetFillColor(kBlue+2);
    h->SetFillStyle(1001);
    h->SetLineColor(kBlue+2);
    h->SetLineWidth(1);
  }
  else 
    std::cout << "WRONG CHOSEN MODEL OPTION" << std::endl;
}

void Nhit(int E, std::string runPeriod,std::string model)
{
  char inputFile[200];
  char hName[100];
  char fDataName[200];
  char text[100];
  char fSimuName[200];
  char plotName[100];

  CaliceStyle();
  int xmin=Xmin(E);
  int xmax=Xmax(E);
  std::vector<float> calibrations=readCalib(Run(E,runPeriod));
  float calib1[3];
  float calib2[3];
  float calib3[3];
  for(unsigned int i=0; i<3; i++){
    calib1[i]=calibrations.at(i);
    calib2[i]=calibrations.at(i+3);
    calib3[i]=calibrations.at(i+6);
  }
  //  sprintf(inputFile,"%s%d%s","/home/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  NhitPion* data=new NhitPion(inputFile);
  data->ShowerBarycenterCut(runPeriod);
  data->Loop();  
  Float_t integral;
  sprintf(hName,"%s%d%s","hNhit_",E,"GeV");
  TH1D *hDat=new TH1D(hName,hName,50,xmin,xmax);
  hDat->GetXaxis()->SetTitle("N_{hit}");
  hDat->GetYaxis()->SetTitle("# events (normalised to unity)");
  hDat->GetYaxis()->SetTitleOffset(1.2);
  hDat->GetXaxis()->SetTitleOffset(0.8);
  hDat->SetMarkerColor(kBlack);
  hDat->SetMarkerStyle(20);
  hDat->SetMarkerSize(0.8);
  hDat->SetLineWidth(2);
  unsigned int n=data->nhit1.size();
  float nhitCalib;
  float coeff1=calib1[1]+calib2[1]+calib3[1];
  float coeff2=calib1[2]+calib2[2]+calib3[2];
  for(unsigned int i=0; i<n; i++){
    nhitCalib=( data->nhit1.at(i)+data->nhit2.at(i)+data->nhit3.at(i) - 
    		coeff1*data->time.at(i) -
    		coeff2*data->time.at(i)*data->time.at(i) );
    hDat->Fill(nhitCalib);
  }
  sprintf(fDataName,"%s%s%s","/home/steen/CBFit/histo/nhit",runPeriod.c_str(),".root");
  if( (E==5&&runPeriod==std::string("AugSep2012")) || (E==10&&runPeriod==std::string("Nov2012")) )
    writeHisto(hDat,fDataName,"RECREATE");
  else 
    writeHisto(hDat,fDataName,"UPDATE");
  hDat->Sumw2();
  integral=hDat->Integral();
  if(integral>0) hDat->Scale(1/integral);


  //sprintf(inputFile,"%s%d%s","/home/steen/sdhcal_analysis/single_pi-_",E,"GeV.root");
  //sprintf(inputFile,"%s%s%s%d%s","/home/steen/resultRootFile/sim_data/",model.c_str(),"/single_pi-_",E,"GeV.root");
  sprintf(inputFile,"%s%s%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/",model.c_str(),"/timecut/single_pi-_",E,"GeV.root");
  std::cout << "inputFile = " << inputFile << std::endl;
  NhitPion* simu=new NhitPion(inputFile);
  simu->Loop();  
  sprintf(hName,"%s%d%s","hNhit_",E,"GeV");
  TH1D *hSim=new TH1D(hName,hName,50,xmin,xmax);
  hSim->GetXaxis()->SetTitle("N_{hit}");
  hSim->GetYaxis()->SetTitle("# events (normalised to unity)");
  hSim->GetYaxis()->SetTitleOffset(1.2);
  hSim->GetXaxis()->SetTitleOffset(0.8);
  SimulationHistoStyleDefinition(hSim,model);
  n=simu->nhit1.size();
  for(unsigned int i=0; i<n; i++){
    hSim->Fill(simu->nhit1.at(i)+simu->nhit2.at(i)+simu->nhit3.at(i));
  }
  sprintf(fSimuName,"%s%s%s","/home/steen/CBFit/histo/nhit",model.c_str(),"_pi-.root");
  if(E==5)
    writeHisto(hSim,fSimuName,"RECREATE");
  else 
    writeHisto(hSim,fSimuName,"UPDATE");
  hSim->Sumw2();
  integral=hSim->Integral();
  if(integral>0) hSim->Scale(1/integral);

  TCanvas *can=new TCanvas();
  can->SetWindowSize(600,600);
  if( hSim->GetBinContent(hSim->GetMaximumBin()) > hDat->GetBinContent(hDat->GetMaximumBin()) )
    hSim->Draw("axis");
  else 
    hDat->Draw("axis");
  hSim->Draw("histsame");
  hDat->Draw("same");
  if( hSim->GetBinContent(hSim->GetMaximumBin()) > hDat->GetBinContent(hDat->GetMaximumBin()) )
    hSim->Draw("axissame");
  else 
    hDat->Draw("axissame");
  TText *tex=new TText();
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  sprintf(text,"%s%d%s","E = ",E," GeV");
  tex->DrawTextNDC(0.13,0.96,text);
  TLegend *leg=new TLegend(0.12,0.67,0.51,0.84);
  leg->SetFillStyle(0);
  std::transform(model.begin(),model.end(),model.begin(),::toupper);
  sprintf(text,"%s%s%s","MC (",model.c_str(),")");
  leg->AddEntry(hSim,text,"f");
  leg->AddEntry(hDat,"Data","p");
  leg->Draw();
  can->Update();
  
  std::transform(model.begin(),model.end(),model.begin(),::tolower);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/nhit_pi-_",E,"GeV_",runPeriod.c_str(),".C");
  can->SaveAs(plotName);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/nhit_pi-_",E,"GeV_",runPeriod.c_str(),".pdf");
  can->SaveAs(plotName);
  //std::cout << E << " " << model.c_str() << " " << (float)simu->nhit1.size()/1000 << " " << sqrt((float)simu->nhit1.size()/1000*(1-(float)simu->nhit1.size()/1000)/1000) << std::endl;
  std::cout << E << " " << model.c_str() << " " << (float)simu->nhit1.size()/20000 << " " << sqrt((float)simu->nhit1.size()/20000*(1-(float)simu->nhit1.size()/20000)/20000) << std::endl;

  // sleep(2);
  // can->Clear();
  // can->Close();
  
}

void ProcessPion()
{
   int energy[]={5,10,15,20,25,30,40,50,60,70,80};
   //int energy[]={10,30,50};
   unsigned int eSize=sizeof(energy)/sizeof(int);
   for(unsigned int i=0; i<eSize; i++){
     Nhit(energy[i],std::string("AugSep2012"),std::string("ftfp_bert_hp"));
     Nhit(energy[i],std::string("AugSep2012"),std::string("qgsp_bert_hp"));
     //Nhit(energy[i],std::string("AugSep2012"),std::string("ftfp_bert"));
     //Nhit(energy[i],std::string("AugSep2012"),std::string("qgsp_bert"));
   }
   //int energyNov[]={10,20,30,40,50,60,70,80};
   //unsigned int eSizeNov=sizeof(energyNov)/sizeof(int);
   //for(unsigned int i=0; i<eSizeNov; i++){
   //  Nhit(energyNov[i],std::string("Nov2012"),std::string("ftfp_bert_hp"));
   //  Nhit(energyNov[i],std::string("Nov2012"),std::string("qgsp_bert_hp"));
   //  Nhit(energyNov[i],std::string("Nov2012"),std::string("ftfp_bert"));
   //  Nhit(energyNov[i],std::string("Nov2012"),std::string("qgsp_bert"));
   //}
}
