#define NclusterPion_cxx
#include "NclusterPion.h"
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

void NclusterPion::ShowerBarycenterCut(std::string runPeriod)
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

void NclusterPion::Loop()
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
      ncluster.push_back(Nclusters);
      time.push_back( (float)spillEventTime*200/pow(10.,9) );
    }      
  }
}

std::vector<double> readCalib(int run)
{
  std::vector<double> calib;
  char inputCalib[200];
  sprintf(inputCalib,"%s%d%s","/home/steen/timeCalib/Nclusters/calib_",run,".txt");
  fstream in;
  in.open(inputCalib);
  if(!in.is_open()){ std::cout << inputCalib << "\t NO SUCH FILE OR DIRECTORY" << std::endl; throw;}
  double coeff[2];
  double coeffError[2];
  while(1){
    if(!in.good()) break;
    in >> coeff[0] >> coeffError[0] >> coeff[1] >> coeffError[1] ;
    calib.push_back(coeff[0]);
    calib.push_back(coeff[1]);
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

void Ncluster(int E, std::string runPeriod,std::string model)
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
  std::vector<double> calibrations=readCalib(Run(E,runPeriod));
  double calib[2];
  calib[0]=calibrations.at(0);
  calib[1]=calibrations.at(1);
  sprintf(inputFile,"%s%d%s","/home/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  NclusterPion* data=new NclusterPion(inputFile);
  data->ShowerBarycenterCut(runPeriod);
  data->Loop();  
  Float_t integral;
  sprintf(hName,"%s%d%s","hNcluster_",E,"GeV");
  TH1D *hDat=new TH1D(hName,hName,Nbin(E),xmin,xmax);
  hDat->GetXaxis()->SetTitle("N_{cluster}");
  hDat->GetXaxis()->SetTitleOffset(0.8);
  hDat->SetMarkerColor(kBlack);
  hDat->SetMarkerStyle(20);
  hDat->SetMarkerSize(0.8);
  hDat->SetLineWidth(2);
  unsigned int n=data->ncluster.size();
  float nclusterCalib;
  float coeff=calib[1];
  for(unsigned int i=0; i<n; i++){
    nclusterCalib=( data->ncluster.at(i) - coeff*data->time.at(i) );
    hDat->Fill(nclusterCalib);
  }
  sprintf(fDataName,"%s%s%s","/home/steen/CBFit/histo/ncluster",runPeriod.c_str(),".root");
  if( (E==5&&runPeriod==std::string("AugSep2012")) || (E==10&&runPeriod==std::string("Nov2012")) )
    writeHisto(hDat,fDataName,"RECREATE");
  else 
    writeHisto(hDat,fDataName,"UPDATE");
  hDat->Sumw2();
  integral=hDat->Integral();
  if(integral>0) hDat->Scale(1/integral);


  //sprintf(inputFile,"%s%s%s%d%s","/home/steen/resultRootFile/sim_data/",model.c_str(),"/single_e-_",E,"GeV.root");
  sprintf(inputFile,"%s%s%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/",model.c_str(),"/timecut/single_pi-_",E,"GeV.root");
  std::cout << "inputFile = " << inputFile << std::endl;
  NclusterPion* simu=new NclusterPion(inputFile);
  simu->Loop();  
  sprintf(hName,"%s%d%s","hNcluster_",E,"GeV");
  TH1D *hSim=new TH1D(hName,hName,Nbin(E),xmin,xmax);
  hSim->GetXaxis()->SetTitle("N_{cluster}");
  hSim->GetXaxis()->SetTitleOffset(0.8);
  SimulationHistoStyleDefinition(hSim,model);
  n=simu->ncluster.size();
  for(unsigned int i=0; i<n; i++){
    hSim->Fill(simu->ncluster.at(i));
  }
  sprintf(fSimuName,"%s%s%s","/home/steen/CBFit/histo/ncluster",model.c_str(),"_pi-.root");
  if(E==5)
    writeHisto(hSim,fSimuName,"RECREATE");
  else 
    writeHisto(hSim,fSimuName,"UPDATE");
  hSim->Sumw2();
  integral=hSim->Integral();
  if(integral>0) hSim->Scale(1/integral);

  TCanvas *can=new TCanvas();
  can->SetWindowSize(600,600);
  if( hSim->GetBinContent( hSim->GetMaximumBin() ) > hDat->GetBinContent( hDat->GetMaximumBin() ) )
    hSim->Draw("axis");
  else hDat->Draw("axis");
  hSim->DrawCopy("samehist");
  hDat->DrawCopy("sameP");
  if( hSim->GetBinContent(hSim->GetMaximumBin()) > hDat->GetBinContent(hDat->GetMaximumBin()) )
    hSim->Draw("axissame");
  else 
    hDat->Draw("axissame");
  TText *tex=new TText();
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.46,0.96,"CALICE Fe-SDHCAL Preliminary");
  sprintf(text,"%s%d%s","E = ",E," GeV");
  tex->DrawTextNDC(0.11,0.96,text);
  TLegend *leg=new TLegend(0.12,0.73,0.51,0.9);
  leg->SetFillStyle(0);
  leg->SetHeader(text);
  std::transform(model.begin(),model.end(),model.begin(),::toupper);
  sprintf(text,"%s%s%s","MC (",model.c_str(),")");
  leg->AddEntry(hSim,text,"f");
  leg->AddEntry(hDat,"Data","p");
  leg->Draw();
  can->Update();
  
  std::transform(model.begin(),model.end(),model.begin(),::tolower);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/ncluster_pi-_",E,"GeV_",runPeriod.c_str(),".C");
  can->SaveAs(plotName);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/ncluster_pi-_",E,"GeV_",runPeriod.c_str(),".pdf");
  can->SaveAs(plotName);
}

void ProcessPion()
{
  int energy[]={5,10,15,20,25,30,40,50,60,70,80};
  //int energy[]={80};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    Ncluster(energy[i],std::string("AugSep2012"),std::string("ftfp_bert_hp"));
    Ncluster(energy[i],std::string("AugSep2012"),std::string("qgsp_bert_hp"));
    //Ncluster(energy[i],std::string("AugSep2012"),std::string("ftfp_bert"));
    //Ncluster(energy[i],std::string("AugSep2012"),std::string("qgsp_bert"));
  }
  int energyNov[]={10,20,30,40,50,60,70,80};
  //int energyNov[]={80};
  unsigned int eSizeNov=sizeof(energyNov)/sizeof(int);
  for(unsigned int i=0; i<eSizeNov; i++){
    //Ncluster(energyNov[i],std::string("Nov2012"),std::string("ftfp_bert_hp"));
    //Ncluster(energyNov[i],std::string("Nov2012"),std::string("qgsp_bert_hp"));
    //Ncluster(energyNov[i],std::string("Nov2012"),std::string("ftfp_bert"));
    //Ncluster(energyNov[i],std::string("Nov2012"),std::string("qgsp_bert"));
  }
}
