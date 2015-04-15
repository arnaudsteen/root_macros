#define NhitElectron_cxx
#include "NhitElectron.h"
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

void NhitElectron::ShowerBarycenterCut(std::string runPeriod)
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

void NhitElectron::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if( (float)spillEventTime*200/pow(10.,9)!=0 && (float)spillEventTime*200/pow(10.,9)<0.3 ) continue;
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
//    if( (float)Nhit/Nlayer>3 &&
//	Begin>=0 &&
//	(float)NInteractinglayer/Nlayer>0.2 &&
//	Single==1 &&
//	Neutral==0 &&
//	Hole==0 &&
//	Begin<5 &&
//	Nlayer<30 &&
//	(float)spillEventTime*200/pow(10.,9)<15 &&
//	CoG[0]>cutBary[0] &&
//	CoG[0]<cutBary[1] &&
//	CoG[2]>cutBary[2] &&
//	CoG[2]<cutBary[3] &&
//	ReconstructedCosTheta>.9&&
//	//F3D/log(Nhit)*CentralRatio>cutFD &&
//	TrackMultiplicity==0){
      nhit1.push_back(Nhit1);
      nhit2.push_back(Nhit2);
      nhit3.push_back(Nhit3);
      time.push_back( (float)spillEventTime*200/pow(10.,9) );
    }      
  }
  std::cout << (float)nhit1.size()/nentries << std::endl;
}

std::vector<double> readCalib(int run)
{
  std::vector<double> calib;
  char inputCalib[200];
  sprintf(inputCalib,"%s%d%s","/home/steen/timeCalib/Nhit/calib_",run,".txt");
  fstream in;
  in.open(inputCalib);
  if(!in.is_open()){ std::cout << inputCalib << "\t NO SUCH FILE OR DIRECTORY" << std::endl; throw;}
  double coeff[5];
  double coeffError[5];
  while(1){
    if(!in.good()) break;
    in >> coeff[0] >> coeffError[0] >> coeff[1] >> coeffError[1] >> coeff[2] >> coeffError[2] >> coeff[3] >> coeffError[3] >> coeff[4] >> coeffError[4] ;
    calib.push_back(coeff[0]);
    calib.push_back(coeff[1]);
    calib.push_back(coeff[2]);
    calib.push_back(coeff[3]);
    calib.push_back(coeff[4]);
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

void Nhit(int E, std::string runPeriod, std::string model)
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
  double calib1[5];
  double calib2[5];
  double calib3[5];
  for(unsigned int i=0; i<5; i++){
    calib1[i]=calibrations.at(i);
    calib2[i]=calibrations.at(i+5);
    calib3[i]=calibrations.at(i+10);
    //std::cout << "calib1[" << i << "] = " <<  calib1[i] << "\t" 
    //	      << "calib2[" << i << "] = " <<  calib2[i] << "\t" 
    //	      << "calib3[" << i << "] = " <<  calib3[i] << std::endl;
  }
  //  sprintf(inputFile,"%s%d%s","/home/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  sprintf(inputFile,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/DHCAL_",Run(E,runPeriod),".root");
  std::cout << "inputFile = " << inputFile << std::endl;
  NhitElectron* data=new NhitElectron(inputFile);
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
  float coeff3=calib1[3]+calib2[3]+calib3[3];
  float coeff4=calib1[4]+calib2[4]+calib3[4];
  for(unsigned int i=0; i<n; i++){
    nhitCalib=( data->nhit1.at(i)+data->nhit2.at(i)+data->nhit3.at(i) - 
    		coeff1*data->time.at(i) -
    		coeff2*data->time.at(i)*data->time.at(i) -
    		coeff3*data->time.at(i)*data->time.at(i)*data->time.at(i) -
    		coeff4*data->time.at(i)*data->time.at(i)*data->time.at(i)*data->time.at(i) );
    hDat->Fill(nhitCalib);
  }
  sprintf(fDataName,"%s%s%s","/home/steen/CBFit/histo/nhit",runPeriod.c_str(),"_electron.root");
  if(E==10)
    writeHisto(hDat,fDataName,"RECREATE");
  else 
    writeHisto(hDat,fDataName,"UPDATE");
  hDat->Sumw2();
  integral=hDat->Integral();
  if(integral>0) hDat->Scale(1/integral);


  //sprintf(inputFile,"%s%s%s%d%s","/home/steen/resultRootFile/sim_data/",model.c_str(),"/single_e-_",E,"GeV.root");
  //  sprintf(inputFile,"%s%d%s","/home/steen/sdhcal_analysis/single_e-_",E,"GeV.root");
  sprintf(inputFile,"%s%s%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/",model.c_str(),"/timecut/single_e-_",E,"GeV.root");
  std::cout << "inputFile = " << inputFile << std::endl;
  NhitElectron* simu=new NhitElectron(inputFile);
  simu->Loop();  
  sprintf(hName,"%s%d%s","hNhit_",E,"GeV");
  TH1D *hSim=new TH1D(hName,hName,50,xmin,xmax);
  hSim->GetXaxis()->SetTitle("N_{hit}");
  hSim->GetXaxis()->SetTitleOffset(0.8);
  hSim->GetYaxis()->SetTitle("# events (normalised to unity)");
  hSim->GetYaxis()->SetTitleOffset(1.2);
  SimulationHistoStyleDefinition(hSim,model);
  n=simu->nhit1.size();
  std::cout << n/20000. << std::endl;
  for(unsigned int i=0; i<n; i++){
    hSim->Fill(simu->nhit1.at(i)+simu->nhit2.at(i)+simu->nhit3.at(i));
  }
  sprintf(fSimuName,"%s%s%s","/home/steen/CBFit/histo/nhit",model.c_str(),"_electron.root");
  if(E==10)
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
  if( hSim->GetBinContent( hSim->GetMaximumBin() ) > hDat->GetBinContent( hDat->GetMaximumBin() ) )
    hSim->Draw("axissame");
  else hDat->Draw("axissame");
  TText *tex=new TText();
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  sprintf(text,"%s%d%s","E = ",E," GeV");
  tex->DrawTextNDC(0.11,0.96,text);
  TLegend *leg=new TLegend(0.12,0.67,0.51,0.84);
  leg->SetFillStyle(0);
  std::transform(model.begin(),model.end(),model.begin(),::toupper);
  sprintf(text,"%s%s%s","MC (",model.c_str(),")");
  leg->AddEntry(hSim,text,"f");
  leg->AddEntry(hDat,"Data","p");
  leg->Draw();
  can->Update();
  
  std::transform(model.begin(),model.end(),model.begin(),::tolower);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/nhit_e-_",E,"GeV_",runPeriod.c_str(),".C");
  can->SaveAs(plotName);
  sprintf(plotName,"%s%s%s%d%s%s%s","./plots/",model.c_str(),"/nhit_e-_",E,"GeV_",runPeriod.c_str(),".pdf");
  can->SaveAs(plotName);

  //sleep(5);
  //can->Clear();
  //can->Close();
}

void ProcessElectron()
{
  int energy[]={10,20,30,40,50};
  unsigned int eSize=sizeof(energy)/sizeof(int);
  for(unsigned int i=0; i<eSize; i++){
    Nhit(energy[i],std::string("AugSep2012"),std::string("ftfp_bert_hp"));
    Nhit(energy[i],std::string("AugSep2012"),std::string("qgsp_bert_hp"));
    //Nhit(energy[i],std::string("AugSep2012"),std::string("ftfp_bert"));
    //Nhit(energy[i],std::string("AugSep2012"),std::string("qgsp_bert"));
  }
}
