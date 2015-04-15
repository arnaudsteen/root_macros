#include <TROOT.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TF1.h>
#include "iostream"
#include "fstream"
#include "algorithm"
#include "map"
#include "stdio.h"
#include "cmath"

void CaliceStyle()
{
  /*CALICE style for figure: use in a ROOT macro like this:*/
  //gROOT->ProcessLine(".L ~/RootStuff/CaliceStyle.C");
  //CaliceStyle();

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetStatFont(42);

  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetCanvasColor(kWhite);  
  gStyle->SetOptStat(0); /*don't show statistics box*/
  gStyle->SetTitleSize(0.05, "xyz"); 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetOptTitle(0);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.05);

  gROOT->ForceStyle();
}
double Polya(double *x, double *par)
{
  double t=1+par[1];
  double q=x[0]/par[0];
  return pow(q*t,par[1])*exp(-q*t);
}

Double_t int_polya(Double_t *x, Double_t *par)
{
  Double_t t = x[0];
  Double_t a = par[0];
  Double_t b = par[1];
  Double_t c = par[2]; // inificiency .. 
  Double_t e = par[3]; // inificiency .. 
  TF1 f("f",Polya,0,1000,3);
  f.SetParameter(0,a);
  f.SetParameter(1,b);
  //  f.SetParameter(2,c);
  
  return e - c*f.Integral(0, t);
}

float charge(int DAQ,int thr)
{
  if(thr==1) return (DAQ-90)/700.0;
  else if(thr==2) return (DAQ-98)/80.0;
  else if(thr==3) return (DAQ-98)/16.3;
  else return 0;
}

struct ThresholdScanStruct
{
  int threshold;
  int daq;
  float efficiency;
  float error_efficiency;
  int nevent;
  float chargeThr;
  int nchamber;
  float multiplicity;
  float error_multiplicity;
  ThresholdScanStruct() : efficiency(0), error_efficiency(0), nevent(0), nchamber(0), multiplicity(0), error_multiplicity(0){}
};

TGraphErrors* DataGraph()
{
  int Run[]={715763,715766,715768,715770,715772,715773,/*715774,*/715775,715776,715777,
	     715778,715779,715780,715781,715782,715783,715784,715785,715786,715787};
  std::vector<int> runVec(Run,Run+sizeof(Run)/sizeof(int));
  fstream in;
  char input[200];
  int DAQ;
  int chamber;
  int nevent;
  int threshold;
  float efficiency;
  float multiplicity;
  std::map<int,ThresholdScanStruct> structMap;
  for(std::vector<int>::iterator run=runVec.begin(); run!=runVec.end(); ++run){
    sprintf(input,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/thrScan/dataThresholdScan-",(*run),".txt");
    in.open(input);
    while ( 1 ){
      if(!in.good())break;
      in >> chamber >> threshold >> DAQ >> nevent >> efficiency >> multiplicity;
      if(chamber==15||chamber==33||chamber==9||chamber==5||chamber==13||chamber==25)continue;
      int key=threshold*1000+DAQ;
      if(structMap.find(key)!=structMap.end()){
	structMap[key].nevent+=nevent;
	structMap[key].efficiency+=efficiency;//*nevent;
	structMap[key].error_efficiency+=efficiency*efficiency;//*nevent;
	structMap[key].nchamber++;
      }
      else{
	structMap[key].threshold=threshold;
	structMap[key].daq=DAQ;
	structMap[key].nevent=nevent;
	structMap[key].efficiency=efficiency;//*nevent;
	structMap[key].error_efficiency=efficiency*efficiency;//*nevent;
	structMap[key].nchamber=1;
      }
      if(chamber==37)break;
    }
    in.close();
  }  
  Int_t index=0;
  TGraphErrors *gr=new TGraphErrors();
  float Charge;
  for(std::map<int,ThresholdScanStruct>::iterator it=structMap.begin(); it!=structMap.end(); ++it){
    it->second.efficiency=it->second.efficiency/it->second.nchamber;//.nevent;
    it->second.error_efficiency=sqrt(it->second.error_efficiency/it->second.nchamber-it->second.efficiency*it->second.efficiency + it->second.efficiency*(1-it->second.efficiency)/it->second.nevent);
    //    it->second.error_efficiency=sqrt(it->second.error_efficiency/it->second.nevent/it->second.nchamber-it->second.efficiency*it->second.efficiency/it->second.nchamber + it->second.efficiency*(1-it->second.efficiency)/it->second.nevent);
    Charge=charge(it->second.daq,it->second.threshold);
    if(it->second.threshold==1&&Charge>0.45)continue;
    if(it->second.threshold==2&&(Charge>4.1||Charge<0.45))continue;
    if(it->second.threshold==3&&(Charge<4.1||Charge>25.))continue;
    gr->SetPoint(index,Charge,it->second.efficiency);
    gr->SetPointError(index,0.,it->second.error_efficiency);
    index++;
  }
  return gr;
}

TGraphErrors* SimGraph()
{
  std::vector<ThresholdScanStruct> structVec;
  ThresholdScanStruct aThresholdScanStruct;
  int nevent;
  float threshold;
  float efficiency;
  float multiplicity;
  float error_multiplicity;
  fstream input;
  input.open("/home/steen/sdhcal_analysis/simThresholdScan.txt");
  while(1){
    if(!input.good())break;
    input >> threshold >> nevent >> efficiency >> multiplicity >> error_multiplicity;
    aThresholdScanStruct.chargeThr=threshold;
    aThresholdScanStruct.nevent=nevent;
    aThresholdScanStruct.efficiency=efficiency;
    aThresholdScanStruct.error_efficiency=sqrt(aThresholdScanStruct.efficiency*(1-aThresholdScanStruct.efficiency)/aThresholdScanStruct.nevent);
    structVec.push_back(aThresholdScanStruct);
  }
  input.close();
  
  structVec.pop_back();

  //for(std::vector<ThresholdScanStruct>::iterator it=structVec.begin(); it!=structVec.end(); ++it)
  //  std::cout << (*it).chargeThr << " " << (*it).nevent << " " << (*it).efficiency << " " << (*it).error_efficiency << std::endl;
  
  Int_t index=0;
  TGraphErrors *gr=new TGraphErrors();
  for(std::vector<ThresholdScanStruct>::iterator it=structVec.begin(); it!=structVec.end(); ++it){
    gr->SetPoint(index,(*it).chargeThr,(*it).efficiency);
    gr->SetPointError(index,0.,(*it).error_efficiency);
    index++;
  }
  return gr;
}

void data_efficiency()
{
  CaliceStyle();
  TCanvas *r1 = new TCanvas();
  r1->Size(800,800);  
  r1->SetWindowSize(800,800);
  r1->SetLogx();
  
  TH1D *he = new TH1D("he"," ",300,0.09,30);
  he->SetMinimum(0.0001);
  he->SetMaximum(1.0);
  he->SetStats(0);
  he->GetXaxis()->SetTitle("Threshold [pC]");
  he->GetXaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetTitle("Efficiency ");
  he->GetYaxis()->SetTitleOffset(1.20);
  he->GetXaxis()->SetTitleOffset(1.);
  he->GetYaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetNdivisions(505);
  he->Draw();

  TGraphErrors* grDat=DataGraph();
  TF1 *func = new TF1("func",int_polya,0.001,30,4);
  func->SetLineWidth(2);
  func->SetLineColor(kBlack);
  func->SetLineStyle(1);
  func->SetParameter(0,4.846);
  func->SetParameter(1,0.7591);
  func->SetParameter(2,0.3758);
  func->SetParameter(3,0.9559);

  grDat->SetMarkerStyle(34);
  grDat->SetMarkerColor(kBlack);
  grDat->SetLineWidth(2);
  grDat->SetLineColor(kBlack);
  grDat->SetMarkerSize(1.2);
  grDat->Draw("psame");
  grDat->Fit(func);
  
  TLatex *tex=new TLatex();
  tex->SetTextSize(0.03);
  tex->SetTextColor(kBlack);
  char myText[200];
  Double_t xPos=0.2;
  Double_t yPos=0.2;
  Double_t delta=0.06;
  tex->DrawLatex(xPos,yPos+4*delta,"Data (SPS H6)");
  sprintf(myText,"%s%.3f%s%.3f%s","#bar{q} = ",func->GetParameter(0)," #pm ",func->GetParError(0)," pC");
  tex->DrawLatex(xPos,yPos+3*delta,myText);
  sprintf(myText,"%s%.3f%s%.3f","#theta = ",func->GetParameter(1)," #pm ",func->GetParError(1));
  tex->DrawLatex(xPos,yPos+2*delta,myText);
  sprintf(myText,"%s%.3f%s%.3f","c = ",func->GetParameter(2)," #pm ",func->GetParError(2));
  tex->DrawLatex(xPos,yPos+1*delta,myText);
  sprintf(myText,"%s%.3f%s%.3f","#varepsilon_{0} = ",func->GetParameter(3)," #pm ",func->GetParError(3));
  tex->DrawLatex(xPos,yPos+0*delta,myText);
  
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.5,0.96,"CALICE Fe-SDHCAL Preliminary");
  r1->SaveAs("./plots/thrScanDat.pdf");
  r1->SaveAs("./plots/thrScanDat.C");

  std::cout << "data chi2 = " << func->GetChisquare()/func->GetNDF() << std::endl;
}

void simu_efficiency()
{
  CaliceStyle();
  TCanvas *r1 = new TCanvas();
  r1->Size(800,800);  
  r1->SetWindowSize(800,800);
  r1->SetLogx();
  
  TH1D *he = new TH1D("he"," ",30,0.09,30);
  he->SetMinimum(0.0001);
  he->SetMaximum(1.0);
  he->SetStats(0);
  he->GetXaxis()->SetTitle("Threshold [pC]");
  he->GetXaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetTitle("Efficiency ");
  he->GetYaxis()->SetTitleOffset(1.20);
  he->GetXaxis()->SetTitleOffset(1.);
  he->GetYaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetNdivisions(505);
  he->Draw();

  TGraphErrors* grSim=SimGraph();
  TF1 *func = new TF1("func",int_polya,0.001,30,4);
  func->SetLineWidth(2);
  func->SetLineColor(kRed-3);
  func->SetLineStyle(1);
  func->SetParameter(0,4.846);
  func->SetParameter(1,0.7591);
  func->SetParameter(2,0.3758);
  func->SetParameter(3,0.9559);
  
  grSim->SetMarkerStyle(21);
  grSim->SetMarkerColor(kRed-3);
  grSim->SetLineWidth(2);
  grSim->SetLineColor(kRed-3);
  grSim->SetMarkerSize(0.8);
  grSim->Draw("psame");
  grSim->Fit(func);

  TLatex *tex=new TLatex();
  tex->SetTextSize(0.03);
  tex->SetTextColor(kRed-3);
  char myText[200];
  Double_t xPos=0.2;
  Double_t yPos=0.2;
  Double_t delta=0.06;
  tex->DrawLatex(xPos,yPos+4*delta,"Simulation (FTFP_BERT_HP)");
  sprintf(myText,"%s%.3f%s%.3f%s","#bar{q} = ",func->GetParameter(0)," #pm ",func->GetParError(0)," pC");
  tex->DrawLatex(xPos,yPos+3*delta,myText);
  sprintf(myText,"%s%.3f%s%.3f","#theta = ",func->GetParameter(1)," #pm ",func->GetParError(1));
  tex->DrawLatex(xPos,yPos+2*delta,myText);
  sprintf(myText,"%s%.3f%s%.3f","c = ",func->GetParameter(2)," #pm ",func->GetParError(2));
  tex->DrawLatex(xPos,yPos+1*delta,myText);
  sprintf(myText,"%s%.3f%s%.3f","#varepsilon_{0} = ",func->GetParameter(3)," #pm ",func->GetParError(3));
  tex->DrawLatex(xPos,yPos+0*delta,myText);
  
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.5,0.96,"CALICE Fe-SDHCAL Preliminary");
  r1->SaveAs("./plots/thrScanSim.pdf");
  r1->SaveAs("./plots/thrScanSim.C");
  std::cout << "simu chi2 = " << func->GetChisquare()/func->GetNDF() << std::endl;
  
}

TGraphErrors* MultiSimGraph()
{
  std::vector<ThresholdScanStruct> structVec;
  ThresholdScanStruct aThresholdScanStruct;
  int nevent;
  float threshold;
  float efficiency;
  float multiplicity;
  float errormultiplicity;
  fstream input;
  input.open("/home/steen/sdhcal_analysis/simThresholdScanForMulti.txt");
  while(1){
    if(!input.good())break;
    input >> threshold >> nevent >> efficiency >> multiplicity >> errormultiplicity;
    aThresholdScanStruct.chargeThr=threshold;
    aThresholdScanStruct.nevent=nevent;
    aThresholdScanStruct.efficiency=efficiency;
    aThresholdScanStruct.multiplicity=multiplicity;
    aThresholdScanStruct.error_multiplicity=errormultiplicity;
    aThresholdScanStruct.error_efficiency=sqrt(aThresholdScanStruct.efficiency*(1-aThresholdScanStruct.efficiency)/aThresholdScanStruct.nevent);
    structVec.push_back(aThresholdScanStruct);
  }
  input.close();
  
  structVec.pop_back();

  Int_t index=0;
  TGraphErrors *gr=new TGraphErrors();
  for(std::vector<ThresholdScanStruct>::iterator it=structVec.begin(); it!=structVec.end(); ++it){
    //if((*it).chargeThr>5.0)break;
    gr->SetPoint(index,(*it).multiplicity,(*it).chargeThr);
    std::cout << (*it).chargeThr << " " << (*it).multiplicity << std::endl;
    gr->SetPointError(index,(*it).error_multiplicity,0.);
    index++;
  }
  return gr;
}

double multiFunc(double *x, double *par)
{
  //return par[0]*exp(x[0]*par[1])+par[2];
  return par[0]/(x[0]*x[0]*x[0])+par[1];
}

void simu_multiplicity()
{
  CaliceStyle();
  TCanvas *r1 = new TCanvas();
  r1->Size(800,800);  
  r1->SetWindowSize(800,800);
  //  r1->SetLogx();
  
  TH1D *he = new TH1D("he"," ",100,1.,5);
  he->SetMinimum(0);
  he->SetMaximum(1.0);
  he->SetStats(0);
  he->GetXaxis()->SetTitle("Threshold [pC]");
  he->GetXaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetTitle("Multiplicity");
  he->GetYaxis()->SetTitleOffset(1.20);
  he->GetXaxis()->SetTitleOffset(1.);
  he->GetYaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetNdivisions(505);
  he->Draw();

  TGraphErrors* grSim=MultiSimGraph();
  TF1 *func = new TF1("func",multiFunc,1.,5,2);
  func->SetLineWidth(2);
  func->SetLineColor(kBlue);
  func->SetLineStyle(1);
  //func->SetParameter(0,1.0);
  //func->SetParameter(1,-1.0);
  //func->SetParameter(2,.05);
  func->SetParameter(0,1.0);
  func->SetParameter(1,1.0);
  
  grSim->SetMarkerStyle(21);
  grSim->SetMarkerColor(kRed-3);
  grSim->SetLineWidth(2);
  grSim->SetLineColor(kRed-3);
  grSim->SetMarkerSize(0.8);
  grSim->Draw("psame");
  grSim->Fit("func","","",1.4,5.);

  TLatex *tex=new TLatex();
  tex->SetTextSize(0.03);
  tex->SetTextColor(kRed-3);
  char myText[200];
  Double_t xPos=0.2;
  Double_t yPos=0.2;
  Double_t delta=0.06;
  tex->DrawLatex(xPos,yPos+4*delta,"Simulation (FTFP_BERT_HP)");
  //sprintf(myText,"%s%.3f%s%.3f%s","#bar{q} = ",func->GetParameter(0)," #pm ",func->GetParError(0)," pC");
  //tex->DrawLatex(xPos,yPos+3*delta,myText);
  //sprintf(myText,"%s%.3f%s%.3f","#theta = ",func->GetParameter(1)," #pm ",func->GetParError(1));
  //tex->DrawLatex(xPos,yPos+2*delta,myText);
  //sprintf(myText,"%s%.3f%s%.3f","c = ",func->GetParameter(2)," #pm ",func->GetParError(2));
  //tex->DrawLatex(xPos,yPos+1*delta,myText);
  //sprintf(myText,"%s%.3f%s%.3f","#varepsilon_{0} = ",func->GetParameter(3)," #pm ",func->GetParError(3));
  //tex->DrawLatex(xPos,yPos+0*delta,myText);
  
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.5,0.96,"CALICE Fe-SDHCAL Preliminary");
  //r1->SaveAs("./plots/thrScanSim.pdf");
  //r1->SaveAs("./plots/thrScanSim.C");
  //std::cout << "simu chi2 = " << func->GetChisquare()/func->GetNDF() << std::endl;
  
}

void efficiency()
{
  CaliceStyle();
  TCanvas *r1 = new TCanvas("r1", "ThrScan");
  r1->Size(600,600);  
  r1->SetWindowSize(600,600);
  r1->SetLogx();
  
  TH1D *he = new TH1D("he"," ",30,0,30);
  he->SetMinimum(0.0001);
  he->SetMaximum(1.0);
  he->SetStats(0);
  he->GetXaxis()->SetTitle("Threshold [pC]");
  he->GetXaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetTitle("Efficiency ");
  he->GetYaxis()->SetTitleOffset(1.20);
  he->GetXaxis()->SetTitleOffset(1.);
  he->GetYaxis()->SetTitleSize(0.04);
  he->GetYaxis()->SetNdivisions(505);
  he->Draw();

  TGraphErrors* grDat=DataGraph();
  TF1 *func = new TF1("func",int_polya,0.001,30,4);
  func->SetLineWidth(2);
  func->SetLineColor(kBlack);
  func->SetLineStyle(1);
  func->SetParameter(0,4.846);
  func->SetParameter(1,0.7591);
  func->SetParameter(2,0.3758);
  func->SetParameter(3,0.9559);

  grDat->SetMarkerStyle(34);
  grDat->SetMarkerColor(kBlack);
  grDat->SetLineWidth(2);
  grDat->SetLineColor(kBlack);
  grDat->SetMarkerSize(1.2);

  TGraphErrors* grSim=SimGraph();
  TF1 *funcSim = new TF1("funcSim",int_polya,0.001,30,4);
  funcSim->SetLineWidth(2);
  funcSim->SetLineColor(kRed-3);
  funcSim->SetLineStyle(1);
  funcSim->SetParameter(0,4.846);
  funcSim->SetParameter(1,0.7591);
  funcSim->SetParameter(2,0.3758);
  funcSim->SetParameter(3,0.9559);
  
  grSim->SetMarkerStyle(21);
  grSim->SetMarkerColor(kRed-3);
  grSim->SetLineWidth(2);
  grSim->SetLineColor(kRed-3);
  grSim->SetMarkerSize(0.8);
  grSim->Draw("psame");
  grSim->Fit(funcSim);
  grDat->Draw("psame");
  grDat->Fit(func);

  //TLegend *leg=new TLegend(0.60,0.74,.85,0.89);
  //leg->SetBorderSize(0);
  //TLegendEntry *entry=leg->AddEntry("grData","DATA","lp");
  //entry->SetLineColor(1);
  //entry->SetLineStyle(1);
  //entry->SetLineWidth(2);
  //entry->SetMarkerColor(kBlack);
  //entry->SetMarkerStyle(34);
  //entry->SetMarkerSize(1.2);
  //entry=leg->AddEntry("gr","Simulation","lp");
  //entry->SetLineColor(kRed-3);
  //entry->SetLineStyle(1);
  //entry->SetLineWidth(2);
  //entry->SetMarkerColor(kRed-3);
  //entry->SetMarkerStyle(20);
  //entry->SetMarkerSize(1.);
  //leg->SetFillStyle(0);
  //leg->SetLineWidth(0);
  //leg->Draw();
  
  TLatex *tex=new TLatex();
  tex->SetTextSize(0.03);
  tex->SetTextColor(kRed-3);
  char myText[200];
  Double_t xPos=0.13;
  Double_t yPos=0.15;
  Double_t delta=0.06;
  // tex->DrawLatex(xPos,yPos+4*delta,"MC fit");
  sprintf(myText,"%s%.2f%s","#bar{Q_{ind}} = ",funcSim->GetParameter(0)," pC");
  tex->DrawLatex(xPos,yPos+3*delta,myText);
  sprintf(myText,"%s%.2f","#theta = ",funcSim->GetParameter(1));
  tex->DrawLatex(xPos,yPos+2*delta,myText);
  sprintf(myText,"%s%.2f","c = ",funcSim->GetParameter(2));
  tex->DrawLatex(xPos,yPos+1*delta,myText);
  sprintf(myText,"%s%.2f","#varepsilon_{0} = ",funcSim->GetParameter(3));
  tex->DrawLatex(xPos,yPos+0*delta,myText);
  //sprintf(myText,"%s%.2f","#chi^{2} = ",funcSim->GetChisquare());
  //tex->DrawLatex(0.02,yPos-4*delta,myText);
  
  xPos=0.7;
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(xPos,yPos+4*delta,"Data fit");
  sprintf(myText,"%s%.2f%s","#bar{Q_{ind}} = ",func->GetParameter(0)," pC");
  tex->DrawLatex(xPos,yPos+3*delta,myText);
  sprintf(myText,"%s%.2f","#theta = ",func->GetParameter(1));
  tex->DrawLatex(xPos,yPos+2*delta,myText);
  sprintf(myText,"%s%.2f","c = ",func->GetParameter(2));
  tex->DrawLatex(xPos,yPos+1*delta,myText);
  sprintf(myText,"%s%.2f","#varepsilon_{0} = ",func->GetParameter(3));
  tex->DrawLatex(xPos,yPos+0*delta,myText);
  
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.44,0.96,"CALICE Fe-SDHCAL Preliminary");

  std::cout << "data chi2 = " << func->GetChisquare()/func->GetNDF() << std::endl;
  std::cout << "sim chi2 = " << funcSim->GetChisquare()/funcSim->GetNDF() << std::endl;
}
