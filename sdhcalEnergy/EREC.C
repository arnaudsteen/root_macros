#include "EREC.h"

std::vector<result> readTXT(std::string fileName)
{
  std::vector<result> vec;
  char input[200];
  sprintf(input,"%s%s","./txtFile/",fileName.c_str());
  fstream in;
  in.open(input);
  int ebeam;
  float erec;
  float erecError;
  float sigma;
  float sigmaError;
  float resol;
  while(1){
    in >> ebeam >> erec >> erecError >> sigma >> sigmaError >> resol;
    result res;
    res.ebeam=ebeam;
    res.erec=erec;
    res.erecError=erecError;
    res.sigma=sigma;
    res.sigmaError=sigmaError;
    res.resol=resol;
    vec.push_back(res);
    if(ebeam==80)break;
  }
  return vec;
}

void DrawEREC()
{
  //=========Macro generated from canvas: r1/Energy
  //=========  (Fri Apr  6 17:14:48 2012) by ROOT version5.28/00g
  TCanvas *r1 = new TCanvas("r1", "Erec (hp phys list)",12,24,550,741);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(1);
  r1->Range(0,0,1,1);
  r1->SetFillColor(0);
  r1->SetBorderMode(0);
  r1->SetBorderSize(0);
  r1->SetTickx(1);
  r1->SetTicky(1);
  r1->SetLeftMargin(0.16);
  r1->SetRightMargin(0.01);
  r1->SetTopMargin(0.0256917);
  r1->SetBottomMargin(0.07692308);
  r1->SetFrameBorderMode();
  
  // ------------>Primitives in pad: r1_1
  TPad *r1_1 = new TPad("r1_1", "Energy_1",0.02,0.37,0.95,0.99);
  r1_1->Draw();
  r1_1->cd();
  r1_1->Range(-19,0.01,95,95);
  r1_1->SetFillColor(0);
  r1_1->SetBorderMode(0);
  r1_1->SetBorderSize(2);
  r1_1->SetTickx(1);
  r1_1->SetTicky(1);
  r1_1->SetLeftMargin(0.16);
  r1_1->SetRightMargin(0.01);
  r1_1->SetTopMargin(0.02);
  r1_1->SetBottomMargin(0);
  r1_1->SetFrameBorderMode(0);
  r1_1->SetFrameBorderMode(0);
   
  TH1D *he = new TH1D("he"," ",85,0,85);
  he->SetMinimum(0.01);
  he->SetMaximum(90);
  he->SetStats(0);
  he->GetXaxis()->SetTitle("E_{beam} [GeV]");
  he->GetXaxis()->SetLabelFont(43);
  he->GetXaxis()->SetLabelSize(0);
  he->GetXaxis()->SetTitleFont(43);
  he->GetXaxis()->SetTitleSize(0); 
  he->GetYaxis()->SetTitle("<E_{reco}>");
  he->GetYaxis()->SetLabelFont(43);
  he->GetYaxis()->SetTitleSize(30);
  he->GetYaxis()->SetLabelSize(20);
  he->GetYaxis()->SetTitleFont(43);
  he->GetYaxis()->SetTitleOffset(1.7);
  he->GetZaxis()->SetLabelFont(42);
  he->GetZaxis()->SetTitleSize(0.05);
  he->GetZaxis()->SetTitleFont(42);
  he->Draw("");
   
  TF1* linearity=new TF1("linearity","x",0,85);
  linearity->SetLineWidth(2);
  linearity->SetLineColor(kBlack);
  linearity->SetLineStyle(7);
  linearity->Draw("same");

  std::vector<result> resultLinear=readTXT(std::string("LinearErec.txt"));
  TGraphErrors *gre = new TGraphErrors(resultLinear.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetLineColor(1);
  gre->SetFillStyle(1);
  gre->SetFillColor(1);
  gre->SetLineWidth(2);
  gre->SetMarkerColor(1);
  gre->SetMarkerStyle(34);
  gre->SetMarkerSize(1.2);
  for(unsigned int i=0; i<resultLinear.size(); i++){
    std::cout << resultLinear.at(i).ebeam << " " << resultLinear.at(i).erec << " " << resultLinear.at(i).erecError << std::endl;
    gre->SetPoint(i,resultLinear.at(i).ebeam,resultLinear.at(i).erec);
    gre->SetPointError(i,0,resultLinear.at(i).erecError);
  }

  TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,0,87.5);
  Graph_Graph3->SetMinimum(0);
  Graph_Graph3->SetMaximum(1193.483);
  Graph_Graph3->SetDirectory(0);
  Graph_Graph3->SetStats(0);
  Graph_Graph3->GetXaxis()->SetLabelFont(42);
  Graph_Graph3->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph3->GetXaxis()->SetTitleFont(42);
  Graph_Graph3->GetYaxis()->SetLabelFont(42);
  Graph_Graph3->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph3->GetYaxis()->SetTitleFont(42);
  Graph_Graph3->GetZaxis()->SetLabelFont(42);
  Graph_Graph3->GetZaxis()->SetTitleSize(0.05);
  Graph_Graph3->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph_Graph3);
   
  gre->Draw("p");
   
  std::vector<result> resultQuadratic=readTXT(std::string("QuadraticErec.txt"));
  gre = new TGraphErrors(resultQuadratic.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetFillColor(kRed-3);
  gre->SetMarkerColor(kRed-3);
  gre->SetLineWidth(2);
  gre->SetMarkerStyle(20);
  for(unsigned int i=0; i<resultQuadratic.size(); i++){
    gre->SetPoint(i,resultQuadratic.at(i).ebeam,resultQuadratic.at(i).erec);
    gre->SetPointError(i,0,resultQuadratic.at(i).erecError);
  }
   
  TH1F *Graph1 = new TH1F("Graph1","Graph",100,0,87.17072);
  Graph1->SetMinimum(2.655724);
  Graph1->SetMaximum(88.56778);
  Graph1->SetDirectory(0);
  Graph1->SetStats(0);
  Graph1->GetXaxis()->SetLabelFont(42);
  Graph1->GetXaxis()->SetTitleSize(0.05);
  Graph1->GetXaxis()->SetTitleFont(42);
  Graph1->GetYaxis()->SetLabelFont(42);
  Graph1->GetYaxis()->SetTitleSize(0.05);
  Graph1->GetYaxis()->SetTitleFont(42);
  Graph1->GetZaxis()->SetLabelFont(42);
  Graph1->GetZaxis()->SetTitleSize(0.05);
  Graph1->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph1);
   
  gre->Draw("p");

  TLegend *leg = new TLegend(0.2,0.7,0.75,0.9,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
   
  TLegendEntry *entry=leg->AddEntry("Graph1","y=E_{beam}","l");
  entry->SetLineColor(kBlack);
  entry->SetLineStyle(7);
  entry->SetLineWidth(2);

  entry=leg->AddEntry("Graph_Graph3","Linear parametrisation","p");
  entry->SetLineColor(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(34);
  entry->SetMarkerSize(1.2);
  
  entry=leg->AddEntry("Graph1","Quadratic parametrisation","p");
  entry->SetLineColor(kRed-3);
  entry->SetLineWidth(kRed-3);
  entry->SetMarkerColor(kRed-3);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.0);

  leg->Draw();

  TText *tex=new TText();
  tex->SetTextSize(0.05);
  tex->SetTextColor(kGray+2);
  tex->DrawTextNDC(0.6,0.05,"Preliminary");
  r1_1->Modified();
  r1->cd();
  
  // ------------>Primitives in pad: r1_2
  TPad *r1_2 = new TPad("r1_2", "Energy_2",0.02,0.0,0.95,0.38);
  r1_2->Draw();
  r1_2->cd();
  r1_2->Range(-19,-0.06545455,95,0.048);
  r1_2->SetFillColor(0);
  r1_2->SetBorderMode(0);
  r1_2->SetBorderSize(2);
  r1_2->SetTickx(1);
  r1_2->SetTicky(1);
  r1_2->SetLeftMargin(0.16);
  r1_2->SetRightMargin(0.01);
  r1_2->SetTopMargin(0.0);
  r1_2->SetBottomMargin(0.23);
  r1_2->SetFrameBorderMode(0);
  r1_2->SetFrameBorderMode(0);
   
  TH1D *hd = new TH1D("hd"," ",85,0,85);
  hd->SetMinimum(-0.12);
  hd->SetMaximum(0.12);
  hd->SetStats(0);
  hd->GetXaxis()->SetTitle("E_{beam} [GeV]");
  hd->GetXaxis()->SetLabelFont(43);
  hd->GetXaxis()->SetLabelSize(20);
  hd->GetXaxis()->SetTitleFont(43);
  hd->GetXaxis()->SetTitleSize(30);
  hd->GetXaxis()->SetTitleOffset(2.);
  hd->GetYaxis()->SetTitle("(#DeltaE)/E_{beam}");
  hd->GetYaxis()->SetLabelFont(43);
  hd->GetYaxis()->SetLabelSize(20);
  hd->GetYaxis()->SetTitleSize(30);
  hd->GetYaxis()->SetTitleOffset(1.7);
  hd->GetYaxis()->SetTitleFont(43);
  hd->GetYaxis()->SetNdivisions(505);
  hd->GetZaxis()->SetLabelFont(42);
  hd->GetZaxis()->SetTitleSize(0.05);
  hd->GetZaxis()->SetTitleFont(42);
  hd->Draw("");
   
  float deltaError;
  float delta;
  gre = new TGraphErrors(resultLinear.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  for(unsigned int i=0; i<resultLinear.size(); i++){
    delta=(resultLinear.at(i).erec-resultLinear.at(i).ebeam)/resultLinear.at(i).ebeam;
    deltaError=sqrt(resultLinear.at(i).erecError/resultLinear.at(i).ebeam/resultLinear.at(i).ebeam);
    gre->SetPoint(i,resultLinear.at(i).ebeam,delta);
    gre->SetPointError(i,0,deltaError);
  }
  gre->SetLineWidth(2);
  gre->SetLineColor(kBlack);
  gre->SetMarkerColor(kBlack);
  gre->SetMarkerSize(1.2);
  gre->SetMarkerStyle(34);
   
  gre->Draw("p");
   
  gre = new TGraphErrors(resultQuadratic.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  for(unsigned int i=0; i<resultQuadratic.size(); i++){
    delta=(resultQuadratic.at(i).erec-resultQuadratic.at(i).ebeam)/resultQuadratic.at(i).ebeam;
    deltaError=sqrt(resultQuadratic.at(i).erecError/resultQuadratic.at(i).ebeam/resultQuadratic.at(i).ebeam);
    gre->SetPoint(i,resultQuadratic.at(i).ebeam,delta);
    gre->SetPointError(i,0,deltaError);
  }
  gre->SetLineWidth(2);
  gre->SetLineColor(kRed-3);
  gre->SetMarkerColor(kRed-3);
  gre->SetMarkerSize(1.0);
  gre->SetMarkerStyle(20);
   
  gre->Draw("p");
   
   
 TF1 *lin1 = new TF1("lin1","0",-0.01,95);
 lin1->SetFillColor(19);
 lin1->SetFillStyle(0);
 lin1->SetLineWidth(1);
 lin1->SetLineStyle(1);
 lin1->SetLineColor(1);
 lin1->GetXaxis()->SetLabelFont(42);
 lin1->GetXaxis()->SetTitleSize(0.05);
 lin1->GetXaxis()->SetTitleFont(42);
 lin1->GetYaxis()->SetLabelFont(42);
 lin1->GetYaxis()->SetTitleSize(0.05);
 lin1->GetYaxis()->SetTitleFont(42);
 lin1->Draw("same");
  
// lin1 = new TF1("lin1","0.03",0.01,95);
// lin1->SetFillColor(1);
// lin1->SetFillStyle(0);
// lin1->SetLineWidth(1);
// lin1->SetLineStyle(2);
// lin1->SetLineColor(17);
// lin1->GetXaxis()->SetLabelFont(42);
// lin1->GetXaxis()->SetTitleSize(0.05);
// lin1->GetXaxis()->SetTitleFont(42);
// lin1->GetYaxis()->SetLabelFont(42);
// lin1->GetYaxis()->SetTitleSize(0.05);
// lin1->GetYaxis()->SetTitleFont(42);
// lin1->Draw("same");
//  
// lin1 = new TF1("lin1","-0.03",0.01,95);
// lin1->SetFillColor(1);
// lin1->SetFillStyle(0);
// lin1->SetLineWidth(1);
// lin1->SetLineStyle(2);
// lin1->SetLineColor(17);
// lin1->GetXaxis()->SetLabelFont(42);
// lin1->GetXaxis()->SetTitleSize(0.05);
// lin1->GetXaxis()->SetTitleFont(42);
// lin1->GetYaxis()->SetLabelFont(42);
// lin1->GetYaxis()->SetTitleSize(0.05);
// lin1->GetYaxis()->SetTitleFont(42);
// lin1->Draw("same");

  r1_2->Modified();
  r1->cd();
  r1->Modified();
  r1->cd();
  r1->SetSelected(r1);
  r1->SaveAs("./plots/EREC.pdf");
}
