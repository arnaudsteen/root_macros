#include "LONGIPION.h"

std::vector<result> readTXT(std::string fileName)
{
  std::vector<result> vec;
  char input[200];
  sprintf(input,"%s%s","/home/steen/root_macros/txtFile/",fileName.c_str());
  std::cout << input << std::endl;
  fstream in;
  in.open(input);
  int ebeam;
  float longiProf;
  float longiProfRMS;
  float longiProfError;
  float longiProfRMSError;
  while(1){
    in >> ebeam >> longiProf >> longiProfError >> longiProfRMS >> longiProfRMSError;
    std::cout << ebeam << " " << longiProf << " " << longiProfError << " " <<  longiProfRMS << " " <<  longiProfRMSError << std::endl;
    result res;
    res.ebeam=ebeam;
    res.longiProf=longiProf;
    res.longiProfRMS=longiProfRMS;
    res.longiProfError=longiProfError;
    res.longiProfRMSError=longiProfRMSError;
    vec.push_back(res);
    if(ebeam==80)break;
  }
  return vec;
}

void DrawLONGI()
{
  //=========Macro generated from canvas: r1/Energy
  //=========  (Fri Apr  6 17:14:48 2012) by ROOT version5.28/00g
  TCanvas *r1 = new TCanvas("r1", "LongiProf",12,24,550,741);
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
  he->SetMaximum(15);
  he->SetStats(0);
  he->GetXaxis()->SetTitle("E_{beam} [GeV]");
  he->GetXaxis()->SetLabelFont(43);
  he->GetXaxis()->SetLabelSize(0);
  he->GetXaxis()->SetTitleFont(43);
  he->GetXaxis()->SetTitleSize(0); 
  he->GetYaxis()->SetTitle("<Z> [layer]");
  he->GetYaxis()->SetLabelFont(43);
  he->GetYaxis()->SetTitleSize(30);
  he->GetYaxis()->SetLabelSize(20);
  he->GetYaxis()->SetTitleFont(43);
  he->GetYaxis()->SetTitleOffset(1.7);
  he->GetZaxis()->SetLabelFont(42);
  he->GetZaxis()->SetTitleSize(0.05);
  he->GetZaxis()->SetTitleFont(42);
  he->Draw("");
   
  std::vector<result> resultData=readTXT(std::string("longiProf_pi-_AugSep2012.txt"));
  TGraphErrors *gre = new TGraphErrors(resultData.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetLineColor(1);
  gre->SetFillStyle(1);
  gre->SetFillColor(1);
  gre->SetLineWidth(2);
  gre->SetMarkerColor(1);
  gre->SetMarkerStyle(34);
  gre->SetMarkerSize(1.2);
  for(unsigned int i=0; i<resultData.size(); i++){
    gre->SetPoint(i,resultData.at(i).ebeam,resultData.at(i).longiProf);
    gre->SetPointError(i,0,resultData.at(i).longiProfError);
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
   
  std::vector<result> resultFTFP=readTXT(std::string("longiProf_pi-_ftfp_bert_hp.txt"));
  gre = new TGraphErrors(resultFTFP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetFillColor(kRed-3);
  gre->SetMarkerColor(kRed-3);
  gre->SetLineWidth(2);
  gre->SetMarkerStyle(20);
  for(unsigned int i=0; i<resultFTFP.size(); i++){
    gre->SetPoint(i,resultFTFP.at(i).ebeam,resultFTFP.at(i).longiProf);
    gre->SetPointError(i,0,resultFTFP.at(i).longiProfError);
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

  std::vector<result> resultQGSP=readTXT(std::string("longiProf_pi-_qgsp_bert_hp.txt"));
  gre = new TGraphErrors(resultQGSP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetFillColor(kBlue-6);
  gre->SetMarkerColor(kBlue-6);
  gre->SetLineWidth(2);
  gre->SetMarkerStyle(25);
  for(unsigned int i=0; i<resultQGSP.size(); i++){
    gre->SetPoint(i,resultQGSP.at(i).ebeam,resultQGSP.at(i).longiProf);
    gre->SetPointError(i,0,resultQGSP.at(i).longiProfError);
  }

  TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,87.5);
  Graph_Graph2->SetMinimum(0);
  Graph_Graph2->SetMaximum(1193.483);
  Graph_Graph2->SetDirectory(0);
  Graph_Graph2->SetStats(0);
  Graph_Graph2->GetXaxis()->SetLabelFont(42);
  Graph_Graph2->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph2->GetXaxis()->SetTitleFont(42);
  Graph_Graph2->GetYaxis()->SetLabelFont(42);
  Graph_Graph2->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph2->GetYaxis()->SetTitleFont(42);
  Graph_Graph2->GetZaxis()->SetLabelFont(42);
  Graph_Graph2->GetZaxis()->SetTitleSize(0.05);
  Graph_Graph2->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph_Graph2);
   
  gre->Draw("p");

  TLegend *leg = new TLegend(0.3,0.4,0.9,0.6,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
   
  TLegendEntry *entry=leg->AddEntry("Graph_Graph3","SDHCAL DATA (H6 CERN SPS)","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(34);
  entry->SetMarkerSize(1.2);

  entry=leg->AddEntry("Graph1","FTFP_BERT_HP","p");
  entry->SetLineColor(kRed-3);
  entry->SetLineStyle(kRed-3);
  entry->SetLineWidth(kRed-3);
  entry->SetMarkerColor(kRed-3);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.0);

  entry=leg->AddEntry("Graph_Graph2","QGSP_BERT_HP","p");
  entry->SetLineColor(kBlue-6);
  entry->SetLineStyle(kBlue-6);
  entry->SetLineWidth(kBlue-6);
  entry->SetMarkerColor(kBlue-6);
  entry->SetMarkerStyle(25);
  entry->SetMarkerSize(0.9);

  leg->Draw();

  TText *tex=new TText();
  tex->SetTextSize(0.05);
  tex->SetTextColor(kGray+2);
  //tex->DrawTextNDC(0.5,0.05,"SDHCAL Preliminary");
  tex->DrawTextNDC(0.3,0.05,"CALICE Fe-SDHCAL Preliminary");
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
  hd->SetMinimum(-0.08);
  hd->SetMaximum(0.08);
  hd->SetStats(0);
  hd->GetXaxis()->SetTitle("E_{beam} [GeV]");
  hd->GetXaxis()->SetLabelFont(43);
  hd->GetXaxis()->SetLabelSize(20);
  hd->GetXaxis()->SetTitleFont(43);
  hd->GetXaxis()->SetTitleSize(30);
  hd->GetXaxis()->SetTitleOffset(2.);
  hd->GetYaxis()->SetTitle("(#DeltaZ)/Z");
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
  gre = new TGraphErrors(resultQGSP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  for(unsigned int i=0; i<resultQGSP.size(); i++){
    delta=(resultQGSP.at(i).longiProf-resultData.at(i).longiProf)/resultData.at(i).longiProf;
    deltaError=1/resultData.at(i).longiProf*
      sqrt(pow(resultQGSP.at(i).longiProfError,2) +
    	   pow(resultQGSP.at(i).longiProf/resultData.at(i).longiProf*resultData.at(i).longiProfError,2));
    gre->SetPoint(i,resultQGSP.at(i).ebeam,delta);
    gre->SetPointError(i,0,deltaError);
  }
  gre->SetLineWidth(2);
  gre->SetLineColor(kBlue-6);
  gre->SetMarkerColor(kBlue-6);
  gre->SetMarkerSize(1.0);
  gre->SetMarkerStyle(25);
   
  gre->Draw("p");
   
  gre = new TGraphErrors(resultFTFP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  for(unsigned int i=0; i<resultFTFP.size(); i++){
    delta=(resultFTFP.at(i).longiProf-resultData.at(i).longiProf)/resultData.at(i).longiProf;
    deltaError=1/resultData.at(i).longiProf*
      sqrt(pow(resultFTFP.at(i).longiProfError,2) +
    	   pow(resultFTFP.at(i).longiProf/resultData.at(i).longiProf*resultData.at(i).longiProfError,2));
    gre->SetPoint(i,resultFTFP.at(i).ebeam,delta);
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
   
  lin1 = new TF1("lin1","0.1",0.01,95);
  lin1->SetFillColor(1);
  lin1->SetFillStyle(0);
  lin1->SetLineWidth(1);
  lin1->SetLineStyle(2);
  lin1->SetLineColor(17);
  lin1->GetXaxis()->SetLabelFont(42);
  lin1->GetXaxis()->SetTitleSize(0.05);
  lin1->GetXaxis()->SetTitleFont(42);
  lin1->GetYaxis()->SetLabelFont(42);
  lin1->GetYaxis()->SetTitleSize(0.05);
  lin1->GetYaxis()->SetTitleFont(42);
  lin1->Draw("same");
   
  lin1 = new TF1("lin1","-0.1",0.01,95);
  lin1->SetFillColor(1);
  lin1->SetFillStyle(0);
  lin1->SetLineWidth(1);
  lin1->SetLineStyle(2);
  lin1->SetLineColor(17);
  lin1->GetXaxis()->SetLabelFont(42);
  lin1->GetXaxis()->SetTitleSize(0.05);
  lin1->GetXaxis()->SetTitleFont(42);
  lin1->GetYaxis()->SetLabelFont(42);
  lin1->GetYaxis()->SetTitleSize(0.05);
  lin1->GetYaxis()->SetTitleFont(42);
  lin1->Draw("same");

  r1_2->Modified();
  r1->cd();
  r1->Modified();
  r1->cd();
  r1->SetSelected(r1);
  r1->SaveAs("../plots/LONGIPION.pdf");
}


void DrawRMS()
{
  //=========Macro generated from canvas: r1/Energy
  //=========  (Fri Apr  6 17:14:48 2012) by ROOT version5.28/00g
  TCanvas *r1 = new TCanvas("rms", "LongiProfRMS",12,24,550,741);
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
  TPad *r1_1 = new TPad("r1_rms", "LongiProfRMS",0.02,0.37,0.95,0.99);
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
  he->SetMaximum(14.0);
  he->SetStats(0);
  he->GetXaxis()->SetTitle("E_{beam} [GeV]");
  he->GetXaxis()->SetLabelFont(43);
  he->GetXaxis()->SetLabelSize(0);
  he->GetXaxis()->SetTitleFont(43);
  he->GetXaxis()->SetTitleSize(0); 
  he->GetYaxis()->SetTitle("#sqrt{<Z^{2}>-<Z>^{2}} [layer]");
  he->GetYaxis()->SetLabelFont(43);
  he->GetYaxis()->SetTitleSize(30);
  he->GetYaxis()->SetLabelSize(20);
  he->GetYaxis()->SetTitleFont(43);
  he->GetYaxis()->SetTitleOffset(1.7);
  he->GetZaxis()->SetLabelFont(42);
  he->GetZaxis()->SetTitleSize(0.05);
  he->GetZaxis()->SetTitleFont(42);
  he->Draw("");
   
  std::vector<result> resultData=readTXT(std::string("longiProf_pi-_AugSep2012.txt"));
  TGraphErrors *gre = new TGraphErrors(resultData.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetLineColor(1);
  gre->SetFillStyle(1);
  gre->SetFillColor(1);
  gre->SetLineWidth(2);
  gre->SetMarkerColor(1);
  gre->SetMarkerStyle(34);
  gre->SetMarkerSize(1.2);
  for(unsigned int i=0; i<resultData.size(); i++){
    gre->SetPoint(i,resultData.at(i).ebeam,resultData.at(i).longiProfRMS);
    gre->SetPointError(i,0,resultData.at(i).longiProfRMSError);
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
   
  std::vector<result> resultFTFP=readTXT(std::string("longiProf_pi-_ftfp_bert_hp.txt"));
  gre = new TGraphErrors(resultFTFP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetFillColor(kRed-3);
  gre->SetMarkerColor(kRed-3);
  gre->SetLineWidth(2);
  gre->SetMarkerStyle(20);
  for(unsigned int i=0; i<resultFTFP.size(); i++){
    gre->SetPoint(i,resultFTFP.at(i).ebeam,resultFTFP.at(i).longiProfRMS);
    gre->SetPointError(i,0,resultFTFP.at(i).longiProfRMSError);
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

  std::vector<result> resultQGSP=readTXT(std::string("longiProf_pi-_qgsp_bert_hp.txt"));
  gre = new TGraphErrors(resultQGSP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetFillColor(kBlue-6);
  gre->SetMarkerColor(kBlue-6);
  gre->SetLineWidth(2);
  gre->SetMarkerStyle(25);
  for(unsigned int i=0; i<resultQGSP.size(); i++){
    gre->SetPoint(i,resultQGSP.at(i).ebeam,resultQGSP.at(i).longiProfRMS);
    gre->SetPointError(i,0,resultQGSP.at(i).longiProfRMSError);
  }

  TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,87.5);
  Graph_Graph2->SetMinimum(0);
  Graph_Graph2->SetMaximum(1193.483);
  Graph_Graph2->SetDirectory(0);
  Graph_Graph2->SetStats(0);
  Graph_Graph2->GetXaxis()->SetLabelFont(42);
  Graph_Graph2->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph2->GetXaxis()->SetTitleFont(42);
  Graph_Graph2->GetYaxis()->SetLabelFont(42);
  Graph_Graph2->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph2->GetYaxis()->SetTitleFont(42);
  Graph_Graph2->GetZaxis()->SetLabelFont(42);
  Graph_Graph2->GetZaxis()->SetTitleSize(0.05);
  Graph_Graph2->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph_Graph2);
   
  gre->Draw("p");

  TLegend *leg = new TLegend(0.3,0.4,0.9,0.6,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
   
  TLegendEntry *entry=leg->AddEntry("Graph_Graph3","SDHCAL DATA (H6 CERN SPS)","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(34);
  entry->SetMarkerSize(1.2);

  entry=leg->AddEntry("Graph1","FTFP_BERT_HP","p");
  entry->SetLineColor(kRed-3);
  entry->SetLineStyle(kRed-3);
  entry->SetLineWidth(kRed-3);
  entry->SetMarkerColor(kRed-3);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.0);

  entry=leg->AddEntry("Graph_Graph2","QGSP_BERT_HP","p");
  entry->SetLineColor(kBlue-6);
  entry->SetLineStyle(kBlue-6);
  entry->SetLineWidth(kBlue-6);
  entry->SetMarkerColor(kBlue-6);
  entry->SetMarkerStyle(25);
  entry->SetMarkerSize(0.9);

  leg->Draw();

  TText *tex=new TText();
  tex->SetTextSize(0.05);
  tex->SetTextColor(kGray+2);
  //tex->DrawTextNDC(0.5,0.05,"SDHCAL Preliminary");
  tex->DrawTextNDC(0.3,0.05,"CALICE Fe-SDHCAL Preliminary");
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
  hd->SetMinimum(-0.08);
  hd->SetMaximum(0.08);
  hd->SetStats(0);
  hd->GetXaxis()->SetTitle("E_{beam} [GeV]");
  hd->GetXaxis()->SetLabelFont(43);
  hd->GetXaxis()->SetLabelSize(20);
  hd->GetXaxis()->SetTitleFont(43);
  hd->GetXaxis()->SetTitleSize(30);
  hd->GetXaxis()->SetTitleOffset(2.);
  hd->GetYaxis()->SetTitle("(#Delta#sigma_{z})/#sigma_{z}");
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
  gre = new TGraphErrors(resultQGSP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  for(unsigned int i=0; i<resultQGSP.size(); i++){
    delta=(resultQGSP.at(i).longiProfRMS-resultData.at(i).longiProfRMS)/resultData.at(i).longiProfRMS;
    deltaError=1/resultData.at(i).longiProfRMS*
      sqrt(pow(resultQGSP.at(i).longiProfRMSError,2) +
    	   pow(resultQGSP.at(i).longiProfRMS/resultData.at(i).longiProfRMS*resultData.at(i).longiProfRMSError,2));
    gre->SetPoint(i,resultQGSP.at(i).ebeam,delta);
    gre->SetPointError(i,0,deltaError);
  }
  gre->SetLineWidth(2);
  gre->SetLineColor(kBlue-6);
  gre->SetMarkerColor(kBlue-6);
  gre->SetMarkerSize(1.0);
  gre->SetMarkerStyle(25);
   
  gre->Draw("p");
   
  gre = new TGraphErrors(resultFTFP.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  for(unsigned int i=0; i<resultFTFP.size(); i++){
    delta=(resultFTFP.at(i).longiProfRMS-resultData.at(i).longiProfRMS)/resultData.at(i).longiProfRMS;
    deltaError=1/resultData.at(i).longiProfRMS*
      sqrt(pow(resultFTFP.at(i).longiProfRMSError,2) +
    	   pow(resultFTFP.at(i).longiProfRMS/resultData.at(i).longiProfRMS*resultData.at(i).longiProfRMSError,2));
    gre->SetPoint(i,resultFTFP.at(i).ebeam,delta);
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
   
  lin1 = new TF1("lin1","0.1",0.01,95);
  lin1->SetFillColor(1);
  lin1->SetFillStyle(0);
  lin1->SetLineWidth(1);
  lin1->SetLineStyle(2);
  lin1->SetLineColor(17);
  lin1->GetXaxis()->SetLabelFont(42);
  lin1->GetXaxis()->SetTitleSize(0.05);
  lin1->GetXaxis()->SetTitleFont(42);
  lin1->GetYaxis()->SetLabelFont(42);
  lin1->GetYaxis()->SetTitleSize(0.05);
  lin1->GetYaxis()->SetTitleFont(42);
  lin1->Draw("same");
   
  lin1 = new TF1("lin1","-0.1",0.01,95);
  lin1->SetFillColor(1);
  lin1->SetFillStyle(0);
  lin1->SetLineWidth(1);
  lin1->SetLineStyle(2);
  lin1->SetLineColor(17);
  lin1->GetXaxis()->SetLabelFont(42);
  lin1->GetXaxis()->SetTitleSize(0.05);
  lin1->GetXaxis()->SetTitleFont(42);
  lin1->GetYaxis()->SetLabelFont(42);
  lin1->GetYaxis()->SetTitleSize(0.05);
  lin1->GetYaxis()->SetTitleFont(42);
  lin1->Draw("same");

  r1_2->Modified();
  r1->cd();
  r1->Modified();
  r1->cd();
  r1->SetSelected(r1);
  r1->SaveAs("../plots/LONGIRMSPION.pdf");
}
