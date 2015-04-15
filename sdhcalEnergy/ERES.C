#include "ERES.h"

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

void DrawERES()
{
  //=========Macro generated from canvas: r1/Energy
  //=========  (Fri Apr  6 17:14:48 2012) by ROOT version5.28/00g
  TCanvas *r2 = new TCanvas("r2", "Resolution",213,34,800,700);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  r2->Range(-17.0625,0.008888889,96.6875,0.2311111);
  r2->SetFillColor(0);
  r2->SetBorderMode(0);
  r2->SetBorderSize(2);
  r2->SetTickx(1);
  r2->SetTicky(1);
  r2->SetLeftMargin(0.15);
  r2->SetRightMargin(0.05);
  r2->SetTopMargin(0.05);
  r2->SetBottomMargin(0.14);
  r2->SetFrameBorderMode(0);
  r2->SetFrameBorderMode(0);

  TH1D *h2 = new TH1D("h2"," ",100,0,100);
  h2->SetMinimum(0.04);
  h2->SetMaximum(0.40);
  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("E_{beam} [GeV]");
  h2->GetXaxis()->SetRange(1,91);
  h2->GetXaxis()->SetLabelFont(43);
  h2->GetXaxis()->SetLabelSize(30);
  h2->GetXaxis()->SetTitleSize(34);
  h2->GetXaxis()->SetTitleFont(43);
  h2->GetYaxis()->SetTitle("#sigma_{reco}/E_{reco}");
  h2->GetYaxis()->SetLabelFont(43);
  h2->GetYaxis()->SetLabelSize(30);
  h2->GetYaxis()->SetTitleSize(34);
  h2->GetYaxis()->SetTitleOffset(1.4);
  h2->GetYaxis()->SetTitleFont(43);
  h2->GetZaxis()->SetLabelFont(43);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetTitleFont(42);
  h2->Draw("");

  std::vector<result> resultLinear=readTXT(std::string("LinearErec.txt"));
  TGraphErrors *gre = new TGraphErrors(resultLinear.size());
  gre->SetName("Graph");
  gre->SetTitle("Graph");
  gre->SetLineWidth(2);
  gre->SetMarkerColor(1);
  gre->SetMarkerStyle(34);
  gre->SetMarkerSize(1.2);
  for(unsigned int i=0; i<resultLinear.size(); i++){
    std::cout << resultLinear.at(i).ebeam << " " << resultLinear.at(i).erec << " " << resultLinear.at(i).erecError << std::endl;
    gre->SetPoint(i,resultLinear.at(i).ebeam,resultLinear.at(i).sigma/resultLinear.at(i).erec);
    gre->SetPointError(i,0,sqrt( pow(resultLinear.at(i).sigmaError/resultLinear.at(i).erec,2) +
				 pow(resultLinear.at(i).sigma*resultLinear.at(i).erecError/(resultLinear.at(i).erec*resultLinear.at(i).erec),2) ));
  }

  TF1 *linFunc=new TF1("qaudFunc",fitfunc,0,90,2);
  linFunc->SetParameters(0.6,0.1);
  //  linFunc->SetParLimits(2,0.0,0.02);
  linFunc->SetLineWidth(2);
  linFunc->SetLineStyle(7);
  linFunc->SetLineColor(kBlack);
  gre->Fit(linFunc);

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
  gre->SetMarkerStyle(20);
  for(unsigned int i=0; i<resultQuadratic.size(); i++){
    std::cout << resultQuadratic.at(i).ebeam << " " << resultQuadratic.at(i).erec << " " << resultQuadratic.at(i).erecError << std::endl;
    gre->SetPoint(i,resultQuadratic.at(i).ebeam,resultQuadratic.at(i).sigma/resultQuadratic.at(i).erec);
    gre->SetPointError(i,0,sqrt( pow(resultQuadratic.at(i).sigmaError/resultQuadratic.at(i).erec,2) +
				 pow(resultQuadratic.at(i).sigma*resultQuadratic.at(i).erecError/(resultQuadratic.at(i).erec*resultQuadratic.at(i).erec),2) ));
  }
  TF1 *quadFunc=new TF1("qaudFunc",fitfunc,0,90,2);
  quadFunc->SetParameters(0.6,0.1);
  //quadFunc->SetParLimits(2,0.0,0.02);
  quadFunc->SetLineWidth(2);
  quadFunc->SetLineStyle(7);
  quadFunc->SetLineColor(kRed-3);
  gre->Fit(quadFunc);

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

  TLegend *leg = new TLegend(0.47,0.707381,0.9032663,0.9122024,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
   
  TLegendEntry *entry=leg->AddEntry("Graph_Graph3","Linear parametrisation","lp");
  entry->SetLineColor(1);
  entry->SetLineWidth(2);
  entry->SetLineStyle(7);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(34);
  entry->SetMarkerSize(1.2);
  
  entry=leg->AddEntry("Graph1","Quadratic parametrisation","lp");
  entry->SetLineColor(kRed-3);
  entry->SetLineWidth(2);
  entry->SetLineStyle(7);
  entry->SetMarkerColor(kRed-3);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.0);

  leg->Draw();

  TText *tex=new TText();
  tex->SetTextSize(0.05);
  tex->SetTextColor(kGray+2);
  tex->DrawTextNDC(0.2,0.2,"Preliminary");
  r2->Modified();
  r2->cd();
  r2->SetSelected(r2);
  r2->SaveAs("./plots/ERES.pdf");
}
