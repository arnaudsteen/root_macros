#define AnaMuons_cxx
#include "AnaMuons.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TText.h>
#include <iostream>

void AnaMuons::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  hTheta=new TH1D("hTheta","",100,0,1);
  hMul=new TH2D("hMul","",100,0,1,10000,0,100);
  hEff=new TH2D("hEff","",100,0,1,100,0,1.1);
  hMulLayer=new TH2D("hMulLayer","",50,0,50,10000,0,100);
  hEffLayer=new TH2D("hEffLayer","",50,0,50,1000,0,2);
  hEff2Layer=new TH2D("hEff2Layer","",50,0,50,1000,0,2);
  hEff3Layer=new TH2D("hEff3Layer","",50,0,50,1000,0,2);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    TVector3 px(-1,0,trackParams[1]);
    TVector3 py(0,-1,trackParams[3]);
    TVector3 p=px.Cross(py);
    if( transversRatio<.05 && chi2Global<100 ){
      hTheta->Fill(p.CosTheta());
      if(p.CosTheta()>0.15){
	hMul->Fill(p.CosTheta(),mulGlobal);
	hEff->Fill(p.CosTheta(),effGlobal);
      }
      if(p.CosTheta()>0.9){//&&trackParams[0]>0&&trackParams[0]<1000&&trackParams[2]>0&&trackParams[2]<1000){
	for(int i=0; i<48; i++){
	  //if(i==44||i==16||i==28)continue;
	  if(Efficiency1[i]>0&&Chi2[i]<10){
	    hMulLayer->Fill(i,Multiplicity[i]);
	    if(mulcount[i]) {
	      mulmap[i]+=Multiplicity[i];
	      mulsquaremap[i]+=Multiplicity[i]*Multiplicity[i];
	      mulcount[i]++;
	    }
	    else{
	      mulmap[i]=Multiplicity[i];
	      mulsquaremap[i]=Multiplicity[i]*Multiplicity[i];
	      mulcount[i]=1;
	    }
	  }
	  if(Efficiency1[i]>-1){
	    hEffLayer->Fill(i,Efficiency1[i]);
	    hEff2Layer->Fill(i,Efficiency2[i]);
	    hEff3Layer->Fill(i,Efficiency3[i]);
	    if(effcount[i]){
	      effmap[i]+=Efficiency1[i];
	      effsquaremap[i]+=Efficiency1[i]*Efficiency1[i];
	      eff2map[i]+=Efficiency2[i];
	      eff2squaremap[i]+=Efficiency2[i]*Efficiency2[i];
	      eff3map[i]+=Efficiency3[i];
	      eff3squaremap[i]+=Efficiency3[i]*Efficiency3[i];
	      effcount[i]++;
	    }
	    else{
	      effmap[i]=Efficiency1[i];
	      effsquaremap[i]=Efficiency1[i]*Efficiency1[i];
	      eff2map[i]=Efficiency2[i];
	      eff2squaremap[i]=Efficiency2[i]*Efficiency2[i];
	      eff3map[i]=Efficiency3[i];
	      eff3squaremap[i]=Efficiency3[i]*Efficiency3[i];
	      effcount[i]=1;
	    }
	  }
	}
      }
    }
  }
}

void ProcessMuons(int run=715747,int energy=30)
{
  char input[200];
  //sprintf(input,"%s%d%s","/home/steen/resultRootFile/tb_data/TDHCAL_",run,".root");
  //sprintf(input,"%s%d%s","/home/steen/sdhcal_analysis/TDHCAL_",run,"_test.root");
  sprintf(input,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/tb_data/TDHCAL_",run,".root");
  AnaMuons *data=new AnaMuons(input);
  data->Loop();
  
  //sprintf(input,"%s%d%s","/lyoserv/home/ilc/steen/resultRootFile/sim_data/single_mu-_",energy,"GeV.root");
  //sprintf(input,"%s%d%s","/home/steen/MarlinReco/trunk/single_mu-_",energy,"GeV_test.root");
  sprintf(input,"%s%d%s","/home/steen/sdhcal_analysis/single_mu-_",energy,"GeV.root");
  //sprintf(input,"%s%d%s","/home/steen/sdhcal_analysis/single_mu-_",energy,"GeV_noCorrection.root");
  AnaMuons *simu=new AnaMuons(input);
  simu->Loop();

  CaliceStyle();
  TLatex *tex=new TLatex();
  tex->SetTextSize(0.04);
  char myText[200];
  TCanvas *cTheta=new TCanvas();
  cTheta->SetWindowSize(600,600);
  cTheta->SetLogy();
  data->hTheta->SetLineColor(kBlack);
  data->hTheta->SetLineWidth(2);
  data->hTheta->GetXaxis()->SetTitle("cos(#theta)");
  data->hTheta->Sumw2();
  data->hTheta->Scale(1/data->hTheta->Integral());
  data->hTheta->Draw("hist");
  sprintf(myText,"%s","Data (SPS H6)");
  tex->DrawTextNDC(0.2,0.8,myText);
  simu->hTheta->SetLineColor(kRed-3);
  simu->hTheta->SetLineWidth(2);
  simu->hTheta->GetXaxis()->SetTitle("cos(#theta)");
  simu->hTheta->Sumw2();
  simu->hTheta->Scale(1/simu->hTheta->Integral());
  simu->hTheta->Draw("histsame");
  sprintf(myText,"%s","Simulation (FTFP_BERT_HP)");
  tex->SetTextColor(kRed-3);
  tex->DrawTextNDC(0.2,0.7,myText);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  cTheta->Update();
  cTheta->SaveAs("theta.pdf");

  TCanvas *cMul=new TCanvas();
  cMul->SetWindowSize(600,600);
  TProfile *dpMul=data->hMul->ProfileX();
  dpMul->SetLineColor(kBlack);
  dpMul->SetLineWidth(2);
  dpMul->SetMarkerStyle(20);
  dpMul->SetMarkerSize(0.8);
  dpMul->SetMarkerColor(kBlack);
  dpMul->GetXaxis()->SetTitle("cos(#theta)");
  dpMul->GetXaxis()->SetTitleSize(0.05);
  dpMul->GetYaxis()->SetTitle("Multiplicity");
  dpMul->GetYaxis()->SetTitleOffset(1.2);
  dpMul->Draw();
  sprintf(myText,"%s","Data (SPS H6)");
  tex->SetTextColor(kBlack);
  tex->DrawTextNDC(0.2,0.4,myText);
  TProfile *spMul=simu->hMul->ProfileX("",1,-1,"");
  spMul->SetLineColor(kRed-3);
  spMul->SetLineWidth(2);
  spMul->SetMarkerStyle(21);
  spMul->SetMarkerSize(0.8);
  spMul->SetMarkerColor(kRed-3);
  spMul->GetXaxis()->SetTitle("cos(#theta)");
  spMul->GetYaxis()->SetTitle("Multiplicity");
  spMul->GetXaxis()->SetTitleSize(0.05);
  spMul->GetYaxis()->SetTitleOffset(1.2);
  spMul->Draw("samep");
  sprintf(myText,"%s","Simulation (FTFP_BERT_HP)");
  tex->SetTextColor(kRed-3);
  tex->DrawTextNDC(0.2,0.3,myText);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  cMul->Update();
  cMul->SaveAs("mul_vs_theta.pdf");

  TCanvas *cEff=new TCanvas();
  cEff->SetWindowSize(600,600);
  TProfile *dpEff=data->hEff->ProfileX();
  dpEff->SetLineColor(kBlack);
  dpEff->SetLineWidth(2);
  dpEff->SetMarkerStyle(20);
  dpEff->SetMarkerSize(0.8);
  dpEff->SetMarkerColor(kBlack);
  dpEff->GetXaxis()->SetTitle("cos(#theta)");
  dpEff->GetYaxis()->SetTitle("Efficiency");
  dpEff->GetYaxis()->SetTitleOffset(1.2);
  dpEff->Draw();
  sprintf(myText,"%s","Data (SPS H6)");
  tex->SetTextColor(kBlack);
  tex->DrawTextNDC(0.2,0.6,myText);
  TProfile *spEff=simu->hEff->ProfileX("",1,-1,"");
  spEff->SetLineColor(kRed-3);
  spEff->SetLineWidth(2);
  spEff->SetMarkerStyle(21);
  spEff->SetMarkerSize(0.8);
  spEff->SetMarkerColor(kRed-3);
  spEff->GetXaxis()->SetTitle("cos(#theta)");
  spEff->GetYaxis()->SetTitle("Multiplicity");
  spEff->GetYaxis()->SetTitleOffset(1.2);
  spEff->Draw("samep");
  sprintf(myText,"%s","Simulation (FTFP_BERT_HP)");
  tex->SetTextColor(kRed-3);
  tex->DrawTextNDC(0.2,0.5,myText);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  cEff->Update();
  cEff->SaveAs("eff_vs_theta.pdf");

  float data_meanMul=0;
  float data_rmsMul=0;
  float data_meanEff=0;
  float data_rmsEff=0;
  float data_meanEff2=0;
  float data_rmsEff2=0;
  float data_meanEff3=0;
  float data_rmsEff3=0;
  float simu_meanMul=0;
  float simu_rmsMul=0;
  float simu_meanEff=0;
  float simu_rmsEff=0;
  float simu_meanEff2=0;
  float simu_rmsEff2=0;
  float simu_meanEff3=0;
  float simu_rmsEff3=0;
  for(int i=0; i<48;i++){
    //if(i==44||i==16||i==28)continue;
    data_meanMul+=data->mulmap[i]/data->mulcount[i];
    data_meanEff+=data->effmap[i]/data->effcount[i];
    data_meanEff2+=data->eff2map[i]/data->effcount[i];
    data_meanEff3+=data->eff3map[i]/data->effcount[i];
    data->mulsquaremap[i]=data->mulmap[i]/data->mulcount[i]*data->mulmap[i]/data->mulcount[i];
    data->effsquaremap[i]=data->effmap[i]/data->effcount[i]*data->effmap[i]/data->effcount[i];
    data->eff2squaremap[i]=data->eff2map[i]/data->effcount[i]*data->eff2map[i]/data->effcount[i];
    data->eff3squaremap[i]=data->eff3map[i]/data->effcount[i]*data->eff3map[i]/data->effcount[i];
    data_rmsMul+=data->mulsquaremap[i];
    data_rmsEff+=data->effsquaremap[i];
    data_rmsEff2+=data->eff2squaremap[i];
    data_rmsEff3+=data->eff3squaremap[i];
 
    simu_meanMul+=simu->mulmap[i]/simu->mulcount[i];
    simu_meanEff+=simu->effmap[i]/simu->effcount[i];
    simu_meanEff2+=simu->eff2map[i]/simu->effcount[i];
    simu_meanEff3+=simu->eff3map[i]/simu->effcount[i];
    simu->mulsquaremap[i]=simu->mulmap[i]/simu->mulcount[i]*simu->mulmap[i]/simu->mulcount[i];
    simu->effsquaremap[i]=simu->effmap[i]/simu->effcount[i]*simu->effmap[i]/simu->effcount[i];
    simu->eff2squaremap[i]=simu->eff2map[i]/simu->effcount[i]*simu->eff2map[i]/simu->effcount[i];
    simu->eff3squaremap[i]=simu->eff3map[i]/simu->effcount[i]*simu->eff3map[i]/simu->effcount[i];
    simu_rmsMul+=simu->mulsquaremap[i];
    simu_rmsEff+=simu->effsquaremap[i];
    simu_rmsEff2+=simu->eff2squaremap[i];
    simu_rmsEff3+=simu->eff3squaremap[i];
  }
  
  data_meanMul=data_meanMul/48;
  data_meanEff=data_meanEff/48;
  data_meanEff2=data_meanEff2/48;
  data_meanEff3=data_meanEff3/48;
  data_rmsMul=sqrt(data_rmsMul/48-data_meanMul*data_meanMul);
  data_rmsEff=sqrt(data_rmsEff/48-data_meanEff*data_meanEff);
  data_rmsEff2=sqrt(data_rmsEff2/48-data_meanEff2*data_meanEff2);
  data_rmsEff3=sqrt(data_rmsEff3/48-data_meanEff3*data_meanEff3);
  std::cout << data_meanMul << " " << data_rmsMul << "\t" 
	    << data_meanEff << " " << data_rmsEff << std::endl;

  simu_meanMul=simu_meanMul/48;
  simu_meanEff=simu_meanEff/48;
  simu_meanEff2=simu_meanEff2/48;
  simu_meanEff3=simu_meanEff3/48;
  simu_rmsMul=sqrt(simu_rmsMul/48-simu_meanMul*simu_meanMul);
  simu_rmsEff=sqrt(simu_rmsEff/48-simu_meanEff*simu_meanEff);
  simu_rmsEff2=sqrt(simu_rmsEff2/48-simu_meanEff2*simu_meanEff2);
  simu_rmsEff3=sqrt(simu_rmsEff3/48-simu_meanEff3*simu_meanEff3);
  std::cout << simu_meanMul << " " << simu_rmsMul << "\t" 
	    << simu_meanEff << " " << simu_rmsEff << std::endl;

  TCanvas *cMulLayer=new TCanvas();
  cMulLayer->SetWindowSize(600,600);
  data->fit_func=new TF1("f","pol0",0,48);
  data->fit_func->SetLineWidth(2);
  data->fit_func->SetLineStyle(7);
  TProfile *dpMulLayer=data->hMulLayer->ProfileX();
  dpMulLayer->SetMarkerStyle(20);
  dpMulLayer->SetMarkerSize(0.8);
  dpMulLayer->SetMarkerColor(kBlack);
  dpMulLayer->SetLineColor(kBlack);
  dpMulLayer->GetXaxis()->SetTitle("Layer");
  dpMulLayer->GetYaxis()->SetTitle("Multiplicity");
  dpMulLayer->GetYaxis()->SetTitleOffset(1.2);
  dpMulLayer->GetYaxis()->SetRangeUser(1.2,2.6);
  dpMulLayer->Draw();
  //  dpMulLayer->Fit(data->fit_func,"n");
  simu->fit_func=new TF1("f","pol0",0,48);
  simu->fit_func->SetLineWidth(2);
  simu->fit_func->SetLineStyle(1);
  simu->fit_func->SetLineColor(kRed-3);
  TProfile *spMulLayer=simu->hMulLayer->ProfileX("",1,-1,"");
  spMulLayer->SetMarkerStyle(21);
  spMulLayer->SetMarkerSize(0.8);
  spMulLayer->SetMarkerColor(kRed-3);
  spMulLayer->SetLineColor(kRed-3);
  spMulLayer->SetLineWidth(1);
  spMulLayer->GetXaxis()->SetTitle("Layer");
  spMulLayer->GetYaxis()->SetTitle("Multiplicity");
  spMulLayer->GetYaxis()->SetTitleOffset(1.2);
  spMulLayer->GetYaxis()->SetRangeUser(1,2.6);
  //  spMulLayer->Fit(simu->fit_func,"n");
  dpMulLayer->Draw();
  data->fit_func->SetParameter(0,data_meanMul);
  spMulLayer->Draw("same");
  simu->fit_func->SetParameter(0,simu_meanMul);
  simu->fit_func->Draw("same");
  data->fit_func->Draw("same");
  sprintf(myText,"%s%.2f%s%.2f%s","#bar{#mu} = ",data_meanMul," #pm ",data_rmsMul," (data)");
  tex->SetTextColor(kBlack);
  tex->DrawLatex(5,1.4,myText);
  sprintf(myText,"%s%.2f%s%.2f%s","#bar{#mu} = ",simu_meanMul," #pm ",simu_rmsMul," (MC)");
  tex->SetTextColor(kRed-3);
  tex->DrawLatex(5,1.3,myText);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  cMulLayer->Update();
  cMulLayer->SaveAs("mulLayer.pdf");

  TCanvas *cEffLayer=new TCanvas();
  cEffLayer->SetWindowSize(600,600);
  data->fit_func=new TF1("f","pol0",0,48);
  data->fit_func->SetLineWidth(2);
  data->fit_func->SetLineStyle(7);
  TProfile *dpEffLayer=data->hEffLayer->ProfileX();
  dpEffLayer->SetMarkerStyle(20);
  dpEffLayer->SetMarkerSize(0.8);
  dpEffLayer->SetMarkerColor(kBlack);
  dpEffLayer->SetLineColor(kBlack);
  dpEffLayer->GetXaxis()->SetTitle("Layer");
  dpEffLayer->GetYaxis()->SetTitle("Efficiency");
  dpEffLayer->GetYaxis()->SetTitleOffset(1.2);
  dpEffLayer->GetYaxis()->SetRangeUser(0.6,1.05);
  dpEffLayer->Fit(data->fit_func);
  //data->fit_func->SetParameter(0,data_meanEff);
  //data->fit_func->Draw("same");
  simu->fit_func=new TF1("f","pol0",0,48);
  simu->fit_func->SetLineWidth(2);
  simu->fit_func->SetLineStyle(1);
  simu->fit_func->SetLineColor(kRed-3);
  TProfile *spEffLayer=simu->hEffLayer->ProfileX("",1,-1,"");
  spEffLayer->SetMarkerStyle(21);
  spEffLayer->SetMarkerSize(0.8);
  spEffLayer->SetMarkerColor(kRed-3);
  spEffLayer->SetLineColor(kRed-3);
  spEffLayer->SetLineWidth(1);
  spEffLayer->GetXaxis()->SetTitle("Layer");
  spEffLayer->GetYaxis()->SetTitle("Efficiency");
  spEffLayer->GetYaxis()->SetTitleOffset(1.2);
  spEffLayer->GetYaxis()->SetRangeUser(0.6,1.05);
  spEffLayer->Draw("same");
  spEffLayer->Fit(simu->fit_func);
  //simu->fit_func->SetParameter(0,simu_meanEff);
  //simu->fit_func->Draw("same");
  dpEffLayer->Draw();
  spEffLayer->Draw("same");
  simu->fit_func->Draw("same");
  data->fit_func->Draw("same");
  sprintf(myText,"%s%.2f%s%.2f%s","#bar{#varepsilon} = ",data->fit_func->GetParameter(0)," #pm ",data_rmsEff," (data)");
  tex->SetTextColor(kBlack);
  tex->DrawLatex(5,0.74666,myText);
  sprintf(myText,"%s%.2f%s%.2f%s","#bar{#varepsilon} = ",simu->fit_func->GetParameter(0)," #pm ",simu_rmsEff," (MC)");
  tex->SetTextColor(kRed-3);
  tex->DrawLatex(5,0.72333,myText);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  cEffLayer->Update();
  cEffLayer->SaveAs("effLayer.pdf"); 


  TCanvas *cEff2Layer=new TCanvas();
  cEff2Layer->SetWindowSize(600,600);
  data->fit_func=new TF1("f","pol0",0,48);
  data->fit_func->SetLineWidth(2);
  data->fit_func->SetLineStyle(7);
  TProfile *dpEff2Layer=data->hEff2Layer->ProfileX();
  dpEff2Layer->SetMarkerStyle(20);
  dpEff2Layer->SetMarkerSize(0.8);
  dpEff2Layer->SetMarkerColor(kBlack);
  dpEff2Layer->SetLineColor(kBlack);
  dpEff2Layer->GetXaxis()->SetTitle("Layer");
  dpEff2Layer->GetYaxis()->SetTitle("Efficiency thr2");
  dpEff2Layer->GetYaxis()->SetTitleOffset(1.2);
  dpEff2Layer->GetYaxis()->SetRangeUser(0.,0.6);
  dpEff2Layer->Fit(data->fit_func);
  //data->fit_func->SetParameter(0,data_meanEff2);
  //data->fit_func->Draw("same");
  simu->fit_func=new TF1("f","pol0",0,48);
  simu->fit_func->SetLineWidth(2);
  simu->fit_func->SetLineStyle(1);
  simu->fit_func->SetLineColor(kRed-3);
  TProfile *spEff2Layer=simu->hEff2Layer->ProfileX("",1,-1,"");
  spEff2Layer->SetMarkerStyle(21);
  spEff2Layer->SetMarkerSize(0.8);
  spEff2Layer->SetMarkerColor(kRed-3);
  spEff2Layer->SetLineColor(kRed-3);
  spEff2Layer->SetLineWidth(1);
  spEff2Layer->GetXaxis()->SetTitle("Layer");
  spEff2Layer->GetYaxis()->SetTitle("Efficiency thr2");
  spEff2Layer->GetYaxis()->SetTitleOffset(1.2);
  spEff2Layer->GetYaxis()->SetRangeUser(0.,0.6);
  spEff2Layer->Draw("same");
  spEff2Layer->Fit(simu->fit_func);
  //simu->fit_func->SetParameter(0,simu_meanEff2);
  //simu->fit_func->Draw("same");
  dpEff2Layer->Draw();
  spEff2Layer->Draw("same");
  simu->fit_func->Draw("same");
  data->fit_func->Draw("same");
  sprintf(myText,"%s%.2f%s%.2f%s","#bar{#varepsilon} = ",data->fit_func->GetParameter(0)," #pm ",data_rmsEff2," (data)");
  tex->SetTextColor(kBlack);
  tex->DrawLatex(5,0.55,myText);
  sprintf(myText,"%s%.2f%s%.2f%s","#bar{#varepsilon} = ",simu->fit_func->GetParameter(0)," #pm ",simu_rmsEff2," (MC)");
  tex->SetTextColor(kRed-3);
  tex->DrawLatex(5,0.5,myText);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  cEff2Layer->Update();
  cEff2Layer->SaveAs("eff2Layer.pdf"); 
  //for(int i=0; i<48; i++){
  //  data->effmap[i]=data->effmap[i]/data->effcount[i];
  //  data->mulmap[i]=data->mulmap[i]/data->mulcount[i];
  //  simu->effmap[i]=simu->effmap[i]/simu->effcount[i];
  //  simu->mulmap[i]=simu->mulmap[i]/simu->mulcount[i];
  //  std::cout << "DATA : data->effmap[" << i << "] = " << data->effmap[i] << "\t data->mulmap[" << i << "] = " << data->mulmap[i] << std::endl;
  //  std::cout << "SIMU : simu->effmap[" << i << "] = " << simu->effmap[i] << "\t simu->mulmap[" << i << "] = " << simu->mulmap[i] << std::endl;
  //}
  TCanvas *cEff3Layer=new TCanvas();
  cEff3Layer->SetWindowSize(600,600);
  data->fit_func=new TF1("f","pol0",0,48);
  data->fit_func->SetLineWidth(2);
  data->fit_func->SetLineStyle(7);
  TProfile *dpEff3Layer=data->hEff3Layer->ProfileX();
  dpEff3Layer->SetMarkerStyle(20);
  dpEff3Layer->SetMarkerSize(0.8);
  dpEff3Layer->SetMarkerColor(kBlack);
  dpEff3Layer->SetLineColor(kBlack);
  dpEff3Layer->GetXaxis()->SetTitle("Layer");
  dpEff3Layer->GetYaxis()->SetTitle("Efficiency thr3");
  dpEff3Layer->GetYaxis()->SetTitleOffset(1.2);
  dpEff3Layer->GetYaxis()->SetRangeUser(0.,0.15);
  dpEff3Layer->Fit(data->fit_func);
  //data->fit_func->SetParameter(0,data_meanEff3);
  //data->fit_func->Draw("same");
  simu->fit_func=new TF1("f","pol0",0,48);
  simu->fit_func->SetLineWidth(2);
  simu->fit_func->SetLineStyle(1);
  simu->fit_func->SetLineColor(kRed-3);
  TProfile *spEff3Layer=simu->hEff3Layer->ProfileX("",1,-1,"");
  spEff3Layer->SetMarkerStyle(21);
  spEff3Layer->SetMarkerSize(0.8);
  spEff3Layer->SetMarkerColor(kRed-3);
  spEff3Layer->SetLineColor(kRed-3);
  spEff3Layer->SetLineWidth(1);
  spEff3Layer->GetXaxis()->SetTitle("Layer");
  spEff3Layer->GetYaxis()->SetTitle("Efficiency thr3");
  spEff3Layer->GetYaxis()->SetTitleOffset(1.2);
  spEff3Layer->GetYaxis()->SetRangeUser(0.,0.15);
  spEff3Layer->Draw("same");
  spEff3Layer->Fit(simu->fit_func);
  //simu->fit_func->SetParameter(0,simu_meanEff3);
  //simu->fit_func->Draw("same");
  dpEff3Layer->Draw();
  spEff3Layer->Draw("same");
  simu->fit_func->Draw("same");
  data->fit_func->Draw("same");
  sprintf(myText,"%s%.3f%s%.3f%s","#bar{#varepsilon} = ",data->fit_func->GetParameter(0)," #pm ",data_rmsEff3," (data)");
  tex->SetTextColor(kBlack);
  tex->DrawLatex(5,0.11,myText);
  sprintf(myText,"%s%.3f%s%.3f%s","#bar{#varepsilon} = ",simu->fit_func->GetParameter(0)," #pm ",simu_rmsEff3," (MC)");
  tex->SetTextColor(kRed-3);
  tex->DrawLatex(5,0.1,myText);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kGray+1);
  tex->DrawTextNDC(0.45,0.96,"CALICE Fe-SDHCAL Preliminary");
  cEff3Layer->Update();
  cEff3Layer->SaveAs("eff3Layer.pdf"); 
}
