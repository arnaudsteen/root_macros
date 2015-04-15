#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TColor.h>
#include <TLegendEntry.h>
#include <TText.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

struct result{
  int ebeam;
  float erec;
  float erecError;
  float sigma;
  float sigmaError;
  float resol;
};

Double_t fitfunc(Double_t* x, Double_t* par)
{
  return sqrt(par[0]*par[0]/x[0] + par[1]*par[1] /*+ par[2]*par[2]/(x[0]*x[0]) + par[2]*par[2]*x[0]*x[0]*/);
}
