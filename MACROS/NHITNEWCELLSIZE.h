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
#include <sstream>
#include <vector>

struct result{
  int ebeam;
  float nhit;
  float nhitError;
  float sigma;
  float sigmaError;
  float resol;
};
