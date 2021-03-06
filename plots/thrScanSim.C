{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Mar 13 18:18:04 2015) by ROOT version5.34/05
   TCanvas *c1 = new TCanvas("c1", "c1",11,52,800,800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-1.342567,-0.1175353,1.625526,1.058818);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogx();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1D *he = new TH1D("he"," ",30,0.09,30);
   he->SetMinimum(0.0001);
   he->SetMaximum(1);
   he->SetStats(0);
   he->GetXaxis()->SetTitle("Threshold [pC]");
   he->GetXaxis()->SetLabelFont(42);
   he->GetXaxis()->SetTitleFont(42);
   he->GetYaxis()->SetTitle("Efficiency ");
   he->GetYaxis()->SetNdivisions(505);
   he->GetYaxis()->SetLabelFont(42);
   he->GetYaxis()->SetTitleOffset(1.2);
   he->GetYaxis()->SetTitleFont(42);
   he->GetZaxis()->SetLabelFont(42);
   he->GetZaxis()->SetTitleSize(0.05);
   he->GetZaxis()->SetTitleFont(42);
   he->Draw("");
   
   TGraphErrors *gre = new TGraphErrors(250);
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#cc3333");
   gre->SetLineColor(ci);
   gre->SetLineWidth(2);

   ci = TColor::GetColor("#cc3333");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(0.8);
   gre->SetPoint(0,0.1,0.953684);
   gre->SetPointError(0,0,0.001571884);
   gre->SetPoint(1,0.1996,0.95167);
   gre->SetPointError(1,0,0.001604);
   gre->SetPoint(2,0.2992,0.947866);
   gre->SetPointError(2,0,0.001662595);
   gre->SetPoint(3,0.3988,0.942776);
   gre->SetPointError(3,0,0.001737185);
   gre->SetPoint(4,0.4984,0.935448);
   gre->SetPointError(4,0,0.001837881);
   gre->SetPoint(5,0.598,0.925882);
   gre->SetPointError(5,0,0.001959261);
   gre->SetPoint(6,0.6976,0.914863);
   gre->SetPointError(6,0,0.002087324);
   gre->SetPoint(7,0.7972,0.903451);
   gre->SetPointError(7,0,0.002208914);
   gre->SetPoint(8,0.8968,0.891089);
   gre->SetPointError(8,0,0.002329963);
   gre->SetPoint(9,0.9964,0.876881);
   gre->SetPointError(9,0,0.002457455);
   gre->SetPoint(10,1.096,0.862281);
   gre->SetPointError(10,0,0.002577354);
   gre->SetPoint(11,1.1956,0.84701);
   gre->SetPointError(11,0,0.002692331);
   gre->SetPoint(12,1.2952,0.830844);
   gre->SetPointError(12,0,0.002803859);
   gre->SetPoint(13,1.3948,0.814566);
   gre->SetPointError(13,0,0.002906769);
   gre->SetPoint(14,1.4944,0.797617);
   gre->SetPointError(14,0,0.003004948);
   gre->SetPoint(15,1.594,0.780668);
   gre->SetPointError(15,0,0.00309483);
   gre->SetPoint(16,1.6936,0.763998);
   gre->SetPointError(16,0,0.003175826);
   gre->SetPoint(17,1.7932,0.749007);
   gre->SetPointError(17,0,0.003242847);
   gre->SetPoint(18,1.8928,0.73217);
   gre->SetPointError(18,0,0.003311984);
   gre->SetPoint(19,1.9924,0.713654);
   gre->SetPointError(19,0,0.003380976);
   gre->SetPoint(20,2.092,0.69609);
   gre->SetPointError(20,0,0.003439995);
   gre->SetPoint(21,2.1916,0.676176);
   gre->SetPointError(21,0,0.00349975);
   gre->SetPoint(22,2.2912,0.659898);
   gre->SetPointError(22,0,0.0035432);
   gre->SetPoint(23,2.3908,0.64362);
   gre->SetPointError(23,0,0.003581988);
   gre->SetPoint(24,2.4904,0.628461);
   gre->SetPointError(24,0,0.003614049);
   gre->SetPoint(25,2.59,0.612743);
   gre->SetPointError(25,0,0.003643271);
   gre->SetPoint(26,2.6896,0.596856);
   gre->SetPointError(26,0,0.003668745);
   gre->SetPoint(27,2.7892,0.581473);
   gre->SetPointError(27,0,0.003689599);
   gre->SetPoint(28,2.8888,0.566091);
   gre->SetPointError(28,0,0.003706766);
   gre->SetPoint(29,2.9884,0.550708);
   gre->SetPointError(29,0,0.003720298);
   gre->SetPoint(30,3.088,0.536164);
   gre->SetPointError(30,0,0.003729784);
   gre->SetPoint(31,3.1876,0.521061);
   gre->SetPointError(31,0,0.00373626);
   gre->SetPoint(32,3.2872,0.505901);
   gre->SetPointError(32,0,0.003739318);
   gre->SetPoint(33,3.3868,0.490966);
   gre->SetPointError(33,0,0.003738968);
   gre->SetPoint(34,3.4864,0.476814);
   gre->SetPointError(34,0,0.003735556);
   gre->SetPoint(35,3.586,0.46378);
   gre->SetPointError(35,0,0.003729754);
   gre->SetPoint(36,3.6856,0.451362);
   gre->SetPointError(36,0,0.003721844);
   gre->SetPoint(37,3.7852,0.438944);
   gre->SetPointError(37,0,0.003711593);
   gre->SetPoint(38,3.8848,0.42591);
   gre->SetPointError(38,0,0.003698295);
   gre->SetPoint(39,3.9844,0.414387);
   gre->SetPointError(39,0,0.003684352);
   gre->SetPoint(40,4.084,0.40141);
   gre->SetPointError(40,0,0.003666161);
   gre->SetPoint(41,4.1836,0.390446);
   gre->SetPointError(41,0,0.003648709);
   gre->SetPoint(42,4.2832,0.378475);
   gre->SetPointError(42,0,0.003627443);
   gre->SetPoint(43,4.3828,0.367959);
   gre->SetPointError(43,0,0.003606824);
   gre->SetPoint(44,4.4824,0.357443);
   gre->SetPointError(44,0,0.003584362);
   gre->SetPoint(45,4.582,0.347094);
   gre->SetPointError(45,0,0.003560423);
   gre->SetPoint(46,4.6816,0.338368);
   gre->SetPointError(46,0,0.003538796);
   gre->SetPoint(47,4.7812,0.328075);
   gre->SetPointError(47,0,0.003511556);
   gre->SetPoint(48,4.8808,0.318566);
   gre->SetPointError(48,0,0.003484691);
   gre->SetPoint(49,4.9804,0.308833);
   gre->SetPointError(49,0,0.003455461);
   gre->SetPoint(50,5.08,0.299994);
   gre->SetPointError(50,0,0.003427361);
   gre->SetPoint(51,5.1796,0.292275);
   gre->SetPointError(51,0,0.003401581);
   gre->SetPoint(52,5.2792,0.284444);
   gre->SetPointError(52,0,0.003374216);
   gre->SetPoint(53,5.3788,0.276556);
   gre->SetPointError(53,0,0.003345389);
   gre->SetPoint(54,5.4784,0.268222);
   gre->SetPointError(54,0,0.00331352);
   gre->SetPoint(55,5.578,0.260279);
   gre->SetPointError(55,0,0.003281755);
   gre->SetPoint(56,5.6776,0.252559);
   gre->SetPointError(56,0,0.003249545);
   gre->SetPoint(57,5.7772,0.24456);
   gre->SetPointError(57,0,0.003214736);
   gre->SetPoint(58,5.8768,0.237512);
   gre->SetPointError(58,0,0.003182819);
   gre->SetPoint(59,5.9764,0.230352);
   gre->SetPointError(59,0,0.00314916);
   gre->SetPoint(60,6.076,0.223863);
   gre->SetPointError(60,0,0.003117547);
   gre->SetPoint(61,6.1756,0.216815);
   gre->SetPointError(61,0,0.003081978);
   gre->SetPoint(62,6.2752,0.210997);
   gre->SetPointError(62,0,0.003051618);
   gre->SetPoint(63,6.3748,0.204676);
   gre->SetPointError(63,0,0.003017576);
   gre->SetPoint(64,6.4744,0.198859);
   gre->SetPointError(64,0,0.002985244);
   gre->SetPoint(65,6.574,0.192929);
   gre->SetPointError(65,0,0.002951259);
   gre->SetPoint(66,6.6736,0.187951);
   gre->SetPointError(66,0,0.002921905);
   gre->SetPoint(67,6.7732,0.182022);
   gre->SetPointError(67,0,0.002885928);
   gre->SetPoint(68,6.8728,0.177882);
   gre->SetPointError(68,0,0.00286013);
   gre->SetPoint(69,6.9724,0.173407);
   gre->SetPointError(69,0,0.0028316);
   gre->SetPoint(70,7.072,0.168932);
   gre->SetPointError(70,0,0.002802379);
   gre->SetPoint(71,7.1716,0.164401);
   gre->SetPointError(71,0,0.002772068);
   gre->SetPoint(72,7.2712,0.159143);
   gre->SetPointError(72,0,0.002735946);
   gre->SetPoint(73,7.3708,0.154556);
   gre->SetPointError(73,0,0.002703573);
   gre->SetPoint(74,7.4704,0.149857);
   gre->SetPointError(74,0,0.002669545);
   gre->SetPoint(75,7.57,0.14555);
   gre->SetPointError(75,0,0.002637559);
   gre->SetPoint(76,7.6696,0.141746);
   gre->SetPointError(76,0,0.002608651);
   gre->SetPoint(77,7.7692,0.137999);
   gre->SetPointError(77,0,0.002579554);
   gre->SetPoint(78,7.8688,0.134251);
   gre->SetPointError(78,0,0.002549808);
   gre->SetPoint(79,7.9684,0.129999);
   gre->SetPointError(79,0,0.002515258);
   gre->SetPoint(80,8.068,0.125916);
   gre->SetPointError(80,0,0.002481245);
   gre->SetPoint(81,8.1676,0.12228);
   gre->SetPointError(81,0,0.002450239);
   gre->SetPoint(82,8.2672,0.118029);
   gre->SetPointError(82,0,0.002413094);
   gre->SetPoint(83,8.3668,0.114505);
   gre->SetPointError(83,0,0.00238154);
   gre->SetPoint(84,8.4664,0.111484);
   gre->SetPointError(84,0,0.002353919);
   gre->SetPoint(85,8.566,0.108128);
   gre->SetPointError(85,0,0.002322593);
   gre->SetPoint(86,8.6656,0.105163);
   gre->SetPointError(86,0,0.002294331);
   gre->SetPoint(87,8.7652,0.101975);
   gre->SetPointError(87,0,0.002263309);
   gre->SetPoint(88,8.8648,0.0986743);
   gre->SetPointError(88,0,0.002230466);
   gre->SetPoint(89,8.9644,0.0960452);
   gre->SetPointError(89,0,0.002203758);
   gre->SetPoint(90,9.064,0.0938636);
   gre->SetPointError(90,0,0.002181213);
   gre->SetPoint(91,9.1636,0.0909549);
   gre->SetPointError(91,0,0.002150594);
   gre->SetPoint(92,9.2632,0.0885495);
   gre->SetPointError(92,0,0.002124772);
   gre->SetPoint(93,9.3628,0.0856967);
   gre->SetPointError(93,0,0.002093533);
   gre->SetPoint(94,9.4624,0.0830676);
   gre->SetPointError(94,0,0.002064131);
   gre->SetPoint(95,9.562,0.0811657);
   gre->SetPointError(95,0,0.002042479);
   gre->SetPoint(96,9.6616,0.0789282);
   gre->SetPointError(96,0,0.00201658);
   gre->SetPoint(97,9.7612,0.0769145);
   gre->SetPointError(97,0,0.001992864);
   gre->SetPoint(98,9.8608,0.074621);
   gre->SetPointError(98,0,0.001965364);
   gre->SetPoint(99,9.9604,0.0724394);
   gre->SetPointError(99,0,0.001938703);
   gre->SetPoint(100,10.06,0.0700901);
   gre->SetPointError(100,0,0.00190942);
   gre->SetPoint(101,10.1596,0.068356);
   gre->SetPointError(101,0,0.001887409);
   gre->SetPoint(102,10.2592,0.0663982);
   gre->SetPointError(102,0,0.001862137);
   gre->SetPoint(103,10.3588,0.0646641);
   gre->SetPointError(103,0,0.001839366);
   gre->SetPoint(104,10.4584,0.0624825);
   gre->SetPointError(104,0,0.001810179);
   gre->SetPoint(105,10.558,0.0606366);
   gre->SetPointError(105,0,0.001784995);
   gre->SetPoint(106,10.6576,0.0596297);
   gre->SetPointError(106,0,0.001771061);
   gre->SetPoint(107,10.7572,0.0577278);
   gre->SetPointError(107,0,0.001744349);
   gre->SetPoint(108,10.8568,0.0561615);
   gre->SetPointError(108,0,0.001721951);
   gre->SetPoint(109,10.9564,0.0547631);
   gre->SetPointError(109,0,0.001701637);
   gre->SetPoint(110,11.056,0.0533647);
   gre->SetPointError(110,0,0.001681013);
   gre->SetPoint(111,11.1556,0.0523578);
   gre->SetPointError(111,0,0.001665964);
   gre->SetPoint(112,11.2552,0.0509593);
   gre->SetPointError(112,0,0.001644776);
   gre->SetPoint(113,11.3548,0.0496728);
   gre->SetPointError(113,0,0.001624982);
   gre->SetPoint(114,11.4544,0.0488896);
   gre->SetPointError(114,0,0.001612785);
   gre->SetPoint(115,11.554,0.0479387);
   gre->SetPointError(115,0,0.001597821);
   gre->SetPoint(116,11.6536,0.0469318);
   gre->SetPointError(116,0,0.001581788);
   gre->SetPoint(117,11.7532,0.0460368);
   gre->SetPointError(117,0,0.001567368);
   gre->SetPoint(118,11.8528,0.0451418);
   gre->SetPointError(118,0,0.001552786);
   gre->SetPoint(119,11.9524,0.0439671);
   gre->SetPointError(119,0,0.001533391);
   gre->SetPoint(120,12.052,0.0425687);
   gre->SetPointError(120,0,0.001509912);
   gre->SetPoint(121,12.1516,0.0415618);
   gre->SetPointError(121,0,0.001492732);
   gre->SetPoint(122,12.2512,0.0406108);
   gre->SetPointError(122,0,0.001476287);
   gre->SetPoint(123,12.3508,0.0393802);
   gre->SetPointError(123,0,0.00145468);
   gre->SetPoint(124,12.4504,0.0380936);
   gre->SetPointError(124,0,0.001431677);
   gre->SetPoint(125,12.55,0.0370868);
   gre->SetPointError(125,0,0.00141337);
   gre->SetPoint(126,12.6496,0.0361918);
   gre->SetPointError(126,0,0.001396861);
   gre->SetPoint(127,12.7492,0.0355205);
   gre->SetPointError(127,0,0.001384327);
   gre->SetPoint(128,12.8488,0.0345136);
   gre->SetPointError(128,0,0.001365277);
   gre->SetPoint(129,12.9484,0.0333389);
   gre->SetPointError(129,0,0.001342658);
   gre->SetPoint(130,13.048,0.0325558);
   gre->SetPointError(130,0,0.001327333);
   gre->SetPoint(131,13.1476,0.0318286);
   gre->SetPointError(131,0,0.001312918);
   gre->SetPoint(132,13.2472,0.0312133);
   gre->SetPointError(132,0,0.001300579);
   gre->SetPoint(133,13.3468,0.030598);
   gre->SetPointError(133,0,0.001288105);
   gre->SetPoint(134,13.4464,0.0299267);
   gre->SetPointError(134,0,0.001274337);
   gre->SetPoint(135,13.546,0.0294233);
   gre->SetPointError(135,0,0.001263902);
   gre->SetPoint(136,13.6456,0.028752);
   gre->SetPointError(136,0,0.001249832);
   gre->SetPoint(137,13.7452,0.0281367);
   gre->SetPointError(137,0,0.001236778);
   gre->SetPoint(138,13.8448,0.0272417);
   gre->SetPointError(138,0,0.001217509);
   gre->SetPoint(139,13.9444,0.0268501);
   gre->SetPointError(139,0,0.00120897);
   gre->SetPoint(140,14.044,0.0262348);
   gre->SetPointError(140,0,0.001195415);
   gre->SetPoint(141,14.1436,0.0257873);
   gre->SetPointError(141,0,0.001185448);
   gre->SetPoint(142,14.2432,0.0250601);
   gre->SetPointError(142,0,0.00116905);
   gre->SetPoint(143,14.3428,0.024277);
   gre->SetPointError(143,0,0.001151101);
   gre->SetPoint(144,14.4424,0.0238295);
   gre->SetPointError(144,0,0.001140704);
   gre->SetPoint(145,14.542,0.0234939);
   gre->SetPointError(145,0,0.001132838);
   gre->SetPoint(146,14.6416,0.0231023);
   gre->SetPointError(146,0,0.001123582);
   gre->SetPoint(147,14.7412,0.0226548);
   gre->SetPointError(147,0,0.001112902);
   gre->SetPoint(148,14.8408,0.0220395);
   gre->SetPointError(148,0,0.00109803);
   gre->SetPoint(149,14.9404,0.0217598);
   gre->SetPointError(149,0,0.001091196);
   gre->SetPoint(150,15.04,0.021592);
   gre->SetPointError(150,0,0.001087074);
   gre->SetPoint(151,15.1396,0.0212004);
   gre->SetPointError(151,0,0.001077387);
   gre->SetPoint(152,15.2392,0.0204173);
   gre->SetPointError(152,0,0.001057724);
   gre->SetPoint(153,15.3388,0.0199698);
   gre->SetPointError(153,0,0.001046307);
   gre->SetPoint(154,15.4384,0.0196901);
   gre->SetPointError(154,0,0.001039102);
   gre->SetPoint(155,15.538,0.0190189);
   gre->SetPointError(155,0,0.001021588);
   gre->SetPoint(156,15.6376,0.0186832);
   gre->SetPointError(156,0,0.001012705);
   gre->SetPoint(157,15.7372,0.0182357);
   gre->SetPointError(157,0,0.001000731);
   gre->SetPoint(158,15.8368,0.0177323);
   gre->SetPointError(158,0,0.0009870749);
   gre->SetPoint(159,15.9364,0.0173967);
   gre->SetPointError(159,0,0.0009778568);
   gre->SetPoint(160,16.036,0.0170051);
   gre->SetPointError(160,0,0.0009669809);
   gre->SetPoint(161,16.1356,0.0168932);
   gre->SetPointError(161,0,0.000963849);
   gre->SetPoint(162,16.2352,0.0163898);
   gre->SetPointError(162,0,0.0009496226);
   gre->SetPoint(163,16.3348,0.0158863);
   gre->SetPointError(163,0,0.0009351616);
   gre->SetPoint(164,16.4344,0.0156626);
   gre->SetPointError(164,0,0.0009286597);
   gre->SetPoint(165,16.534,0.0154388);
   gre->SetPointError(165,0,0.0009221059);
   gre->SetPoint(166,16.6336,0.0151591);
   gre->SetPointError(166,0,0.0009138447);
   gre->SetPoint(167,16.7332,0.0149913);
   gre->SetPointError(167,0,0.0009088503);
   gre->SetPoint(168,16.8328,0.0145998);
   gre->SetPointError(168,0,0.0008970826);
   gre->SetPoint(169,16.9324,0.0140404);
   gre->SetPointError(169,0,0.0008799784);
   gre->SetPoint(170,17.032,0.0137048);
   gre->SetPointError(170,0,0.0008695459);
   gre->SetPoint(171,17.1316,0.013481);
   gre->SetPointError(171,0,0.0008625146);
   gre->SetPoint(172,17.2312,0.0133691);
   gre->SetPointError(172,0,0.0008589762);
   gre->SetPoint(173,17.3308,0.0130894);
   gre->SetPointError(173,0,0.0008500636);
   gre->SetPoint(174,17.4304,0.0128098);
   gre->SetPointError(174,0,0.0008410548);
   gre->SetPoint(175,17.53,0.0127538);
   gre->SetPointError(175,0,0.0008392381);
   gre->SetPoint(176,17.6296,0.0123623);
   gre->SetPointError(176,0,0.0008264206);
   gre->SetPoint(177,17.7292,0.0119148);
   gre->SetPointError(177,0,0.0008115089);
   gre->SetPoint(178,17.8288,0.011691);
   gre->SetPointError(178,0,0.0008039423);
   gre->SetPoint(179,17.9284,0.0115232);
   gre->SetPointError(179,0,0.0007982198);
   gre->SetPoint(180,18.028,0.0112994);
   gre->SetPointError(180,0,0.0007905198);
   gre->SetPoint(181,18.1276,0.0110757);
   gre->SetPointError(181,0,0.0007827441);
   gre->SetPoint(182,18.2272,0.0109638);
   gre->SetPointError(182,0,0.000778824);
   gre->SetPoint(183,18.3268,0.0107401);
   gre->SetPointError(183,0,0.0007709248);
   gre->SetPoint(184,18.4264,0.0106841);
   gre->SetPointError(184,0,0.0007689342);
   gre->SetPoint(185,18.526,0.0104044);
   gre->SetPointError(185,0,0.0007589096);
   gre->SetPoint(186,18.6256,0.0101807);
   gre->SetPointError(186,0,0.0007507917);
   gre->SetPoint(187,18.7252,0.0100129);
   gre->SetPointError(187,0,0.0007446418);
   gre->SetPoint(188,18.8248,0.00973318);
   gre->SetPointError(188,0,0.0007342706);
   gre->SetPoint(189,18.9244,0.00967724);
   gre->SetPointError(189,0,0.0007321782);
   gre->SetPoint(190,19.024,0.00950943);
   gre->SetPointError(190,0,0.0007258637);
   gre->SetPoint(191,19.1236,0.00934161);
   gre->SetPointError(191,0,0.0007194912);
   gre->SetPoint(192,19.2232,0.00922974);
   gre->SetPointError(192,0,0.0007152105);
   gre->SetPoint(193,19.3228,0.00889411);
   gre->SetPointError(193,0,0.0007022051);
   gre->SetPoint(194,19.4224,0.00889411);
   gre->SetPointError(194,0,0.0007022051);
   gre->SetPoint(195,19.522,0.00883817);
   gre->SetPointError(195,0,0.000700013);
   gre->SetPoint(196,19.6216,0.00861442);
   gre->SetPointError(196,0,0.0006911734);
   gre->SetPoint(197,19.7212,0.00855848);
   gre->SetPointError(197,0,0.000688945);
   gre->SetPoint(198,19.8208,0.00833473);
   gre->SetPointError(198,0,0.0006799563);
   gre->SetPoint(199,19.9204,0.00816692);
   gre->SetPointError(199,0,0.0006731334);
   gre->SetPoint(200,20.02,0.00794317);
   gre->SetPointError(200,0,0.0006639233);
   gre->SetPoint(201,20.1196,0.00788723);
   gre->SetPointError(201,0,0.0006615999);
   gre->SetPoint(202,20.2192,0.00783129);
   gre->SetPointError(202,0,0.0006592682);
   gre->SetPoint(203,20.3188,0.00771942);
   gre->SetPointError(203,0,0.0006545793);
   gre->SetPoint(204,20.4184,0.0075516);
   gre->SetPointError(204,0,0.0006474797);
   gre->SetPoint(205,20.518,0.00743973);
   gre->SetPointError(205,0,0.0006427021);
   gre->SetPoint(206,20.6176,0.00727191);
   gre->SetPointError(206,0,0.0006354657);
   gre->SetPoint(207,20.7172,0.00716004);
   gre->SetPointError(207,0,0.0006305943);
   gre->SetPoint(208,20.8168,0.00699222);
   gre->SetPointError(208,0,0.0006232131);
   gre->SetPoint(209,20.9164,0.00693629);
   gre->SetPointError(209,0,0.000620733);
   gre->SetPoint(210,21.016,0.00688035);
   gre->SetPointError(210,0,0.0006182423);
   gre->SetPoint(211,21.1156,0.00682441);
   gre->SetPointError(211,0,0.0006157412);
   gre->SetPoint(212,21.2152,0.00660066);
   gre->SetPointError(212,0,0.0006056313);
   gre->SetPoint(213,21.3148,0.00660066);
   gre->SetPointError(213,0,0.0006056313);
   gre->SetPoint(214,21.4144,0.00648878);
   gre->SetPointError(214,0,0.0006005105);
   gre->SetPoint(215,21.514,0.00637691);
   gre->SetPointError(215,0,0.0005953449);
   gre->SetPoint(216,21.6136,0.00632097);
   gre->SetPointError(216,0,0.0005927446);
   gre->SetPoint(217,21.7132,0.0062091);
   gre->SetPointError(217,0,0.000587509);
   gre->SetPoint(218,21.8128,0.00604128);
   gre->SetPointError(218,0,0.0005795639);
   gre->SetPoint(219,21.9124,0.00604128);
   gre->SetPointError(219,0,0.0005795639);
   gre->SetPoint(220,22.012,0.00592941);
   gre->SetPointError(220,0,0.0005742051);
   gre->SetPoint(221,22.1116,0.00581753);
   gre->SetPointError(221,0,0.000568794);
   gre->SetPoint(222,22.2112,0.00553784);
   gre->SetPointError(222,0,0.0005550307);
   gre->SetPoint(223,22.3108,0.0054819);
   gre->SetPointError(223,0,0.0005522358);
   gre->SetPoint(224,22.4104,0.0054819);
   gre->SetPointError(224,0,0.0005522358);
   gre->SetPoint(225,22.51,0.00542597);
   gre->SetPointError(225,0,0.0005494269);
   gre->SetPoint(226,22.6096,0.00531409);
   gre->SetPointError(226,0,0.0005437636);
   gre->SetPoint(227,22.7092,0.00531409);
   gre->SetPointError(227,0,0.0005437636);
   gre->SetPoint(228,22.8088,0.00520222);
   gre->SetPointError(228,0,0.0005380398);
   gre->SetPoint(229,22.9084,0.00520222);
   gre->SetPointError(229,0,0.0005380398);
   gre->SetPoint(230,23.008,0.00514628);
   gre->SetPointError(230,0,0.0005351543);
   gre->SetPoint(231,23.1076,0.00509034);
   gre->SetPointError(231,0,0.0005322527);
   gre->SetPoint(232,23.2072,0.00509034);
   gre->SetPointError(232,0,0.0005322527);
   gre->SetPoint(233,23.3068,0.0050344);
   gre->SetPointError(233,0,0.0005293349);
   gre->SetPoint(234,23.4064,0.00486659);
   gre->SetPointError(234,0,0.000520482);
   gre->SetPoint(235,23.506,0.00486659);
   gre->SetPointError(235,0,0.000520482);
   gre->SetPoint(236,23.6056,0.00464284);
   gre->SetPointError(236,0,0.0005084333);
   gre->SetPoint(237,23.7052,0.00436315);
   gre->SetPointError(237,0,0.0004929504);
   gre->SetPoint(238,23.8048,0.00425127);
   gre->SetPointError(238,0,0.0004866166);
   gre->SetPoint(239,23.9044,0.0041394);
   gre->SetPointError(239,0,0.0004801983);
   gre->SetPoint(240,24.004,0.00397158);
   gre->SetPointError(240,0,0.0004704031);
   gre->SetPoint(241,24.1036,0.00380377);
   gre->SetPointError(241,0,0.0004603968);
   gre->SetPoint(242,24.2032,0.00380377);
   gre->SetPointError(242,0,0.0004603968);
   gre->SetPoint(243,24.3028,0.00374783);
   gre->SetPointError(243,0,0.0004570116);
   gre->SetPoint(244,24.4024,0.00363596);
   gre->SetPointError(244,0,0.0004501645);
   gre->SetPoint(245,24.502,0.00358002);
   gre->SetPointError(245,0,0.0004467007);
   gre->SetPoint(246,24.6016,0.00346814);
   gre->SetPointError(246,0,0.00043969);
   gre->SetPoint(247,24.7012,0.00346814);
   gre->SetPointError(247,0,0.00043969);
   gre->SetPoint(248,24.8008,0.00335627);
   gre->SetPointError(248,0,0.0004325647);
   gre->SetPoint(249,24.9004,0.00335627);
   gre->SetPointError(249,0,0.0004325647);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",250,0.09,27.38044);
   Graph_Graph1->SetMinimum(0);
   Graph_Graph1->SetMaximum(1.050489);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1);
   
   
   TF1 *func = new TF1("*func",0.09,27.38044,4);
    //The original function :  had originally been created by:
    //TF1 *func = new TF1("func",,0.09,27.38044,4);
   func->SetRange(0.09,27.38044);
   func->SetName("func");
   func->SetTitle("");
   func->SetSavedPoint(0,0.9716154);
   func->SetSavedPoint(1,0.9489539);
   func->SetSavedPoint(2,0.9156439);
   func->SetSavedPoint(3,0.8765235);
   func->SetSavedPoint(4,0.8341549);
   func->SetSavedPoint(5,0.7901777);
   func->SetSavedPoint(6,0.745714);
   func->SetSavedPoint(7,0.7015542);
   func->SetSavedPoint(8,0.6582603);
   func->SetSavedPoint(9,0.61623);
   func->SetSavedPoint(10,0.5757394);
   func->SetSavedPoint(11,0.5369731);
   func->SetSavedPoint(12,0.5000467);
   func->SetSavedPoint(13,0.4650229);
   func->SetSavedPoint(14,0.4319239);
   func->SetSavedPoint(15,0.4007416);
   func->SetSavedPoint(16,0.3714444);
   func->SetSavedPoint(17,0.343984);
   func->SetSavedPoint(18,0.3182992);
   func->SetSavedPoint(19,0.2943203);
   func->SetSavedPoint(20,0.2719716);
   func->SetSavedPoint(21,0.2511738);
   func->SetSavedPoint(22,0.2318459);
   func->SetSavedPoint(23,0.2139063);
   func->SetSavedPoint(24,0.1972744);
   func->SetSavedPoint(25,0.181871);
   func->SetSavedPoint(26,0.1676192);
   func->SetSavedPoint(27,0.1544447);
   func->SetSavedPoint(28,0.1422759);
   func->SetSavedPoint(29,0.1310449);
   func->SetSavedPoint(30,0.1206868);
   func->SetSavedPoint(31,0.1111401);
   func->SetSavedPoint(32,0.1023467);
   func->SetSavedPoint(33,0.09425207);
   func->SetSavedPoint(34,0.08680466);
   func->SetSavedPoint(35,0.07995634);
   func->SetSavedPoint(36,0.07366198);
   func->SetSavedPoint(37,0.06787947);
   func->SetSavedPoint(38,0.06256952);
   func->SetSavedPoint(39,0.05769552);
   func->SetSavedPoint(40,0.05322346);
   func->SetSavedPoint(41,0.04912172);
   func->SetSavedPoint(42,0.04536099);
   func->SetSavedPoint(43,0.04191407);
   func->SetSavedPoint(44,0.03875583);
   func->SetSavedPoint(45,0.03586296);
   func->SetSavedPoint(46,0.03321396);
   func->SetSavedPoint(47,0.03078896);
   func->SetSavedPoint(48,0.02856962);
   func->SetSavedPoint(49,0.02653902);
   func->SetSavedPoint(50,0.02468158);
   func->SetSavedPoint(51,0.02298294);
   func->SetSavedPoint(52,0.02142987);
   func->SetSavedPoint(53,0.02001023);
   func->SetSavedPoint(54,0.01871281);
   func->SetSavedPoint(55,0.01752735);
   func->SetSavedPoint(56,0.01644439);
   func->SetSavedPoint(57,0.01545526);
   func->SetSavedPoint(58,0.014552);
   func->SetSavedPoint(59,0.0137273);
   func->SetSavedPoint(60,0.01297446);
   func->SetSavedPoint(61,0.01228732);
   func->SetSavedPoint(62,0.01166027);
   func->SetSavedPoint(63,0.01108812);
   func->SetSavedPoint(64,0.01056616);
   func->SetSavedPoint(65,0.01009004);
   func->SetSavedPoint(66,0.009655812);
   func->SetSavedPoint(67,0.009259835);
   func->SetSavedPoint(68,0.008898791);
   func->SetSavedPoint(69,0.00856964);
   func->SetSavedPoint(70,0.008269602);
   func->SetSavedPoint(71,0.007996135);
   func->SetSavedPoint(72,0.007746917);
   func->SetSavedPoint(73,0.007519822);
   func->SetSavedPoint(74,0.007312911);
   func->SetSavedPoint(75,0.00712441);
   func->SetSavedPoint(76,0.006952699);
   func->SetSavedPoint(77,0.006796298);
   func->SetSavedPoint(78,0.006653858);
   func->SetSavedPoint(79,0.006524144);
   func->SetSavedPoint(80,0.006406032);
   func->SetSavedPoint(81,0.006298492);
   func->SetSavedPoint(82,0.006200589);
   func->SetSavedPoint(83,0.006111466);
   func->SetSavedPoint(84,0.006030342);
   func->SetSavedPoint(85,0.005956507);
   func->SetSavedPoint(86,0.005889311);
   func->SetSavedPoint(87,0.005828161);
   func->SetSavedPoint(88,0.005772519);
   func->SetSavedPoint(89,0.005721891);
   func->SetSavedPoint(90,0.00567583);
   func->SetSavedPoint(91,0.005633927);
   func->SetSavedPoint(92,0.005595809);
   func->SetSavedPoint(93,0.005561137);
   func->SetSavedPoint(94,0.005529601);
   func->SetSavedPoint(95,0.005500919);
   func->SetSavedPoint(96,0.005474836);
   func->SetSavedPoint(97,0.005451117);
   func->SetSavedPoint(98,0.005429549);
   func->SetSavedPoint(99,0.005409938);
   func->SetSavedPoint(100,0.005392108);
   func->SetSavedPoint(101,0.09);
   func->SetSavedPoint(102,27.38044);
   func->SetFillColor(19);
   func->SetFillStyle(0);

   ci = TColor::GetColor("#cc3333");
   func->SetLineColor(ci);
   func->SetLineWidth(2);
   func->SetChisquare(1684.302);
   func->SetNDF(246);
   func->GetXaxis()->SetLabelFont(42);
   func->GetXaxis()->SetTitleSize(0.05);
   func->GetXaxis()->SetTitleFont(42);
   func->GetYaxis()->SetLabelFont(42);
   func->GetYaxis()->SetTitleSize(0.05);
   func->GetYaxis()->SetTitleFont(42);
   func->SetParameter(0,4.300539);
   func->SetParError(0,0.0047171);
   func->SetParLimits(0,0,0);
   func->SetParameter(1,0.5945082);
   func->SetParError(1,0.005720835);
   func->SetParLimits(1,0,0);
   func->SetParameter(2,0.4025105);
   func->SetParError(2,0.0008844036);
   func->SetParLimits(2,0,0);
   func->SetParameter(3,0.974564);
   func->SetParError(3,0.0008091988);
   func->SetParLimits(3,0,0);
   gre->GetListOfFunctions()->Add(func);
   gre->Draw("p");
   TLatex *   tex = new TLatex(0.2,0.44,"Simulation (FTFP_BERT_HP)");

   ci = TColor::GetColor("#cc3333");
   tex->SetTextColor(ci);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.38,"#bar{q} = 4.301 #pm 0.005 pC");

   ci = TColor::GetColor("#cc3333");
   tex->SetTextColor(ci);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.32,"#theta = 0.595 #pm 0.006");

   ci = TColor::GetColor("#cc3333");
   tex->SetTextColor(ci);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.26,"c = 0.403 #pm 0.001");

   ci = TColor::GetColor("#cc3333");
   tex->SetTextColor(ci);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.2,"#varepsilon_{0} = 0.975 #pm 0.001");

   ci = TColor::GetColor("#cc3333");
   tex->SetTextColor(ci);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   TText *text = new TText(0.5,0.96,"CALICE Fe-SDHCAL Preliminary");
   text->SetNDC();

   ci = TColor::GetColor("#999999");
   text->SetTextColor(ci);
   text->SetTextSize(0.03);
   text->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
