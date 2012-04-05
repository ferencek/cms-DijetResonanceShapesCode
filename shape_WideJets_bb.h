#ifndef _SHAPE_H
#define _SHAPE_H


  double mass;

  double y500[50] = {8.63919753581, 10.565972805, 15.8584427014, 13.4988155067, 15.6897206157, 13.6339574978, 12.1723733246, 16.7235977799, 23.0183547288, 18.4156662226, 21.8455603942, 38.1367510408, 38.0371696576, 41.2867505252, 50.2177226916, 68.2909747213, 63.1335769072, 69.8968175948, 89.1973632947, 88.4260510132, 91.1276186407, 122.839684539, 155.224847145, 132.173345581, 161.11053402, 201.87949045, 217.452048026, 268.326398447, 306.03839083, 342.558861978, 361.987751886, 421.028459638, 485.285583228, 470.217274211, 450.81449578, 426.89573741, 356.90013735, 327.751857892, 187.827385634, 197.071182325, 129.494068407, 90.7300399467, 84.6154854074, 79.0503593907, 54.4894461781, 48.8232118934, 46.2270486951, 36.7662843913, 37.1843713671, 33.1825486422};

  double y700[50] = {4.87256503105, 12.0859398991, 6.45331510901, 12.8367160782, 10.0482997596, 19.5265097618, 16.8166934922, 21.6265918091, 15.1718264446, 23.9173628315, 29.1609811559, 37.949563995, 41.8084382415, 39.243013896, 55.68784298, 51.1770440713, 62.1421202049, 90.8551491424, 84.0130588189, 93.1299896315, 102.571446821, 120.706339769, 141.007545128, 170.406401426, 178.988948338, 190.113676913, 239.263770439, 256.137735859, 321.573653355, 372.373888567, 371.331002295, 441.644325517, 508.169467874, 526.534256823, 559.911106236, 441.517638572, 381.228109136, 255.97462558, 193.455719225, 134.076120652, 120.053920105, 101.667740799, 64.8055178151, 58.9952218458, 46.2457710505, 47.4511483535, 51.3627624363, 56.2819570452, 34.3656224534, 32.5063744634};

  double y1200[50] = {12.9906328619, 9.13238433748, 11.7810270786, 18.4756633416, 12.2236349955, 14.7309914529, 12.9405116439, 21.853356868, 23.1107827127, 19.1289537698, 27.6639861166, 34.5522132292, 36.9457847625, 46.2358465493, 63.4012108594, 66.6386065111, 56.5460076928, 58.6696620286, 95.6907039955, 89.0968650505, 109.256530195, 113.215220161, 150.09122254, 150.045135938, 192.775283292, 179.6151058, 252.015843101, 262.350732401, 358.446893059, 322.846842118, 391.926952429, 435.593522035, 482.525183231, 557.586941376, 552.898277983, 396.242357917, 353.43811655, 219.625860795, 120.786776736, 93.0482418835, 94.9597862512, 60.7640557736, 41.2592954934, 40.1192034334, 27.6479288116, 28.1032529548, 22.6308020055, 20.8612934723, 8.9441806376, 22.0063240528};

  double y2000[50] = {12.7808403969, 12.8971474171, 11.5245826095, 6.91264951974, 14.8460571915, 17.4230606034, 18.5579680726, 18.5090816095, 15.7750904933, 26.8413795754, 22.7342077866, 18.3849611878, 39.0352316797, 45.8532707542, 49.8993217498, 58.49915988, 83.1395265386, 78.4664730802, 81.6897537783, 116.465699479, 117.171752185, 114.109611511, 154.041480169, 166.541259013, 191.229624756, 231.573601, 245.325736843, 267.113146611, 326.885295145, 362.045482546, 407.557469837, 498.852247849, 540.949125223, 591.002313294, 569.244985722, 474.275295042, 300.85782446, 188.892152779, 101.665556587, 63.5028879121, 65.1184010878, 38.1909905002, 23.4359780401, 18.595880717, 18.7766828686, 7.86708504707, 10.1879217923, 7.43163228035, 2.83843743801, 10.4736208394};

  double y3500[50] = {30.9753884152, 39.0601946563, 39.7852644101, 45.6292080581, 37.1264160648, 44.0496849343, 39.1401996166, 34.3554473296, 53.4538888633, 46.8280930817, 64.4330922738, 75.5920766294, 64.6524626687, 70.7992405742, 71.9176153392, 78.6727952734, 78.5335683376, 97.8712078407, 97.6181919128, 126.511010461, 117.267872103, 133.547833122, 163.465247363, 138.051942274, 162.94176434, 203.978575662, 224.355625473, 239.260954745, 258.298511922, 321.47675211, 344.282890864, 401.985454746, 448.687736742, 553.469047569, 484.413155623, 381.286750063, 261.988334619, 98.2706677392, 69.3809216097, 28.5237169564, 21.9603617564, 15.0764331967, 4.39932066947, 6.48114935309, 1.44293212891, 4.52877554297, 1.42406487465, 3.37568724155, 1.42406487465, 1.44293212891};


  double bincenter[50]={
  0.31,  0.33,  0.35,  0.37,  0.39,  0.41,  0.43,  0.45,  0.47,  0.49,
  0.51,  0.53,  0.55,  0.57,  0.59,  0.61,  0.63,  0.65,  0.67,  0.69,
  0.71,  0.73,  0.75,  0.77,  0.79,  0.81,  0.83,  0.85,  0.87,  0.89,
  0.91,  0.93,  0.95,  0.97,  0.99,  1.01,  1.03,  1.05,  1.07,  1.09,
  1.11,  1.13,  1.15,  1.17,  1.19,  1.21,  1.23,  1.25,  1.27,  1.29};

 
  std::vector<double> v;

  double mqstar[5] = {500., 700., 1200., 2000., 3500.};

  const int nMassBins = 103;
  double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,
  354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
  1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,
  4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 
  10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};

  double binwidth[nMassBins+1] = {2, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 23, 24, 26, 26, 29, 29, 32, 33, 34, 36, 37, 39, 41, 43, 
  44, 47, 48, 50, 52, 54, 56, 58, 60, 63, 65, 67, 70, 72, 75, 77, 80, 83, 86, 89, 92, 95, 99, 101, 106, 108, 113, 116, 120, 124, 128, 132, 137, 142, 146,
  150, 156, 161, 166, 172, 177, 183, 189, 195, 202, 208, 214, 222, 229, 236, 244, 252, 260, 269, 277, 286, 295, 305, 315, 324, 335, 346, 358, 368, 381, 
  392, 406, 418, 432, 445, 460, 268};

  double massnew[nMassBins+1];
  double FoundQstarBinnedProb[nMassBins+1];
  
  const unsigned int n = 50;
  double x[n];
  double y[n];
  double f[n];


#endif
