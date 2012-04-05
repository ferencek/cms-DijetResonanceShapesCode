#ifndef _SHAPE_H
#define _SHAPE_H


  double mass;

  double y500[50] = {233.211758025, 212.548281379, 263.885323055, 324.417027481, 347.310794719, 408.370811597, 376.727311507, 427.641124882, 533.253159232, 533.836596705, 559.337332994, 597.888621502, 596.879688606, 674.302453399, 658.223142825, 678.604302846, 799.648474477, 805.342265993, 762.322929464, 849.008808747, 774.866985507, 856.11215394, 858.172066428, 906.238379888, 843.950899035, 942.659356073, 941.483239956, 1015.89918748, 1113.92783407, 1111.5946821, 1157.58392455, 1106.67341589, 1114.09747763, 982.956133716, 902.219825879, 720.155996904, 550.018648639, 415.70027151, 295.955113284, 210.731094986, 175.564183146, 164.415960468, 113.092347749, 145.295570694, 96.1697834805, 90.5854704678, 89.4562695399, 80.0072520003, 63.0916513205, 69.5978983566};

  double y700[50] = {211.502931848, 263.335296229, 261.746750064, 292.920953311, 389.231716692, 357.175609514, 408.146611266, 423.999271676, 464.881930381, 524.316823184, 520.130386598, 600.242903315, 665.553674519, 669.215208597, 767.623337559, 786.49999854, 781.85998065, 843.836332209, 792.82130041, 868.618666396, 998.897435173, 877.159534521, 981.180457115, 978.671715006, 971.868330687, 1054.44797949, 1010.78806446, 1116.810585, 1220.83093685, 1236.19607235, 1404.95131401, 1363.07260329, 1332.77024482, 1252.04691617, 1051.81602374, 737.812952295, 602.457645461, 393.112122811, 309.434168994, 175.615613289, 156.192336954, 141.740721002, 117.039136484, 74.6855634376, 87.3555637822, 80.6283562705, 72.2604909614, 75.804269731, 67.0507851914, 59.1718202755};

  double y1200[50] = {180.170186952, 202.789023325, 239.048505716, 254.021308474, 316.554121025, 301.613063104, 364.056341395, 405.156977169, 467.479268178, 509.085287139, 537.886336975, 567.881705053, 602.029779598, 711.240822352, 762.740588993, 796.507143006, 808.930614613, 872.489572667, 918.554247141, 1029.34733713, 988.129097499, 994.219447531, 1064.28926568, 1059.96848109, 1111.82244322, 1103.531253, 1165.41521338, 1226.05544438, 1351.16697152, 1446.16990396, 1664.83496569, 1805.08093908, 1852.48803541, 1657.62050958, 1361.01401091, 898.186476037, 570.171217456, 283.55006057, 149.856788531, 133.016292103, 89.9263125211, 81.1358371824, 65.3749836758, 63.7764514387, 37.3461385965, 41.148951076, 32.7470789477, 29.2953353152, 36.6849231198, 28.6370839104};

  double y2000[50] = {140.225067213, 204.109986119, 215.905848622, 222.213761687, 275.366276339, 294.398044772, 335.092498943, 382.477839962, 438.606709041, 500.89656169, 540.159552321, 592.213620394, 670.797617301, 704.10590107, 806.132229552, 864.3581824, 956.255439706, 1126.43333878, 1018.72338944, 1132.42822852, 1103.31461114, 1228.75633851, 1264.39920662, 1208.55884646, 1408.53231827, 1320.37929665, 1426.84418418, 1494.54814389, 1674.00850032, 1912.47915522, 2078.94624668, 2282.54342113, 2467.78347927, 2294.34993497, 1627.20944949, 992.507549487, 430.283912428, 243.876592621, 105.123512462, 80.7714792937, 64.4637066722, 46.5737551972, 36.4496404827, 26.8272076175, 18.5252668411, 13.378109403, 17.9482120425, 6.07790538669, 12.966632843, 6.57923346758};

  double y3500[50] = {217.726469308, 279.296559572, 281.775968991, 295.664727964, 343.269526131, 340.044178762, 409.288326569, 451.883740649, 472.456743404, 566.444614425, 554.948992305, 688.53123495, 693.137019448, 703.599098094, 875.459185563, 995.909688801, 1017.83691815, 1007.86736631, 1141.39145879, 1139.5323437, 1195.40881892, 1222.61752961, 1278.21331722, 1218.15810487, 1218.75464276, 1381.01361821, 1456.97796965, 1500.88029896, 1705.45737582, 1890.64534738, 2102.9630746, 2504.71355588, 2670.53920499, 2441.02352202, 1618.28657776, 836.472976632, 377.27449318, 112.07400246, 52.887611635, 30.1016468182, 18.3147185743, 20.8489170745, 2.80900895596, 4.73424816132, 1.02418291569, 0.0, 0.0, 1.44680285454, 1.23776292801, 0.0};
  

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
