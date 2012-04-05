#ifndef _SHAPE_H
#define _SHAPE_H


  double mass;

  double y500[50] = {8.63919753581, 11.3024380207, 15.8584427014, 15.3931927532, 17.0557974428, 14.858171761, 17.9613996372, 21.5163119137, 26.3670294508, 23.7777045444, 26.3572352752, 47.6879247278, 41.4561919123, 49.217860058, 61.2632548138, 91.3929892182, 90.6517746374, 87.9752013907, 138.153524689, 124.196289197, 117.552624218, 167.419196606, 172.642690532, 161.944241747, 202.279491924, 210.527216725, 241.494089544, 293.21346505, 341.832266249, 338.006466977, 377.318334289, 442.313541926, 489.474844709, 470.700277135, 425.668305278, 419.484482966, 339.110264316, 292.332306169, 158.530987613, 163.468288973, 95.3030296713, 61.0215042606, 54.6944052354, 47.2335919216, 28.6930571347, 30.215424329, 23.3655640036, 20.698700048, 15.2426378801, 22.7470517978};

  double y700[50] = {5.30353105068, 12.2813392729, 10.9052353501, 14.3920379207, 10.8687372208, 26.4897435904, 17.2165754512, 22.785209544, 20.843367748, 25.0811594427, 40.6859205887, 38.4258036166, 54.776133962, 57.7176492736, 66.0663515627, 72.8019715548, 84.9957882389, 112.526206881, 110.711182684, 113.925595231, 133.939157799, 158.271149963, 159.5062062, 202.948050253, 213.213083275, 228.221999682, 269.526773848, 293.10071484, 357.181743279, 404.432040185, 417.358013541, 472.322200306, 497.956234403, 510.210262008, 542.64866209, 410.97986199, 346.779346444, 234.192795858, 146.922793701, 90.5631642863, 80.6817563921, 61.3078062758, 42.0409829393, 45.0610121638, 21.3483301774, 20.9388991222, 18.7363274917, 23.666512996, 17.9521396533, 18.5226578936};

  double y1200[50] = {14.207210511, 8.15552312881, 14.2504320145, 14.7009337917, 15.5740203038, 15.5625260472, 18.4270207062, 24.358046338, 33.9819531664, 20.2020184249, 43.4378397167, 42.3228456751, 47.8512445837, 69.3685995787, 82.1399324536, 96.9782082438, 82.8012153953, 94.3188958988, 108.29862804, 126.447411753, 150.927859798, 159.746058494, 210.112534352, 201.112885259, 234.602033161, 213.243600763, 270.985036314, 291.113485664, 378.627058163, 319.586711064, 397.41969502, 454.880339563, 483.578367487, 502.094692156, 511.530486748, 366.040802203, 267.957338125, 160.552777275, 77.474365972, 57.6219993457, 58.7627855986, 30.8686512336, 23.6835200861, 14.7572008073, 18.2461525872, 11.0766398385, 8.6384467259, 9.91433282942, 4.52274981141, 6.52066453546};

  double y2000[50] = {9.78640592098, 9.65582057089, 13.2640735954, 12.4051274136, 16.7760623842, 21.1519175842, 17.1040615365, 26.2191997841, 14.7841992229, 22.1603811458, 32.6788486093, 29.8307448924, 53.9310202822, 66.3737024888, 76.4306181371, 91.4748538882, 106.007001124, 100.385467723, 114.039964564, 153.278117329, 160.529664889, 186.6939255, 195.588241287, 220.711158134, 208.989359252, 277.699347094, 267.648323931, 301.381285705, 338.290989913, 368.71228347, 458.313553549, 468.627392203, 522.980567284, 516.582283422, 473.969481334, 388.802338563, 238.098995477, 115.105356, 57.3208309263, 34.5291408151, 36.0743995756, 11.7084975913, 8.61343426257, 6.05464789271, 11.724133648, 3.80965352058, 8.33860701323, 2.92728567123, 0.0, 0.756138391793};

  double y3500[50] = {44.2140443921, 37.524687551, 43.8335413784, 41.9905421361, 37.4054850712, 48.796015501, 44.3169848621, 42.6579943299, 47.193351537, 73.536796391, 57.8307792991, 78.6060414463, 74.8252939954, 77.5588680729, 86.2953493968, 103.320554852, 117.523922265, 121.119435824, 115.211651981, 130.993466415, 158.854361698, 139.204948321, 188.817327447, 172.676110879, 193.754187562, 218.513780721, 249.009037435, 273.524151325, 270.786866978, 350.494017251, 353.00248801, 395.889403977, 397.03063333, 471.857343473, 434.288255081, 298.682841852, 196.469229616, 65.3868975416, 34.2963696569, 11.4417192116, 4.66061362624, 6.83380997181, 1.77627918124, 0.0, 1.44293212891, 2.88815414906, 1.42406487465, 0.0, 0.0, 0.0};

  
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
