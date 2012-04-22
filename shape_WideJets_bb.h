#ifndef _SHAPE_H
#define _SHAPE_H


  double mass;

  double y500[50] = {8.9629959017, 15.1863757372, 16.4676121548, 16.8919432536, 9.83839890361, 16.5812256485, 11.1234261766, 28.1816038638, 30.4479567036, 24.9473285526, 38.1039339676, 35.9700872451, 38.4450912029, 46.5983887538, 57.8077797517, 58.0990860388, 78.3471239433, 96.2605822682, 78.4604656771, 102.900998041, 120.82212472, 140.167786814, 142.728572167, 171.65849185, 187.428627387, 226.220513269, 278.066167995, 316.861188792, 338.196054481, 365.856857538, 433.775320277, 461.475252271, 449.375938922, 461.21319963, 450.331651114, 387.561248802, 286.096328795, 219.141819932, 151.782094531, 153.382153027, 86.2425823137, 78.9578000456, 63.1831390187, 66.1058953106, 54.0795645118, 42.4938624948, 26.2898144946, 36.782093741, 23.1353883818, 28.7367328703};

  double y700[50] = {13.5842507556, 6.60070860386, 9.51364234835, 15.1638499573, 15.3558292016, 14.9859536514, 17.5406283364, 24.9717584178, 21.9094005004, 31.1813206598, 43.1168310791, 28.619670026, 40.1697228476, 48.8104340509, 57.3820520788, 64.6501423195, 74.0698261932, 95.6070931107, 66.7010582089, 106.79005526, 128.79342106, 133.085696436, 171.014418006, 174.060829043, 171.691687465, 240.145861581, 273.507452443, 305.285815492, 352.032564208, 363.496983647, 440.122034959, 475.072352543, 541.418798819, 504.398861334, 490.376683816, 392.224528395, 272.648879863, 208.509712309, 163.690590709, 107.840759352, 86.1678214297, 89.9285783619, 86.4167347029, 40.5184970573, 48.3115272298, 59.7361123189, 36.234590061, 48.6927760392, 27.1894153059, 25.9378989562};

  double y1200[50] = {7.40415696055, 13.6237627938, 21.3547618985, 3.93783252686, 8.40813491493, 14.8061902672, 18.6486998647, 22.7997717634, 24.6103823557, 28.1040418297, 26.8368044272, 34.354634881, 42.2222889736, 54.1759063229, 50.2696295753, 63.6785546616, 60.4427244216, 79.2931593433, 90.4447763339, 105.208130024, 96.8683386222, 129.427052639, 154.616738573, 198.001427285, 186.570745215, 227.617765747, 235.197378457, 274.81920062, 367.533909068, 337.053860277, 406.701433077, 481.082525857, 518.345989525, 553.85468404, 487.292065375, 405.047087282, 273.25550089, 173.583573952, 95.724605374, 85.5588411912, 78.2041481212, 56.4919239432, 44.6363881677, 35.6263326406, 30.2225242183, 19.7913957387, 27.7121640071, 15.886581026, 18.2824142501, 13.9796235189};

  double y2000[50] = {12.3564771637, 7.61293031275, 11.4465908185, 13.6975083724, 7.78323234618, 24.0764616877, 19.3026077598, 15.6831256598, 21.1944480464, 20.8450015858, 17.886527881, 27.5248542204, 35.8881525546, 41.7146070302, 66.4188623801, 57.9036692306, 73.6079054475, 67.0823541656, 102.684784293, 103.557308845, 118.150619961, 132.432967693, 166.275766775, 176.828221351, 194.344825447, 245.70009511, 236.751501337, 301.834296584, 300.822007976, 376.608742446, 454.405628867, 503.249250568, 576.242831126, 556.758251339, 561.775507443, 425.648019038, 286.329419829, 148.185408421, 88.3227125034, 63.7023587376, 52.8501829281, 29.0822808295, 29.309288986, 15.6611463428, 14.6866334155, 10.9531324953, 6.94325608015, 9.03804254532, 7.15610802174, 4.97130734473};

  double y3500[50] = {29.0780562088, 36.663995564, 43.8590622321, 49.8205624819, 35.1285423785, 37.7298346311, 45.2202432081, 38.5342534631, 47.5543106943, 46.6267103925, 70.1802324876, 79.3539003879, 61.7172632813, 68.2611512616, 68.3795707375, 84.7630704567, 73.2316328958, 102.346149355, 93.2704570815, 128.502058558, 119.273966081, 126.162745699, 177.18279551, 142.955737859, 171.58130423, 203.147706695, 221.921080522, 246.82237272, 286.935500532, 302.450479209, 364.814042158, 408.612908632, 455.934937552, 562.366674766, 476.981642127, 380.366965123, 204.650544383, 85.7731666043, 64.6106219366, 27.5580112934, 17.1640530974, 14.2591787055, 9.29173914343, 1.46489519626, 1.44293212891, 4.52877554297, 1.47236061096, 2.89642548561, 1.44293212891, 1.42406487465};

  
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
