void SetVDWRADIUS(double *vdw) 
{ 
      double vdw_radii_file[104]; 
      int i;
      vdw_radii_file[0]  = 0.0;
      vdw_radii_file[1] = 2.72687;
      vdw_radii_file[2] = 2.23177;
      vdw_radii_file[3] = 2.31586;
      vdw_radii_file[4] = 2.59365;
      vdw_radii_file[5] = 3.85788;
      vdw_radii_file[6] = 3.63867;
      vdw_radii_file[7] = 3.4582;
      vdw_radii_file[8] = 3.30702;
      vdw_radii_file[9] = 3.17852;
      vdw_radii_file[10] = 3.06419;
      vdw_radii_file[11] = 2.81853;
      vdw_radii_file[12] = 2.85443;
      vdw_radii_file[13] = 4.25094;
      vdw_radii_file[14] = 4.05819;
      vdw_radii_file[15] = 3.91835;
      vdw_radii_file[16] = 3.81252;
      vdw_radii_file[17] = 3.72937;
      vdw_radii_file[18] = 3.65473;
      vdw_radii_file[19] = 3.60182;
      vdw_radii_file[20] = 3.21159;
      vdw_radii_file[21] = 3.11332;
      vdw_radii_file[22] = 2.99994;
      vdw_radii_file[23] = 2.97065;
      vdw_radii_file[24] = 2.85632;
      vdw_radii_file[25] = 2.79774;
      vdw_radii_file[26] = 2.75144;
      vdw_radii_file[27] = 2.71365;
      vdw_radii_file[28] = 2.67774;
      vdw_radii_file[29] = 3.3023;
      vdw_radii_file[30] = 2.61066;
      vdw_radii_file[31] = 4.14133;
      vdw_radii_file[32] = 4.04401;
      vdw_radii_file[33] = 3.99677;
      vdw_radii_file[34] = 3.97315;
      vdw_radii_file[35] = 3.95803;
      vdw_radii_file[36] = 3.91268;
      vdw_radii_file[37] = 3.88717;
      vdw_radii_file[38] = 3.44025;
      vdw_radii_file[39] = 3.16057;
      vdw_radii_file[40] = 2.95175;
      vdw_radii_file[41] = 2.99049;
      vdw_radii_file[42] = 2.88372;
      vdw_radii_file[43] = 2.8327;
      vdw_radii_file[44] = 2.79963;
      vdw_radii_file[45] = 2.7675;
      vdw_radii_file[46] = 2.73916;
      vdw_radii_file[47] = 2.97443;
      vdw_radii_file[48] = 2.69097;
      vdw_radii_file[49] = 4.21692;
      vdw_radii_file[50] = 4.14984;
      vdw_radii_file[51] = 4.17629;
      vdw_radii_file[52] = 4.22354;
      vdw_radii_file[53] = 4.25188;
      vdw_radii_file[54] = 4.16118;
      vdw_radii_file[55] = 4.26795;
      vdw_radii_file[56] = 3.49883;
      vdw_radii_file[57] = 3.32781;
      vdw_radii_file[58] = 3.35993;
      vdw_radii_file[59] = 3.40718;
      vdw_radii_file[60] = 3.37789;
      vdw_radii_file[61] = 3.35143;
      vdw_radii_file[62] = 3.32592;
      vdw_radii_file[63] = 3.30041;
      vdw_radii_file[64] = 3.1823;
      vdw_radii_file[65] = 3.26072;
      vdw_radii_file[66] = 3.23899;
      vdw_radii_file[67] = 3.22104;
      vdw_radii_file[68] = 3.20403;
      vdw_radii_file[69] = 3.18797;
      vdw_radii_file[70] = 3.17002;
      vdw_radii_file[71] = 3.4393;
      vdw_radii_file[72] = 2.96781;
      vdw_radii_file[73] = 2.99522;
      vdw_radii_file[74] = 2.89978;
      vdw_radii_file[75] = 2.79113;
      vdw_radii_file[76] = 2.94797;
      vdw_radii_file[77] = 2.68341;
      vdw_radii_file[78] = 2.60215;
      vdw_radii_file[79] = 3.11143;
      vdw_radii_file[80] = 2.55585;
      vdw_radii_file[81] = 4.10732;
      vdw_radii_file[82] = 4.06008;
      vdw_radii_file[83] = 4.12905;
      vdw_radii_file[84] = 4.44936;
      vdw_radii_file[85] = 4.4881;
      vdw_radii_file[86] = 4.50227;
      vdw_radii_file[87] = 4.62983;
      vdw_radii_file[88] = 3.47426;
      vdw_radii_file[89] = 3.28623;
      vdw_radii_file[90] = 3.20875;
      vdw_radii_file[91] = 3.23521;
      vdw_radii_file[92] = 3.20781;
      vdw_radii_file[93] = 3.23521;
      vdw_radii_file[94] = 3.23521;
      vdw_radii_file[95] = 3.19458;
      vdw_radii_file[96] = 3.14261;
      vdw_radii_file[97] = 3.1549;
      vdw_radii_file[98] = 3.13033;
      vdw_radii_file[99] = 3.1171;
      vdw_radii_file[100] = 3.10482;
      vdw_radii_file[101] = 3.09348;
      vdw_radii_file[102] = 3.06892;
      vdw_radii_file[103] = 3.05758;
      for(i = 1; i <= 103; i++) vdw[i-1] = vdw_radii_file[i];   
      return; 
}     