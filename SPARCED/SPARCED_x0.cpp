#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "SPARCED_k.h"

namespace amici {
namespace model_SPARCED {

void x0_SPARCED(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1897.3920000000001;
    x0[1] = 296.6386;
    x0[2] = 6.2455280000000002;
    x0[3] = 205.63849999999999;
    x0[4] = 2.2307359999999998;
    x0[7] = 6.2455280000000002;
    x0[8] = 6.2455280000000002;
    x0[9] = 6.2455280000000002;
    x0[10] = 6.2455280000000002;
    x0[11] = 6.2455280000000002;
    x0[12] = 6.2455280000000002;
    x0[13] = 6.2455280000000002;
    x0[14] = 6.2455280000000002;
    x0[15] = 6.2455280000000002;
    x0[16] = 6.2455280000000002;
    x0[17] = 6.2455280000000002;
    x0[18] = 6.2455280000000002;
    x0[19] = 6.2455280000000002;
    x0[20] = 6.2455280000000002;
    x0[21] = 6.2455280000000002;
    x0[22] = 6.2455280000000002;
    x0[23] = 6.2455280000000002;
    x0[24] = 6.2455280000000002;
    x0[25] = 6.2455280000000002;
    x0[26] = 6.2455280000000002;
    x0[27] = 0.82180839999999999;
    x0[28] = 45.50432;
    x0[29] = 120.2569;
    x0[30] = 0.0002431635;
    x0[31] = 6.0277599999999999e-6;
    x0[35] = 0.027216089999999998;
    x0[37] = 1.503781;
    x0[38] = 2.4003570000000001;
    x0[39] = 7996.759;
    x0[40] = 0.5355972;
    x0[41] = 0.0069708620000000004;
    x0[42] = 1119.2159999999999;
    x0[43] = 0.23333590000000001;
    x0[44] = 0.58917489999999995;
    x0[45] = 3.8541919999999998;
    x0[46] = 0.7878773;
    x0[47] = 0.25528509999999999;
    x0[48] = 0.88940549999999996;
    x0[49] = 2.746016;
    x0[50] = 11.40893;
    x0[51] = 1.1226860000000001;
    x0[52] = 141.28649999999999;
    x0[53] = 49.187809999999999;
    x0[54] = 16.023430000000001;
    x0[55] = 55.298499999999997;
    x0[56] = 4.0293700000000001;
    x0[57] = 14.0382;
    x0[58] = 5.2308079999999997;
    x0[59] = 1.2831440000000001;
    x0[60] = 20.44209;
    x0[61] = 0.5978504;
    x0[62] = 108.8043;
    x0[63] = 11.601430000000001;
    x0[64] = 4.2046130000000002;
    x0[65] = 41.624369999999999;
    x0[66] = 3.4434999999999998;
    x0[67] = 3.9009480000000001;
    x0[68] = 0.68424439999999997;
    x0[69] = 71.167820000000006;
    x0[70] = 0.73011530000000002;
    x0[71] = 9.4029919999999994;
    x0[72] = 26.51906;
    x0[73] = 16.740469999999998;
    x0[74] = 0.083924310000000002;
    x0[75] = 0.1497549;
    x0[76] = 8950.1010000000006;
    x0[77] = 0.2997245;
    x0[78] = 20.871130000000001;
    x0[79] = 0.48264800000000002;
    x0[80] = 0.2292498;
    x0[81] = 0.2620149;
    x0[82] = 0.017137139999999999;
    x0[83] = 2.2099300000000001e-17;
    x0[84] = 8.0312350000000006;
    x0[85] = 2.2366280000000001e-16;
    x0[86] = 7.563593e-20;
    x0[87] = 0.52592799999999995;
    x0[88] = 1.054758e-19;
    x0[89] = 26.546610000000001;
    x0[90] = 6.3426979999999997e-21;
    x0[91] = 0.00036742090000000001;
    x0[92] = 0.84152179999999999;
    x0[93] = 0.00095053659999999997;
    x0[94] = 80.800669999999997;
    x0[95] = 0.0046882759999999999;
    x0[96] = 0.00013715259999999999;
    x0[97] = 47.695749999999997;
    x0[98] = 2.0664309999999999e-8;
    x0[99] = 0.00038833779999999999;
    x0[100] = 9.769605e-8;
    x0[101] = 5.477824;
    x0[102] = 4.702948e-5;
    x0[103] = 640.14700000000005;
    x0[104] = 2.7487399999999999e-6;
    x0[105] = 1.097548;
    x0[106] = 50.476860000000002;
    x0[107] = 5.3176369999999999e-5;
    x0[108] = 1.1610359999999999e-6;
    x0[109] = 12.5002;
    x0[110] = 4.0903119999999999e-5;
    x0[111] = 88.260329999999996;
    x0[112] = 3.23704e-8;
    x0[113] = 5.753417e-5;
    x0[114] = 0.00077206970000000003;
    x0[115] = 1.783563;
    x0[116] = 0.0038809560000000001;
    x0[117] = 1.8847699999999999e-5;
    x0[118] = 9.4741580000000006e-5;
    x0[119] = 1.1232129999999999e-8;
    x0[120] = 5.6460459999999998e-8;
    x0[121] = 2259.9989999999998;
    x0[122] = 8.0187850000000002e-8;
    x0[123] = 8.0187839999999998e-6;
    x0[124] = 15543.35;
    x0[125] = 0.00071657670000000004;
    x0[126] = 0.0071638719999999999;
    x0[127] = 1474.624;
    x0[128] = 6.7982860000000007e-5;
    x0[129] = 0.00054288779999999998;
    x0[130] = 0.00050140459999999999;
    x0[131] = 0.52001459999999999;
    x0[132] = 4.1182180000000001e-7;
    x0[133] = 0.001145514;
    x0[134] = 4.5279230000000004;
    x0[135] = 0.00041517330000000001;
    x0[136] = 0.00047952100000000001;
    x0[137] = 2.842428e-5;
    x0[138] = 0.013983219999999999;
    x0[139] = 0.0033506970000000001;
    x0[140] = 2.0201730000000002;
    x0[141] = 0.094630580000000006;
    x0[142] = 0.029054199999999999;
    x0[143] = 0.005548758;
    x0[144] = 3.3338209999999999;
    x0[145] = 1.0235749999999999;
    x0[146] = 0.1943735;
    x0[147] = 0.060049329999999998;
    x0[148] = 1.67411e-8;
    x0[149] = 2.115529;
    x0[150] = 0.0006063485;
    x0[151] = 0.10747429999999999;
    x0[152] = 3.3490030000000002e-5;
    x0[153] = 0.90534539999999997;
    x0[154] = 0.00095558060000000002;
    x0[155] = 9.8216909999999998e-17;
    x0[156] = 9.3593389999999997e-15;
    x0[157] = 2.4719279999999999e-15;
    x0[158] = 2.8596060000000001e-19;
    x0[159] = 4.1326780000000001e-16;
    x0[160] = 1.2566920000000001e-18;
    x0[161] = 8.9529110000000004e-15;
    x0[162] = 50.034570000000002;
    x0[163] = 5.051704;
    x0[164] = 3.7981889999999998;
    x0[165] = 0.38350119999999999;
    x0[166] = 1.739962;
    x0[167] = 8.2153650000000004e-18;
    x0[168] = 8.2959949999999996e-19;
    x0[169] = 4.6520009999999997e-11;
    x0[170] = 8.5733449999999998;
    x0[171] = 0.2296793;
    x0[172] = 0.18671570000000001;
    x0[173] = 0.02181495;
    x0[174] = 0.029796429999999999;
    x0[175] = 2.5034200000000002;
    x0[176] = 0.19003790000000001;
    x0[177] = 0.087056960000000003;
    x0[178] = 4.1104640000000002e-19;
    x0[179] = 0.014426090000000001;
    x0[180] = 0.0066086239999999996;
    x0[181] = 3.1203200000000001e-20;
    x0[182] = 1.429426e-20;
    x0[184] = 0.073501769999999994;
    x0[185] = 3.4862510000000001e-5;
    x0[186] = 0.72690980000000005;
    x0[187] = 0.30894440000000001;
    x0[188] = 1.835326e-17;
    x0[189] = 1.5671890000000001e-16;
    x0[190] = 1.749132e-28;
    x0[191] = 7.6561050000000004e-14;
    x0[192] = 6.5595120000000004e-20;
    x0[193] = 3.0151970000000002e-16;
    x0[194] = 2.001079e-16;
    x0[195] = 4.8481019999999997e-27;
    x0[196] = 5.0127749999999999e-16;
    x0[197] = 9.0915519999999999e-17;
    x0[198] = 2.6910819999999999e-27;
    x0[199] = 5.3963790000000001e-16;
    x0[200] = 6.2846840000000003e-15;
    x0[201] = 2.7234920000000002e-16;
    x0[202] = 7.1921509999999998e-27;
    x0[203] = 4.6761090000000003e-27;
    x0[204] = 6.4426780000000001e-28;
    x0[205] = 8.0920089999999997e-39;
    x0[206] = 7.1206799999999994e-27;
    x0[207] = 1.973005e-27;
    x0[208] = 1.6563599999999999e-39;
    x0[209] = 5.6272080000000001e-13;
    x0[210] = 4.5255890000000001e-24;
    x0[211] = 4.9382530000000004e-16;
    x0[212] = 5.1490360000000001e-19;
    x0[213] = 6.8660789999999996e-17;
    x0[215] = 4.7361539999999997e-32;
    x0[216] = 8.5715829999999998e-43;
    x0[217] = 9.3115590000000004e-32;
    x0[218] = 2.3505299999999999e-42;
    x0[220] = 1.0254130000000001e-26;
    x0[222] = 8.1137640000000001e-30;
    x0[223] = 7.6065730000000003e-38;
    x0[224] = 1.0170390000000001e-27;
    x0[225] = 0.505247;
    x0[226] = 0.038354109999999997;
    x0[227] = 8.2959340000000004e-20;
    x0[228] = 1.9918730000000001e-16;
    x0[229] = 4.8440140000000001e-27;
    x0[230] = 4.9505909999999997e-16;
    x0[232] = 9.0494540000000001e-17;
    x0[233] = 4.7142239999999999e-32;
    x0[234] = 2.7120489999999999e-27;
    x0[235] = 8.5503610000000008e-43;
    x0[236] = 5.392754e-16;
    x0[237] = 7.1873800000000004e-27;
    x0[238] = 6.255583e-15;
    x0[239] = 4.7203380000000002e-27;
    x0[240] = 2.3572630000000001e-42;
    x0[241] = 6.4437450000000004e-28;
    x0[242] = 8.087255e-39;
    x0[243] = 7.0899959999999995e-27;
    x0[244] = 1.9906209999999998e-27;
    x0[245] = 1.6737189999999999e-39;
    x0[247] = 5.6236669999999997e-13;
    x0[248] = 1.023324e-26;
    x0[250] = 4.5249230000000003e-24;
    x0[251] = 8.1380809999999995e-30;
    x0[252] = 4.9350309999999997e-16;
    x0[253] = 5.1460629999999999e-19;
    x0[254] = 6.8599729999999995e-17;
    x0[255] = 7.606579e-38;
    x0[256] = 1.016605e-27;
    x0[257] = 3.751139e-19;
    x0[258] = 4.9603169999999997e-17;
    x0[259] = 5.8846069999999997e-38;
    x0[260] = 7.7438050000000004e-28;
    x0[261] = 3.301344e-18;
    x0[262] = 1.1917929999999999e-28;
    x0[263] = 8.2699860000000006e-18;
    x0[265] = 1.499908e-18;
    x0[266] = 7.8136220000000005e-34;
    x0[267] = 4.452381e-29;
    x0[268] = 1.424342e-44;
    x0[269] = 8.9028090000000003e-18;
    x0[270] = 1.1865429999999999e-28;
    x0[271] = 1.036831e-16;
    x0[272] = 7.7001619999999995e-29;
    x0[273] = 3.8678260000000002e-44;
    x0[274] = 1.059008e-29;
    x0[275] = 1.3299419999999999e-40;
    x0[276] = 1.1703610000000001e-28;
    x0[277] = 3.2482060000000001e-29;
    x0[278] = 2.7311169999999998e-41;
    x0[280] = 9.2836429999999998e-15;
    x0[281] = 1.7483040000000001e-28;
    x0[283] = 8.9105440000000003e-26;
    x0[284] = 2.2704410000000002e-31;
    x0[285] = 8.1470039999999995e-18;
    x0[286] = 8.4947369999999993e-21;
    x0[287] = 1.132748e-18;
    x0[288] = 1.4799230000000001e-39;
    x0[289] = 1.8432749999999999e-29;
    x0[290] = 7.1412450000000003e-16;
    x0[291] = 4.8721860000000002e-26;
    x0[292] = 2.4682290000000001e-15;
    x0[294] = 3.2444700000000002e-16;
    x0[295] = 1.6901749999999999e-31;
    x0[296] = 1.017803e-26;
    x0[297] = 3.1072159999999999e-42;
    x0[298] = 2.08357e-16;
    x0[299] = 2.7769309999999999e-27;
    x0[300] = 2.2427760000000001e-14;
    x0[301] = 1.861872e-27;
    x0[302] = 9.1536030000000004e-43;
    x0[303] = 2.4796199999999998e-28;
    x0[304] = 3.1077699999999997e-39;
    x0[305] = 2.5401929999999999e-26;
    x0[306] = 7.6796660000000002e-28;
    x0[307] = 6.4449959999999999e-40;
    x0[309] = 2.1728669999999999e-13;
    x0[310] = 1.9845970000000001e-26;
    x0[312] = 3.248321e-24;
    x0[313] = 6.3355100000000001e-30;
    x0[314] = 1.9067090000000001e-16;
    x0[315] = 8.0908680000000002e-19;
    x0[316] = 1.5405820000000001e-16;
    x0[317] = 3.50231e-36;
    x0[318] = 3.3028869999999999e-25;
    x0[319] = 7.1488509999999996e-16;
    x0[320] = 4.8752370000000002e-26;
    x0[321] = 2.472361e-15;
    x0[323] = 3.2479249999999999e-16;
    x0[324] = 1.6919750000000001e-31;
    x0[325] = 1.0188230000000001e-26;
    x0[326] = 3.1104990000000001e-42;
    x0[327] = 2.0866689999999999e-16;
    x0[328] = 2.7810600000000001e-27;
    x0[329] = 2.245164e-14;
    x0[330] = 1.8646450000000001e-27;
    x0[331] = 9.1672420000000006e-43;
    x0[332] = 2.4833169999999999e-28;
    x0[333] = 3.1124050000000001e-39;
    x0[334] = 2.5429090000000001e-26;
    x0[335] = 7.6911059999999997e-28;
    x0[336] = 6.4545859999999998e-40;
    x0[338] = 2.1760990000000001e-13;
    x0[339] = 1.9875630000000001e-26;
    x0[341] = 3.2519679999999999e-24;
    x0[342] = 6.3425369999999995e-30;
    x0[343] = 1.909544e-16;
    x0[344] = 8.0915789999999996e-19;
    x0[345] = 1.540766e-16;
    x0[346] = 3.5023199999999998e-36;
    x0[347] = 3.3029430000000002e-25;
    x0[348] = 6.0032509999999995e-19;
    x0[349] = 1.1431129999999999e-16;
    x0[350] = 2.6161840000000001e-36;
    x0[351] = 2.4717060000000001e-25;
    x0[352] = 1.9171270000000001e-19;
    x0[353] = 1.0443319999999999e-29;
    x0[354] = 5.2942420000000004e-19;
    x0[356] = 7.7421049999999999e-20;
    x0[357] = 4.0331740000000002e-35;
    x0[358] = 3.2250999999999999e-30;
    x0[359] = 1.188276e-45;
    x0[360] = 4.0370089999999998e-19;
    x0[361] = 6.1491000000000001e-30;
    x0[362] = 5.3518570000000002e-18;
    x0[363] = 5.7093869999999999e-30;
    x0[364] = 2.5160389999999999e-45;
    x0[365] = 6.8881279999999996e-31;
    x0[366] = 9.5289749999999999e-42;
    x0[367] = 8.3531650000000004e-30;
    x0[368] = 2.6474409999999999e-30;
    x0[369] = 2.154496e-42;
    x0[371] = 1.202825e-16;
    x0[372] = 2.0966600000000001e-30;
    x0[373] = 2.6900650000000002e-34;
    x0[374] = 4.8414519999999999e-28;
    x0[375] = 8.7110060000000005e-33;
    x0[376] = 4.2221400000000001e-19;
    x0[377] = 2.0054640000000001e-22;
    x0[378] = 3.7126929999999998e-20;
    x0[379] = 3.4915419999999999e-41;
    x0[380] = 6.1698600000000001e-31;
    x0[381] = 6.8802970000000002e-19;
    x0[382] = 5.1094170000000004e-29;
    x0[383] = 2.643136e-18;
    x0[385] = 2.7785900000000002e-19;
    x0[386] = 1.4474790000000001e-34;
    x0[387] = 1.171755e-29;
    x0[388] = 3.6593180000000001e-45;
    x0[389] = 1.562084e-19;
    x0[390] = 2.3793220000000001e-30;
    x0[391] = 1.9207309999999999e-17;
    x0[392] = 1.8209479999999998e-30;
    x0[393] = 8.8718e-46;
    x0[394] = 2.6574019999999999e-31;
    x0[395] = 3.6612920000000001e-42;
    x0[396] = 2.998567e-29;
    x0[397] = 7.4112799999999999e-31;
    x0[398] = 8.2832510000000002e-43;
    x0[400] = 4.6543749999999997e-17;
    x0[401] = 1.398964e-29;
    x0[403] = 3.4496420000000002e-28;
    x0[404] = 5.5154690000000001e-33;
    x0[405] = 1.633701e-19;
    x0[406] = 3.2100319999999999e-22;
    x0[407] = 8.5573609999999998e-20;
    x0[408] = 1.494436e-39;
    x0[409] = 1.832008e-28;
    x0[410] = 3.4579240000000001e-18;
    x0[411] = 8.5126420000000002e-29;
    x0[412] = 8.5943079999999997e-18;
    x0[414] = 1.1782500000000001e-18;
    x0[415] = 6.1379780000000003e-34;
    x0[416] = 6.0363309999999998e-29;
    x0[417] = 1.2898819999999999e-44;
    x0[418] = 7.0218239999999997e-18;
    x0[419] = 9.3585790000000002e-29;
    x0[420] = 8.1448460000000003e-17;
    x0[421] = 4.839926e-29;
    x0[422] = 2.2827170000000001e-44;
    x0[423] = 8.3778279999999997e-30;
    x0[424] = 1.055775e-40;
    x0[425] = 9.6396720000000005e-29;
    x0[426] = 2.2687149999999999e-29;
    x0[427] = 1.49543e-41;
    x0[429] = 9.7644399999999998e-15;
    x0[430] = 1.707595e-28;
    x0[431] = 4.3063720000000001e-32;
    x0[432] = 9.8488699999999995e-26;
    x0[433] = 1.5502779999999999e-31;
    x0[434] = 1.2853130000000001e-17;
    x0[435] = 2.4399390000000001e-21;
    x0[436] = 7.5283880000000002e-19;
    x0[437] = 3.8320920000000003e-40;
    x0[438] = 1.1809609999999999e-29;
    x0[439] = 2.282943e-17;
    x0[440] = 2.9688629999999999e-28;
    x0[441] = 2.8370079999999997e-17;
    x0[443] = 2.5929600000000001e-17;
    x0[444] = 1.3507770000000001e-32;
    x0[445] = 6.0546250000000003e-28;
    x0[446] = 1.043771e-43;
    x0[447] = 1.8550480000000001e-16;
    x0[448] = 2.0603159999999998e-27;
    x0[449] = 1.7924249999999999e-15;
    x0[450] = 2.1827689999999999e-27;
    x0[451] = 8.0341290000000003e-43;
    x0[452] = 1.1843479999999999e-28;
    x0[453] = 9.2972529999999998e-40;
    x0[454] = 8.1257160000000009e-28;
    x0[455] = 9.1695909999999993e-28;
    x0[456] = 2.956571e-40;
    x0[458] = 9.6811029999999997e-14;
    x0[459] = 1.688809e-27;
    x0[460] = 2.8590709999999999e-31;
    x0[461] = 7.9692779999999996e-25;
    x0[462] = 3.4219879999999998e-30;
    x0[463] = 2.2654979999999999e-16;
    x0[464] = 1.121071e-19;
    x0[465] = 2.0471890000000001e-17;
    x0[466] = 1.7945149999999999e-38;
    x0[467] = 3.2364109999999999e-28;
    x0[468] = 1.9610130000000001e-17;
    x0[469] = 2.402515e-28;
    x0[470] = 2.4369460000000001e-17;
    x0[472] = 2.227313e-17;
    x0[473] = 1.160297e-32;
    x0[474] = 2.670145e-28;
    x0[475] = 8.4188699999999996e-44;
    x0[476] = 1.592826e-16;
    x0[477] = 1.7690780000000001e-27;
    x0[478] = 1.539666e-15;
    x0[479] = 1.415068e-27;
    x0[480] = 7.9620030000000003e-43;
    x0[481] = 9.5159669999999996e-29;
    x0[482] = 7.9621859999999996e-40;
    x0[483] = 6.9797410000000004e-28;
    x0[484] = 5.6001480000000004e-28;
    x0[485] = 2.4860960000000002e-40;
    x0[487] = 8.3055409999999995e-14;
    x0[488] = 1.4262689999999999e-27;
    x0[489] = 1.7420069999999999e-31;
    x0[490] = 6.697302e-25;
    x0[491] = 2.3047679999999998e-30;
    x0[492] = 1.9435999999999999e-16;
    x0[493] = 9.6882330000000003e-20;
    x0[494] = 1.76917e-17;
    x0[495] = 1.520179e-38;
    x0[496] = 2.7618489999999999e-28;
    x0[497] = 1.0986579999999999e-19;
    x0[498] = 5.459551e-30;
    x0[499] = 4.218002e-19;
    x0[501] = 4.4369030000000002e-20;
    x0[502] = 2.3113609999999999e-35;
    x0[503] = 1.182781e-30;
    x0[504] = 5.8535620000000001e-46;
    x0[505] = 2.4946750000000001e-20;
    x0[506] = 3.799818e-31;
    x0[507] = 3.0670580000000001e-18;
    x0[508] = 3.5838770000000001e-31;
    x0[509] = 1.5415669999999999e-46;
    x0[510] = 4.286908e-32;
    x0[511] = 5.845726e-43;
    x0[512] = 4.9768069999999999e-30;
    x0[513] = 1.205966e-31;
    x0[514] = 1.3228180000000001e-43;
    x0[516] = 7.4331170000000003e-18;
    x0[517] = 2.6251269999999998e-31;
    x0[518] = 3.4365210000000001e-31;
    x0[519] = 5.4652889999999997e-29;
    x0[520] = 1.1082870000000001e-33;
    x0[521] = 2.609048e-20;
    x0[522] = 5.1264769999999998e-23;
    x0[523] = 1.3666250000000001e-20;
    x0[524] = 2.380338e-40;
    x0[525] = 2.4518450000000001e-29;
    x0[526] = 1.546283e-18;
    x0[527] = 8.4915849999999995e-29;
    x0[528] = 4.270138e-18;
    x0[530] = 6.2444929999999997e-19;
    x0[531] = 3.2530070000000002e-34;
    x0[532] = 2.6016969999999999e-29;
    x0[533] = 9.607913e-45;
    x0[534] = 3.2561009999999999e-18;
    x0[535] = 4.959636e-29;
    x0[536] = 4.3166079999999998e-17;
    x0[537] = 4.6201400000000001e-29;
    x0[538] = 2.032558e-44;
    x0[539] = 5.5556689999999997e-30;
    x0[540] = 7.6859169999999999e-41;
    x0[541] = 6.7375180000000001e-29;
    x0[542] = 2.1448369999999999e-29;
    x0[543] = 1.7378219999999999e-41;
    x0[545] = 9.7015369999999992e-16;
    x0[546] = 1.6898899999999999e-29;
    x0[547] = 2.24406e-33;
    x0[548] = 3.9049620000000002e-27;
    x0[549] = 7.0486790000000004e-32;
    x0[550] = 3.4054210000000001e-18;
    x0[551] = 1.617529e-21;
    x0[552] = 2.9945139999999998e-19;
    x0[553] = 2.8041539999999999e-40;
    x0[554] = 4.9812290000000001e-30;
    x0[555] = 3.1429799999999999e-17;
    x0[556] = 7.7415140000000001e-28;
    x0[557] = 7.8115480000000001e-17;
    x0[559] = 1.070936e-17;
    x0[560] = 5.5789380000000001e-33;
    x0[561] = 5.6027350000000004e-28;
    x0[562] = 1.1805909999999999e-43;
    x0[563] = 6.3823009999999995e-17;
    x0[564] = 8.5062319999999991e-28;
    x0[565] = 7.403023e-16;
    x0[566] = 4.433553e-28;
    x0[567] = 2.0857820000000001e-43;
    x0[568] = 7.6142269999999997e-29;
    x0[569] = 9.5974539999999994e-40;
    x0[570] = 8.7806459999999999e-28;
    x0[571] = 2.0871660000000002e-28;
    x0[572] = 1.361205e-40;
    x0[574] = 8.8751750000000001e-14;
    x0[575] = 1.5521369999999999e-27;
    x0[576] = 3.9174280000000001e-31;
    x0[577] = 8.9530590000000007e-25;
    x0[578] = 1.448168e-30;
    x0[579] = 1.168257e-16;
    x0[580] = 2.217625e-20;
    x0[581] = 6.842444e-18;
    x0[582] = 3.48335e-39;
    x0[583] = 1.0736209999999999e-28;
    x0[584] = 2.2599149999999998e-16;
    x0[585] = 2.9400450000000001e-27;
    x0[586] = 2.8083899999999998e-16;
    x0[588] = 2.5668050000000002e-16;
    x0[589] = 1.337152e-31;
    x0[590] = 6.0140199999999999e-27;
    x0[591] = 1.033682e-42;
    x0[592] = 1.836341e-15;
    x0[593] = 2.0395389999999999e-26;
    x0[594] = 1.7743440000000001e-14;
    x0[595] = 2.16464e-26;
    x0[596] = 7.9525939999999998e-42;
    x0[597] = 1.172933e-27;
    x0[598] = 9.2036610000000003e-39;
    x0[599] = 8.0437520000000003e-27;
    x0[600] = 9.0932489999999994e-27;
    x0[601] = 2.927296e-39;
    x0[603] = 9.5835299999999997e-13;
    x0[604] = 1.6718530000000001e-26;
    x0[605] = 2.8336290000000001e-30;
    x0[606] = 7.8900919999999994e-24;
    x0[607] = 3.3883080000000001e-29;
    x0[608] = 2.2426649999999999e-15;
    x0[609] = 1.1097660000000001e-18;
    x0[610] = 2.0265460000000001e-16;
    x0[611] = 1.776466e-37;
    x0[612] = 3.2038329999999999e-27;
    x0[613] = 1.4048600000000001e-44;
    x0[615] = 1.7458150000000001e-44;
    x0[616] = 4.5512739999999998e-55;
    x0[617] = 1.5956360000000001e-44;
    x0[618] = 8.2546380000000002e-60;
    x0[621] = 4.3957959999999999e-32;
    x0[622] = 1.624054e-54;
    x0[623] = 1.413445e-42;
    x0[624] = 6.8614000000000002e-54;
    x0[632] = 6.8349050000000004e-41;
    x0[633] = 7.5897460000000002e-54;
    x0[634] = 9.1436910000000001e-57;
    x0[635] = 9.3754219999999999e-51;
    x0[636] = 1.1155350000000001e-54;
    x0[637] = 2.4815500000000001e-42;
    x0[638] = 2.3418069999999999e-32;
    x0[639] = 5.6069079999999996e-31;
    x0[642] = 0.9354751;
    x0[643] = 0.01808599;
    x0[644] = 5.2994450000000004;
    x0[645] = 60.845419999999997;
    x0[646] = 0.2780629;
    x0[647] = 0.2804741;
    x0[648] = 2.2140840000000002;
    x0[649] = 0.023694280000000002;
    x0[650] = 8.9622150000000005;
    x0[651] = 2.0104549999999999;
    x0[652] = 0.33867150000000001;
    x0[653] = 1.72326;
    x0[654] = 2.742829;
    x0[655] = 0.18619630000000001;
    x0[656] = 9.3214210000000003e-25;
    x0[657] = 3.0692659999999998e-27;
    x0[658] = 8.3132619999999992e-12;
    x0[659] = 9.6769239999999996;
    x0[660] = 8.0446799999999996e-13;
    x0[661] = 2.8007040000000001;
    x0[662] = 161.31270000000001;
    x0[663] = 0.7812905;
    x0[665] = 3.269177;
    x0[666] = 3.2405629999999999;
    x0[667] = 0.90510190000000001;
    x0[668] = 0.41723700000000002;
    x0[669] = 0.1053176;
    x0[670] = 125.1738;
    x0[671] = 16.265820000000001;
    x0[672] = 21.132400000000001;
    x0[673] = 0.0019851539999999998;
    x0[674] = 0.40725460000000002;
    x0[675] = 7.3282290000000003;
    x0[676] = 69.527760000000001;
    x0[677] = 0.70209129999999997;
    x0[678] = 1.781693;
    x0[679] = 0.0020186779999999999;
    x0[680] = 0.015870550000000001;
    x0[681] = 2.9838840000000002;
    x0[682] = 0.2469073;
    x0[683] = 2.871016;
    x0[684] = 0.070547899999999997;
    x0[685] = 45.439999999999998;
    x0[686] = 35.489719999999998;
    x0[687] = 17.758019999999998;
    x0[688] = 179.8177;
    x0[689] = 1.10864e-12;
    x0[690] = 9998.2700000000004;
    x0[691] = 1.6072439999999999;
    x0[692] = 2.3156970000000001;
    x0[693] = 0.73510200000000003;
    x0[694] = 48.16328;
    x0[695] = 0.038927400000000001;
    x0[696] = 0.38929639999999999;
    x0[697] = 11.50741;
    x0[698] = 0.1849392;
    x0[699] = 0.063594929999999994;
    x0[700] = 0.16712009999999999;
    x0[701] = 1.18038;
    x0[702] = 1.839429;
    x0[703] = 0.0039578649999999996;
    x0[704] = 40.765639999999998;
    x0[705] = 1.4418679999999999;
    x0[706] = 0.056440299999999999;
    x0[707] = 0.34843020000000002;
    x0[708] = 0.087125610000000006;
    x0[709] = 0.1965498;
    x0[710] = 20.153189999999999;
    x0[711] = 1.675284e-12;
    x0[712] = 1.581535e-9;
    x0[713] = 1039.0599999999999;
    x0[715] = 220.58029999999999;
    x0[716] = 22.083390000000001;
    x0[717] = 11.10791;
    x0[718] = 0.3011412;
    x0[719] = 0.51166619999999996;
    x0[720] = 105.30889999999999;
    x0[722] = 0.29566920000000002;
    x0[723] = 18.415199999999999;
    x0[724] = 1.0388189999999999;
    x0[725] = 0.651945;
    x0[726] = 1.100209;
    x0[727] = 0.00056015419999999995;
    x0[728] = 0.0028078769999999999;
    x0[729] = 0.032723229999999999;
    x0[730] = 1.179728e-13;
    x0[731] = 0.0036469359999999999;
    x0[733] = 0.21971640000000001;
    x0[734] = 0.028551239999999999;
    x0[735] = 2.3306909999999998;
    x0[736] = 0.23333709999999999;
    x0[737] = 0.70209589999999999;
    x0[738] = 3.268933e-6;
    x0[739] = 1.344692e-5;
    x0[740] = 0.0048319330000000001;
    x0[741] = 0.019874079999999999;
    x0[742] = 1.582469e-12;
    x0[743] = 0.0033835359999999999;
    x0[744] = 0.00045316369999999999;
    x0[745] = 0.0003899277;
    x0[746] = 0.01442712;
    x0[747] = 0.00067821120000000003;
    x0[748] = 0.001935113;
    x0[750] = 0.0011035400000000001;
    x0[751] = 0.0031171219999999999;
    x0[752] = 0.00053287539999999996;
    x0[753] = 0.00034030049999999999;
    x0[754] = 3.9844269999999999e-5;
    x0[756] = 0.00010565600000000001;
    x0[757] = 0.1163975;
    x0[758] = 0.001354809;
    x0[759] = 0.001541695;
    x0[760] = 0.0002003369;
    x0[761] = 235.65119999999999;
    x0[762] = 15.316789999999999;
    x0[763] = 0.081914589999999995;
    x0[764] = 0.51267560000000001;
    x0[765] = 0.066620079999999998;
    x0[766] = 9.0304939999999991;
    x0[771] = 37.688139999999997;
    x0[772] = 88.811070000000001;
    x0[773] = 0.0025299839999999999;
    x0[774] = 0.015812400000000001;
    x0[775] = 0.0006324961;
    x0[776] = 0.00094874409999999999;
    x0[777] = 0.0006324961;
    x0[778] = 0.005376217;
    x0[779] = 0.005376217;
    x0[780] = 0.005376217;
    x0[781] = 0.005376217;
    x0[782] = 0.0018974879999999999;
    x0[783] = 0.00031624799999999998;
    x0[784] = 0.0006324961;
    x0[785] = 0.005376217;
    x0[786] = 0.005376217;
    x0[787] = 0.005376217;
    x0[788] = 0.005376217;
    x0[789] = 0.005376217;
    x0[790] = 0.005376217;
    x0[791] = 0.005376217;
    x0[792] = 0.005376217;
    x0[793] = 0.005376217;
    x0[794] = 0.005376217;
    x0[795] = 0.005376217;
    x0[796] = 0.005376217;
    x0[797] = 0.005376217;
    x0[798] = 0.0037949759999999998;
    x0[799] = 0.013914910000000001;
    x0[800] = 0.0006324961;
    x0[801] = 0.0034787279999999999;
    x0[802] = 0.002213736;
    x0[804] = 0.00031624799999999998;
    x0[805] = 0.0041112240000000001;
    x0[806] = 0.0031624800000000001;
    x0[807] = 0.0025299839999999999;
    x0[808] = 0.00158124;
    x0[809] = 0.01138493;
    x0[810] = 0.0041112240000000001;
    x0[811] = 0.0072737050000000001;
    x0[812] = 0.0085386969999999996;
    x0[813] = 0.00158124;
    x0[814] = 0.0066412090000000003;
    x0[815] = 0.0047437199999999999;
    x0[817] = 0.00094874409999999999;
    x0[818] = 0.0066412090000000003;
    x0[819] = 0.00094874409999999999;
    x0[820] = 0.02909482;
    x0[821] = 0.0028462320000000002;
    x0[822] = 0.0006324961;
    x0[823] = 0.0006324961;
    x0[824] = 0.0041112240000000001;
    x0[825] = 0.00031624799999999998;
    x0[826] = 0.0012649919999999999;
    x0[827] = 0.00094874409999999999;
    x0[830] = 0.0025299839999999999;
    x0[831] = 0.0012649919999999999;
    x0[832] = 0.005376217;
    x0[835] = 0.0056924649999999999;
    x0[837] = 0.00031624799999999998;
    x0[838] = 0.00031624799999999998;
    x0[840] = 0.0006324961;
    x0[841] = 0.00094874409999999999;
    x0[842] = 0.00158124;
    x0[843] = 0.00031624799999999998;
    x0[844] = 0.00031624799999999998;
    x0[845] = 0.0006324961;
    x0[846] = 0.0006324961;
    x0[849] = 0.0006324961;
    x0[850] = 0.00031624799999999998;
    x0[851] = 0.0028462320000000002;
    x0[852] = 0.00031624799999999998;
    x0[853] = 0.00158124;
    x0[855] = 0.00031624799999999998;
    x0[856] = 0.004427472;
    x0[857] = 0.0025299839999999999;
    x0[858] = 0.00094874409999999999;
    x0[859] = 0.0012649919999999999;
    x0[860] = 0.0018974879999999999;
    x0[861] = 0.00031624799999999998;
    x0[862] = 0.0006324961;
    x0[863] = 0.00094874409999999999;
    x0[864] = 0.00031624799999999998;
    x0[865] = 0.00031624799999999998;
    x0[866] = 0.0006324961;
    x0[867] = 0.0012649919999999999;
    x0[868] = 0.00031624799999999998;
    x0[869] = 0.00094874409999999999;
    x0[870] = 0.00094874409999999999;
    x0[871] = 0.002213736;
    x0[872] = 0.0006324961;
    x0[873] = 0.01043619;
    x0[874] = 0.0072737050000000001;
    x0[875] = 0.00031624799999999998;
    x0[876] = 0.00094874409999999999;
    x0[877] = 0.00094874409999999999;
    x0[878] = 0.0012649919999999999;
    x0[879] = 0.00031624799999999998;
    x0[880] = 0.00094874409999999999;
    x0[881] = 0.0012649919999999999;
    x0[882] = 0.00031624799999999998;
    x0[883] = 0.0006324961;
    x0[884] = 0.00031624799999999998;
    x0[886] = 0.00094874409999999999;
    x0[887] = 0.01170118;
    x0[888] = 0.0018974879999999999;
    x0[889] = 0.0006324961;
    x0[890] = 0.0012649919999999999;
    x0[891] = 0.0098036890000000008;
    x0[892] = 0.00031624799999999998;
    x0[893] = 0.00094874409999999999;
    x0[894] = 0.0006324961;
    x0[895] = 0.012017430000000001;
    x0[896] = 0.00094874409999999999;
    x0[898] = 0.002213736;
    x0[899] = 0.00031624799999999998;
    x0[900] = 0.0006324961;
    x0[903] = 0.012017430000000001;
    x0[904] = 0.00094874409999999999;
    x0[905] = 0.0006324961;
    x0[908] = 0.00094874409999999999;
    x0[909] = 0.00158124;
    x0[910] = 0.0006324961;
    x0[911] = 0.00158124;
    x0[912] = 0.0018974879999999999;
}

} // namespace model_SPARCED
} // namespace amici
