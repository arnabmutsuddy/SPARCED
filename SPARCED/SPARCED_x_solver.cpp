#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "x_rdata.h"

namespace amici {
namespace model_SPARCED {

void x_solver_SPARCED(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Ribosome;
    x_solver[1] = p53inac;
    x_solver[2] = p53ac;
    x_solver[3] = MDM2;
    x_solver[4] = Wip1;
    x_solver[5] = ATMP;
    x_solver[6] = ATRac;
    x_solver[7] = MDM2product1;
    x_solver[8] = MDM2product2;
    x_solver[9] = MDM2product3;
    x_solver[10] = MDM2product4;
    x_solver[11] = MDM2product5;
    x_solver[12] = MDM2product6;
    x_solver[13] = MDM2product7;
    x_solver[14] = MDM2product8;
    x_solver[15] = MDM2product9;
    x_solver[16] = MDM2pro;
    x_solver[17] = Wip1product1;
    x_solver[18] = Wip1product2;
    x_solver[19] = Wip1product3;
    x_solver[20] = Wip1product4;
    x_solver[21] = Wip1product5;
    x_solver[22] = Wip1product6;
    x_solver[23] = Wip1product7;
    x_solver[24] = Wip1product8;
    x_solver[25] = Wip1product9;
    x_solver[26] = Wip1pro;
    x_solver[27] = BRCA2;
    x_solver[28] = MSH6;
    x_solver[29] = MGMT;
    x_solver[30] = damageDSB;
    x_solver[31] = damageSSB;
    x_solver[32] = ppAKT_MDM2;
    x_solver[33] = pMDM2;
    x_solver[34] = ARF;
    x_solver[35] = MDM4;
    x_solver[36] = p53ac_MDM4;
    x_solver[37] = ATMinac;
    x_solver[38] = ATRinac;
    x_solver[39] = pRB;
    x_solver[40] = pRBp;
    x_solver[41] = pRBpp;
    x_solver[42] = E2F;
    x_solver[43] = Cd;
    x_solver[44] = Cd_Cdk46;
    x_solver[45] = Cd_Cdk46_p27;
    x_solver[46] = Ce;
    x_solver[47] = Ce_Cdk2;
    x_solver[48] = Skp2;
    x_solver[49] = Ce_Cdk2_p27;
    x_solver[50] = Pe;
    x_solver[51] = Pai;
    x_solver[52] = Pei;
    x_solver[53] = Pbi;
    x_solver[54] = Ca;
    x_solver[55] = Ca_Cdk2;
    x_solver[56] = Ca_Cdk2_p27;
    x_solver[57] = p27;
    x_solver[58] = Cdh1i;
    x_solver[59] = Cdh1a;
    x_solver[60] = E2Fp;
    x_solver[61] = p27p;
    x_solver[62] = Pa;
    x_solver[63] = Cb;
    x_solver[64] = Cb_Cdk1;
    x_solver[65] = Cdc20i;
    x_solver[66] = Cdc20a;
    x_solver[67] = Pb;
    x_solver[68] = Wee1;
    x_solver[69] = Wee1p;
    x_solver[70] = Cb_Cdk1_p27;
    x_solver[71] = Chk1;
    x_solver[72] = pRB_E2F;
    x_solver[73] = pRBp_E2F;
    x_solver[74] = p21;
    x_solver[75] = Cd_Cdk46_p21;
    x_solver[76] = Ce_Cdk2_p21;
    x_solver[77] = Ca_Cdk2_p21;
    x_solver[78] = Cb_Cdk1_p21;
    x_solver[79] = L;
    x_solver[80] = R;
    x_solver[81] = L_R;
    x_solver[82] = Ractive;
    x_solver[83] = flip;
    x_solver[84] = Ractive_flip;
    x_solver[85] = pC8;
    x_solver[86] = Ractive_pC8;
    x_solver[87] = C8;
    x_solver[88] = Bar;
    x_solver[89] = C8_Bar;
    x_solver[90] = pC3;
    x_solver[91] = C8_pC3;
    x_solver[92] = C3;
    x_solver[93] = pC6;
    x_solver[94] = C3_pC6;
    x_solver[95] = C6;
    x_solver[96] = C6_C8;
    x_solver[97] = XIAP;
    x_solver[98] = C3_XIAP;
    x_solver[99] = PARP;
    x_solver[100] = C3_PARP;
    x_solver[101] = cPARP;
    x_solver[102] = Bid;
    x_solver[103] = C8_Bid;
    x_solver[104] = tBid;
    x_solver[105] = Bcl2c;
    x_solver[106] = tBid_Bcl2c;
    x_solver[107] = Bax;
    x_solver[108] = tBid_Bax;
    x_solver[109] = Baxactive;
    x_solver[110] = Baxm;
    x_solver[111] = Bcl2;
    x_solver[112] = Baxm_Bcl2;
    x_solver[113] = Bax2;
    x_solver[114] = Bax2_Bcl2;
    x_solver[115] = Bax4;
    x_solver[116] = Bax4_Bcl2;
    x_solver[117] = M;
    x_solver[118] = Bax4_M;
    x_solver[119] = Mactive;
    x_solver[120] = CytoCm;
    x_solver[121] = Mactive_CytoCm;
    x_solver[122] = CytoCr;
    x_solver[123] = Smacm;
    x_solver[124] = Mactive_Smacm;
    x_solver[125] = Smacr;
    x_solver[126] = CytoC;
    x_solver[127] = Apaf;
    x_solver[128] = CytoC_Apaf;
    x_solver[129] = Apafactive;
    x_solver[130] = pC9;
    x_solver[131] = Apop;
    x_solver[132] = Apop_C3;
    x_solver[133] = Smac;
    x_solver[134] = Apop_XIAP;
    x_solver[135] = Smac_XIAP;
    x_solver[136] = C3_Ub;
    x_solver[137] = BAD;
    x_solver[138] = PUMA;
    x_solver[139] = NOXA;
    x_solver[140] = Bcl2c_BAD;
    x_solver[141] = Bcl2c_PUMA;
    x_solver[142] = Bcl2c_NOXA;
    x_solver[143] = BIM;
    x_solver[144] = BIM_Bax;
    x_solver[145] = Bcl2c_BIM;
    x_solver[146] = ppERK_BIM;
    x_solver[147] = pBIM;
    x_solver[148] = ppAKT_BAD;
    x_solver[149] = pBAD;
    x_solver[150] = ppERK_BAD;
    x_solver[151] = E;
    x_solver[152] = H;
    x_solver[153] = HGF;
    x_solver[154] = P;
    x_solver[155] = F;
    x_solver[156] = I;
    x_solver[157] = INS;
    x_solver[158] = EGFR;
    x_solver[159] = pE1;
    x_solver[160] = E2;
    x_solver[161] = pE2;
    x_solver[162] = E3;
    x_solver[163] = E4;
    x_solver[164] = pE4;
    x_solver[165] = Ev3;
    x_solver[166] = Met;
    x_solver[167] = Pr;
    x_solver[168] = Fr;
    x_solver[169] = Ir;
    x_solver[170] = Isr;
    x_solver[171] = E1E1;
    x_solver[172] = E1E2;
    x_solver[173] = E1E3;
    x_solver[174] = E1E4;
    x_solver[175] = E2E2;
    x_solver[176] = E2E3;
    x_solver[177] = E2E4;
    x_solver[178] = E3E4;
    x_solver[179] = E4E4;
    x_solver[180] = Met_Met;
    x_solver[181] = FrFr;
    x_solver[182] = IrIr;
    x_solver[183] = Isr_Isr;
    x_solver[184] = EE1;
    x_solver[185] = HE3;
    x_solver[186] = HE4;
    x_solver[187] = HGF_Met;
    x_solver[188] = PPr;
    x_solver[189] = FFr;
    x_solver[190] = EE1E2;
    x_solver[191] = EE1Ev3;
    x_solver[192] = EE1E1;
    x_solver[193] = EE1E3;
    x_solver[194] = EE1E4;
    x_solver[195] = E2HE3;
    x_solver[196] = E1HE3;
    x_solver[197] = HE3E3;
    x_solver[198] = HE3Ev3;
    x_solver[199] = HE3E4;
    x_solver[200] = E2HE4;
    x_solver[201] = HE4Ev3;
    x_solver[202] = E1HE4;
    x_solver[203] = E3HE4;
    x_solver[204] = HE4E4;
    x_solver[205] = HGF_Met_Met;
    x_solver[206] = PPrPr;
    x_solver[207] = FFrFr;
    x_solver[208] = IIrIr;
    x_solver[209] = INS_Isr_Isr;
    x_solver[210] = EE1EE1;
    x_solver[211] = EE1HE3;
    x_solver[212] = EE1HE4;
    x_solver[213] = HE3HE3;
    x_solver[214] = HE3HE4;
    x_solver[215] = HE4HE4;
    x_solver[216] = HGF_Met_HGF_Met;
    x_solver[217] = PPrPPr;
    x_solver[218] = FFrFFr;
    x_solver[219] = IIrIrI;
    x_solver[220] = INS_Isr_Isr_INS;
    x_solver[221] = E1_ppERK;
    x_solver[222] = E2_ppERK;
    x_solver[223] = E4_ppERK;
    x_solver[224] = pEE1E2;
    x_solver[225] = pEE1Ev3;
    x_solver[226] = pEE1E1;
    x_solver[227] = pEE1EE1;
    x_solver[228] = pEE1E3;
    x_solver[229] = pEE1HE3;
    x_solver[230] = pEE1E4;
    x_solver[231] = pEE1HE4;
    x_solver[232] = pE2HE3;
    x_solver[233] = pHE3Ev3;
    x_solver[234] = pE1HE3;
    x_solver[235] = pHE3E4;
    x_solver[236] = pHE3HE4;
    x_solver[237] = pE2HE4;
    x_solver[238] = pHE4Ev3;
    x_solver[239] = pE1HE4;
    x_solver[240] = pE3HE4;
    x_solver[241] = pHE4E4;
    x_solver[242] = pHE4HE4;
    x_solver[243] = pHGF_Met_Met;
    x_solver[244] = pHGF_Met_HGF_Met;
    x_solver[245] = pPPrPPr;
    x_solver[246] = pPPrPr;
    x_solver[247] = pFFrFFr;
    x_solver[248] = pFFrFr;
    x_solver[249] = pIIrIr;
    x_solver[250] = pINS_Isr_Isr;
    x_solver[251] = pIIrIrI;
    x_solver[252] = pINS_Isr_Isr_INS;
    x_solver[253] = pIIrIr_IRS;
    x_solver[254] = pINS_Isr_Isr_IRS;
    x_solver[255] = pIIrIrI_IRS;
    x_solver[256] = pINS_Isr_Isr_INS_IRS;
    x_solver[257] = Sp_EE1E2;
    x_solver[258] = Sp_EE1Ev3;
    x_solver[259] = Sp_EE1E1;
    x_solver[260] = Sp_EE1EE1;
    x_solver[261] = Sp_EE1E3;
    x_solver[262] = Sp_EE1HE3;
    x_solver[263] = Sp_EE1E4;
    x_solver[264] = Sp_EE1HE4;
    x_solver[265] = Sp_E2HE3;
    x_solver[266] = Sp_HE3Ev3;
    x_solver[267] = Sp_E1HE3;
    x_solver[268] = Sp_HE3E4;
    x_solver[269] = Sp_HE3HE4;
    x_solver[270] = Sp_E2HE4;
    x_solver[271] = Sp_HE4Ev3;
    x_solver[272] = Sp_E1HE4;
    x_solver[273] = Sp_E3HE4;
    x_solver[274] = Sp_HE4E4;
    x_solver[275] = Sp_HE4HE4;
    x_solver[276] = Sp_HGF_Met_Met;
    x_solver[277] = Sp_HGF_Met_HGF_Met;
    x_solver[278] = Sp_PPrPPr;
    x_solver[279] = Sp_PPrPr;
    x_solver[280] = Sp_FFrFFr;
    x_solver[281] = Sp_FFrFr;
    x_solver[282] = Sp_IIrIr;
    x_solver[283] = Sp_INS_Isr_Isr;
    x_solver[284] = Sp_IIrIrI;
    x_solver[285] = Sp_INS_Isr_Isr_INS;
    x_solver[286] = EE1E2int;
    x_solver[287] = EE1Ev3int;
    x_solver[288] = EE1E1int;
    x_solver[289] = EE1EE1int;
    x_solver[290] = EE1E3int;
    x_solver[291] = EE1HE3int;
    x_solver[292] = EE1E4int;
    x_solver[293] = EE1HE4int;
    x_solver[294] = E2HE3int;
    x_solver[295] = HE3Ev3int;
    x_solver[296] = E1HE3int;
    x_solver[297] = HE3E4int;
    x_solver[298] = HE3HE4int;
    x_solver[299] = E2HE4int;
    x_solver[300] = HE4Ev3int;
    x_solver[301] = E1HE4int;
    x_solver[302] = E3HE4int;
    x_solver[303] = HE4E4int;
    x_solver[304] = HE4HE4int;
    x_solver[305] = HGF_Met_Metint;
    x_solver[306] = HGF_Met_HGF_Metint;
    x_solver[307] = PPrPPrint;
    x_solver[308] = PPrPrint;
    x_solver[309] = FFrFFrint;
    x_solver[310] = FFrFrint;
    x_solver[311] = IIrIr_int;
    x_solver[312] = INS_Isr_Isr_int;
    x_solver[313] = IIrIrI_int;
    x_solver[314] = INS_Isr_Isr_INS_int;
    x_solver[315] = pEE1E2int;
    x_solver[316] = pEE1Ev3int;
    x_solver[317] = pEE1E1int;
    x_solver[318] = pEE1EE1int;
    x_solver[319] = pEE1E3int;
    x_solver[320] = pEE1HE3int;
    x_solver[321] = pEE1E4int;
    x_solver[322] = pEE1HE4int;
    x_solver[323] = pE2HE3int;
    x_solver[324] = pHE3Ev3int;
    x_solver[325] = pE1HE3int;
    x_solver[326] = pHE3E4int;
    x_solver[327] = pHE3HE4int;
    x_solver[328] = pE2HE4int;
    x_solver[329] = pHE4Ev3int;
    x_solver[330] = pE1HE4int;
    x_solver[331] = pE3HE4int;
    x_solver[332] = pHE4E4int;
    x_solver[333] = pHE4HE4int;
    x_solver[334] = pHGF_Met_Metint;
    x_solver[335] = pHGF_Met_HGF_Metint;
    x_solver[336] = pPPrPPrint;
    x_solver[337] = pPPrPrint;
    x_solver[338] = pFFrFFrint;
    x_solver[339] = pFFrFrint;
    x_solver[340] = pIIrIr_int;
    x_solver[341] = pINS_Isr_Isr_int;
    x_solver[342] = pIIrIrI_int;
    x_solver[343] = pINS_Isr_Isr_INS_int;
    x_solver[344] = pIIrIr_int_IRS;
    x_solver[345] = pINS_Isr_Isr_int_IRS;
    x_solver[346] = pIIrIrI_int_IRS;
    x_solver[347] = pINS_Isr_Isr_INS_int_IRS;
    x_solver[348] = pEE1E2_G2_SOS;
    x_solver[349] = pEE1Ev3_G2_SOS;
    x_solver[350] = pEE1E1_G2_SOS;
    x_solver[351] = pEE1EE1_G2_SOS;
    x_solver[352] = pEE1E3_G2_SOS;
    x_solver[353] = pEE1HE3_G2_SOS;
    x_solver[354] = pEE1E4_G2_SOS;
    x_solver[355] = pEE1HE4_G2_SOS;
    x_solver[356] = pE2HE3_G2_SOS;
    x_solver[357] = pHE3Ev3_G2_SOS;
    x_solver[358] = pE1HE3_G2_SOS;
    x_solver[359] = pHE3E4_G2_SOS;
    x_solver[360] = pHE3HE4_G2_SOS;
    x_solver[361] = pE2HE4_G2_SOS;
    x_solver[362] = pHE4Ev3_G2_SOS;
    x_solver[363] = pE1HE4_G2_SOS;
    x_solver[364] = pE3HE4_G2_SOS;
    x_solver[365] = pHE4E4_G2_SOS;
    x_solver[366] = pHE4HE4_G2_SOS;
    x_solver[367] = pHGF_Met_Met_G2_SOS;
    x_solver[368] = pHGF_Met_HGF_Met_G2_SOS;
    x_solver[369] = pPPrPPr_G2_SOS;
    x_solver[370] = pPPrPr_G2_SOS;
    x_solver[371] = pFFrFFr_G2_SOS;
    x_solver[372] = pFFrFr_G2_SOS;
    x_solver[373] = pIIrIr_IRS_G2_SOS;
    x_solver[374] = pINS_Isr_Isr_IRS_G2_SOS;
    x_solver[375] = pIIrIrI_IRS_G2_SOS;
    x_solver[376] = pINS_Isr_Isr_INS_IRS_G2_SOS;
    x_solver[377] = pEE1E2int_G2_SOS;
    x_solver[378] = pEE1Ev3int_G2_SOS;
    x_solver[379] = pEE1E1int_G2_SOS;
    x_solver[380] = pEE1EE1int_G2_SOS;
    x_solver[381] = pEE1E3int_G2_SOS;
    x_solver[382] = pEE1HE3int_G2_SOS;
    x_solver[383] = pEE1E4int_G2_SOS;
    x_solver[384] = pEE1HE4int_G2_SOS;
    x_solver[385] = pE2HE3int_G2_SOS;
    x_solver[386] = pHE3Ev3int_G2_SOS;
    x_solver[387] = pE1HE3int_G2_SOS;
    x_solver[388] = pHE3E4int_G2_SOS;
    x_solver[389] = pHE3HE4int_G2_SOS;
    x_solver[390] = pE2HE4int_G2_SOS;
    x_solver[391] = pHE4Ev3int_G2_SOS;
    x_solver[392] = pE1HE4int_G2_SOS;
    x_solver[393] = pE3HE4int_G2_SOS;
    x_solver[394] = pHE4E4int_G2_SOS;
    x_solver[395] = pHE4HE4int_G2_SOS;
    x_solver[396] = pHGF_Met_Metint_G2_SOS;
    x_solver[397] = pHGF_Met_HGF_Metint_G2_SOS;
    x_solver[398] = pPPrPPrint_G2_SOS;
    x_solver[399] = pPPrPrint_G2_SOS;
    x_solver[400] = pFFrFFrint_G2_SOS;
    x_solver[401] = pFFrFrint_G2_SOS;
    x_solver[402] = pIIrIr_int_IRS_G2_SOS;
    x_solver[403] = pINS_Isr_Isr_int_IRS_G2_SOS;
    x_solver[404] = pIIrIrI_int_IRS_G2_SOS;
    x_solver[405] = pINS_Isr_Isr_INS_int_IRS_G2_SOS;
    x_solver[406] = pEE1E2_PLCg;
    x_solver[407] = pEE1Ev3_PLCg;
    x_solver[408] = pEE1E1_PLCg;
    x_solver[409] = pEE1EE1_PLCg;
    x_solver[410] = pEE1E3_PLCg;
    x_solver[411] = pEE1HE3_PLCg;
    x_solver[412] = pEE1E4_PLCg;
    x_solver[413] = pEE1HE4_PLCg;
    x_solver[414] = pE2HE3_PLCg;
    x_solver[415] = pHE3Ev3_PLCg;
    x_solver[416] = pE1HE3_PLCg;
    x_solver[417] = pHE3E4_PLCg;
    x_solver[418] = pHE3HE4_PLCg;
    x_solver[419] = pE2HE4_PLCg;
    x_solver[420] = pHE4Ev3_PLCg;
    x_solver[421] = pE1HE4_PLCg;
    x_solver[422] = pE3HE4_PLCg;
    x_solver[423] = pHE4E4_PLCg;
    x_solver[424] = pHE4HE4_PLCg;
    x_solver[425] = pHGF_Met_Met_PLCg;
    x_solver[426] = pHGF_Met_HGF_Met_PLCg;
    x_solver[427] = pPPrPPr_PLCg;
    x_solver[428] = pPPrPr_PLCg;
    x_solver[429] = pFFrFFr_PLCg;
    x_solver[430] = pFFrFr_PLCg;
    x_solver[431] = pIIrIr_IRS_PLCg;
    x_solver[432] = pINS_Isr_Isr_IRS_PLCg;
    x_solver[433] = pIIrIrI_IRS_PLCg;
    x_solver[434] = pINS_Isr_Isr_INS_IRS_PLCg;
    x_solver[435] = pEE1E2_PI3K1;
    x_solver[436] = pEE1Ev3_PI3K1;
    x_solver[437] = pEE1E1_PI3K1;
    x_solver[438] = pEE1EE1_PI3K1;
    x_solver[439] = pEE1E3_PI3K1;
    x_solver[440] = pEE1HE3_PI3K1;
    x_solver[441] = pEE1E4_PI3K1;
    x_solver[442] = pEE1HE4_PI3K1;
    x_solver[443] = pE2HE3_PI3K1;
    x_solver[444] = pHE3Ev3_PI3K1;
    x_solver[445] = pE1HE3_PI3K1;
    x_solver[446] = pHE3E4_PI3K1;
    x_solver[447] = pHE3HE4_PI3K1;
    x_solver[448] = pE2HE4_PI3K1;
    x_solver[449] = pHE4Ev3_PI3K1;
    x_solver[450] = pE1HE4_PI3K1;
    x_solver[451] = pE3HE4_PI3K1;
    x_solver[452] = pHE4E4_PI3K1;
    x_solver[453] = pHE4HE4_PI3K1;
    x_solver[454] = pHGF_Met_Met_PI3K1;
    x_solver[455] = pHGF_Met_HGF_Met_PI3K1;
    x_solver[456] = pPPrPPr_PI3K1;
    x_solver[457] = pPPrPr_PI3K1;
    x_solver[458] = pFFrFFr_PI3K1;
    x_solver[459] = pFFrFr_PI3K1;
    x_solver[460] = pIIrIr_IRS_PI3K1;
    x_solver[461] = pINS_Isr_Isr_IRS_PI3K1;
    x_solver[462] = pIIrIrI_IRS_PI3K1;
    x_solver[463] = pINS_Isr_Isr_INS_IRS_PI3K1;
    x_solver[464] = pEE1E2_PI3K2;
    x_solver[465] = pEE1Ev3_PI3K2;
    x_solver[466] = pEE1E1_PI3K2;
    x_solver[467] = pEE1EE1_PI3K2;
    x_solver[468] = pEE1E3_PI3K2;
    x_solver[469] = pEE1HE3_PI3K2;
    x_solver[470] = pEE1E4_PI3K2;
    x_solver[471] = pEE1HE4_PI3K2;
    x_solver[472] = pE2HE3_PI3K2;
    x_solver[473] = pHE3Ev3_PI3K2;
    x_solver[474] = pE1HE3_PI3K2;
    x_solver[475] = pHE3E4_PI3K2;
    x_solver[476] = pHE3HE4_PI3K2;
    x_solver[477] = pE2HE4_PI3K2;
    x_solver[478] = pHE4Ev3_PI3K2;
    x_solver[479] = pE1HE4_PI3K2;
    x_solver[480] = pE3HE4_PI3K2;
    x_solver[481] = pHE4E4_PI3K2;
    x_solver[482] = pHE4HE4_PI3K2;
    x_solver[483] = pHGF_Met_Met_PI3K2;
    x_solver[484] = pHGF_Met_HGF_Met_PI3K2;
    x_solver[485] = pPPrPPr_PI3K2;
    x_solver[486] = pPPrPr_PI3K2;
    x_solver[487] = pFFrFFr_PI3K2;
    x_solver[488] = pFFrFr_PI3K2;
    x_solver[489] = pIIrIr_IRS_PI3K2;
    x_solver[490] = pINS_Isr_Isr_IRS_PI3K2;
    x_solver[491] = pIIrIrI_IRS_PI3K2;
    x_solver[492] = pINS_Isr_Isr_INS_IRS_PI3K2;
    x_solver[493] = pEE1E2int_G2_SOS_RasD;
    x_solver[494] = pEE1Ev3int_G2_SOS_RasD;
    x_solver[495] = pEE1E1int_G2_SOS_RasD;
    x_solver[496] = pEE1EE1int_G2_SOS_RasD;
    x_solver[497] = pEE1E3int_G2_SOS_RasD;
    x_solver[498] = pEE1HE3int_G2_SOS_RasD;
    x_solver[499] = pEE1E4int_G2_SOS_RasD;
    x_solver[500] = pEE1HE4int_G2_SOS_RasD;
    x_solver[501] = pE2HE3int_G2_SOS_RasD;
    x_solver[502] = pHE3Ev3int_G2_SOS_RasD;
    x_solver[503] = pE1HE3int_G2_SOS_RasD;
    x_solver[504] = pHE3E4int_G2_SOS_RasD;
    x_solver[505] = pHE3HE4int_G2_SOS_RasD;
    x_solver[506] = pE2HE4int_G2_SOS_RasD;
    x_solver[507] = pHE4Ev3int_G2_SOS_RasD;
    x_solver[508] = pE1HE4int_G2_SOS_RasD;
    x_solver[509] = pE3HE4int_G2_SOS_RasD;
    x_solver[510] = pHE4E4int_G2_SOS_RasD;
    x_solver[511] = pHE4HE4int_G2_SOS_RasD;
    x_solver[512] = pHGF_Met_Metint_G2_SOS_RasD;
    x_solver[513] = pHGF_Met_HGF_Metint_G2_SOS_RasD;
    x_solver[514] = pPPrPPrint_G2_SOS_RasD;
    x_solver[515] = pPPrPrint_G2_SOS_RasD;
    x_solver[516] = pFFrFFrint_G2_SOS_RasD;
    x_solver[517] = pFFrFrint_G2_SOS_RasD;
    x_solver[518] = pIIrIr_int_IRS_G2_SOS_RasD;
    x_solver[519] = pINS_Isr_Isr_int_IRS_G2_SOS_RasD;
    x_solver[520] = pIIrIrI_int_IRS_G2_SOS_RasD;
    x_solver[521] = pINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD;
    x_solver[522] = pEE1E2_G2_SOS_RasD;
    x_solver[523] = pEE1Ev3_G2_SOS_RasD;
    x_solver[524] = pEE1E1_G2_SOS_RasD;
    x_solver[525] = pEE1EE1_G2_SOS_RasD;
    x_solver[526] = pEE1E3_G2_SOS_RasD;
    x_solver[527] = pEE1HE3_G2_SOS_RasD;
    x_solver[528] = pEE1E4_G2_SOS_RasD;
    x_solver[529] = pEE1HE4_G2_SOS_RasD;
    x_solver[530] = pE2HE3_G2_SOS_RasD;
    x_solver[531] = pHE3Ev3_G2_SOS_RasD;
    x_solver[532] = pE1HE3_G2_SOS_RasD;
    x_solver[533] = pHE3E4_G2_SOS_RasD;
    x_solver[534] = pHE3HE4_G2_SOS_RasD;
    x_solver[535] = pE2HE4_G2_SOS_RasD;
    x_solver[536] = pHE4Ev3_G2_SOS_RasD;
    x_solver[537] = pE1HE4_G2_SOS_RasD;
    x_solver[538] = pE3HE4_G2_SOS_RasD;
    x_solver[539] = pHE4E4_G2_SOS_RasD;
    x_solver[540] = pHE4HE4_G2_SOS_RasD;
    x_solver[541] = pHGF_Met_Met_G2_SOS_RasD;
    x_solver[542] = pHGF_Met_HGF_Met_G2_SOS_RasD;
    x_solver[543] = pPPrPPr_G2_SOS_RasD;
    x_solver[544] = pPPrPr_G2_SOS_RasD;
    x_solver[545] = pFFrFFr_G2_SOS_RasD;
    x_solver[546] = pFFrFr_G2_SOS_RasD;
    x_solver[547] = pIIrIr_IRS_G2_SOS_RasD;
    x_solver[548] = pINS_Isr_Isr_IRS_G2_SOS_RasD;
    x_solver[549] = pIIrIrI_IRS_G2_SOS_RasD;
    x_solver[550] = pINS_Isr_Isr_INS_IRS_G2_SOS_RasD;
    x_solver[551] = pEE1E2_PLCg_PIP2;
    x_solver[552] = pEE1Ev3_PLCg_PIP2;
    x_solver[553] = pEE1E1_PLCg_PIP2;
    x_solver[554] = pEE1EE1_PLCg_PIP2;
    x_solver[555] = pEE1E3_PLCg_PIP2;
    x_solver[556] = pEE1HE3_PLCg_PIP2;
    x_solver[557] = pEE1E4_PLCg_PIP2;
    x_solver[558] = pEE1HE4_PLCg_PIP2;
    x_solver[559] = pE2HE3_PLCg_PIP2;
    x_solver[560] = pHE3Ev3_PLCg_PIP2;
    x_solver[561] = pE1HE3_PLCg_PIP2;
    x_solver[562] = pHE3E4_PLCg_PIP2;
    x_solver[563] = pHE3HE4_PLCg_PIP2;
    x_solver[564] = pE2HE4_PLCg_PIP2;
    x_solver[565] = pHE4Ev3_PLCg_PIP2;
    x_solver[566] = pE1HE4_PLCg_PIP2;
    x_solver[567] = pE3HE4_PLCg_PIP2;
    x_solver[568] = pHE4E4_PLCg_PIP2;
    x_solver[569] = pHE4HE4_PLCg_PIP2;
    x_solver[570] = pHGF_Met_Met_PLCg_PIP2;
    x_solver[571] = pHGF_Met_HGF_Met_PLCg_PIP2;
    x_solver[572] = pPPrPPr_PLCg_PIP2;
    x_solver[573] = pPPrPr_PLCg_PIP2;
    x_solver[574] = pFFrFFr_PLCg_PIP2;
    x_solver[575] = pFFrFr_PLCg_PIP2;
    x_solver[576] = pIIrIr_IRS_PLCg_PIP2;
    x_solver[577] = pINS_Isr_Isr_IRS_PLCg_PIP2;
    x_solver[578] = pIIrIrI_IRS_PLCg_PIP2;
    x_solver[579] = pINS_Isr_Isr_INS_IRS_PLCg_PIP2;
    x_solver[580] = pEE1E2_PI3K1_PIP2;
    x_solver[581] = pEE1Ev3_PI3K1_PIP2;
    x_solver[582] = pEE1E1_PI3K1_PIP2;
    x_solver[583] = pEE1EE1_PI3K1_PIP2;
    x_solver[584] = pEE1E3_PI3K1_PIP2;
    x_solver[585] = pEE1HE3_PI3K1_PIP2;
    x_solver[586] = pEE1E4_PI3K1_PIP2;
    x_solver[587] = pEE1HE4_PI3K1_PIP2;
    x_solver[588] = pE2HE3_PI3K1_PIP2;
    x_solver[589] = pHE3Ev3_PI3K1_PIP2;
    x_solver[590] = pE1HE3_PI3K1_PIP2;
    x_solver[591] = pHE3E4_PI3K1_PIP2;
    x_solver[592] = pHE3HE4_PI3K1_PIP2;
    x_solver[593] = pE2HE4_PI3K1_PIP2;
    x_solver[594] = pHE4Ev3_PI3K1_PIP2;
    x_solver[595] = pE1HE4_PI3K1_PIP2;
    x_solver[596] = pE3HE4_PI3K1_PIP2;
    x_solver[597] = pHE4E4_PI3K1_PIP2;
    x_solver[598] = pHE4HE4_PI3K1_PIP2;
    x_solver[599] = pHGF_Met_Met_PI3K1_PIP2;
    x_solver[600] = pHGF_Met_HGF_Met_PI3K1_PIP2;
    x_solver[601] = pPPrPPr_PI3K1_PIP2;
    x_solver[602] = pPPrPr_PI3K1_PIP2;
    x_solver[603] = pFFrFFr_PI3K1_PIP2;
    x_solver[604] = pFFrFr_PI3K1_PIP2;
    x_solver[605] = pIIrIr_IRS_PI3K1_PIP2;
    x_solver[606] = pINS_Isr_Isr_IRS_PI3K1_PIP2;
    x_solver[607] = pIIrIrI_IRS_PI3K1_PIP2;
    x_solver[608] = pINS_Isr_Isr_INS_IRS_PI3K1_PIP2;
    x_solver[609] = pEE1E2_PI3K2_PIP;
    x_solver[610] = pEE1Ev3_PI3K2_PIP;
    x_solver[611] = pEE1E1_PI3K2_PIP;
    x_solver[612] = pEE1EE1_PI3K2_PIP;
    x_solver[613] = pEE1E3_PI3K2_PIP;
    x_solver[614] = pEE1HE3_PI3K2_PIP;
    x_solver[615] = pEE1E4_PI3K2_PIP;
    x_solver[616] = pEE1HE4_PI3K2_PIP;
    x_solver[617] = pE2HE3_PI3K2_PIP;
    x_solver[618] = pHE3Ev3_PI3K2_PIP;
    x_solver[619] = pE1HE3_PI3K2_PIP;
    x_solver[620] = pHE3E4_PI3K2_PIP;
    x_solver[621] = pHE3HE4_PI3K2_PIP;
    x_solver[622] = pE2HE4_PI3K2_PIP;
    x_solver[623] = pHE4Ev3_PI3K2_PIP;
    x_solver[624] = pE1HE4_PI3K2_PIP;
    x_solver[625] = pE3HE4_PI3K2_PIP;
    x_solver[626] = pHE4E4_PI3K2_PIP;
    x_solver[627] = pHE4HE4_PI3K2_PIP;
    x_solver[628] = pHGF_Met_Met_PI3K2_PIP;
    x_solver[629] = pHGF_Met_HGF_Met_PI3K2_PIP;
    x_solver[630] = pPPrPPr_PI3K2_PIP;
    x_solver[631] = pPPrPr_PI3K2_PIP;
    x_solver[632] = pFFrFFr_PI3K2_PIP;
    x_solver[633] = pFFrFr_PI3K2_PIP;
    x_solver[634] = pIIrIr_IRS_PI3K2_PIP;
    x_solver[635] = pINS_Isr_Isr_IRS_PI3K2_PIP;
    x_solver[636] = pIIrIrI_IRS_PI3K2_PIP;
    x_solver[637] = pINS_Isr_Isr_INS_IRS_PI3K2_PIP;
    x_solver[638] = IRS;
    x_solver[639] = Sp;
    x_solver[640] = Cbl;
    x_solver[641] = G2;
    x_solver[642] = G2_SOS;
    x_solver[643] = G2_pSOS;
    x_solver[644] = PLCg;
    x_solver[645] = PI3KC1;
    x_solver[646] = PI3KR1;
    x_solver[647] = PI3K1;
    x_solver[648] = pPI3K1;
    x_solver[649] = PI3K2;
    x_solver[650] = mTORC1;
    x_solver[651] = mTORC1active;
    x_solver[652] = PIP;
    x_solver[653] = PI3P;
    x_solver[654] = DAG;
    x_solver[655] = GRP;
    x_solver[656] = DAG_GRP;
    x_solver[657] = RasT;
    x_solver[658] = RasD;
    x_solver[659] = NF1;
    x_solver[660] = pNF1;
    x_solver[661] = pCRaf;
    x_solver[662] = CRaf;
    x_solver[663] = RasT_CRaf;
    x_solver[664] = BRaf;
    x_solver[665] = RasT_CRaf_BRaf;
    x_solver[666] = MEK;
    x_solver[667] = pMEK;
    x_solver[668] = ppMEK;
    x_solver[669] = MKP3;
    x_solver[670] = ERKnuc;
    x_solver[671] = ppERKnuc;
    x_solver[672] = RSK;
    x_solver[673] = pRSK;
    x_solver[674] = pRSKnuc;
    x_solver[675] = MKP1;
    x_solver[676] = pMKP1;
    x_solver[677] = cFos;
    x_solver[678] = pcFos;
    x_solver[679] = cJun;
    x_solver[680] = pcFos_cJun;
    x_solver[681] = cMyc;
    x_solver[682] = bCATENINnuc;
    x_solver[683] = bCATENIN;
    x_solver[684] = pbCATENIN;
    x_solver[685] = IP3;
    x_solver[686] = PIP2;
    x_solver[687] = PIP3;
    x_solver[688] = PTEN;
    x_solver[689] = PIP3_AKT;
    x_solver[690] = AKT;
    x_solver[691] = pAKT;
    x_solver[692] = ppAKT;
    x_solver[693] = PDK1;
    x_solver[694] = PIP3_PDK1;
    x_solver[695] = PIP3_pAKT;
    x_solver[696] = Rictor;
    x_solver[697] = mTOR;
    x_solver[698] = mTORC2;
    x_solver[699] = PIP3_ppAKT;
    x_solver[700] = GSK3b;
    x_solver[701] = pGSK3b;
    x_solver[702] = TSC1;
    x_solver[703] = TSC2;
    x_solver[704] = pTSC2;
    x_solver[705] = TSC;
    x_solver[706] = PKC;
    x_solver[707] = DAG_PKC;
    x_solver[708] = pRKIP;
    x_solver[709] = RKIP;
    x_solver[710] = RKIP_CRaf;
    x_solver[711] = ERK;
    x_solver[712] = pERK;
    x_solver[713] = ppERK;
    x_solver[714] = FOXO;
    x_solver[715] = pFOXO;
    x_solver[716] = RhebD;
    x_solver[717] = RhebT;
    x_solver[718] = Raptor;
    x_solver[719] = S6K;
    x_solver[720] = pS6K;
    x_solver[721] = EIF4EBP1;
    x_solver[722] = pEIF4EBP1;
    x_solver[723] = SOS;
    x_solver[724] = G2_SOS_ppERK;
    x_solver[725] = CRaf_ppERK;
    x_solver[726] = RasD_DAG_GRP;
    x_solver[727] = RasT_NF1;
    x_solver[728] = NF1_ppERK;
    x_solver[729] = MEK_RasT_CRaf_BRaf;
    x_solver[730] = pMEK_RasT_CRaf_BRaf;
    x_solver[731] = ERK_ppMEK;
    x_solver[732] = pERK_ppMEK;
    x_solver[733] = RSK_ppERK;
    x_solver[734] = pRSKnuc_MKP1;
    x_solver[735] = ppERKnuc_MKP1;
    x_solver[736] = cFos_pRSKnuc;
    x_solver[737] = cFos_ppERKnuc;
    x_solver[738] = RKIP_DAG_PKC;
    x_solver[739] = PIP3_PTEN;
    x_solver[740] = PIP3_AKT_PIP3_PDK1;
    x_solver[741] = PIP3_pAKT_mTORC2;
    x_solver[742] = GSK3b_ppAKT;
    x_solver[743] = TSC2_ppAKT;
    x_solver[744] = TSC2_ppERK;
    x_solver[745] = RhebT_TSC;
    x_solver[746] = EIF4EBP1_mTORC1active;
    x_solver[747] = S6K_mTORC1active;
    x_solver[748] = FOXO_ppAKT;
    x_solver[749] = PI3K1_mTORC1active;
    x_solver[750] = pERK_MKP3;
    x_solver[751] = ppERK_MKP3;
    x_solver[752] = ppERKnuc_pMKP1;
    x_solver[753] = RasT_BRaf;
    x_solver[754] = RasT_BRaf_BRaf;
    x_solver[755] = MEK_RasT_BRaf_BRaf;
    x_solver[756] = pMEK_RasT_BRaf_BRaf;
    x_solver[757] = EIF4E;
    x_solver[758] = EIF4EBP1_EIF4E;
    x_solver[759] = RasT_CRaf_CRaf;
    x_solver[760] = MEK_RasT_CRaf_CRaf;
    x_solver[761] = pMEK_RasT_CRaf_CRaf;
    x_solver[762] = FOXOnuc;
    x_solver[763] = MEKi;
    x_solver[764] = MEKi_ppMEK;
    x_solver[765] = AKTi;
    x_solver[766] = AKTi_AKT;
    x_solver[767] = mT;
    x_solver[768] = EIF4E_mT;
    x_solver[769] = Cdk1;
    x_solver[770] = Cdk2;
    x_solver[771] = Cdk46;
    x_solver[772] = pRBpp_E2F;
    x_solver[773] = Cd_Cdk46_pRB;
    x_solver[774] = Cd_Cdk46_pRB_E2F;
    x_solver[775] = Ce_Cdk2_pRBp;
    x_solver[776] = Ce_Cdk2_pRBp_E2F;
    x_solver[777] = p18;
}

} // namespace amici
} // namespace model_SPARCED