#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "tcl.h"
#include "x.h"

namespace amici {
namespace model_SPARCED {

void x_rdata_SPARCED(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = Ribosome;
    x_rdata[1] = p53inac;
    x_rdata[2] = p53ac;
    x_rdata[3] = MDM2;
    x_rdata[4] = Wip1;
    x_rdata[5] = ATMP;
    x_rdata[6] = ATRac;
    x_rdata[7] = MDM2product1;
    x_rdata[8] = MDM2product2;
    x_rdata[9] = MDM2product3;
    x_rdata[10] = MDM2product4;
    x_rdata[11] = MDM2product5;
    x_rdata[12] = MDM2product6;
    x_rdata[13] = MDM2product7;
    x_rdata[14] = MDM2product8;
    x_rdata[15] = MDM2product9;
    x_rdata[16] = MDM2pro;
    x_rdata[17] = Wip1product1;
    x_rdata[18] = Wip1product2;
    x_rdata[19] = Wip1product3;
    x_rdata[20] = Wip1product4;
    x_rdata[21] = Wip1product5;
    x_rdata[22] = Wip1product6;
    x_rdata[23] = Wip1product7;
    x_rdata[24] = Wip1product8;
    x_rdata[25] = Wip1product9;
    x_rdata[26] = Wip1pro;
    x_rdata[27] = BRCA2;
    x_rdata[28] = MSH6;
    x_rdata[29] = MGMT;
    x_rdata[30] = damageDSB;
    x_rdata[31] = damageSSB;
    x_rdata[32] = ppAKT_MDM2;
    x_rdata[33] = pMDM2;
    x_rdata[34] = ARF;
    x_rdata[35] = MDM4;
    x_rdata[36] = p53ac_MDM4;
    x_rdata[37] = ATMinac;
    x_rdata[38] = ATRinac;
    x_rdata[39] = pRB;
    x_rdata[40] = pRBp;
    x_rdata[41] = pRBpp;
    x_rdata[42] = E2F;
    x_rdata[43] = Cd;
    x_rdata[44] = Cd_Cdk46;
    x_rdata[45] = Cd_Cdk46_p27;
    x_rdata[46] = Ce;
    x_rdata[47] = Ce_Cdk2;
    x_rdata[48] = Skp2;
    x_rdata[49] = Ce_Cdk2_p27;
    x_rdata[50] = Pe;
    x_rdata[51] = Pai;
    x_rdata[52] = Pei;
    x_rdata[53] = Pbi;
    x_rdata[54] = Ca;
    x_rdata[55] = Ca_Cdk2;
    x_rdata[56] = Ca_Cdk2_p27;
    x_rdata[57] = p27;
    x_rdata[58] = Cdh1i;
    x_rdata[59] = Cdh1a;
    x_rdata[60] = E2Fp;
    x_rdata[61] = p27p;
    x_rdata[62] = Pa;
    x_rdata[63] = Cb;
    x_rdata[64] = Cb_Cdk1;
    x_rdata[65] = Cdc20i;
    x_rdata[66] = Cdc20a;
    x_rdata[67] = Pb;
    x_rdata[68] = Wee1;
    x_rdata[69] = Wee1p;
    x_rdata[70] = Cb_Cdk1_p27;
    x_rdata[71] = Chk1;
    x_rdata[72] = pRB_E2F;
    x_rdata[73] = pRBp_E2F;
    x_rdata[74] = p21;
    x_rdata[75] = Cd_Cdk46_p21;
    x_rdata[76] = Ce_Cdk2_p21;
    x_rdata[77] = Ca_Cdk2_p21;
    x_rdata[78] = Cb_Cdk1_p21;
    x_rdata[79] = L;
    x_rdata[80] = R;
    x_rdata[81] = L_R;
    x_rdata[82] = Ractive;
    x_rdata[83] = flip;
    x_rdata[84] = Ractive_flip;
    x_rdata[85] = pC8;
    x_rdata[86] = Ractive_pC8;
    x_rdata[87] = C8;
    x_rdata[88] = Bar;
    x_rdata[89] = C8_Bar;
    x_rdata[90] = pC3;
    x_rdata[91] = C8_pC3;
    x_rdata[92] = C3;
    x_rdata[93] = pC6;
    x_rdata[94] = C3_pC6;
    x_rdata[95] = C6;
    x_rdata[96] = C6_C8;
    x_rdata[97] = XIAP;
    x_rdata[98] = C3_XIAP;
    x_rdata[99] = PARP;
    x_rdata[100] = C3_PARP;
    x_rdata[101] = cPARP;
    x_rdata[102] = Bid;
    x_rdata[103] = C8_Bid;
    x_rdata[104] = tBid;
    x_rdata[105] = Bcl2c;
    x_rdata[106] = tBid_Bcl2c;
    x_rdata[107] = Bax;
    x_rdata[108] = tBid_Bax;
    x_rdata[109] = Baxactive;
    x_rdata[110] = Baxm;
    x_rdata[111] = Bcl2;
    x_rdata[112] = Baxm_Bcl2;
    x_rdata[113] = Bax2;
    x_rdata[114] = Bax2_Bcl2;
    x_rdata[115] = Bax4;
    x_rdata[116] = Bax4_Bcl2;
    x_rdata[117] = M;
    x_rdata[118] = Bax4_M;
    x_rdata[119] = Mactive;
    x_rdata[120] = CytoCm;
    x_rdata[121] = Mactive_CytoCm;
    x_rdata[122] = CytoCr;
    x_rdata[123] = Smacm;
    x_rdata[124] = Mactive_Smacm;
    x_rdata[125] = Smacr;
    x_rdata[126] = CytoC;
    x_rdata[127] = Apaf;
    x_rdata[128] = CytoC_Apaf;
    x_rdata[129] = Apafactive;
    x_rdata[130] = pC9;
    x_rdata[131] = Apop;
    x_rdata[132] = Apop_C3;
    x_rdata[133] = Smac;
    x_rdata[134] = Apop_XIAP;
    x_rdata[135] = Smac_XIAP;
    x_rdata[136] = C3_Ub;
    x_rdata[137] = BAD;
    x_rdata[138] = PUMA;
    x_rdata[139] = NOXA;
    x_rdata[140] = Bcl2c_BAD;
    x_rdata[141] = Bcl2c_PUMA;
    x_rdata[142] = Bcl2c_NOXA;
    x_rdata[143] = BIM;
    x_rdata[144] = BIM_Bax;
    x_rdata[145] = Bcl2c_BIM;
    x_rdata[146] = ppERK_BIM;
    x_rdata[147] = pBIM;
    x_rdata[148] = ppAKT_BAD;
    x_rdata[149] = pBAD;
    x_rdata[150] = ppERK_BAD;
    x_rdata[151] = E;
    x_rdata[152] = H;
    x_rdata[153] = HGF;
    x_rdata[154] = P;
    x_rdata[155] = F;
    x_rdata[156] = I;
    x_rdata[157] = INS;
    x_rdata[158] = EGFR;
    x_rdata[159] = pE1;
    x_rdata[160] = E2;
    x_rdata[161] = pE2;
    x_rdata[162] = E3;
    x_rdata[163] = E4;
    x_rdata[164] = pE4;
    x_rdata[165] = Ev3;
    x_rdata[166] = Met;
    x_rdata[167] = Pr;
    x_rdata[168] = Fr;
    x_rdata[169] = Ir;
    x_rdata[170] = Isr;
    x_rdata[171] = E1E1;
    x_rdata[172] = E1E2;
    x_rdata[173] = E1E3;
    x_rdata[174] = E1E4;
    x_rdata[175] = E2E2;
    x_rdata[176] = E2E3;
    x_rdata[177] = E2E4;
    x_rdata[178] = E3E4;
    x_rdata[179] = E4E4;
    x_rdata[180] = Met_Met;
    x_rdata[181] = FrFr;
    x_rdata[182] = IrIr;
    x_rdata[183] = Isr_Isr;
    x_rdata[184] = EE1;
    x_rdata[185] = HE3;
    x_rdata[186] = HE4;
    x_rdata[187] = HGF_Met;
    x_rdata[188] = PPr;
    x_rdata[189] = FFr;
    x_rdata[190] = EE1E2;
    x_rdata[191] = EE1Ev3;
    x_rdata[192] = EE1E1;
    x_rdata[193] = EE1E3;
    x_rdata[194] = EE1E4;
    x_rdata[195] = E2HE3;
    x_rdata[196] = E1HE3;
    x_rdata[197] = HE3E3;
    x_rdata[198] = HE3Ev3;
    x_rdata[199] = HE3E4;
    x_rdata[200] = E2HE4;
    x_rdata[201] = HE4Ev3;
    x_rdata[202] = E1HE4;
    x_rdata[203] = E3HE4;
    x_rdata[204] = HE4E4;
    x_rdata[205] = HGF_Met_Met;
    x_rdata[206] = PPrPr;
    x_rdata[207] = FFrFr;
    x_rdata[208] = IIrIr;
    x_rdata[209] = INS_Isr_Isr;
    x_rdata[210] = EE1EE1;
    x_rdata[211] = EE1HE3;
    x_rdata[212] = EE1HE4;
    x_rdata[213] = HE3HE3;
    x_rdata[214] = HE3HE4;
    x_rdata[215] = HE4HE4;
    x_rdata[216] = HGF_Met_HGF_Met;
    x_rdata[217] = PPrPPr;
    x_rdata[218] = FFrFFr;
    x_rdata[219] = IIrIrI;
    x_rdata[220] = INS_Isr_Isr_INS;
    x_rdata[221] = E1_ppERK;
    x_rdata[222] = E2_ppERK;
    x_rdata[223] = E4_ppERK;
    x_rdata[224] = pEE1E2;
    x_rdata[225] = pEE1Ev3;
    x_rdata[226] = pEE1E1;
    x_rdata[227] = pEE1EE1;
    x_rdata[228] = pEE1E3;
    x_rdata[229] = pEE1HE3;
    x_rdata[230] = pEE1E4;
    x_rdata[231] = pEE1HE4;
    x_rdata[232] = pE2HE3;
    x_rdata[233] = pHE3Ev3;
    x_rdata[234] = pE1HE3;
    x_rdata[235] = pHE3E4;
    x_rdata[236] = pHE3HE4;
    x_rdata[237] = pE2HE4;
    x_rdata[238] = pHE4Ev3;
    x_rdata[239] = pE1HE4;
    x_rdata[240] = pE3HE4;
    x_rdata[241] = pHE4E4;
    x_rdata[242] = pHE4HE4;
    x_rdata[243] = pHGF_Met_Met;
    x_rdata[244] = pHGF_Met_HGF_Met;
    x_rdata[245] = pPPrPPr;
    x_rdata[246] = pPPrPr;
    x_rdata[247] = pFFrFFr;
    x_rdata[248] = pFFrFr;
    x_rdata[249] = pIIrIr;
    x_rdata[250] = pINS_Isr_Isr;
    x_rdata[251] = pIIrIrI;
    x_rdata[252] = pINS_Isr_Isr_INS;
    x_rdata[253] = pIIrIr_IRS;
    x_rdata[254] = pINS_Isr_Isr_IRS;
    x_rdata[255] = pIIrIrI_IRS;
    x_rdata[256] = pINS_Isr_Isr_INS_IRS;
    x_rdata[257] = Sp_EE1E2;
    x_rdata[258] = Sp_EE1Ev3;
    x_rdata[259] = Sp_EE1E1;
    x_rdata[260] = Sp_EE1EE1;
    x_rdata[261] = Sp_EE1E3;
    x_rdata[262] = Sp_EE1HE3;
    x_rdata[263] = Sp_EE1E4;
    x_rdata[264] = Sp_EE1HE4;
    x_rdata[265] = Sp_E2HE3;
    x_rdata[266] = Sp_HE3Ev3;
    x_rdata[267] = Sp_E1HE3;
    x_rdata[268] = Sp_HE3E4;
    x_rdata[269] = Sp_HE3HE4;
    x_rdata[270] = Sp_E2HE4;
    x_rdata[271] = Sp_HE4Ev3;
    x_rdata[272] = Sp_E1HE4;
    x_rdata[273] = Sp_E3HE4;
    x_rdata[274] = Sp_HE4E4;
    x_rdata[275] = Sp_HE4HE4;
    x_rdata[276] = Sp_HGF_Met_Met;
    x_rdata[277] = Sp_HGF_Met_HGF_Met;
    x_rdata[278] = Sp_PPrPPr;
    x_rdata[279] = Sp_PPrPr;
    x_rdata[280] = Sp_FFrFFr;
    x_rdata[281] = Sp_FFrFr;
    x_rdata[282] = Sp_IIrIr;
    x_rdata[283] = Sp_INS_Isr_Isr;
    x_rdata[284] = Sp_IIrIrI;
    x_rdata[285] = Sp_INS_Isr_Isr_INS;
    x_rdata[286] = EE1E2int;
    x_rdata[287] = EE1Ev3int;
    x_rdata[288] = EE1E1int;
    x_rdata[289] = EE1EE1int;
    x_rdata[290] = EE1E3int;
    x_rdata[291] = EE1HE3int;
    x_rdata[292] = EE1E4int;
    x_rdata[293] = EE1HE4int;
    x_rdata[294] = E2HE3int;
    x_rdata[295] = HE3Ev3int;
    x_rdata[296] = E1HE3int;
    x_rdata[297] = HE3E4int;
    x_rdata[298] = HE3HE4int;
    x_rdata[299] = E2HE4int;
    x_rdata[300] = HE4Ev3int;
    x_rdata[301] = E1HE4int;
    x_rdata[302] = E3HE4int;
    x_rdata[303] = HE4E4int;
    x_rdata[304] = HE4HE4int;
    x_rdata[305] = HGF_Met_Metint;
    x_rdata[306] = HGF_Met_HGF_Metint;
    x_rdata[307] = PPrPPrint;
    x_rdata[308] = PPrPrint;
    x_rdata[309] = FFrFFrint;
    x_rdata[310] = FFrFrint;
    x_rdata[311] = IIrIr_int;
    x_rdata[312] = INS_Isr_Isr_int;
    x_rdata[313] = IIrIrI_int;
    x_rdata[314] = INS_Isr_Isr_INS_int;
    x_rdata[315] = pEE1E2int;
    x_rdata[316] = pEE1Ev3int;
    x_rdata[317] = pEE1E1int;
    x_rdata[318] = pEE1EE1int;
    x_rdata[319] = pEE1E3int;
    x_rdata[320] = pEE1HE3int;
    x_rdata[321] = pEE1E4int;
    x_rdata[322] = pEE1HE4int;
    x_rdata[323] = pE2HE3int;
    x_rdata[324] = pHE3Ev3int;
    x_rdata[325] = pE1HE3int;
    x_rdata[326] = pHE3E4int;
    x_rdata[327] = pHE3HE4int;
    x_rdata[328] = pE2HE4int;
    x_rdata[329] = pHE4Ev3int;
    x_rdata[330] = pE1HE4int;
    x_rdata[331] = pE3HE4int;
    x_rdata[332] = pHE4E4int;
    x_rdata[333] = pHE4HE4int;
    x_rdata[334] = pHGF_Met_Metint;
    x_rdata[335] = pHGF_Met_HGF_Metint;
    x_rdata[336] = pPPrPPrint;
    x_rdata[337] = pPPrPrint;
    x_rdata[338] = pFFrFFrint;
    x_rdata[339] = pFFrFrint;
    x_rdata[340] = pIIrIr_int;
    x_rdata[341] = pINS_Isr_Isr_int;
    x_rdata[342] = pIIrIrI_int;
    x_rdata[343] = pINS_Isr_Isr_INS_int;
    x_rdata[344] = pIIrIr_int_IRS;
    x_rdata[345] = pINS_Isr_Isr_int_IRS;
    x_rdata[346] = pIIrIrI_int_IRS;
    x_rdata[347] = pINS_Isr_Isr_INS_int_IRS;
    x_rdata[348] = pEE1E2_G2_SOS;
    x_rdata[349] = pEE1Ev3_G2_SOS;
    x_rdata[350] = pEE1E1_G2_SOS;
    x_rdata[351] = pEE1EE1_G2_SOS;
    x_rdata[352] = pEE1E3_G2_SOS;
    x_rdata[353] = pEE1HE3_G2_SOS;
    x_rdata[354] = pEE1E4_G2_SOS;
    x_rdata[355] = pEE1HE4_G2_SOS;
    x_rdata[356] = pE2HE3_G2_SOS;
    x_rdata[357] = pHE3Ev3_G2_SOS;
    x_rdata[358] = pE1HE3_G2_SOS;
    x_rdata[359] = pHE3E4_G2_SOS;
    x_rdata[360] = pHE3HE4_G2_SOS;
    x_rdata[361] = pE2HE4_G2_SOS;
    x_rdata[362] = pHE4Ev3_G2_SOS;
    x_rdata[363] = pE1HE4_G2_SOS;
    x_rdata[364] = pE3HE4_G2_SOS;
    x_rdata[365] = pHE4E4_G2_SOS;
    x_rdata[366] = pHE4HE4_G2_SOS;
    x_rdata[367] = pHGF_Met_Met_G2_SOS;
    x_rdata[368] = pHGF_Met_HGF_Met_G2_SOS;
    x_rdata[369] = pPPrPPr_G2_SOS;
    x_rdata[370] = pPPrPr_G2_SOS;
    x_rdata[371] = pFFrFFr_G2_SOS;
    x_rdata[372] = pFFrFr_G2_SOS;
    x_rdata[373] = pIIrIr_IRS_G2_SOS;
    x_rdata[374] = pINS_Isr_Isr_IRS_G2_SOS;
    x_rdata[375] = pIIrIrI_IRS_G2_SOS;
    x_rdata[376] = pINS_Isr_Isr_INS_IRS_G2_SOS;
    x_rdata[377] = pEE1E2int_G2_SOS;
    x_rdata[378] = pEE1Ev3int_G2_SOS;
    x_rdata[379] = pEE1E1int_G2_SOS;
    x_rdata[380] = pEE1EE1int_G2_SOS;
    x_rdata[381] = pEE1E3int_G2_SOS;
    x_rdata[382] = pEE1HE3int_G2_SOS;
    x_rdata[383] = pEE1E4int_G2_SOS;
    x_rdata[384] = pEE1HE4int_G2_SOS;
    x_rdata[385] = pE2HE3int_G2_SOS;
    x_rdata[386] = pHE3Ev3int_G2_SOS;
    x_rdata[387] = pE1HE3int_G2_SOS;
    x_rdata[388] = pHE3E4int_G2_SOS;
    x_rdata[389] = pHE3HE4int_G2_SOS;
    x_rdata[390] = pE2HE4int_G2_SOS;
    x_rdata[391] = pHE4Ev3int_G2_SOS;
    x_rdata[392] = pE1HE4int_G2_SOS;
    x_rdata[393] = pE3HE4int_G2_SOS;
    x_rdata[394] = pHE4E4int_G2_SOS;
    x_rdata[395] = pHE4HE4int_G2_SOS;
    x_rdata[396] = pHGF_Met_Metint_G2_SOS;
    x_rdata[397] = pHGF_Met_HGF_Metint_G2_SOS;
    x_rdata[398] = pPPrPPrint_G2_SOS;
    x_rdata[399] = pPPrPrint_G2_SOS;
    x_rdata[400] = pFFrFFrint_G2_SOS;
    x_rdata[401] = pFFrFrint_G2_SOS;
    x_rdata[402] = pIIrIr_int_IRS_G2_SOS;
    x_rdata[403] = pINS_Isr_Isr_int_IRS_G2_SOS;
    x_rdata[404] = pIIrIrI_int_IRS_G2_SOS;
    x_rdata[405] = pINS_Isr_Isr_INS_int_IRS_G2_SOS;
    x_rdata[406] = pEE1E2_PLCg;
    x_rdata[407] = pEE1Ev3_PLCg;
    x_rdata[408] = pEE1E1_PLCg;
    x_rdata[409] = pEE1EE1_PLCg;
    x_rdata[410] = pEE1E3_PLCg;
    x_rdata[411] = pEE1HE3_PLCg;
    x_rdata[412] = pEE1E4_PLCg;
    x_rdata[413] = pEE1HE4_PLCg;
    x_rdata[414] = pE2HE3_PLCg;
    x_rdata[415] = pHE3Ev3_PLCg;
    x_rdata[416] = pE1HE3_PLCg;
    x_rdata[417] = pHE3E4_PLCg;
    x_rdata[418] = pHE3HE4_PLCg;
    x_rdata[419] = pE2HE4_PLCg;
    x_rdata[420] = pHE4Ev3_PLCg;
    x_rdata[421] = pE1HE4_PLCg;
    x_rdata[422] = pE3HE4_PLCg;
    x_rdata[423] = pHE4E4_PLCg;
    x_rdata[424] = pHE4HE4_PLCg;
    x_rdata[425] = pHGF_Met_Met_PLCg;
    x_rdata[426] = pHGF_Met_HGF_Met_PLCg;
    x_rdata[427] = pPPrPPr_PLCg;
    x_rdata[428] = pPPrPr_PLCg;
    x_rdata[429] = pFFrFFr_PLCg;
    x_rdata[430] = pFFrFr_PLCg;
    x_rdata[431] = pIIrIr_IRS_PLCg;
    x_rdata[432] = pINS_Isr_Isr_IRS_PLCg;
    x_rdata[433] = pIIrIrI_IRS_PLCg;
    x_rdata[434] = pINS_Isr_Isr_INS_IRS_PLCg;
    x_rdata[435] = pEE1E2_PI3K1;
    x_rdata[436] = pEE1Ev3_PI3K1;
    x_rdata[437] = pEE1E1_PI3K1;
    x_rdata[438] = pEE1EE1_PI3K1;
    x_rdata[439] = pEE1E3_PI3K1;
    x_rdata[440] = pEE1HE3_PI3K1;
    x_rdata[441] = pEE1E4_PI3K1;
    x_rdata[442] = pEE1HE4_PI3K1;
    x_rdata[443] = pE2HE3_PI3K1;
    x_rdata[444] = pHE3Ev3_PI3K1;
    x_rdata[445] = pE1HE3_PI3K1;
    x_rdata[446] = pHE3E4_PI3K1;
    x_rdata[447] = pHE3HE4_PI3K1;
    x_rdata[448] = pE2HE4_PI3K1;
    x_rdata[449] = pHE4Ev3_PI3K1;
    x_rdata[450] = pE1HE4_PI3K1;
    x_rdata[451] = pE3HE4_PI3K1;
    x_rdata[452] = pHE4E4_PI3K1;
    x_rdata[453] = pHE4HE4_PI3K1;
    x_rdata[454] = pHGF_Met_Met_PI3K1;
    x_rdata[455] = pHGF_Met_HGF_Met_PI3K1;
    x_rdata[456] = pPPrPPr_PI3K1;
    x_rdata[457] = pPPrPr_PI3K1;
    x_rdata[458] = pFFrFFr_PI3K1;
    x_rdata[459] = pFFrFr_PI3K1;
    x_rdata[460] = pIIrIr_IRS_PI3K1;
    x_rdata[461] = pINS_Isr_Isr_IRS_PI3K1;
    x_rdata[462] = pIIrIrI_IRS_PI3K1;
    x_rdata[463] = pINS_Isr_Isr_INS_IRS_PI3K1;
    x_rdata[464] = pEE1E2_PI3K2;
    x_rdata[465] = pEE1Ev3_PI3K2;
    x_rdata[466] = pEE1E1_PI3K2;
    x_rdata[467] = pEE1EE1_PI3K2;
    x_rdata[468] = pEE1E3_PI3K2;
    x_rdata[469] = pEE1HE3_PI3K2;
    x_rdata[470] = pEE1E4_PI3K2;
    x_rdata[471] = pEE1HE4_PI3K2;
    x_rdata[472] = pE2HE3_PI3K2;
    x_rdata[473] = pHE3Ev3_PI3K2;
    x_rdata[474] = pE1HE3_PI3K2;
    x_rdata[475] = pHE3E4_PI3K2;
    x_rdata[476] = pHE3HE4_PI3K2;
    x_rdata[477] = pE2HE4_PI3K2;
    x_rdata[478] = pHE4Ev3_PI3K2;
    x_rdata[479] = pE1HE4_PI3K2;
    x_rdata[480] = pE3HE4_PI3K2;
    x_rdata[481] = pHE4E4_PI3K2;
    x_rdata[482] = pHE4HE4_PI3K2;
    x_rdata[483] = pHGF_Met_Met_PI3K2;
    x_rdata[484] = pHGF_Met_HGF_Met_PI3K2;
    x_rdata[485] = pPPrPPr_PI3K2;
    x_rdata[486] = pPPrPr_PI3K2;
    x_rdata[487] = pFFrFFr_PI3K2;
    x_rdata[488] = pFFrFr_PI3K2;
    x_rdata[489] = pIIrIr_IRS_PI3K2;
    x_rdata[490] = pINS_Isr_Isr_IRS_PI3K2;
    x_rdata[491] = pIIrIrI_IRS_PI3K2;
    x_rdata[492] = pINS_Isr_Isr_INS_IRS_PI3K2;
    x_rdata[493] = pEE1E2int_G2_SOS_RasD;
    x_rdata[494] = pEE1Ev3int_G2_SOS_RasD;
    x_rdata[495] = pEE1E1int_G2_SOS_RasD;
    x_rdata[496] = pEE1EE1int_G2_SOS_RasD;
    x_rdata[497] = pEE1E3int_G2_SOS_RasD;
    x_rdata[498] = pEE1HE3int_G2_SOS_RasD;
    x_rdata[499] = pEE1E4int_G2_SOS_RasD;
    x_rdata[500] = pEE1HE4int_G2_SOS_RasD;
    x_rdata[501] = pE2HE3int_G2_SOS_RasD;
    x_rdata[502] = pHE3Ev3int_G2_SOS_RasD;
    x_rdata[503] = pE1HE3int_G2_SOS_RasD;
    x_rdata[504] = pHE3E4int_G2_SOS_RasD;
    x_rdata[505] = pHE3HE4int_G2_SOS_RasD;
    x_rdata[506] = pE2HE4int_G2_SOS_RasD;
    x_rdata[507] = pHE4Ev3int_G2_SOS_RasD;
    x_rdata[508] = pE1HE4int_G2_SOS_RasD;
    x_rdata[509] = pE3HE4int_G2_SOS_RasD;
    x_rdata[510] = pHE4E4int_G2_SOS_RasD;
    x_rdata[511] = pHE4HE4int_G2_SOS_RasD;
    x_rdata[512] = pHGF_Met_Metint_G2_SOS_RasD;
    x_rdata[513] = pHGF_Met_HGF_Metint_G2_SOS_RasD;
    x_rdata[514] = pPPrPPrint_G2_SOS_RasD;
    x_rdata[515] = pPPrPrint_G2_SOS_RasD;
    x_rdata[516] = pFFrFFrint_G2_SOS_RasD;
    x_rdata[517] = pFFrFrint_G2_SOS_RasD;
    x_rdata[518] = pIIrIr_int_IRS_G2_SOS_RasD;
    x_rdata[519] = pINS_Isr_Isr_int_IRS_G2_SOS_RasD;
    x_rdata[520] = pIIrIrI_int_IRS_G2_SOS_RasD;
    x_rdata[521] = pINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD;
    x_rdata[522] = pEE1E2_G2_SOS_RasD;
    x_rdata[523] = pEE1Ev3_G2_SOS_RasD;
    x_rdata[524] = pEE1E1_G2_SOS_RasD;
    x_rdata[525] = pEE1EE1_G2_SOS_RasD;
    x_rdata[526] = pEE1E3_G2_SOS_RasD;
    x_rdata[527] = pEE1HE3_G2_SOS_RasD;
    x_rdata[528] = pEE1E4_G2_SOS_RasD;
    x_rdata[529] = pEE1HE4_G2_SOS_RasD;
    x_rdata[530] = pE2HE3_G2_SOS_RasD;
    x_rdata[531] = pHE3Ev3_G2_SOS_RasD;
    x_rdata[532] = pE1HE3_G2_SOS_RasD;
    x_rdata[533] = pHE3E4_G2_SOS_RasD;
    x_rdata[534] = pHE3HE4_G2_SOS_RasD;
    x_rdata[535] = pE2HE4_G2_SOS_RasD;
    x_rdata[536] = pHE4Ev3_G2_SOS_RasD;
    x_rdata[537] = pE1HE4_G2_SOS_RasD;
    x_rdata[538] = pE3HE4_G2_SOS_RasD;
    x_rdata[539] = pHE4E4_G2_SOS_RasD;
    x_rdata[540] = pHE4HE4_G2_SOS_RasD;
    x_rdata[541] = pHGF_Met_Met_G2_SOS_RasD;
    x_rdata[542] = pHGF_Met_HGF_Met_G2_SOS_RasD;
    x_rdata[543] = pPPrPPr_G2_SOS_RasD;
    x_rdata[544] = pPPrPr_G2_SOS_RasD;
    x_rdata[545] = pFFrFFr_G2_SOS_RasD;
    x_rdata[546] = pFFrFr_G2_SOS_RasD;
    x_rdata[547] = pIIrIr_IRS_G2_SOS_RasD;
    x_rdata[548] = pINS_Isr_Isr_IRS_G2_SOS_RasD;
    x_rdata[549] = pIIrIrI_IRS_G2_SOS_RasD;
    x_rdata[550] = pINS_Isr_Isr_INS_IRS_G2_SOS_RasD;
    x_rdata[551] = pEE1E2_PLCg_PIP2;
    x_rdata[552] = pEE1Ev3_PLCg_PIP2;
    x_rdata[553] = pEE1E1_PLCg_PIP2;
    x_rdata[554] = pEE1EE1_PLCg_PIP2;
    x_rdata[555] = pEE1E3_PLCg_PIP2;
    x_rdata[556] = pEE1HE3_PLCg_PIP2;
    x_rdata[557] = pEE1E4_PLCg_PIP2;
    x_rdata[558] = pEE1HE4_PLCg_PIP2;
    x_rdata[559] = pE2HE3_PLCg_PIP2;
    x_rdata[560] = pHE3Ev3_PLCg_PIP2;
    x_rdata[561] = pE1HE3_PLCg_PIP2;
    x_rdata[562] = pHE3E4_PLCg_PIP2;
    x_rdata[563] = pHE3HE4_PLCg_PIP2;
    x_rdata[564] = pE2HE4_PLCg_PIP2;
    x_rdata[565] = pHE4Ev3_PLCg_PIP2;
    x_rdata[566] = pE1HE4_PLCg_PIP2;
    x_rdata[567] = pE3HE4_PLCg_PIP2;
    x_rdata[568] = pHE4E4_PLCg_PIP2;
    x_rdata[569] = pHE4HE4_PLCg_PIP2;
    x_rdata[570] = pHGF_Met_Met_PLCg_PIP2;
    x_rdata[571] = pHGF_Met_HGF_Met_PLCg_PIP2;
    x_rdata[572] = pPPrPPr_PLCg_PIP2;
    x_rdata[573] = pPPrPr_PLCg_PIP2;
    x_rdata[574] = pFFrFFr_PLCg_PIP2;
    x_rdata[575] = pFFrFr_PLCg_PIP2;
    x_rdata[576] = pIIrIr_IRS_PLCg_PIP2;
    x_rdata[577] = pINS_Isr_Isr_IRS_PLCg_PIP2;
    x_rdata[578] = pIIrIrI_IRS_PLCg_PIP2;
    x_rdata[579] = pINS_Isr_Isr_INS_IRS_PLCg_PIP2;
    x_rdata[580] = pEE1E2_PI3K1_PIP2;
    x_rdata[581] = pEE1Ev3_PI3K1_PIP2;
    x_rdata[582] = pEE1E1_PI3K1_PIP2;
    x_rdata[583] = pEE1EE1_PI3K1_PIP2;
    x_rdata[584] = pEE1E3_PI3K1_PIP2;
    x_rdata[585] = pEE1HE3_PI3K1_PIP2;
    x_rdata[586] = pEE1E4_PI3K1_PIP2;
    x_rdata[587] = pEE1HE4_PI3K1_PIP2;
    x_rdata[588] = pE2HE3_PI3K1_PIP2;
    x_rdata[589] = pHE3Ev3_PI3K1_PIP2;
    x_rdata[590] = pE1HE3_PI3K1_PIP2;
    x_rdata[591] = pHE3E4_PI3K1_PIP2;
    x_rdata[592] = pHE3HE4_PI3K1_PIP2;
    x_rdata[593] = pE2HE4_PI3K1_PIP2;
    x_rdata[594] = pHE4Ev3_PI3K1_PIP2;
    x_rdata[595] = pE1HE4_PI3K1_PIP2;
    x_rdata[596] = pE3HE4_PI3K1_PIP2;
    x_rdata[597] = pHE4E4_PI3K1_PIP2;
    x_rdata[598] = pHE4HE4_PI3K1_PIP2;
    x_rdata[599] = pHGF_Met_Met_PI3K1_PIP2;
    x_rdata[600] = pHGF_Met_HGF_Met_PI3K1_PIP2;
    x_rdata[601] = pPPrPPr_PI3K1_PIP2;
    x_rdata[602] = pPPrPr_PI3K1_PIP2;
    x_rdata[603] = pFFrFFr_PI3K1_PIP2;
    x_rdata[604] = pFFrFr_PI3K1_PIP2;
    x_rdata[605] = pIIrIr_IRS_PI3K1_PIP2;
    x_rdata[606] = pINS_Isr_Isr_IRS_PI3K1_PIP2;
    x_rdata[607] = pIIrIrI_IRS_PI3K1_PIP2;
    x_rdata[608] = pINS_Isr_Isr_INS_IRS_PI3K1_PIP2;
    x_rdata[609] = pEE1E2_PI3K2_PIP;
    x_rdata[610] = pEE1Ev3_PI3K2_PIP;
    x_rdata[611] = pEE1E1_PI3K2_PIP;
    x_rdata[612] = pEE1EE1_PI3K2_PIP;
    x_rdata[613] = pEE1E3_PI3K2_PIP;
    x_rdata[614] = pEE1HE3_PI3K2_PIP;
    x_rdata[615] = pEE1E4_PI3K2_PIP;
    x_rdata[616] = pEE1HE4_PI3K2_PIP;
    x_rdata[617] = pE2HE3_PI3K2_PIP;
    x_rdata[618] = pHE3Ev3_PI3K2_PIP;
    x_rdata[619] = pE1HE3_PI3K2_PIP;
    x_rdata[620] = pHE3E4_PI3K2_PIP;
    x_rdata[621] = pHE3HE4_PI3K2_PIP;
    x_rdata[622] = pE2HE4_PI3K2_PIP;
    x_rdata[623] = pHE4Ev3_PI3K2_PIP;
    x_rdata[624] = pE1HE4_PI3K2_PIP;
    x_rdata[625] = pE3HE4_PI3K2_PIP;
    x_rdata[626] = pHE4E4_PI3K2_PIP;
    x_rdata[627] = pHE4HE4_PI3K2_PIP;
    x_rdata[628] = pHGF_Met_Met_PI3K2_PIP;
    x_rdata[629] = pHGF_Met_HGF_Met_PI3K2_PIP;
    x_rdata[630] = pPPrPPr_PI3K2_PIP;
    x_rdata[631] = pPPrPr_PI3K2_PIP;
    x_rdata[632] = pFFrFFr_PI3K2_PIP;
    x_rdata[633] = pFFrFr_PI3K2_PIP;
    x_rdata[634] = pIIrIr_IRS_PI3K2_PIP;
    x_rdata[635] = pINS_Isr_Isr_IRS_PI3K2_PIP;
    x_rdata[636] = pIIrIrI_IRS_PI3K2_PIP;
    x_rdata[637] = pINS_Isr_Isr_INS_IRS_PI3K2_PIP;
    x_rdata[638] = IRS;
    x_rdata[639] = Sp;
    x_rdata[640] = Cbl;
    x_rdata[641] = G2;
    x_rdata[642] = G2_SOS;
    x_rdata[643] = G2_pSOS;
    x_rdata[644] = PLCg;
    x_rdata[645] = PI3KC1;
    x_rdata[646] = PI3KR1;
    x_rdata[647] = PI3K1;
    x_rdata[648] = pPI3K1;
    x_rdata[649] = PI3K2;
    x_rdata[650] = mTORC1;
    x_rdata[651] = mTORC1active;
    x_rdata[652] = PIP;
    x_rdata[653] = PI3P;
    x_rdata[654] = DAG;
    x_rdata[655] = GRP;
    x_rdata[656] = DAG_GRP;
    x_rdata[657] = RasT;
    x_rdata[658] = RasD;
    x_rdata[659] = NF1;
    x_rdata[660] = pNF1;
    x_rdata[661] = pCRaf;
    x_rdata[662] = CRaf;
    x_rdata[663] = RasT_CRaf;
    x_rdata[664] = BRaf;
    x_rdata[665] = RasT_CRaf_BRaf;
    x_rdata[666] = MEK;
    x_rdata[667] = pMEK;
    x_rdata[668] = ppMEK;
    x_rdata[669] = MKP3;
    x_rdata[670] = ERKnuc;
    x_rdata[671] = ppERKnuc;
    x_rdata[672] = RSK;
    x_rdata[673] = pRSK;
    x_rdata[674] = pRSKnuc;
    x_rdata[675] = MKP1;
    x_rdata[676] = pMKP1;
    x_rdata[677] = cFos;
    x_rdata[678] = pcFos;
    x_rdata[679] = cJun;
    x_rdata[680] = pcFos_cJun;
    x_rdata[681] = cMyc;
    x_rdata[682] = bCATENINnuc;
    x_rdata[683] = bCATENIN;
    x_rdata[684] = pbCATENIN;
    x_rdata[685] = IP3;
    x_rdata[686] = PIP2;
    x_rdata[687] = PIP3;
    x_rdata[688] = PTEN;
    x_rdata[689] = PIP3_AKT;
    x_rdata[690] = AKT;
    x_rdata[691] = pAKT;
    x_rdata[692] = ppAKT;
    x_rdata[693] = PDK1;
    x_rdata[694] = PIP3_PDK1;
    x_rdata[695] = PIP3_pAKT;
    x_rdata[696] = Rictor;
    x_rdata[697] = mTOR;
    x_rdata[698] = mTORC2;
    x_rdata[699] = PIP3_ppAKT;
    x_rdata[700] = GSK3b;
    x_rdata[701] = pGSK3b;
    x_rdata[702] = TSC1;
    x_rdata[703] = TSC2;
    x_rdata[704] = pTSC2;
    x_rdata[705] = TSC;
    x_rdata[706] = PKC;
    x_rdata[707] = DAG_PKC;
    x_rdata[708] = pRKIP;
    x_rdata[709] = RKIP;
    x_rdata[710] = RKIP_CRaf;
    x_rdata[711] = ERK;
    x_rdata[712] = pERK;
    x_rdata[713] = ppERK;
    x_rdata[714] = FOXO;
    x_rdata[715] = pFOXO;
    x_rdata[716] = RhebD;
    x_rdata[717] = RhebT;
    x_rdata[718] = Raptor;
    x_rdata[719] = S6K;
    x_rdata[720] = pS6K;
    x_rdata[721] = EIF4EBP1;
    x_rdata[722] = pEIF4EBP1;
    x_rdata[723] = SOS;
    x_rdata[724] = G2_SOS_ppERK;
    x_rdata[725] = CRaf_ppERK;
    x_rdata[726] = RasD_DAG_GRP;
    x_rdata[727] = RasT_NF1;
    x_rdata[728] = NF1_ppERK;
    x_rdata[729] = MEK_RasT_CRaf_BRaf;
    x_rdata[730] = pMEK_RasT_CRaf_BRaf;
    x_rdata[731] = ERK_ppMEK;
    x_rdata[732] = pERK_ppMEK;
    x_rdata[733] = RSK_ppERK;
    x_rdata[734] = pRSKnuc_MKP1;
    x_rdata[735] = ppERKnuc_MKP1;
    x_rdata[736] = cFos_pRSKnuc;
    x_rdata[737] = cFos_ppERKnuc;
    x_rdata[738] = RKIP_DAG_PKC;
    x_rdata[739] = PIP3_PTEN;
    x_rdata[740] = PIP3_AKT_PIP3_PDK1;
    x_rdata[741] = PIP3_pAKT_mTORC2;
    x_rdata[742] = GSK3b_ppAKT;
    x_rdata[743] = TSC2_ppAKT;
    x_rdata[744] = TSC2_ppERK;
    x_rdata[745] = RhebT_TSC;
    x_rdata[746] = EIF4EBP1_mTORC1active;
    x_rdata[747] = S6K_mTORC1active;
    x_rdata[748] = FOXO_ppAKT;
    x_rdata[749] = PI3K1_mTORC1active;
    x_rdata[750] = pERK_MKP3;
    x_rdata[751] = ppERK_MKP3;
    x_rdata[752] = ppERKnuc_pMKP1;
    x_rdata[753] = RasT_BRaf;
    x_rdata[754] = RasT_BRaf_BRaf;
    x_rdata[755] = MEK_RasT_BRaf_BRaf;
    x_rdata[756] = pMEK_RasT_BRaf_BRaf;
    x_rdata[757] = EIF4E;
    x_rdata[758] = EIF4EBP1_EIF4E;
    x_rdata[759] = RasT_CRaf_CRaf;
    x_rdata[760] = MEK_RasT_CRaf_CRaf;
    x_rdata[761] = pMEK_RasT_CRaf_CRaf;
    x_rdata[762] = FOXOnuc;
    x_rdata[763] = MEKi;
    x_rdata[764] = MEKi_ppMEK;
    x_rdata[765] = AKTi;
    x_rdata[766] = AKTi_AKT;
    x_rdata[767] = mT;
    x_rdata[768] = EIF4E_mT;
    x_rdata[769] = Cdk1;
    x_rdata[770] = Cdk2;
    x_rdata[771] = Cdk46;
    x_rdata[772] = pRBpp_E2F;
    x_rdata[773] = Cd_Cdk46_pRB;
    x_rdata[774] = Cd_Cdk46_pRB_E2F;
    x_rdata[775] = Ce_Cdk2_pRBp;
    x_rdata[776] = Ce_Cdk2_pRBp_E2F;
    x_rdata[777] = p18;
    x_rdata[778] = p19;
    x_rdata[779] = p57;
    x_rdata[780] = Cd_Cdk46_p18;
    x_rdata[781] = Cd_Cdk46_p19;
    x_rdata[782] = Cd_Cdk46_p57;
    x_rdata[783] = Ce_Cdk2_p57;
    x_rdata[784] = Ca_Cdk2_p57;
    x_rdata[785] = Cb_Cdk1_p57;
}

} // namespace amici
} // namespace model_SPARCED