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
    x_rdata[45] = Md;
    x_rdata[46] = Cd_Cdk46_p27;
    x_rdata[47] = Ce;
    x_rdata[48] = Ce_Cdk2;
    x_rdata[49] = Me;
    x_rdata[50] = Skp2;
    x_rdata[51] = Ce_Cdk2_p27;
    x_rdata[52] = Pe;
    x_rdata[53] = Pai;
    x_rdata[54] = Pei;
    x_rdata[55] = Pbi;
    x_rdata[56] = Ca;
    x_rdata[57] = Ca_Cdk2;
    x_rdata[58] = Ma;
    x_rdata[59] = Ca_Cdk2_p27;
    x_rdata[60] = p27;
    x_rdata[61] = Cdh1i;
    x_rdata[62] = Cdh1a;
    x_rdata[63] = E2Fp;
    x_rdata[64] = p27p;
    x_rdata[65] = Pa;
    x_rdata[66] = Cb;
    x_rdata[67] = Cb_Cdk1;
    x_rdata[68] = Mb;
    x_rdata[69] = Cdc20i;
    x_rdata[70] = Cdc20a;
    x_rdata[71] = Pb;
    x_rdata[72] = Wee1;
    x_rdata[73] = Wee1p;
    x_rdata[74] = Cb_Cdk1_p27;
    x_rdata[75] = Chk1;
    x_rdata[76] = pRB_E2F;
    x_rdata[77] = pRBp_E2F;
    x_rdata[78] = p21;
    x_rdata[79] = Cd_Cdk46_p21;
    x_rdata[80] = Ce_Cdk2_p21;
    x_rdata[81] = Ca_Cdk2_p21;
    x_rdata[82] = Cb_Cdk1_p21;
    x_rdata[83] = L;
    x_rdata[84] = R;
    x_rdata[85] = L_R;
    x_rdata[86] = Ractive;
    x_rdata[87] = flip;
    x_rdata[88] = Ractive_flip;
    x_rdata[89] = pC8;
    x_rdata[90] = Ractive_pC8;
    x_rdata[91] = C8;
    x_rdata[92] = Bar;
    x_rdata[93] = C8_Bar;
    x_rdata[94] = pC3;
    x_rdata[95] = C8_pC3;
    x_rdata[96] = C3;
    x_rdata[97] = pC6;
    x_rdata[98] = C3_pC6;
    x_rdata[99] = C6;
    x_rdata[100] = C6_C8;
    x_rdata[101] = XIAP;
    x_rdata[102] = C3_XIAP;
    x_rdata[103] = PARP;
    x_rdata[104] = C3_PARP;
    x_rdata[105] = cPARP;
    x_rdata[106] = Bid;
    x_rdata[107] = C8_Bid;
    x_rdata[108] = tBid;
    x_rdata[109] = Bcl2c;
    x_rdata[110] = tBid_Bcl2c;
    x_rdata[111] = Bax;
    x_rdata[112] = tBid_Bax;
    x_rdata[113] = Baxactive;
    x_rdata[114] = Baxm;
    x_rdata[115] = Bcl2;
    x_rdata[116] = Baxm_Bcl2;
    x_rdata[117] = Bax2;
    x_rdata[118] = Bax2_Bcl2;
    x_rdata[119] = Bax4;
    x_rdata[120] = Bax4_Bcl2;
    x_rdata[121] = M;
    x_rdata[122] = Bax4_M;
    x_rdata[123] = Mactive;
    x_rdata[124] = CytoCm;
    x_rdata[125] = Mactive_CytoCm;
    x_rdata[126] = CytoCr;
    x_rdata[127] = Smacm;
    x_rdata[128] = Mactive_Smacm;
    x_rdata[129] = Smacr;
    x_rdata[130] = CytoC;
    x_rdata[131] = Apaf;
    x_rdata[132] = CytoC_Apaf;
    x_rdata[133] = Apafactive;
    x_rdata[134] = pC9;
    x_rdata[135] = Apop;
    x_rdata[136] = Apop_C3;
    x_rdata[137] = Smac;
    x_rdata[138] = Apop_XIAP;
    x_rdata[139] = Smac_XIAP;
    x_rdata[140] = C3_Ub;
    x_rdata[141] = BAD;
    x_rdata[142] = PUMA;
    x_rdata[143] = NOXA;
    x_rdata[144] = Bcl2c_BAD;
    x_rdata[145] = Bcl2c_PUMA;
    x_rdata[146] = Bcl2c_NOXA;
    x_rdata[147] = BIM;
    x_rdata[148] = BIM_Bax;
    x_rdata[149] = Bcl2c_BIM;
    x_rdata[150] = ppERK_BIM;
    x_rdata[151] = pBIM;
    x_rdata[152] = ppAKT_BAD;
    x_rdata[153] = pBAD;
    x_rdata[154] = ppERK_BAD;
    x_rdata[155] = E;
    x_rdata[156] = H;
    x_rdata[157] = HGF;
    x_rdata[158] = P;
    x_rdata[159] = F;
    x_rdata[160] = I;
    x_rdata[161] = INS;
    x_rdata[162] = EGFR;
    x_rdata[163] = pE1;
    x_rdata[164] = E2;
    x_rdata[165] = pE2;
    x_rdata[166] = E3;
    x_rdata[167] = E4;
    x_rdata[168] = pE4;
    x_rdata[169] = Ev3;
    x_rdata[170] = Met;
    x_rdata[171] = Pr;
    x_rdata[172] = Fr;
    x_rdata[173] = Ir;
    x_rdata[174] = Isr;
    x_rdata[175] = E1E1;
    x_rdata[176] = E1E2;
    x_rdata[177] = E1E3;
    x_rdata[178] = E1E4;
    x_rdata[179] = E2E2;
    x_rdata[180] = E2E3;
    x_rdata[181] = E2E4;
    x_rdata[182] = E3E4;
    x_rdata[183] = E4E4;
    x_rdata[184] = Met_Met;
    x_rdata[185] = FrFr;
    x_rdata[186] = IrIr;
    x_rdata[187] = Isr_Isr;
    x_rdata[188] = EE1;
    x_rdata[189] = HE3;
    x_rdata[190] = HE4;
    x_rdata[191] = HGF_Met;
    x_rdata[192] = PPr;
    x_rdata[193] = FFr;
    x_rdata[194] = EE1E2;
    x_rdata[195] = EE1Ev3;
    x_rdata[196] = EE1E1;
    x_rdata[197] = EE1E3;
    x_rdata[198] = EE1E4;
    x_rdata[199] = E2HE3;
    x_rdata[200] = E1HE3;
    x_rdata[201] = HE3E3;
    x_rdata[202] = HE3Ev3;
    x_rdata[203] = HE3E4;
    x_rdata[204] = E2HE4;
    x_rdata[205] = HE4Ev3;
    x_rdata[206] = E1HE4;
    x_rdata[207] = E3HE4;
    x_rdata[208] = HE4E4;
    x_rdata[209] = HGF_Met_Met;
    x_rdata[210] = PPrPr;
    x_rdata[211] = FFrFr;
    x_rdata[212] = IIrIr;
    x_rdata[213] = INS_Isr_Isr;
    x_rdata[214] = EE1EE1;
    x_rdata[215] = EE1HE3;
    x_rdata[216] = EE1HE4;
    x_rdata[217] = HE3HE3;
    x_rdata[218] = HE3HE4;
    x_rdata[219] = HE4HE4;
    x_rdata[220] = HGF_Met_HGF_Met;
    x_rdata[221] = PPrPPr;
    x_rdata[222] = FFrFFr;
    x_rdata[223] = IIrIrI;
    x_rdata[224] = INS_Isr_Isr_INS;
    x_rdata[225] = E1_ppERK;
    x_rdata[226] = E2_ppERK;
    x_rdata[227] = E4_ppERK;
    x_rdata[228] = pEE1E2;
    x_rdata[229] = pEE1Ev3;
    x_rdata[230] = pEE1E1;
    x_rdata[231] = pEE1EE1;
    x_rdata[232] = pEE1E3;
    x_rdata[233] = pEE1HE3;
    x_rdata[234] = pEE1E4;
    x_rdata[235] = pEE1HE4;
    x_rdata[236] = pE2HE3;
    x_rdata[237] = pHE3Ev3;
    x_rdata[238] = pE1HE3;
    x_rdata[239] = pHE3E4;
    x_rdata[240] = pHE3HE4;
    x_rdata[241] = pE2HE4;
    x_rdata[242] = pHE4Ev3;
    x_rdata[243] = pE1HE4;
    x_rdata[244] = pE3HE4;
    x_rdata[245] = pHE4E4;
    x_rdata[246] = pHE4HE4;
    x_rdata[247] = pHGF_Met_Met;
    x_rdata[248] = pHGF_Met_HGF_Met;
    x_rdata[249] = pPPrPPr;
    x_rdata[250] = pPPrPr;
    x_rdata[251] = pFFrFFr;
    x_rdata[252] = pFFrFr;
    x_rdata[253] = pIIrIr;
    x_rdata[254] = pINS_Isr_Isr;
    x_rdata[255] = pIIrIrI;
    x_rdata[256] = pINS_Isr_Isr_INS;
    x_rdata[257] = pIIrIr_IRS;
    x_rdata[258] = pINS_Isr_Isr_IRS;
    x_rdata[259] = pIIrIrI_IRS;
    x_rdata[260] = pINS_Isr_Isr_INS_IRS;
    x_rdata[261] = Sp_EE1E2;
    x_rdata[262] = Sp_EE1Ev3;
    x_rdata[263] = Sp_EE1E1;
    x_rdata[264] = Sp_EE1EE1;
    x_rdata[265] = Sp_EE1E3;
    x_rdata[266] = Sp_EE1HE3;
    x_rdata[267] = Sp_EE1E4;
    x_rdata[268] = Sp_EE1HE4;
    x_rdata[269] = Sp_E2HE3;
    x_rdata[270] = Sp_HE3Ev3;
    x_rdata[271] = Sp_E1HE3;
    x_rdata[272] = Sp_HE3E4;
    x_rdata[273] = Sp_HE3HE4;
    x_rdata[274] = Sp_E2HE4;
    x_rdata[275] = Sp_HE4Ev3;
    x_rdata[276] = Sp_E1HE4;
    x_rdata[277] = Sp_E3HE4;
    x_rdata[278] = Sp_HE4E4;
    x_rdata[279] = Sp_HE4HE4;
    x_rdata[280] = Sp_HGF_Met_Met;
    x_rdata[281] = Sp_HGF_Met_HGF_Met;
    x_rdata[282] = Sp_PPrPPr;
    x_rdata[283] = Sp_PPrPr;
    x_rdata[284] = Sp_FFrFFr;
    x_rdata[285] = Sp_FFrFr;
    x_rdata[286] = Sp_IIrIr;
    x_rdata[287] = Sp_INS_Isr_Isr;
    x_rdata[288] = Sp_IIrIrI;
    x_rdata[289] = Sp_INS_Isr_Isr_INS;
    x_rdata[290] = EE1E2int;
    x_rdata[291] = EE1Ev3int;
    x_rdata[292] = EE1E1int;
    x_rdata[293] = EE1EE1int;
    x_rdata[294] = EE1E3int;
    x_rdata[295] = EE1HE3int;
    x_rdata[296] = EE1E4int;
    x_rdata[297] = EE1HE4int;
    x_rdata[298] = E2HE3int;
    x_rdata[299] = HE3Ev3int;
    x_rdata[300] = E1HE3int;
    x_rdata[301] = HE3E4int;
    x_rdata[302] = HE3HE4int;
    x_rdata[303] = E2HE4int;
    x_rdata[304] = HE4Ev3int;
    x_rdata[305] = E1HE4int;
    x_rdata[306] = E3HE4int;
    x_rdata[307] = HE4E4int;
    x_rdata[308] = HE4HE4int;
    x_rdata[309] = HGF_Met_Metint;
    x_rdata[310] = HGF_Met_HGF_Metint;
    x_rdata[311] = PPrPPrint;
    x_rdata[312] = PPrPrint;
    x_rdata[313] = FFrFFrint;
    x_rdata[314] = FFrFrint;
    x_rdata[315] = IIrIr_int;
    x_rdata[316] = INS_Isr_Isr_int;
    x_rdata[317] = IIrIrI_int;
    x_rdata[318] = INS_Isr_Isr_INS_int;
    x_rdata[319] = pEE1E2int;
    x_rdata[320] = pEE1Ev3int;
    x_rdata[321] = pEE1E1int;
    x_rdata[322] = pEE1EE1int;
    x_rdata[323] = pEE1E3int;
    x_rdata[324] = pEE1HE3int;
    x_rdata[325] = pEE1E4int;
    x_rdata[326] = pEE1HE4int;
    x_rdata[327] = pE2HE3int;
    x_rdata[328] = pHE3Ev3int;
    x_rdata[329] = pE1HE3int;
    x_rdata[330] = pHE3E4int;
    x_rdata[331] = pHE3HE4int;
    x_rdata[332] = pE2HE4int;
    x_rdata[333] = pHE4Ev3int;
    x_rdata[334] = pE1HE4int;
    x_rdata[335] = pE3HE4int;
    x_rdata[336] = pHE4E4int;
    x_rdata[337] = pHE4HE4int;
    x_rdata[338] = pHGF_Met_Metint;
    x_rdata[339] = pHGF_Met_HGF_Metint;
    x_rdata[340] = pPPrPPrint;
    x_rdata[341] = pPPrPrint;
    x_rdata[342] = pFFrFFrint;
    x_rdata[343] = pFFrFrint;
    x_rdata[344] = pIIrIr_int;
    x_rdata[345] = pINS_Isr_Isr_int;
    x_rdata[346] = pIIrIrI_int;
    x_rdata[347] = pINS_Isr_Isr_INS_int;
    x_rdata[348] = pIIrIr_int_IRS;
    x_rdata[349] = pINS_Isr_Isr_int_IRS;
    x_rdata[350] = pIIrIrI_int_IRS;
    x_rdata[351] = pINS_Isr_Isr_INS_int_IRS;
    x_rdata[352] = pEE1E2_G2_SOS;
    x_rdata[353] = pEE1Ev3_G2_SOS;
    x_rdata[354] = pEE1E1_G2_SOS;
    x_rdata[355] = pEE1EE1_G2_SOS;
    x_rdata[356] = pEE1E3_G2_SOS;
    x_rdata[357] = pEE1HE3_G2_SOS;
    x_rdata[358] = pEE1E4_G2_SOS;
    x_rdata[359] = pEE1HE4_G2_SOS;
    x_rdata[360] = pE2HE3_G2_SOS;
    x_rdata[361] = pHE3Ev3_G2_SOS;
    x_rdata[362] = pE1HE3_G2_SOS;
    x_rdata[363] = pHE3E4_G2_SOS;
    x_rdata[364] = pHE3HE4_G2_SOS;
    x_rdata[365] = pE2HE4_G2_SOS;
    x_rdata[366] = pHE4Ev3_G2_SOS;
    x_rdata[367] = pE1HE4_G2_SOS;
    x_rdata[368] = pE3HE4_G2_SOS;
    x_rdata[369] = pHE4E4_G2_SOS;
    x_rdata[370] = pHE4HE4_G2_SOS;
    x_rdata[371] = pHGF_Met_Met_G2_SOS;
    x_rdata[372] = pHGF_Met_HGF_Met_G2_SOS;
    x_rdata[373] = pPPrPPr_G2_SOS;
    x_rdata[374] = pPPrPr_G2_SOS;
    x_rdata[375] = pFFrFFr_G2_SOS;
    x_rdata[376] = pFFrFr_G2_SOS;
    x_rdata[377] = pIIrIr_IRS_G2_SOS;
    x_rdata[378] = pINS_Isr_Isr_IRS_G2_SOS;
    x_rdata[379] = pIIrIrI_IRS_G2_SOS;
    x_rdata[380] = pINS_Isr_Isr_INS_IRS_G2_SOS;
    x_rdata[381] = pEE1E2int_G2_SOS;
    x_rdata[382] = pEE1Ev3int_G2_SOS;
    x_rdata[383] = pEE1E1int_G2_SOS;
    x_rdata[384] = pEE1EE1int_G2_SOS;
    x_rdata[385] = pEE1E3int_G2_SOS;
    x_rdata[386] = pEE1HE3int_G2_SOS;
    x_rdata[387] = pEE1E4int_G2_SOS;
    x_rdata[388] = pEE1HE4int_G2_SOS;
    x_rdata[389] = pE2HE3int_G2_SOS;
    x_rdata[390] = pHE3Ev3int_G2_SOS;
    x_rdata[391] = pE1HE3int_G2_SOS;
    x_rdata[392] = pHE3E4int_G2_SOS;
    x_rdata[393] = pHE3HE4int_G2_SOS;
    x_rdata[394] = pE2HE4int_G2_SOS;
    x_rdata[395] = pHE4Ev3int_G2_SOS;
    x_rdata[396] = pE1HE4int_G2_SOS;
    x_rdata[397] = pE3HE4int_G2_SOS;
    x_rdata[398] = pHE4E4int_G2_SOS;
    x_rdata[399] = pHE4HE4int_G2_SOS;
    x_rdata[400] = pHGF_Met_Metint_G2_SOS;
    x_rdata[401] = pHGF_Met_HGF_Metint_G2_SOS;
    x_rdata[402] = pPPrPPrint_G2_SOS;
    x_rdata[403] = pPPrPrint_G2_SOS;
    x_rdata[404] = pFFrFFrint_G2_SOS;
    x_rdata[405] = pFFrFrint_G2_SOS;
    x_rdata[406] = pIIrIr_int_IRS_G2_SOS;
    x_rdata[407] = pINS_Isr_Isr_int_IRS_G2_SOS;
    x_rdata[408] = pIIrIrI_int_IRS_G2_SOS;
    x_rdata[409] = pINS_Isr_Isr_INS_int_IRS_G2_SOS;
    x_rdata[410] = pEE1E2_PLCg;
    x_rdata[411] = pEE1Ev3_PLCg;
    x_rdata[412] = pEE1E1_PLCg;
    x_rdata[413] = pEE1EE1_PLCg;
    x_rdata[414] = pEE1E3_PLCg;
    x_rdata[415] = pEE1HE3_PLCg;
    x_rdata[416] = pEE1E4_PLCg;
    x_rdata[417] = pEE1HE4_PLCg;
    x_rdata[418] = pE2HE3_PLCg;
    x_rdata[419] = pHE3Ev3_PLCg;
    x_rdata[420] = pE1HE3_PLCg;
    x_rdata[421] = pHE3E4_PLCg;
    x_rdata[422] = pHE3HE4_PLCg;
    x_rdata[423] = pE2HE4_PLCg;
    x_rdata[424] = pHE4Ev3_PLCg;
    x_rdata[425] = pE1HE4_PLCg;
    x_rdata[426] = pE3HE4_PLCg;
    x_rdata[427] = pHE4E4_PLCg;
    x_rdata[428] = pHE4HE4_PLCg;
    x_rdata[429] = pHGF_Met_Met_PLCg;
    x_rdata[430] = pHGF_Met_HGF_Met_PLCg;
    x_rdata[431] = pPPrPPr_PLCg;
    x_rdata[432] = pPPrPr_PLCg;
    x_rdata[433] = pFFrFFr_PLCg;
    x_rdata[434] = pFFrFr_PLCg;
    x_rdata[435] = pIIrIr_IRS_PLCg;
    x_rdata[436] = pINS_Isr_Isr_IRS_PLCg;
    x_rdata[437] = pIIrIrI_IRS_PLCg;
    x_rdata[438] = pINS_Isr_Isr_INS_IRS_PLCg;
    x_rdata[439] = pEE1E2_PI3K1;
    x_rdata[440] = pEE1Ev3_PI3K1;
    x_rdata[441] = pEE1E1_PI3K1;
    x_rdata[442] = pEE1EE1_PI3K1;
    x_rdata[443] = pEE1E3_PI3K1;
    x_rdata[444] = pEE1HE3_PI3K1;
    x_rdata[445] = pEE1E4_PI3K1;
    x_rdata[446] = pEE1HE4_PI3K1;
    x_rdata[447] = pE2HE3_PI3K1;
    x_rdata[448] = pHE3Ev3_PI3K1;
    x_rdata[449] = pE1HE3_PI3K1;
    x_rdata[450] = pHE3E4_PI3K1;
    x_rdata[451] = pHE3HE4_PI3K1;
    x_rdata[452] = pE2HE4_PI3K1;
    x_rdata[453] = pHE4Ev3_PI3K1;
    x_rdata[454] = pE1HE4_PI3K1;
    x_rdata[455] = pE3HE4_PI3K1;
    x_rdata[456] = pHE4E4_PI3K1;
    x_rdata[457] = pHE4HE4_PI3K1;
    x_rdata[458] = pHGF_Met_Met_PI3K1;
    x_rdata[459] = pHGF_Met_HGF_Met_PI3K1;
    x_rdata[460] = pPPrPPr_PI3K1;
    x_rdata[461] = pPPrPr_PI3K1;
    x_rdata[462] = pFFrFFr_PI3K1;
    x_rdata[463] = pFFrFr_PI3K1;
    x_rdata[464] = pIIrIr_IRS_PI3K1;
    x_rdata[465] = pINS_Isr_Isr_IRS_PI3K1;
    x_rdata[466] = pIIrIrI_IRS_PI3K1;
    x_rdata[467] = pINS_Isr_Isr_INS_IRS_PI3K1;
    x_rdata[468] = pEE1E2_PI3K2;
    x_rdata[469] = pEE1Ev3_PI3K2;
    x_rdata[470] = pEE1E1_PI3K2;
    x_rdata[471] = pEE1EE1_PI3K2;
    x_rdata[472] = pEE1E3_PI3K2;
    x_rdata[473] = pEE1HE3_PI3K2;
    x_rdata[474] = pEE1E4_PI3K2;
    x_rdata[475] = pEE1HE4_PI3K2;
    x_rdata[476] = pE2HE3_PI3K2;
    x_rdata[477] = pHE3Ev3_PI3K2;
    x_rdata[478] = pE1HE3_PI3K2;
    x_rdata[479] = pHE3E4_PI3K2;
    x_rdata[480] = pHE3HE4_PI3K2;
    x_rdata[481] = pE2HE4_PI3K2;
    x_rdata[482] = pHE4Ev3_PI3K2;
    x_rdata[483] = pE1HE4_PI3K2;
    x_rdata[484] = pE3HE4_PI3K2;
    x_rdata[485] = pHE4E4_PI3K2;
    x_rdata[486] = pHE4HE4_PI3K2;
    x_rdata[487] = pHGF_Met_Met_PI3K2;
    x_rdata[488] = pHGF_Met_HGF_Met_PI3K2;
    x_rdata[489] = pPPrPPr_PI3K2;
    x_rdata[490] = pPPrPr_PI3K2;
    x_rdata[491] = pFFrFFr_PI3K2;
    x_rdata[492] = pFFrFr_PI3K2;
    x_rdata[493] = pIIrIr_IRS_PI3K2;
    x_rdata[494] = pINS_Isr_Isr_IRS_PI3K2;
    x_rdata[495] = pIIrIrI_IRS_PI3K2;
    x_rdata[496] = pINS_Isr_Isr_INS_IRS_PI3K2;
    x_rdata[497] = pEE1E2int_G2_SOS_RasD;
    x_rdata[498] = pEE1Ev3int_G2_SOS_RasD;
    x_rdata[499] = pEE1E1int_G2_SOS_RasD;
    x_rdata[500] = pEE1EE1int_G2_SOS_RasD;
    x_rdata[501] = pEE1E3int_G2_SOS_RasD;
    x_rdata[502] = pEE1HE3int_G2_SOS_RasD;
    x_rdata[503] = pEE1E4int_G2_SOS_RasD;
    x_rdata[504] = pEE1HE4int_G2_SOS_RasD;
    x_rdata[505] = pE2HE3int_G2_SOS_RasD;
    x_rdata[506] = pHE3Ev3int_G2_SOS_RasD;
    x_rdata[507] = pE1HE3int_G2_SOS_RasD;
    x_rdata[508] = pHE3E4int_G2_SOS_RasD;
    x_rdata[509] = pHE3HE4int_G2_SOS_RasD;
    x_rdata[510] = pE2HE4int_G2_SOS_RasD;
    x_rdata[511] = pHE4Ev3int_G2_SOS_RasD;
    x_rdata[512] = pE1HE4int_G2_SOS_RasD;
    x_rdata[513] = pE3HE4int_G2_SOS_RasD;
    x_rdata[514] = pHE4E4int_G2_SOS_RasD;
    x_rdata[515] = pHE4HE4int_G2_SOS_RasD;
    x_rdata[516] = pHGF_Met_Metint_G2_SOS_RasD;
    x_rdata[517] = pHGF_Met_HGF_Metint_G2_SOS_RasD;
    x_rdata[518] = pPPrPPrint_G2_SOS_RasD;
    x_rdata[519] = pPPrPrint_G2_SOS_RasD;
    x_rdata[520] = pFFrFFrint_G2_SOS_RasD;
    x_rdata[521] = pFFrFrint_G2_SOS_RasD;
    x_rdata[522] = pIIrIr_int_IRS_G2_SOS_RasD;
    x_rdata[523] = pINS_Isr_Isr_int_IRS_G2_SOS_RasD;
    x_rdata[524] = pIIrIrI_int_IRS_G2_SOS_RasD;
    x_rdata[525] = pINS_Isr_Isr_INS_int_IRS_G2_SOS_RasD;
    x_rdata[526] = pEE1E2_G2_SOS_RasD;
    x_rdata[527] = pEE1Ev3_G2_SOS_RasD;
    x_rdata[528] = pEE1E1_G2_SOS_RasD;
    x_rdata[529] = pEE1EE1_G2_SOS_RasD;
    x_rdata[530] = pEE1E3_G2_SOS_RasD;
    x_rdata[531] = pEE1HE3_G2_SOS_RasD;
    x_rdata[532] = pEE1E4_G2_SOS_RasD;
    x_rdata[533] = pEE1HE4_G2_SOS_RasD;
    x_rdata[534] = pE2HE3_G2_SOS_RasD;
    x_rdata[535] = pHE3Ev3_G2_SOS_RasD;
    x_rdata[536] = pE1HE3_G2_SOS_RasD;
    x_rdata[537] = pHE3E4_G2_SOS_RasD;
    x_rdata[538] = pHE3HE4_G2_SOS_RasD;
    x_rdata[539] = pE2HE4_G2_SOS_RasD;
    x_rdata[540] = pHE4Ev3_G2_SOS_RasD;
    x_rdata[541] = pE1HE4_G2_SOS_RasD;
    x_rdata[542] = pE3HE4_G2_SOS_RasD;
    x_rdata[543] = pHE4E4_G2_SOS_RasD;
    x_rdata[544] = pHE4HE4_G2_SOS_RasD;
    x_rdata[545] = pHGF_Met_Met_G2_SOS_RasD;
    x_rdata[546] = pHGF_Met_HGF_Met_G2_SOS_RasD;
    x_rdata[547] = pPPrPPr_G2_SOS_RasD;
    x_rdata[548] = pPPrPr_G2_SOS_RasD;
    x_rdata[549] = pFFrFFr_G2_SOS_RasD;
    x_rdata[550] = pFFrFr_G2_SOS_RasD;
    x_rdata[551] = pIIrIr_IRS_G2_SOS_RasD;
    x_rdata[552] = pINS_Isr_Isr_IRS_G2_SOS_RasD;
    x_rdata[553] = pIIrIrI_IRS_G2_SOS_RasD;
    x_rdata[554] = pINS_Isr_Isr_INS_IRS_G2_SOS_RasD;
    x_rdata[555] = pEE1E2_PLCg_PIP2;
    x_rdata[556] = pEE1Ev3_PLCg_PIP2;
    x_rdata[557] = pEE1E1_PLCg_PIP2;
    x_rdata[558] = pEE1EE1_PLCg_PIP2;
    x_rdata[559] = pEE1E3_PLCg_PIP2;
    x_rdata[560] = pEE1HE3_PLCg_PIP2;
    x_rdata[561] = pEE1E4_PLCg_PIP2;
    x_rdata[562] = pEE1HE4_PLCg_PIP2;
    x_rdata[563] = pE2HE3_PLCg_PIP2;
    x_rdata[564] = pHE3Ev3_PLCg_PIP2;
    x_rdata[565] = pE1HE3_PLCg_PIP2;
    x_rdata[566] = pHE3E4_PLCg_PIP2;
    x_rdata[567] = pHE3HE4_PLCg_PIP2;
    x_rdata[568] = pE2HE4_PLCg_PIP2;
    x_rdata[569] = pHE4Ev3_PLCg_PIP2;
    x_rdata[570] = pE1HE4_PLCg_PIP2;
    x_rdata[571] = pE3HE4_PLCg_PIP2;
    x_rdata[572] = pHE4E4_PLCg_PIP2;
    x_rdata[573] = pHE4HE4_PLCg_PIP2;
    x_rdata[574] = pHGF_Met_Met_PLCg_PIP2;
    x_rdata[575] = pHGF_Met_HGF_Met_PLCg_PIP2;
    x_rdata[576] = pPPrPPr_PLCg_PIP2;
    x_rdata[577] = pPPrPr_PLCg_PIP2;
    x_rdata[578] = pFFrFFr_PLCg_PIP2;
    x_rdata[579] = pFFrFr_PLCg_PIP2;
    x_rdata[580] = pIIrIr_IRS_PLCg_PIP2;
    x_rdata[581] = pINS_Isr_Isr_IRS_PLCg_PIP2;
    x_rdata[582] = pIIrIrI_IRS_PLCg_PIP2;
    x_rdata[583] = pINS_Isr_Isr_INS_IRS_PLCg_PIP2;
    x_rdata[584] = pEE1E2_PI3K1_PIP2;
    x_rdata[585] = pEE1Ev3_PI3K1_PIP2;
    x_rdata[586] = pEE1E1_PI3K1_PIP2;
    x_rdata[587] = pEE1EE1_PI3K1_PIP2;
    x_rdata[588] = pEE1E3_PI3K1_PIP2;
    x_rdata[589] = pEE1HE3_PI3K1_PIP2;
    x_rdata[590] = pEE1E4_PI3K1_PIP2;
    x_rdata[591] = pEE1HE4_PI3K1_PIP2;
    x_rdata[592] = pE2HE3_PI3K1_PIP2;
    x_rdata[593] = pHE3Ev3_PI3K1_PIP2;
    x_rdata[594] = pE1HE3_PI3K1_PIP2;
    x_rdata[595] = pHE3E4_PI3K1_PIP2;
    x_rdata[596] = pHE3HE4_PI3K1_PIP2;
    x_rdata[597] = pE2HE4_PI3K1_PIP2;
    x_rdata[598] = pHE4Ev3_PI3K1_PIP2;
    x_rdata[599] = pE1HE4_PI3K1_PIP2;
    x_rdata[600] = pE3HE4_PI3K1_PIP2;
    x_rdata[601] = pHE4E4_PI3K1_PIP2;
    x_rdata[602] = pHE4HE4_PI3K1_PIP2;
    x_rdata[603] = pHGF_Met_Met_PI3K1_PIP2;
    x_rdata[604] = pHGF_Met_HGF_Met_PI3K1_PIP2;
    x_rdata[605] = pPPrPPr_PI3K1_PIP2;
    x_rdata[606] = pPPrPr_PI3K1_PIP2;
    x_rdata[607] = pFFrFFr_PI3K1_PIP2;
    x_rdata[608] = pFFrFr_PI3K1_PIP2;
    x_rdata[609] = pIIrIr_IRS_PI3K1_PIP2;
    x_rdata[610] = pINS_Isr_Isr_IRS_PI3K1_PIP2;
    x_rdata[611] = pIIrIrI_IRS_PI3K1_PIP2;
    x_rdata[612] = pINS_Isr_Isr_INS_IRS_PI3K1_PIP2;
    x_rdata[613] = pEE1E2_PI3K2_PIP;
    x_rdata[614] = pEE1Ev3_PI3K2_PIP;
    x_rdata[615] = pEE1E1_PI3K2_PIP;
    x_rdata[616] = pEE1EE1_PI3K2_PIP;
    x_rdata[617] = pEE1E3_PI3K2_PIP;
    x_rdata[618] = pEE1HE3_PI3K2_PIP;
    x_rdata[619] = pEE1E4_PI3K2_PIP;
    x_rdata[620] = pEE1HE4_PI3K2_PIP;
    x_rdata[621] = pE2HE3_PI3K2_PIP;
    x_rdata[622] = pHE3Ev3_PI3K2_PIP;
    x_rdata[623] = pE1HE3_PI3K2_PIP;
    x_rdata[624] = pHE3E4_PI3K2_PIP;
    x_rdata[625] = pHE3HE4_PI3K2_PIP;
    x_rdata[626] = pE2HE4_PI3K2_PIP;
    x_rdata[627] = pHE4Ev3_PI3K2_PIP;
    x_rdata[628] = pE1HE4_PI3K2_PIP;
    x_rdata[629] = pE3HE4_PI3K2_PIP;
    x_rdata[630] = pHE4E4_PI3K2_PIP;
    x_rdata[631] = pHE4HE4_PI3K2_PIP;
    x_rdata[632] = pHGF_Met_Met_PI3K2_PIP;
    x_rdata[633] = pHGF_Met_HGF_Met_PI3K2_PIP;
    x_rdata[634] = pPPrPPr_PI3K2_PIP;
    x_rdata[635] = pPPrPr_PI3K2_PIP;
    x_rdata[636] = pFFrFFr_PI3K2_PIP;
    x_rdata[637] = pFFrFr_PI3K2_PIP;
    x_rdata[638] = pIIrIr_IRS_PI3K2_PIP;
    x_rdata[639] = pINS_Isr_Isr_IRS_PI3K2_PIP;
    x_rdata[640] = pIIrIrI_IRS_PI3K2_PIP;
    x_rdata[641] = pINS_Isr_Isr_INS_IRS_PI3K2_PIP;
    x_rdata[642] = IRS;
    x_rdata[643] = Sp;
    x_rdata[644] = Cbl;
    x_rdata[645] = G2;
    x_rdata[646] = G2_SOS;
    x_rdata[647] = G2_pSOS;
    x_rdata[648] = PLCg;
    x_rdata[649] = PI3KC1;
    x_rdata[650] = PI3KR1;
    x_rdata[651] = PI3K1;
    x_rdata[652] = pPI3K1;
    x_rdata[653] = PI3K2;
    x_rdata[654] = mTORC1;
    x_rdata[655] = mTORC1active;
    x_rdata[656] = PIP;
    x_rdata[657] = PI3P;
    x_rdata[658] = DAG;
    x_rdata[659] = GRP;
    x_rdata[660] = DAG_GRP;
    x_rdata[661] = RasT;
    x_rdata[662] = RasD;
    x_rdata[663] = NF1;
    x_rdata[664] = pNF1;
    x_rdata[665] = pCRaf;
    x_rdata[666] = CRaf;
    x_rdata[667] = RasT_CRaf;
    x_rdata[668] = BRaf;
    x_rdata[669] = RasT_CRaf_BRaf;
    x_rdata[670] = MEK;
    x_rdata[671] = pMEK;
    x_rdata[672] = ppMEK;
    x_rdata[673] = MKP3;
    x_rdata[674] = ERKnuc;
    x_rdata[675] = ppERKnuc;
    x_rdata[676] = RSK;
    x_rdata[677] = pRSK;
    x_rdata[678] = pRSKnuc;
    x_rdata[679] = MKP1;
    x_rdata[680] = pMKP1;
    x_rdata[681] = cFos;
    x_rdata[682] = pcFos;
    x_rdata[683] = cJun;
    x_rdata[684] = pcFos_cJun;
    x_rdata[685] = cMyc;
    x_rdata[686] = bCATENINnuc;
    x_rdata[687] = bCATENIN;
    x_rdata[688] = pbCATENIN;
    x_rdata[689] = IP3;
    x_rdata[690] = PIP2;
    x_rdata[691] = PIP3;
    x_rdata[692] = PTEN;
    x_rdata[693] = PIP3_AKT;
    x_rdata[694] = AKT;
    x_rdata[695] = pAKT;
    x_rdata[696] = ppAKT;
    x_rdata[697] = PDK1;
    x_rdata[698] = PIP3_PDK1;
    x_rdata[699] = PIP3_pAKT;
    x_rdata[700] = Rictor;
    x_rdata[701] = mTOR;
    x_rdata[702] = mTORC2;
    x_rdata[703] = PIP3_ppAKT;
    x_rdata[704] = GSK3b;
    x_rdata[705] = pGSK3b;
    x_rdata[706] = TSC1;
    x_rdata[707] = TSC2;
    x_rdata[708] = pTSC2;
    x_rdata[709] = TSC;
    x_rdata[710] = PKC;
    x_rdata[711] = DAG_PKC;
    x_rdata[712] = pRKIP;
    x_rdata[713] = RKIP;
    x_rdata[714] = RKIP_CRaf;
    x_rdata[715] = ERK;
    x_rdata[716] = pERK;
    x_rdata[717] = ppERK;
    x_rdata[718] = FOXO;
    x_rdata[719] = pFOXO;
    x_rdata[720] = RhebD;
    x_rdata[721] = RhebT;
    x_rdata[722] = Raptor;
    x_rdata[723] = S6K;
    x_rdata[724] = pS6K;
    x_rdata[725] = EIF4EBP1;
    x_rdata[726] = pEIF4EBP1;
    x_rdata[727] = SOS;
    x_rdata[728] = G2_SOS_ppERK;
    x_rdata[729] = CRaf_ppERK;
    x_rdata[730] = RasD_DAG_GRP;
    x_rdata[731] = RasT_NF1;
    x_rdata[732] = NF1_ppERK;
    x_rdata[733] = MEK_RasT_CRaf_BRaf;
    x_rdata[734] = pMEK_RasT_CRaf_BRaf;
    x_rdata[735] = ERK_ppMEK;
    x_rdata[736] = pERK_ppMEK;
    x_rdata[737] = RSK_ppERK;
    x_rdata[738] = pRSKnuc_MKP1;
    x_rdata[739] = ppERKnuc_MKP1;
    x_rdata[740] = cFos_pRSKnuc;
    x_rdata[741] = cFos_ppERKnuc;
    x_rdata[742] = RKIP_DAG_PKC;
    x_rdata[743] = PIP3_PTEN;
    x_rdata[744] = PIP3_AKT_PIP3_PDK1;
    x_rdata[745] = PIP3_pAKT_mTORC2;
    x_rdata[746] = GSK3b_ppAKT;
    x_rdata[747] = TSC2_ppAKT;
    x_rdata[748] = TSC2_ppERK;
    x_rdata[749] = RhebT_TSC;
    x_rdata[750] = EIF4EBP1_mTORC1active;
    x_rdata[751] = S6K_mTORC1active;
    x_rdata[752] = FOXO_ppAKT;
    x_rdata[753] = PI3K1_mTORC1active;
    x_rdata[754] = pERK_MKP3;
    x_rdata[755] = ppERK_MKP3;
    x_rdata[756] = ppERKnuc_pMKP1;
    x_rdata[757] = RasT_BRaf;
    x_rdata[758] = RasT_BRaf_BRaf;
    x_rdata[759] = MEK_RasT_BRaf_BRaf;
    x_rdata[760] = pMEK_RasT_BRaf_BRaf;
    x_rdata[761] = EIF4E;
    x_rdata[762] = EIF4EBP1_EIF4E;
    x_rdata[763] = RasT_CRaf_CRaf;
    x_rdata[764] = MEK_RasT_CRaf_CRaf;
    x_rdata[765] = pMEK_RasT_CRaf_CRaf;
    x_rdata[766] = FOXOnuc;
    x_rdata[767] = MEKi;
    x_rdata[768] = MEKi_ppMEK;
    x_rdata[769] = AKTi;
    x_rdata[770] = AKTi_AKT;
    x_rdata[771] = mT;
    x_rdata[772] = EIF4E_mT;
    x_rdata[773] = Cdk1;
    x_rdata[774] = Cdk2;
    x_rdata[775] = Cdk46;
    x_rdata[776] = tcl_pRBpp_E2F;
    x_rdata[777] = tcl_Cd_Cdk46_pRB;
    x_rdata[778] = tcl_Cd_Cdk46_pRB_E2F;
    x_rdata[779] = tcl_Ce_Cdk2_pRBp;
    x_rdata[780] = tcl_Ce_Cdk2_pRBp_E2F;
    x_rdata[781] = tcl_Cd_Cdk46_pRBp;
    x_rdata[782] = tcl_Cd_Cdk46_pRBp_E2F;
    x_rdata[783] = tcl_Ce_Cdk2_pRBpp;
    x_rdata[784] = tcl_Ce_Cdk2_pRBpp_E2F;
}

} // namespace amici
} // namespace model_SPARCED