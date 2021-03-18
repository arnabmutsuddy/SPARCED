#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 2443> dxdotdw_colptrs_SPARCED_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 268, 271, 272, 273, 274, 275, 276, 278, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 339, 342, 343, 344, 345, 347, 350, 353, 356, 359, 362, 365, 368, 371, 374, 377, 380, 383, 386, 389, 392, 395, 398, 401, 404, 407, 410, 413, 416, 419, 422, 425, 428, 431, 434, 437, 440, 443, 446, 449, 452, 455, 458, 461, 464, 467, 470, 473, 476, 479, 482, 485, 488, 491, 494, 497, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 540, 543, 545, 548, 551, 554, 557, 560, 563, 566, 569, 572, 575, 578, 581, 584, 587, 590, 593, 596, 599, 602, 605, 608, 611, 614, 617, 620, 623, 626, 629, 632, 635, 637, 639, 642, 645, 647, 649, 652, 655, 657, 659, 662, 665, 668, 671, 674, 677, 680, 683, 686, 689, 692, 694, 696, 699, 702, 705, 708, 711, 714, 717, 720, 722, 724, 727, 730, 733, 736, 739, 742, 745, 748, 751, 754, 757, 760, 763, 766, 769, 771, 773, 775, 777, 779, 781, 784, 787, 790, 793, 796, 799, 802, 805, 807, 809, 812, 815, 818, 821, 824, 827, 830, 833, 836, 839, 842, 845, 848, 851, 854, 857, 860, 863, 866, 869, 872, 875, 878, 881, 884, 887, 889, 891, 894, 897, 900, 903, 906, 909, 912, 915, 918, 921, 924, 927, 930, 933, 936, 939, 942, 945, 948, 951, 954, 957, 960, 963, 966, 969, 971, 973, 976, 978, 981, 983, 986, 988, 991, 994, 997, 1000, 1003, 1006, 1008, 1010, 1013, 1016, 1018, 1020, 1023, 1026, 1029, 1032, 1034, 1036, 1039, 1042, 1045, 1048, 1050, 1052, 1054, 1056, 1059, 1062, 1065, 1068, 1071, 1074, 1076, 1078, 1081, 1084, 1087, 1090, 1093, 1096, 1098, 1100, 1103, 1106, 1109, 1112, 1115, 1118, 1120, 1122, 1125, 1128, 1131, 1134, 1137, 1140, 1143, 1146, 1149, 1152, 1155, 1158, 1161, 1164, 1167, 1170, 1173, 1176, 1179, 1182, 1184, 1186, 1189, 1192, 1194, 1196, 1199, 1202, 1205, 1208, 1211, 1214, 1217, 1220, 1223, 1226, 1228, 1230, 1232, 1234, 1236, 1238, 1240, 1242, 1244, 1246, 1248, 1250, 1252, 1254, 1256, 1258, 1260, 1262, 1264, 1266, 1268, 1270, 1272, 1274, 1276, 1278, 1280, 1282, 1284, 1286, 1288, 1290, 1292, 1294, 1296, 1298, 1300, 1302, 1304, 1306, 1308, 1310, 1312, 1314, 1316, 1318, 1320, 1322, 1324, 1326, 1328, 1330, 1332, 1334, 1336, 1338, 1340, 1342, 1345, 1348, 1351, 1354, 1357, 1360, 1363, 1366, 1369, 1372, 1375, 1378, 1381, 1384, 1387, 1390, 1393, 1396, 1399, 1402, 1405, 1408, 1411, 1414, 1417, 1420, 1423, 1426, 1429, 1432, 1435, 1438, 1441, 1444, 1447, 1450, 1453, 1456, 1459, 1462, 1465, 1468, 1471, 1474, 1477, 1480, 1483, 1486, 1489, 1492, 1495, 1498, 1501, 1504, 1507, 1510, 1513, 1516, 1518, 1520, 1522, 1524, 1526, 1528, 1530, 1532, 1534, 1536, 1538, 1540, 1542, 1544, 1546, 1548, 1550, 1552, 1554, 1556, 1558, 1560, 1562, 1564, 1566, 1568, 1570, 1572, 1574, 1576, 1578, 1580, 1582, 1584, 1586, 1588, 1590, 1592, 1594, 1596, 1598, 1600, 1602, 1604, 1606, 1608, 1610, 1612, 1614, 1616, 1618, 1620, 1622, 1624, 1626, 1628, 1630, 1632, 1634, 1636, 1638, 1640, 1642, 1644, 1646, 1648, 1650, 1652, 1654, 1656, 1658, 1660, 1662, 1664, 1666, 1668, 1670, 1672, 1674, 1676, 1678, 1680, 1682, 1684, 1686, 1688, 1690, 1692, 1694, 1696, 1698, 1700, 1702, 1704, 1706, 1708, 1710, 1712, 1714, 1716, 1718, 1720, 1722, 1724, 1726, 1728, 1730, 1732, 1734, 1736, 1738, 1740, 1742, 1744, 1746, 1748, 1751, 1754, 1757, 1760, 1763, 1766, 1769, 1772, 1775, 1778, 1781, 1784, 1787, 1790, 1793, 1796, 1799, 1802, 1805, 1808, 1811, 1814, 1817, 1820, 1823, 1826, 1829, 1832, 1835, 1838, 1841, 1844, 1847, 1850, 1853, 1856, 1859, 1862, 1865, 1868, 1871, 1874, 1877, 1880, 1883, 1886, 1889, 1892, 1895, 1898, 1901, 1904, 1907, 1910, 1913, 1916, 1919, 1922, 1925, 1928, 1931, 1934, 1937, 1940, 1943, 1946, 1949, 1952, 1955, 1958, 1961, 1964, 1967, 1970, 1973, 1976, 1979, 1982, 1985, 1988, 1991, 1994, 1997, 2000, 2003, 2006, 2009, 2012, 2015, 2018, 2021, 2024, 2027, 2030, 2033, 2036, 2039, 2042, 2045, 2048, 2051, 2054, 2057, 2060, 2063, 2066, 2069, 2072, 2075, 2078, 2081, 2084, 2087, 2090, 2093, 2096, 2099, 2102, 2105, 2108, 2111, 2114, 2117, 2120, 2123, 2126, 2129, 2132, 2135, 2138, 2141, 2144, 2147, 2150, 2153, 2156, 2159, 2162, 2165, 2168, 2171, 2174, 2177, 2180, 2183, 2186, 2189, 2192, 2195, 2198, 2201, 2204, 2207, 2210, 2213, 2216, 2219, 2222, 2225, 2228, 2231, 2234, 2237, 2240, 2243, 2246, 2249, 2252, 2255, 2258, 2261, 2264, 2267, 2270, 2273, 2276, 2279, 2282, 2285, 2288, 2291, 2294, 2297, 2300, 2303, 2306, 2309, 2312, 2315, 2318, 2321, 2324, 2327, 2330, 2333, 2336, 2339, 2342, 2345, 2348, 2351, 2354, 2357, 2360, 2363, 2366, 2369, 2372, 2375, 2378, 2381, 2384, 2387, 2390, 2393, 2396, 2399, 2402, 2405, 2408, 2411, 2414, 2417, 2420, 2423, 2426, 2429, 2432, 2435, 2438, 2441, 2444, 2447, 2450, 2453, 2456, 2459, 2462, 2465, 2468, 2471, 2474, 2477, 2480, 2483, 2486, 2489, 2492, 2495, 2498, 2501, 2504, 2507, 2510, 2513, 2516, 2519, 2522, 2525, 2528, 2531, 2534, 2537, 2540, 2543, 2546, 2549, 2552, 2555, 2558, 2561, 2564, 2567, 2570, 2573, 2576, 2579, 2582, 2585, 2588, 2591, 2594, 2597, 2600, 2603, 2606, 2609, 2612, 2615, 2618, 2621, 2624, 2627, 2630, 2633, 2636, 2639, 2642, 2645, 2648, 2651, 2654, 2657, 2660, 2663, 2666, 2669, 2672, 2675, 2678, 2681, 2684, 2687, 2690, 2693, 2696, 2699, 2702, 2705, 2708, 2711, 2714, 2717, 2720, 2723, 2726, 2729, 2732, 2735, 2738, 2741, 2744, 2747, 2750, 2753, 2756, 2759, 2762, 2765, 2768, 2771, 2774, 2777, 2780, 2783, 2786, 2789, 2792, 2795, 2798, 2801, 2804, 2807, 2810, 2813, 2816, 2819, 2822, 2825, 2828, 2831, 2834, 2837, 2840, 2843, 2846, 2849, 2852, 2855, 2858, 2861, 2864, 2867, 2870, 2873, 2876, 2879, 2882, 2885, 2888, 2891, 2894, 2897, 2900, 2903, 2906, 2909, 2912, 2915, 2918, 2921, 2924, 2927, 2930, 2933, 2936, 2939, 2942, 2945, 2948, 2951, 2954, 2957, 2960, 2963, 2966, 2969, 2972, 2975, 2978, 2981, 2984, 2987, 2990, 2993, 2996, 2999, 3002, 3005, 3008, 3011, 3014, 3017, 3020, 3023, 3026, 3029, 3032, 3035, 3038, 3041, 3044, 3047, 3050, 3053, 3056, 3059, 3062, 3065, 3068, 3071, 3074, 3077, 3080, 3083, 3086, 3089, 3092, 3095, 3098, 3101, 3104, 3107, 3110, 3113, 3116, 3119, 3122, 3125, 3128, 3131, 3134, 3137, 3140, 3143, 3146, 3149, 3152, 3155, 3158, 3161, 3164, 3167, 3170, 3173, 3176, 3179, 3182, 3185, 3188, 3191, 3194, 3197, 3200, 3203, 3206, 3209, 3212, 3215, 3218, 3221, 3224, 3227, 3230, 3233, 3236, 3239, 3242, 3245, 3248, 3251, 3254, 3257, 3260, 3263, 3266, 3269, 3272, 3275, 3278, 3281, 3284, 3287, 3290, 3293, 3296, 3299, 3302, 3305, 3308, 3311, 3314, 3318, 3322, 3326, 3330, 3334, 3338, 3342, 3346, 3350, 3354, 3358, 3362, 3366, 3370, 3374, 3378, 3382, 3386, 3390, 3394, 3398, 3402, 3406, 3410, 3414, 3418, 3422, 3426, 3430, 3433, 3436, 3439, 3442, 3445, 3448, 3451, 3454, 3457, 3460, 3463, 3466, 3469, 3472, 3475, 3478, 3481, 3484, 3487, 3490, 3493, 3496, 3499, 3502, 3505, 3508, 3511, 3514, 3517, 3520, 3523, 3526, 3529, 3532, 3535, 3538, 3541, 3544, 3547, 3550, 3553, 3556, 3559, 3562, 3565, 3568, 3571, 3574, 3577, 3580, 3583, 3586, 3589, 3592, 3595, 3598, 3601, 3604, 3607, 3610, 3613, 3616, 3619, 3622, 3625, 3628, 3631, 3634, 3637, 3640, 3643, 3646, 3649, 3652, 3655, 3658, 3661, 3664, 3667, 3670, 3673, 3676, 3679, 3682, 3685, 3688, 3691, 3694, 3697, 3700, 3703, 3706, 3709, 3712, 3715, 3718, 3721, 3724, 3727, 3730, 3733, 3736, 3739, 3742, 3745, 3748, 3751, 3754, 3757, 3760, 3763, 3766, 3769, 3772, 3775, 3778, 3781, 3784, 3787, 3790, 3793, 3796, 3799, 3802, 3805, 3808, 3811, 3814, 3817, 3820, 3823, 3826, 3829, 3832, 3835, 3838, 3841, 3844, 3847, 3850, 3853, 3856, 3859, 3862, 3865, 3868, 3871, 3874, 3877, 3880, 3883, 3886, 3889, 3892, 3895, 3898, 3901, 3904, 3907, 3910, 3913, 3916, 3919, 3922, 3925, 3928, 3931, 3934, 3937, 3940, 3943, 3946, 3949, 3952, 3955, 3958, 3961, 3964, 3967, 3970, 3973, 3976, 3979, 3982, 3985, 3988, 3991, 3994, 3997, 4000, 4003, 4006, 4009, 4012, 4014, 4016, 4018, 4020, 4023, 4026, 4029, 4032, 4035, 4038, 4041, 4044, 4047, 4050, 4053, 4056, 4059, 4062, 4065, 4068, 4071, 4074, 4077, 4080, 4082, 4084, 4087, 4090, 4093, 4096, 4099, 4102, 4105, 4108, 4111, 4113, 4115, 4117, 4119, 4122, 4125, 4128, 4130, 4133, 4136, 4139, 4141, 4144, 4147, 4150, 4153, 4156, 4159, 4162, 4165, 4168, 4170, 4173, 4176, 4179, 4181, 4184, 4187, 4190, 4193, 4196, 4198, 4200, 4202, 4205, 4208, 4211, 4213, 4215, 4218, 4221, 4224, 4227, 4230, 4233, 4236, 4239, 4242, 4245, 4248, 4251, 4254, 4257, 4260, 4263, 4266, 4269, 4272, 4274, 4277, 4280, 4283, 4286, 4289, 4292, 4295, 4298, 4301, 4304, 4307, 4310, 4313, 4316, 4318, 4321, 4324, 4327, 4329, 4332, 4335, 4338, 4341, 4344, 4346, 4348, 4350, 4352, 4354, 4357, 4360, 4363, 4366, 4369, 4372, 4375, 4378, 4380, 4382, 4385, 4388, 4391, 4394, 4397, 4399, 4402, 4405, 4408, 4410, 4413, 4416, 4419, 4421, 4424, 4427, 4430, 4432, 4435, 4438, 4440, 4443, 4446, 4449, 4452, 4453, 4454, 4456, 4458, 4460, 4462, 4464, 4466, 4468, 4470, 4472, 4474, 4476, 4478, 4481, 4484, 4487, 4490, 4493, 4496, 4499, 4501, 4504, 4507, 4510, 4512, 4515, 4518, 4521, 4524, 4527, 4530, 4532, 4533, 4534, 4535, 4536, 4537, 4538, 4539, 4540, 4541, 4542, 4543, 4544, 4545, 4546, 4547, 4548, 4549, 4550, 4551, 4552, 4553, 4554, 4555, 4556, 4557, 4558, 4559, 4560, 4561, 4562, 4563, 4564, 4565, 4566, 4567, 4568, 4569, 4570, 4571, 4572, 4573, 4574, 4575, 4576, 4577, 4578, 4579, 4580, 4581, 4582, 4583, 4584, 4585, 4586, 4587, 4588, 4589, 4590, 4591, 4592, 4593, 4594, 4595, 4596, 4597, 4598, 4599, 4600, 4601, 4602, 4603, 4604, 4605, 4606, 4607, 4608, 4609, 4610, 4611, 4612, 4613, 4614, 4615, 4616, 4617, 4618, 4619, 4620, 4621, 4622, 4623, 4624, 4625, 4626, 4627, 4628, 4629, 4630, 4631, 4632, 4633, 4634, 4635, 4636, 4637, 4638, 4639, 4640, 4641, 4642, 4643, 4644, 4645, 4646, 4647, 4648, 4649, 4650, 4651, 4652, 4653, 4654, 4655, 4656, 4657, 4658, 4659, 4660, 4661, 4662, 4663, 4664, 4665, 4666, 4667, 4668, 4669, 4670, 4671, 4672, 4673, 4674, 4675, 4676, 4677, 4678, 4679, 4680, 4681, 4682, 4683, 4684, 4685, 4686, 4687, 4688, 4689, 4690, 4691, 4692, 4693, 4694, 4695, 4696, 4697, 4698, 4699, 4700, 4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716, 4717, 4718, 4719, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 4734, 4735, 4736, 4737, 4738, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 4755, 4756, 4757, 4758, 4759, 4760, 4761, 4762, 4763, 4764, 4765, 4766, 4767, 4768, 4769, 4770, 4771, 4772, 4773, 4774, 4775, 4776, 4777, 4778, 4779, 4780, 4781, 4782, 4783, 4784, 4785, 4786, 4787, 4788, 4789, 4790, 4791, 4792, 4793, 4794, 4795, 4796, 4797, 4798, 4799, 4800, 4801, 4802, 4803, 4804, 4805, 4806, 4807, 4808, 4809, 4810, 4811, 4812, 4813, 4814, 4815, 4816, 4817, 4818, 4819, 4820, 4821, 4822, 4823, 4824, 4825, 4826, 4827, 4828, 4829, 4830, 4831, 4832, 4833, 4834, 4835, 4836, 4837, 4838, 4839, 4840, 4841, 4842, 4843, 4844, 4845, 4846, 4847, 4848, 4849, 4850, 4851, 4852, 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 4899, 4900, 4901, 4902, 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4937, 4938, 4939, 4940, 4941, 4942, 4943, 4944, 4945, 4946, 4947, 4948, 4949, 4950, 4951, 4952, 4953, 4954, 4955, 4956, 4957, 4958, 4959, 4960, 4961, 4962, 4963, 4964, 4965, 4966, 4967, 4968, 4969, 4970, 4971, 4972, 4973, 4974, 4975, 4976, 4977, 4978, 4979, 4980, 4981, 4982, 4983, 4984, 4985, 4986, 4987, 4988, 4989, 4990, 4991, 4992, 4993, 4994, 4995, 4996, 4997, 4998, 4999, 5000, 5001, 5002, 5003, 5004, 5005, 5006, 5007, 5008, 5009, 5010, 5011, 5012, 5013, 5014, 5015, 5016, 5017, 5018, 5019, 5020, 5021, 5022, 5023, 5024, 5025, 5026, 5027, 5028, 5029, 5030, 5031, 5032, 5033, 5034, 5035, 5036, 5037, 5038, 5039, 5040, 5041, 5042, 5043, 5044, 5045, 5046, 5047, 5048, 5049, 5050, 5051, 5052, 5053, 5054, 5055, 5056, 5057, 5058, 5059, 5060, 5061, 5062, 5063, 5064, 5065, 5066, 5067, 5068, 5069, 5070, 5071, 5072, 5073, 5074, 5075, 5076, 5077, 5078, 5079, 5080, 5081, 5082, 5083, 5084, 5085, 5086, 5087, 5088, 5089, 5090, 5091, 5092, 5093, 5094, 5095, 5096, 5097, 5098, 5099, 5100, 5101, 5102, 5103, 5104, 5105, 5106, 5107, 5108, 5109, 5110, 5111, 5112, 5113, 5114, 5115, 5116, 5117, 5118, 5119, 5120, 5121, 5122, 5123, 5124, 5125, 5126, 5127, 5128, 5129, 5130, 5131, 5132, 5133, 5134, 5135, 5136, 5137
};

void dxdotdw_colptrs_SPARCED(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED