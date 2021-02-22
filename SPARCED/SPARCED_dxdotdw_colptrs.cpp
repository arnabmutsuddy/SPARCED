#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 2470> dxdotdw_colptrs_SPARCED_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 251, 254, 255, 256, 257, 258, 259, 261, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 322, 325, 326, 327, 328, 330, 333, 336, 339, 342, 345, 348, 351, 354, 357, 360, 363, 366, 369, 372, 375, 378, 381, 384, 387, 390, 393, 394, 396, 398, 401, 404, 405, 407, 409, 410, 412, 414, 417, 420, 421, 422, 423, 424, 426, 428, 429, 430, 431, 433, 435, 436, 437, 439, 441, 444, 447, 449, 451, 452, 453, 455, 457, 458, 459, 460, 461, 463, 465, 466, 467, 468, 470, 472, 473, 474, 476, 478, 481, 484, 485, 487, 489, 490, 491, 492, 494, 496, 497, 498, 499, 501, 503, 504, 505, 506, 507, 508, 509, 510, 511, 514, 517, 520, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 559, 562, 564, 567, 570, 573, 576, 579, 582, 585, 588, 591, 594, 597, 600, 603, 606, 609, 612, 615, 618, 621, 624, 627, 630, 633, 636, 639, 642, 645, 648, 651, 654, 656, 658, 661, 664, 666, 668, 671, 674, 676, 678, 681, 684, 687, 690, 693, 696, 699, 702, 705, 708, 711, 713, 715, 718, 721, 724, 727, 730, 733, 736, 739, 741, 743, 746, 749, 752, 755, 758, 761, 764, 767, 770, 773, 776, 779, 782, 785, 788, 790, 792, 794, 796, 798, 800, 803, 806, 809, 812, 815, 818, 821, 824, 826, 828, 831, 834, 837, 840, 843, 846, 849, 852, 855, 858, 861, 864, 867, 870, 873, 876, 879, 882, 885, 888, 891, 894, 897, 900, 903, 906, 908, 910, 913, 916, 919, 922, 925, 928, 931, 934, 937, 940, 943, 946, 949, 952, 955, 958, 961, 964, 967, 970, 973, 976, 979, 982, 985, 988, 990, 992, 995, 997, 1000, 1002, 1005, 1007, 1010, 1013, 1016, 1019, 1022, 1025, 1027, 1029, 1032, 1035, 1037, 1039, 1042, 1045, 1048, 1051, 1053, 1055, 1058, 1061, 1064, 1067, 1069, 1071, 1073, 1075, 1078, 1081, 1084, 1087, 1090, 1093, 1095, 1097, 1100, 1103, 1106, 1109, 1112, 1115, 1117, 1119, 1122, 1125, 1128, 1131, 1134, 1137, 1139, 1141, 1144, 1147, 1150, 1153, 1156, 1159, 1162, 1165, 1168, 1171, 1174, 1177, 1180, 1183, 1186, 1189, 1192, 1195, 1198, 1201, 1203, 1205, 1208, 1211, 1213, 1215, 1218, 1221, 1224, 1227, 1230, 1233, 1236, 1239, 1242, 1245, 1247, 1249, 1251, 1253, 1255, 1257, 1259, 1261, 1263, 1265, 1267, 1269, 1271, 1273, 1275, 1277, 1279, 1281, 1283, 1285, 1287, 1289, 1291, 1293, 1295, 1297, 1299, 1301, 1303, 1305, 1307, 1309, 1311, 1313, 1315, 1317, 1319, 1321, 1323, 1325, 1327, 1329, 1331, 1333, 1335, 1337, 1339, 1341, 1343, 1345, 1347, 1349, 1351, 1353, 1355, 1357, 1359, 1361, 1364, 1367, 1370, 1373, 1376, 1379, 1382, 1385, 1388, 1391, 1394, 1397, 1400, 1403, 1406, 1409, 1412, 1415, 1418, 1421, 1424, 1427, 1430, 1433, 1436, 1439, 1442, 1445, 1448, 1451, 1454, 1457, 1460, 1463, 1466, 1469, 1472, 1475, 1478, 1481, 1484, 1487, 1490, 1493, 1496, 1499, 1502, 1505, 1508, 1511, 1514, 1517, 1520, 1523, 1526, 1529, 1532, 1535, 1537, 1539, 1541, 1543, 1545, 1547, 1549, 1551, 1553, 1555, 1557, 1559, 1561, 1563, 1565, 1567, 1569, 1571, 1573, 1575, 1577, 1579, 1581, 1583, 1585, 1587, 1589, 1591, 1593, 1595, 1597, 1599, 1601, 1603, 1605, 1607, 1609, 1611, 1613, 1615, 1617, 1619, 1621, 1623, 1625, 1627, 1629, 1631, 1633, 1635, 1637, 1639, 1641, 1643, 1645, 1647, 1649, 1651, 1653, 1655, 1657, 1659, 1661, 1663, 1665, 1667, 1669, 1671, 1673, 1675, 1677, 1679, 1681, 1683, 1685, 1687, 1689, 1691, 1693, 1695, 1697, 1699, 1701, 1703, 1705, 1707, 1709, 1711, 1713, 1715, 1717, 1719, 1721, 1723, 1725, 1727, 1729, 1731, 1733, 1735, 1737, 1739, 1741, 1743, 1745, 1747, 1749, 1751, 1753, 1755, 1757, 1759, 1761, 1763, 1765, 1767, 1770, 1773, 1776, 1779, 1782, 1785, 1788, 1791, 1794, 1797, 1800, 1803, 1806, 1809, 1812, 1815, 1818, 1821, 1824, 1827, 1830, 1833, 1836, 1839, 1842, 1845, 1848, 1851, 1854, 1857, 1860, 1863, 1866, 1869, 1872, 1875, 1878, 1881, 1884, 1887, 1890, 1893, 1896, 1899, 1902, 1905, 1908, 1911, 1914, 1917, 1920, 1923, 1926, 1929, 1932, 1935, 1938, 1941, 1944, 1947, 1950, 1953, 1956, 1959, 1962, 1965, 1968, 1971, 1974, 1977, 1980, 1983, 1986, 1989, 1992, 1995, 1998, 2001, 2004, 2007, 2010, 2013, 2016, 2019, 2022, 2025, 2028, 2031, 2034, 2037, 2040, 2043, 2046, 2049, 2052, 2055, 2058, 2061, 2064, 2067, 2070, 2073, 2076, 2079, 2082, 2085, 2088, 2091, 2094, 2097, 2100, 2103, 2106, 2109, 2112, 2115, 2118, 2121, 2124, 2127, 2130, 2133, 2136, 2139, 2142, 2145, 2148, 2151, 2154, 2157, 2160, 2163, 2166, 2169, 2172, 2175, 2178, 2181, 2184, 2187, 2190, 2193, 2196, 2199, 2202, 2205, 2208, 2211, 2214, 2217, 2220, 2223, 2226, 2229, 2232, 2235, 2238, 2241, 2244, 2247, 2250, 2253, 2256, 2259, 2262, 2265, 2268, 2271, 2274, 2277, 2280, 2283, 2286, 2289, 2292, 2295, 2298, 2301, 2304, 2307, 2310, 2313, 2316, 2319, 2322, 2325, 2328, 2331, 2334, 2337, 2340, 2343, 2346, 2349, 2352, 2355, 2358, 2361, 2364, 2367, 2370, 2373, 2376, 2379, 2382, 2385, 2388, 2391, 2394, 2397, 2400, 2403, 2406, 2409, 2412, 2415, 2418, 2421, 2424, 2427, 2430, 2433, 2436, 2439, 2442, 2445, 2448, 2451, 2454, 2457, 2460, 2463, 2466, 2469, 2472, 2475, 2478, 2481, 2484, 2487, 2490, 2493, 2496, 2499, 2502, 2505, 2508, 2511, 2514, 2517, 2520, 2523, 2526, 2529, 2532, 2535, 2538, 2541, 2544, 2547, 2550, 2553, 2556, 2559, 2562, 2565, 2568, 2571, 2574, 2577, 2580, 2583, 2586, 2589, 2592, 2595, 2598, 2601, 2604, 2607, 2610, 2613, 2616, 2619, 2622, 2625, 2628, 2631, 2634, 2637, 2640, 2643, 2646, 2649, 2652, 2655, 2658, 2661, 2664, 2667, 2670, 2673, 2676, 2679, 2682, 2685, 2688, 2691, 2694, 2697, 2700, 2703, 2706, 2709, 2712, 2715, 2718, 2721, 2724, 2727, 2730, 2733, 2736, 2739, 2742, 2745, 2748, 2751, 2754, 2757, 2760, 2763, 2766, 2769, 2772, 2775, 2778, 2781, 2784, 2787, 2790, 2793, 2796, 2799, 2802, 2805, 2808, 2811, 2814, 2817, 2820, 2823, 2826, 2829, 2832, 2835, 2838, 2841, 2844, 2847, 2850, 2853, 2856, 2859, 2862, 2865, 2868, 2871, 2874, 2877, 2880, 2883, 2886, 2889, 2892, 2895, 2898, 2901, 2904, 2907, 2910, 2913, 2916, 2919, 2922, 2925, 2928, 2931, 2934, 2937, 2940, 2943, 2946, 2949, 2952, 2955, 2958, 2961, 2964, 2967, 2970, 2973, 2976, 2979, 2982, 2985, 2988, 2991, 2994, 2997, 3000, 3003, 3006, 3009, 3012, 3015, 3018, 3021, 3024, 3027, 3030, 3033, 3036, 3039, 3042, 3045, 3048, 3051, 3054, 3057, 3060, 3063, 3066, 3069, 3072, 3075, 3078, 3081, 3084, 3087, 3090, 3093, 3096, 3099, 3102, 3105, 3108, 3111, 3114, 3117, 3120, 3123, 3126, 3129, 3132, 3135, 3138, 3141, 3144, 3147, 3150, 3153, 3156, 3159, 3162, 3165, 3168, 3171, 3174, 3177, 3180, 3183, 3186, 3189, 3192, 3195, 3198, 3201, 3204, 3207, 3210, 3213, 3216, 3219, 3222, 3225, 3228, 3231, 3234, 3237, 3240, 3243, 3246, 3249, 3252, 3255, 3258, 3261, 3264, 3267, 3270, 3273, 3276, 3279, 3282, 3285, 3288, 3291, 3294, 3297, 3300, 3303, 3306, 3309, 3312, 3315, 3318, 3321, 3324, 3327, 3330, 3333, 3337, 3341, 3345, 3349, 3353, 3357, 3361, 3365, 3369, 3373, 3377, 3381, 3385, 3389, 3393, 3397, 3401, 3405, 3409, 3413, 3417, 3421, 3425, 3429, 3433, 3437, 3441, 3445, 3449, 3452, 3455, 3458, 3461, 3464, 3467, 3470, 3473, 3476, 3479, 3482, 3485, 3488, 3491, 3494, 3497, 3500, 3503, 3506, 3509, 3512, 3515, 3518, 3521, 3524, 3527, 3530, 3533, 3536, 3539, 3542, 3545, 3548, 3551, 3554, 3557, 3560, 3563, 3566, 3569, 3572, 3575, 3578, 3581, 3584, 3587, 3590, 3593, 3596, 3599, 3602, 3605, 3608, 3611, 3614, 3617, 3620, 3623, 3626, 3629, 3632, 3635, 3638, 3641, 3644, 3647, 3650, 3653, 3656, 3659, 3662, 3665, 3668, 3671, 3674, 3677, 3680, 3683, 3686, 3689, 3692, 3695, 3698, 3701, 3704, 3707, 3710, 3713, 3716, 3719, 3722, 3725, 3728, 3731, 3734, 3737, 3740, 3743, 3746, 3749, 3752, 3755, 3758, 3761, 3764, 3767, 3770, 3773, 3776, 3779, 3782, 3785, 3788, 3791, 3794, 3797, 3800, 3803, 3806, 3809, 3812, 3815, 3818, 3821, 3824, 3827, 3830, 3833, 3836, 3839, 3842, 3845, 3848, 3851, 3854, 3857, 3860, 3863, 3866, 3869, 3872, 3875, 3878, 3881, 3884, 3887, 3890, 3893, 3896, 3899, 3902, 3905, 3908, 3911, 3914, 3917, 3920, 3923, 3926, 3929, 3932, 3935, 3938, 3941, 3944, 3947, 3950, 3953, 3956, 3959, 3962, 3965, 3968, 3971, 3974, 3977, 3980, 3983, 3986, 3989, 3992, 3995, 3998, 4001, 4004, 4007, 4010, 4013, 4016, 4019, 4022, 4025, 4028, 4031, 4033, 4035, 4037, 4039, 4042, 4045, 4048, 4051, 4054, 4057, 4060, 4063, 4066, 4069, 4072, 4075, 4078, 4081, 4084, 4087, 4090, 4093, 4096, 4099, 4101, 4103, 4106, 4109, 4112, 4115, 4118, 4121, 4124, 4127, 4130, 4132, 4134, 4136, 4138, 4141, 4144, 4147, 4149, 4152, 4155, 4158, 4160, 4163, 4166, 4169, 4172, 4175, 4178, 4181, 4184, 4187, 4189, 4192, 4195, 4198, 4200, 4203, 4206, 4209, 4212, 4215, 4217, 4219, 4221, 4224, 4227, 4230, 4232, 4234, 4237, 4240, 4243, 4246, 4249, 4252, 4255, 4258, 4261, 4264, 4267, 4270, 4273, 4276, 4279, 4282, 4285, 4288, 4291, 4293, 4296, 4299, 4302, 4305, 4308, 4311, 4314, 4317, 4320, 4323, 4326, 4329, 4332, 4335, 4337, 4340, 4343, 4346, 4348, 4351, 4354, 4357, 4360, 4363, 4365, 4367, 4369, 4371, 4373, 4376, 4379, 4382, 4385, 4388, 4391, 4394, 4397, 4399, 4401, 4404, 4407, 4410, 4413, 4416, 4418, 4421, 4424, 4427, 4429, 4432, 4435, 4438, 4440, 4443, 4446, 4449, 4451, 4454, 4457, 4459, 4462, 4465, 4468, 4471, 4472, 4473, 4475, 4477, 4479, 4481, 4483, 4485, 4487, 4489, 4491, 4493, 4495, 4497, 4500, 4503, 4506, 4509, 4512, 4515, 4518, 4520, 4523, 4526, 4529, 4531, 4534, 4537, 4540, 4543, 4546, 4549, 4551, 4552, 4553, 4554, 4555, 4556, 4557, 4558, 4559, 4560, 4561, 4562, 4563, 4564, 4565, 4566, 4567, 4568, 4569, 4570, 4571, 4572, 4573, 4574, 4575, 4576, 4577, 4578, 4579, 4580, 4581, 4582, 4583, 4584, 4585, 4586, 4587, 4588, 4589, 4590, 4591, 4592, 4593, 4594, 4595, 4596, 4597, 4598, 4599, 4600, 4601, 4602, 4603, 4604, 4605, 4606, 4607, 4608, 4609, 4610, 4611, 4612, 4613, 4614, 4615, 4616, 4617, 4618, 4619, 4620, 4621, 4622, 4623, 4624, 4625, 4626, 4627, 4628, 4629, 4630, 4631, 4632, 4633, 4634, 4635, 4636, 4637, 4638, 4639, 4640, 4641, 4642, 4643, 4644, 4645, 4646, 4647, 4648, 4649, 4650, 4651, 4652, 4653, 4654, 4655, 4656, 4657, 4658, 4659, 4660, 4661, 4662, 4663, 4664, 4665, 4666, 4667, 4668, 4669, 4670, 4671, 4672, 4673, 4674, 4675, 4676, 4677, 4678, 4679, 4680, 4681, 4682, 4683, 4684, 4685, 4686, 4687, 4688, 4689, 4690, 4691, 4692, 4693, 4694, 4695, 4696, 4697, 4698, 4699, 4700, 4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716, 4717, 4718, 4719, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 4734, 4735, 4736, 4737, 4738, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 4755, 4756, 4757, 4758, 4759, 4760, 4761, 4762, 4763, 4764, 4765, 4766, 4767, 4768, 4769, 4770, 4771, 4772, 4773, 4774, 4775, 4776, 4777, 4778, 4779, 4780, 4781, 4782, 4783, 4784, 4785, 4786, 4787, 4788, 4789, 4790, 4791, 4792, 4793, 4794, 4795, 4796, 4797, 4798, 4799, 4800, 4801, 4802, 4803, 4804, 4805, 4806, 4807, 4808, 4809, 4810, 4811, 4812, 4813, 4814, 4815, 4816, 4817, 4818, 4819, 4820, 4821, 4822, 4823, 4824, 4825, 4826, 4827, 4828, 4829, 4830, 4831, 4832, 4833, 4834, 4835, 4836, 4837, 4838, 4839, 4840, 4841, 4842, 4843, 4844, 4845, 4846, 4847, 4848, 4849, 4850, 4851, 4852, 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 4899, 4900, 4901, 4902, 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4937, 4938, 4939, 4940, 4941, 4942, 4943, 4944, 4945, 4946, 4947, 4948, 4949, 4950, 4951, 4952, 4953, 4954, 4955, 4956, 4957, 4958, 4959, 4960, 4961, 4962, 4963, 4964, 4965, 4966, 4967, 4968, 4969, 4970, 4971, 4972, 4973, 4974, 4975, 4976, 4977, 4978, 4979, 4980, 4981, 4982, 4983, 4984, 4985, 4986, 4987, 4988, 4989, 4990, 4991, 4992, 4993, 4994, 4995, 4996, 4997, 4998, 4999, 5000, 5001, 5002, 5003, 5004, 5005, 5006, 5007, 5008, 5009, 5010, 5011, 5012, 5013, 5014, 5015, 5016, 5017, 5018, 5019, 5020, 5021, 5022, 5023, 5024, 5025, 5026, 5027, 5028, 5029, 5030, 5031, 5032, 5033, 5034, 5035, 5036, 5037, 5038, 5039, 5040, 5041, 5042, 5043, 5044, 5045, 5046, 5047, 5048, 5049, 5050, 5051, 5052, 5053, 5054, 5055, 5056, 5057, 5058, 5059, 5060, 5061, 5062, 5063, 5064, 5065, 5066, 5067, 5068, 5069, 5070, 5071, 5072, 5073, 5074, 5075, 5076, 5077, 5078, 5079, 5080, 5081, 5082, 5083, 5084, 5085, 5086, 5087, 5088, 5089, 5090, 5091, 5092, 5093, 5094, 5095, 5096, 5097, 5098, 5099, 5100, 5101, 5102, 5103, 5104, 5105, 5106, 5107, 5108, 5109, 5110, 5111, 5112, 5113, 5114, 5115, 5116, 5117, 5118, 5119, 5120, 5121, 5122, 5123, 5124, 5125, 5126, 5127, 5128, 5129, 5130, 5131, 5132, 5133, 5134, 5135, 5136, 5137, 5138, 5139, 5140, 5141, 5142, 5143, 5144, 5145, 5146, 5147, 5148, 5149, 5150, 5151, 5152, 5153, 5154, 5155, 5156
};

void dxdotdw_colptrs_SPARCED(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED