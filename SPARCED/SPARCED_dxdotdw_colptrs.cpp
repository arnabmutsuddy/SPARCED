#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 2400> dxdotdw_colptrs_SPARCED_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 257, 260, 261, 262, 263, 264, 265, 267, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 328, 331, 332, 333, 334, 336, 339, 342, 345, 348, 351, 354, 357, 360, 363, 366, 369, 372, 375, 378, 381, 384, 387, 390, 393, 396, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 437, 440, 442, 445, 448, 451, 454, 457, 460, 463, 466, 469, 472, 475, 478, 481, 484, 487, 490, 493, 496, 499, 502, 505, 508, 511, 514, 517, 520, 523, 526, 529, 532, 534, 536, 539, 542, 544, 546, 549, 552, 554, 556, 559, 562, 565, 568, 571, 574, 577, 580, 583, 586, 589, 591, 593, 596, 599, 602, 605, 608, 611, 614, 617, 619, 621, 624, 627, 630, 633, 636, 639, 642, 645, 648, 651, 654, 657, 660, 663, 666, 668, 670, 672, 674, 676, 678, 681, 684, 687, 690, 693, 696, 699, 702, 704, 706, 709, 712, 715, 718, 721, 724, 727, 730, 733, 736, 739, 742, 745, 748, 751, 754, 757, 760, 763, 766, 769, 772, 775, 778, 781, 784, 786, 788, 791, 794, 797, 800, 803, 806, 809, 812, 815, 818, 821, 824, 827, 830, 833, 836, 839, 842, 845, 848, 851, 854, 857, 860, 863, 866, 868, 870, 873, 875, 878, 880, 883, 885, 888, 891, 894, 897, 900, 903, 905, 907, 910, 913, 915, 917, 920, 923, 926, 929, 931, 933, 936, 939, 942, 945, 947, 949, 951, 953, 956, 959, 962, 965, 968, 971, 973, 975, 978, 981, 984, 987, 990, 993, 995, 997, 1000, 1003, 1006, 1009, 1012, 1015, 1017, 1019, 1022, 1025, 1028, 1031, 1034, 1037, 1040, 1043, 1046, 1049, 1052, 1055, 1058, 1061, 1064, 1067, 1070, 1073, 1076, 1079, 1081, 1083, 1086, 1089, 1091, 1093, 1096, 1099, 1102, 1105, 1108, 1111, 1114, 1117, 1120, 1123, 1125, 1127, 1129, 1131, 1133, 1135, 1137, 1139, 1141, 1143, 1145, 1147, 1149, 1151, 1153, 1155, 1157, 1159, 1161, 1163, 1165, 1167, 1169, 1171, 1173, 1175, 1177, 1179, 1181, 1183, 1185, 1187, 1189, 1191, 1193, 1195, 1197, 1199, 1201, 1203, 1205, 1207, 1209, 1211, 1213, 1215, 1217, 1219, 1221, 1223, 1225, 1227, 1229, 1231, 1233, 1235, 1237, 1239, 1242, 1245, 1248, 1251, 1254, 1257, 1260, 1263, 1266, 1269, 1272, 1275, 1278, 1281, 1284, 1287, 1290, 1293, 1296, 1299, 1302, 1305, 1308, 1311, 1314, 1317, 1320, 1323, 1326, 1329, 1332, 1335, 1338, 1341, 1344, 1347, 1350, 1353, 1356, 1359, 1362, 1365, 1368, 1371, 1374, 1377, 1380, 1383, 1386, 1389, 1392, 1395, 1398, 1401, 1404, 1407, 1410, 1413, 1415, 1417, 1419, 1421, 1423, 1425, 1427, 1429, 1431, 1433, 1435, 1437, 1439, 1441, 1443, 1445, 1447, 1449, 1451, 1453, 1455, 1457, 1459, 1461, 1463, 1465, 1467, 1469, 1471, 1473, 1475, 1477, 1479, 1481, 1483, 1485, 1487, 1489, 1491, 1493, 1495, 1497, 1499, 1501, 1503, 1505, 1507, 1509, 1511, 1513, 1515, 1517, 1519, 1521, 1523, 1525, 1527, 1529, 1531, 1533, 1535, 1537, 1539, 1541, 1543, 1545, 1547, 1549, 1551, 1553, 1555, 1557, 1559, 1561, 1563, 1565, 1567, 1569, 1571, 1573, 1575, 1577, 1579, 1581, 1583, 1585, 1587, 1589, 1591, 1593, 1595, 1597, 1599, 1601, 1603, 1605, 1607, 1609, 1611, 1613, 1615, 1617, 1619, 1621, 1623, 1625, 1627, 1629, 1631, 1633, 1635, 1637, 1639, 1641, 1643, 1645, 1648, 1651, 1654, 1657, 1660, 1663, 1666, 1669, 1672, 1675, 1678, 1681, 1684, 1687, 1690, 1693, 1696, 1699, 1702, 1705, 1708, 1711, 1714, 1717, 1720, 1723, 1726, 1729, 1732, 1735, 1738, 1741, 1744, 1747, 1750, 1753, 1756, 1759, 1762, 1765, 1768, 1771, 1774, 1777, 1780, 1783, 1786, 1789, 1792, 1795, 1798, 1801, 1804, 1807, 1810, 1813, 1816, 1819, 1822, 1825, 1828, 1831, 1834, 1837, 1840, 1843, 1846, 1849, 1852, 1855, 1858, 1861, 1864, 1867, 1870, 1873, 1876, 1879, 1882, 1885, 1888, 1891, 1894, 1897, 1900, 1903, 1906, 1909, 1912, 1915, 1918, 1921, 1924, 1927, 1930, 1933, 1936, 1939, 1942, 1945, 1948, 1951, 1954, 1957, 1960, 1963, 1966, 1969, 1972, 1975, 1978, 1981, 1984, 1987, 1990, 1993, 1996, 1999, 2002, 2005, 2008, 2011, 2014, 2017, 2020, 2023, 2026, 2029, 2032, 2035, 2038, 2041, 2044, 2047, 2050, 2053, 2056, 2059, 2062, 2065, 2068, 2071, 2074, 2077, 2080, 2083, 2086, 2089, 2092, 2095, 2098, 2101, 2104, 2107, 2110, 2113, 2116, 2119, 2122, 2125, 2128, 2131, 2134, 2137, 2140, 2143, 2146, 2149, 2152, 2155, 2158, 2161, 2164, 2167, 2170, 2173, 2176, 2179, 2182, 2185, 2188, 2191, 2194, 2197, 2200, 2203, 2206, 2209, 2212, 2215, 2218, 2221, 2224, 2227, 2230, 2233, 2236, 2239, 2242, 2245, 2248, 2251, 2254, 2257, 2260, 2263, 2266, 2269, 2272, 2275, 2278, 2281, 2284, 2287, 2290, 2293, 2296, 2299, 2302, 2305, 2308, 2311, 2314, 2317, 2320, 2323, 2326, 2329, 2332, 2335, 2338, 2341, 2344, 2347, 2350, 2353, 2356, 2359, 2362, 2365, 2368, 2371, 2374, 2377, 2380, 2383, 2386, 2389, 2392, 2395, 2398, 2401, 2404, 2407, 2410, 2413, 2416, 2419, 2422, 2425, 2428, 2431, 2434, 2437, 2440, 2443, 2446, 2449, 2452, 2455, 2458, 2461, 2464, 2467, 2470, 2473, 2476, 2479, 2482, 2485, 2488, 2491, 2494, 2497, 2500, 2503, 2506, 2509, 2512, 2515, 2518, 2521, 2524, 2527, 2530, 2533, 2536, 2539, 2542, 2545, 2548, 2551, 2554, 2557, 2560, 2563, 2566, 2569, 2572, 2575, 2578, 2581, 2584, 2587, 2590, 2593, 2596, 2599, 2602, 2605, 2608, 2611, 2614, 2617, 2620, 2623, 2626, 2629, 2632, 2635, 2638, 2641, 2644, 2647, 2650, 2653, 2656, 2659, 2662, 2665, 2668, 2671, 2674, 2677, 2680, 2683, 2686, 2689, 2692, 2695, 2698, 2701, 2704, 2707, 2710, 2713, 2716, 2719, 2722, 2725, 2728, 2731, 2734, 2737, 2740, 2743, 2746, 2749, 2752, 2755, 2758, 2761, 2764, 2767, 2770, 2773, 2776, 2779, 2782, 2785, 2788, 2791, 2794, 2797, 2800, 2803, 2806, 2809, 2812, 2815, 2818, 2821, 2824, 2827, 2830, 2833, 2836, 2839, 2842, 2845, 2848, 2851, 2854, 2857, 2860, 2863, 2866, 2869, 2872, 2875, 2878, 2881, 2884, 2887, 2890, 2893, 2896, 2899, 2902, 2905, 2908, 2911, 2914, 2917, 2920, 2923, 2926, 2929, 2932, 2935, 2938, 2941, 2944, 2947, 2950, 2953, 2956, 2959, 2962, 2965, 2968, 2971, 2974, 2977, 2980, 2983, 2986, 2989, 2992, 2995, 2998, 3001, 3004, 3007, 3010, 3013, 3016, 3019, 3022, 3025, 3028, 3031, 3034, 3037, 3040, 3043, 3046, 3049, 3052, 3055, 3058, 3061, 3064, 3067, 3070, 3073, 3076, 3079, 3082, 3085, 3088, 3091, 3094, 3097, 3100, 3103, 3106, 3109, 3112, 3115, 3118, 3121, 3124, 3127, 3130, 3133, 3136, 3139, 3142, 3145, 3148, 3151, 3154, 3157, 3160, 3163, 3166, 3169, 3172, 3175, 3178, 3181, 3184, 3187, 3190, 3193, 3196, 3199, 3202, 3205, 3208, 3211, 3215, 3219, 3223, 3227, 3231, 3235, 3239, 3243, 3247, 3251, 3255, 3259, 3263, 3267, 3271, 3275, 3279, 3283, 3287, 3291, 3295, 3299, 3303, 3307, 3311, 3315, 3319, 3323, 3327, 3330, 3333, 3336, 3339, 3342, 3345, 3348, 3351, 3354, 3357, 3360, 3363, 3366, 3369, 3372, 3375, 3378, 3381, 3384, 3387, 3390, 3393, 3396, 3399, 3402, 3405, 3408, 3411, 3414, 3417, 3420, 3423, 3426, 3429, 3432, 3435, 3438, 3441, 3444, 3447, 3450, 3453, 3456, 3459, 3462, 3465, 3468, 3471, 3474, 3477, 3480, 3483, 3486, 3489, 3492, 3495, 3498, 3501, 3504, 3507, 3510, 3513, 3516, 3519, 3522, 3525, 3528, 3531, 3534, 3537, 3540, 3543, 3546, 3549, 3552, 3555, 3558, 3561, 3564, 3567, 3570, 3573, 3576, 3579, 3582, 3585, 3588, 3591, 3594, 3597, 3600, 3603, 3606, 3609, 3612, 3615, 3618, 3621, 3624, 3627, 3630, 3633, 3636, 3639, 3642, 3645, 3648, 3651, 3654, 3657, 3660, 3663, 3666, 3669, 3672, 3675, 3678, 3681, 3684, 3687, 3690, 3693, 3696, 3699, 3702, 3705, 3708, 3711, 3714, 3717, 3720, 3723, 3726, 3729, 3732, 3735, 3738, 3741, 3744, 3747, 3750, 3753, 3756, 3759, 3762, 3765, 3768, 3771, 3774, 3777, 3780, 3783, 3786, 3789, 3792, 3795, 3798, 3801, 3804, 3807, 3810, 3813, 3816, 3819, 3822, 3825, 3828, 3831, 3834, 3837, 3840, 3843, 3846, 3849, 3852, 3855, 3858, 3861, 3864, 3867, 3870, 3873, 3876, 3879, 3882, 3885, 3888, 3891, 3894, 3897, 3900, 3903, 3906, 3909, 3911, 3913, 3915, 3917, 3920, 3923, 3926, 3929, 3932, 3935, 3938, 3941, 3944, 3947, 3950, 3953, 3956, 3959, 3962, 3965, 3968, 3971, 3974, 3977, 3979, 3981, 3984, 3987, 3990, 3993, 3996, 3999, 4002, 4005, 4008, 4010, 4012, 4014, 4016, 4019, 4022, 4025, 4027, 4030, 4033, 4036, 4038, 4041, 4044, 4047, 4050, 4053, 4056, 4059, 4062, 4065, 4067, 4070, 4073, 4076, 4078, 4081, 4084, 4087, 4090, 4093, 4095, 4097, 4099, 4102, 4105, 4108, 4110, 4112, 4115, 4118, 4121, 4124, 4127, 4130, 4133, 4136, 4139, 4142, 4145, 4148, 4151, 4154, 4157, 4160, 4163, 4166, 4169, 4171, 4174, 4177, 4180, 4183, 4186, 4189, 4192, 4195, 4198, 4201, 4204, 4207, 4210, 4213, 4215, 4218, 4221, 4224, 4226, 4229, 4232, 4235, 4238, 4241, 4243, 4245, 4247, 4249, 4251, 4254, 4257, 4260, 4263, 4266, 4269, 4272, 4275, 4277, 4279, 4282, 4285, 4288, 4291, 4294, 4296, 4299, 4302, 4305, 4307, 4310, 4313, 4316, 4318, 4321, 4324, 4327, 4329, 4332, 4335, 4337, 4340, 4343, 4346, 4349, 4350, 4351, 4353, 4355, 4357, 4359, 4361, 4363, 4365, 4367, 4369, 4371, 4373, 4375, 4378, 4381, 4384, 4387, 4390, 4393, 4396, 4398, 4401, 4404, 4407, 4409, 4412, 4415, 4418, 4421, 4424, 4427, 4429, 4430, 4431, 4432, 4433, 4434, 4435, 4436, 4437, 4438, 4439, 4440, 4441, 4442, 4443, 4444, 4445, 4446, 4447, 4448, 4449, 4450, 4451, 4452, 4453, 4454, 4455, 4456, 4457, 4458, 4459, 4460, 4461, 4462, 4463, 4464, 4465, 4466, 4467, 4468, 4469, 4470, 4471, 4472, 4473, 4474, 4475, 4476, 4477, 4478, 4479, 4480, 4481, 4482, 4483, 4484, 4485, 4486, 4487, 4488, 4489, 4490, 4491, 4492, 4493, 4494, 4495, 4496, 4497, 4498, 4499, 4500, 4501, 4502, 4503, 4504, 4505, 4506, 4507, 4508, 4509, 4510, 4511, 4512, 4513, 4514, 4515, 4516, 4517, 4518, 4519, 4520, 4521, 4522, 4523, 4524, 4525, 4526, 4527, 4528, 4529, 4530, 4531, 4532, 4533, 4534, 4535, 4536, 4537, 4538, 4539, 4540, 4541, 4542, 4543, 4544, 4545, 4546, 4547, 4548, 4549, 4550, 4551, 4552, 4553, 4554, 4555, 4556, 4557, 4558, 4559, 4560, 4561, 4562, 4563, 4564, 4565, 4566, 4567, 4568, 4569, 4570, 4571, 4572, 4573, 4574, 4575, 4576, 4577, 4578, 4579, 4580, 4581, 4582, 4583, 4584, 4585, 4586, 4587, 4588, 4589, 4590, 4591, 4592, 4593, 4594, 4595, 4596, 4597, 4598, 4599, 4600, 4601, 4602, 4603, 4604, 4605, 4606, 4607, 4608, 4609, 4610, 4611, 4612, 4613, 4614, 4615, 4616, 4617, 4618, 4619, 4620, 4621, 4622, 4623, 4624, 4625, 4626, 4627, 4628, 4629, 4630, 4631, 4632, 4633, 4634, 4635, 4636, 4637, 4638, 4639, 4640, 4641, 4642, 4643, 4644, 4645, 4646, 4647, 4648, 4649, 4650, 4651, 4652, 4653, 4654, 4655, 4656, 4657, 4658, 4659, 4660, 4661, 4662, 4663, 4664, 4665, 4666, 4667, 4668, 4669, 4670, 4671, 4672, 4673, 4674, 4675, 4676, 4677, 4678, 4679, 4680, 4681, 4682, 4683, 4684, 4685, 4686, 4687, 4688, 4689, 4690, 4691, 4692, 4693, 4694, 4695, 4696, 4697, 4698, 4699, 4700, 4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710, 4711, 4712, 4713, 4714, 4715, 4716, 4717, 4718, 4719, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4729, 4730, 4731, 4732, 4733, 4734, 4735, 4736, 4737, 4738, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4747, 4748, 4749, 4750, 4751, 4752, 4753, 4754, 4755, 4756, 4757, 4758, 4759, 4760, 4761, 4762, 4763, 4764, 4765, 4766, 4767, 4768, 4769, 4770, 4771, 4772, 4773, 4774, 4775, 4776, 4777, 4778, 4779, 4780, 4781, 4782, 4783, 4784, 4785, 4786, 4787, 4788, 4789, 4790, 4791, 4792, 4793, 4794, 4795, 4796, 4797, 4798, 4799, 4800, 4801, 4802, 4803, 4804, 4805, 4806, 4807, 4808, 4809, 4810, 4811, 4812, 4813, 4814, 4815, 4816, 4817, 4818, 4819, 4820, 4821, 4822, 4823, 4824, 4825, 4826, 4827, 4828, 4829, 4830, 4831, 4832, 4833, 4834, 4835, 4836, 4837, 4838, 4839, 4840, 4841, 4842, 4843, 4844, 4845, 4846, 4847, 4848, 4849, 4850, 4851, 4852, 4853, 4854, 4855, 4856, 4857, 4858, 4859, 4860, 4861, 4862, 4863, 4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4872, 4873, 4874, 4875, 4876, 4877, 4878, 4879, 4880, 4881, 4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4893, 4894, 4895, 4896, 4897, 4898, 4899, 4900, 4901, 4902, 4903, 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922, 4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4937, 4938, 4939, 4940, 4941, 4942, 4943, 4944, 4945, 4946, 4947, 4948, 4949, 4950, 4951, 4952, 4953, 4954, 4955, 4956, 4957, 4958, 4959, 4960, 4961, 4962, 4963, 4964, 4965, 4966, 4967, 4968, 4969, 4970, 4971, 4972, 4973, 4974, 4975, 4976, 4977, 4978, 4979, 4980, 4981, 4982, 4983, 4984, 4985, 4986, 4987, 4988, 4989, 4990, 4991, 4992, 4993, 4994, 4995, 4996, 4997, 4998, 4999, 5000, 5001, 5002, 5003, 5004, 5005, 5006, 5007, 5008, 5009, 5010, 5011, 5012, 5013, 5014, 5015, 5016, 5017, 5018, 5019, 5020, 5021, 5022, 5023, 5024, 5025, 5026, 5027, 5028, 5029, 5030, 5031, 5032, 5033, 5034
};

void dxdotdw_colptrs_SPARCED(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED