#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 782> dwdx_colptrs_SPARCED_ = {
    0, 1, 5, 10, 18, 22, 26, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 81, 83, 85, 88, 90, 91, 92, 93, 100, 106, 106, 113, 117, 123, 127, 129, 135, 144, 153, 158, 162, 166, 169, 172, 175, 180, 186, 195, 199, 206, 209, 214, 215, 219, 223, 228, 233, 242, 245, 250, 254, 260, 263, 266, 272, 275, 278, 284, 286, 290, 294, 297, 299, 301, 304, 307, 309, 311, 315, 318, 322, 324, 326, 329, 332, 336, 338, 341, 343, 346, 350, 353, 355, 358, 359, 361, 364, 367, 374, 376, 379, 382, 384, 388, 393, 395, 399, 401, 405, 407, 408, 410, 413, 415, 417, 420, 422, 424, 427, 430, 432, 435, 437, 439, 443, 446, 449, 451, 453, 454, 458, 460, 462, 464, 466, 468, 472, 475, 477, 480, 482, 485, 487, 490, 499, 514, 518, 520, 524, 527, 530, 540, 542, 551, 553, 561, 571, 573, 577, 581, 583, 587, 589, 591, 594, 597, 601, 605, 607, 610, 613, 616, 619, 622, 625, 628, 631, 641, 651, 661, 665, 668, 672, 677, 681, 687, 693, 699, 704, 710, 713, 717, 722, 727, 731, 737, 743, 749, 755, 759, 765, 770, 775, 780, 786, 792, 795, 801, 806, 811, 816, 821, 825, 829, 832, 835, 838, 845, 852, 859, 866, 873, 880, 887, 894, 901, 908, 915, 922, 929, 936, 943, 950, 957, 964, 971, 978, 985, 992, 999, 1006, 1013, 1017, 1021, 1025, 1029, 1035, 1041, 1047, 1053, 1055, 1057, 1059, 1061, 1063, 1065, 1067, 1069, 1071, 1073, 1075, 1077, 1079, 1081, 1083, 1085, 1087, 1089, 1091, 1093, 1095, 1097, 1099, 1101, 1103, 1105, 1107, 1109, 1111, 1114, 1117, 1120, 1123, 1126, 1129, 1132, 1135, 1138, 1141, 1144, 1147, 1150, 1153, 1156, 1159, 1162, 1165, 1168, 1171, 1174, 1177, 1180, 1183, 1186, 1189, 1192, 1195, 1198, 1201, 1204, 1207, 1210, 1213, 1216, 1219, 1222, 1225, 1228, 1231, 1234, 1237, 1240, 1243, 1246, 1249, 1252, 1255, 1258, 1261, 1264, 1267, 1270, 1273, 1276, 1279, 1282, 1285, 1288, 1291, 1294, 1297, 1300, 1303, 1306, 1309, 1312, 1315, 1318, 1321, 1324, 1327, 1330, 1333, 1336, 1339, 1342, 1345, 1348, 1351, 1354, 1357, 1360, 1363, 1366, 1369, 1372, 1375, 1378, 1381, 1384, 1387, 1390, 1393, 1396, 1399, 1402, 1405, 1408, 1411, 1414, 1417, 1420, 1423, 1426, 1429, 1432, 1435, 1438, 1441, 1444, 1447, 1450, 1453, 1456, 1459, 1462, 1465, 1468, 1471, 1474, 1477, 1480, 1483, 1486, 1489, 1492, 1495, 1498, 1501, 1504, 1507, 1510, 1513, 1516, 1519, 1522, 1525, 1528, 1531, 1534, 1537, 1540, 1543, 1546, 1549, 1552, 1555, 1558, 1561, 1564, 1567, 1570, 1573, 1576, 1579, 1582, 1585, 1588, 1591, 1594, 1597, 1600, 1603, 1606, 1609, 1612, 1615, 1618, 1621, 1624, 1627, 1630, 1633, 1636, 1639, 1642, 1645, 1648, 1651, 1654, 1657, 1660, 1663, 1666, 1669, 1672, 1675, 1678, 1681, 1684, 1687, 1690, 1693, 1696, 1699, 1702, 1705, 1708, 1711, 1714, 1717, 1720, 1723, 1726, 1729, 1732, 1735, 1738, 1741, 1744, 1747, 1750, 1753, 1756, 1759, 1762, 1765, 1768, 1771, 1774, 1777, 1780, 1783, 1786, 1789, 1792, 1795, 1798, 1801, 1804, 1807, 1810, 1813, 1816, 1819, 1822, 1825, 1828, 1831, 1834, 1837, 1840, 1843, 1846, 1849, 1852, 1855, 1858, 1861, 1864, 1867, 1870, 1873, 1876, 1879, 1882, 1885, 1888, 1891, 1894, 1897, 1900, 1903, 1906, 1909, 1912, 1915, 1918, 1921, 1924, 1927, 1930, 1933, 1936, 1939, 1942, 1945, 1948, 1951, 1954, 1957, 1960, 1963, 1966, 1969, 1972, 1975, 1978, 1981, 1984, 1987, 1990, 1993, 1996, 1999, 2002, 2005, 2008, 2011, 2014, 2017, 2020, 2023, 2026, 2029, 2032, 2035, 2038, 2041, 2044, 2047, 2050, 2053, 2056, 2059, 2062, 2065, 2068, 2071, 2074, 2077, 2080, 2083, 2086, 2089, 2092, 2095, 2098, 2101, 2104, 2107, 2110, 2113, 2116, 2119, 2122, 2125, 2128, 2131, 2134, 2137, 2140, 2143, 2146, 2149, 2152, 2155, 2158, 2161, 2164, 2167, 2176, 2206, 2236, 2238, 2299, 2301, 2331, 2333, 2335, 2367, 2369, 2399, 2402, 2407, 2437, 2439, 2442, 2444, 2447, 2452, 2513, 2516, 2518, 2520, 2524, 2528, 2530, 2534, 2538, 2543, 2548, 2550, 2553, 2558, 2560, 2563, 2567, 2570, 2573, 2576, 2579, 2581, 2583, 2584, 2586, 2589, 2592, 2593, 2653, 2659, 2661, 2664, 2667, 2669, 2677, 2679, 2682, 2685, 2687, 2690, 2693, 2696, 2699, 2701, 2703, 2707, 2709, 2713, 2715, 2718, 2720, 2723, 2725, 2727, 2731, 2745, 2748, 2750, 2751, 2753, 2755, 2757, 2760, 2763, 2765, 2767, 2770, 2773, 2776, 2779, 2782, 2785, 2788, 2791, 2794, 2797, 2800, 2804, 2807, 2810, 2813, 2816, 2819, 2822, 2825, 2828, 2831, 2832, 2835, 2838, 2841, 2844, 2847, 2850, 2853, 2857, 2861, 2864, 2867, 2993, 2995, 2999, 3002, 3005, 3007, 3008, 3010, 3011, 3013, 3015, 3017, 3019, 3023, 3026, 3028, 3031, 3034, 3037, 3040
};

void dwdx_colptrs_SPARCED(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED