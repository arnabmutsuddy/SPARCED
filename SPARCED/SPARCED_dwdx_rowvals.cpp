#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 3049> dwdx_rowvals_SPARCED_ = {
    1, 143, 254, 255, 257, 256, 258, 272, 292, 313, 144, 254, 258, 261, 262, 263, 311, 1864, 145, 256, 265, 267, 255, 261, 267, 268, 255, 262, 270, 408, 271, 274, 273, 276, 275, 278, 277, 280, 279, 282, 281, 284, 283, 286, 285, 288, 287, 290, 3, 289, 291, 294, 293, 296, 295, 298, 297, 300, 299, 302, 301, 304, 303, 306, 305, 308, 307, 310, 4, 309, 146, 315, 147, 316, 148, 317, 266, 315, 269, 316, 317, 1865, 1866, 312, 1867, 149, 311, 312, 150, 313, 314, 151, 152, 11, 12, 13, 14, 15, 20, 21, 153, 319, 323, 11, 12, 13, 14, 15, 20, 21, 327, 329, 330, 11, 12, 13, 14, 15, 20, 21, 154, 319, 327, 332, 334, 155, 321, 337, 339, 340, 423, 322, 323, 325, 337, 338, 342, 424, 431, 337, 341, 343, 414, 432, 337, 344, 433, 156, 346, 348, 410, 425, 346, 347, 350, 362, 426, 434, 318, 346, 349, 351, 358, 362, 370, 416, 435, 157, 348, 354, 355, 372, 346, 352, 362, 436, 350, 357, 360, 437, 158, 381, 382, 159, 358, 359, 160, 400, 401, 161, 362, 364, 365, 427, 346, 362, 363, 367, 428, 438, 318, 332, 346, 362, 366, 368, 374, 381, 418, 439, 346, 362, 369, 440, 162, 343, 351, 368, 370, 391, 413, 375, 376, 441, 163, 354, 374, 378, 387, 333, 335, 442, 371, 372, 373, 443, 367, 380, 383, 444, 164, 385, 387, 388, 429, 385, 386, 390, 430, 445, 374, 385, 389, 391, 394, 400, 404, 420, 446, 165, 394, 396, 364, 387, 395, 397, 447, 390, 399, 402, 448, 166, 349, 366, 389, 404, 406, 405, 407, 449, 385, 392, 450, 167, 357, 380, 399, 408, 409, 320, 325, 451, 328, 452, 168, 414, 416, 418, 420, 422, 337, 415, 453, 346, 362, 417, 454, 346, 362, 419, 455, 385, 421, 456, 169, 466, 170, 466, 467, 468, 1879, 469, 471, 1880, 171, 469, 470, 1881, 172, 471, 482, 552, 472, 473, 1882, 474, 476, 491, 1883, 173, 474, 475, 1884, 174, 476, 527, 477, 478, 1885, 479, 485, 488, 1886, 175, 479, 480, 481, 1887, 482, 1888, 483, 484, 1889, 176, 485, 532, 534, 486, 487, 1890, 177, 488, 489, 490, 1891, 1892, 178, 491, 492, 493, 1893, 494, 496, 1894, 179, 494, 536, 538, 540, 545, 550, 495, 1895, 180, 496, 542, 497, 498, 1896, 499, 1897, 500, 501, 503, 1898, 501, 505, 509, 551, 1899, 502, 1900, 504, 505, 507, 1901, 506, 1902, 508, 509, 511, 1903, 510, 1904, 511, 512, 513, 514, 517, 547, 181, 514, 515, 516, 520, 549, 1905, 182, 517, 518, 519, 530, 548, 1906, 521, 522, 1907, 183, 522, 523, 524, 1908, 525, 1909, 184, 525, 526, 527, 532, 1910, 528, 529, 1911, 531, 534, 1912, 533, 1913, 535, 1914, 1915, 185, 536, 1872, 1875, 186, 538, 187, 540, 537, 1916, 539, 1917, 541, 1918, 188, 542, 545, 1868, 543, 544, 1919, 546, 1920, 1869, 1870, 1921, 1871, 1922, 1873, 1874, 1923, 1878, 1924, 1876, 1877, 1925, 189, 553, 559, 583, 605, 675, 677, 679, 683, 190, 565, 573, 577, 587, 595, 599, 609, 615, 681, 685, 687, 689, 691, 693, 191, 625, 629, 697, 192, 633, 193, 639, 644, 701, 194, 703, 705, 195, 707, 709, 196, 553, 557, 581, 603, 619, 657, 659, 661, 663, 620, 1926, 197, 555, 579, 601, 621, 659, 665, 667, 669, 622, 1927, 198, 563, 577, 585, 607, 661, 667, 671, 199, 571, 593, 599, 613, 623, 663, 669, 671, 673, 624, 1928, 200, 569, 591, 611, 201, 625, 627, 695, 202, 633, 203, 639, 646, 699, 204, 647, 205, 649, 658, 675, 1929, 660, 677, 1930, 662, 679, 681, 1931, 664, 683, 685, 1932, 666, 1933, 668, 687, 1934, 670, 689, 1935, 672, 691, 1936, 674, 693, 1937, 696, 697, 1938, 700, 701, 1939, 648, 703, 1940, 650, 707, 1941, 554, 555, 557, 562, 563, 568, 569, 571, 576, 1942, 568, 578, 579, 581, 585, 590, 591, 593, 598, 1943, 576, 598, 600, 601, 603, 607, 611, 613, 618, 1944, 626, 627, 632, 1945, 634, 635, 1946, 640, 641, 646, 1947, 556, 678, 711, 769, 1948, 570, 712, 770, 1949, 558, 559, 676, 713, 771, 1950, 564, 565, 680, 715, 773, 1951, 572, 573, 684, 717, 775, 1952, 580, 688, 719, 777, 1953, 582, 583, 682, 721, 779, 1954, 586, 587, 1955, 592, 720, 778, 1956, 594, 595, 722, 780, 1957, 602, 690, 724, 782, 1958, 612, 725, 783, 1959, 604, 605, 686, 726, 784, 1960, 608, 609, 692, 727, 785, 1961, 614, 615, 694, 728, 786, 1962, 628, 629, 698, 730, 788, 1963, 638, 733, 791, 1964, 644, 645, 702, 735, 793, 1965, 704, 705, 736, 794, 1966, 708, 709, 737, 795, 1967, 560, 561, 714, 772, 1968, 566, 567, 584, 716, 774, 1969, 574, 575, 606, 718, 776, 1970, 588, 589, 1971, 596, 597, 610, 723, 781, 1972, 616, 617, 729, 787, 1973, 630, 631, 731, 789, 1974, 636, 637, 732, 790, 1975, 642, 643, 734, 792, 1976, 706, 738, 796, 1977, 710, 739, 797, 1978, 651, 652, 1979, 653, 654, 1980, 655, 656, 1981, 740, 827, 1001, 1059, 1117, 1175, 1982, 741, 828, 1002, 1060, 1118, 1176, 1983, 742, 829, 1003, 1061, 1119, 1177, 1984, 743, 830, 1004, 1062, 1120, 1178, 1985, 744, 831, 1005, 1063, 1121, 1179, 1986, 745, 832, 1006, 1064, 1122, 1180, 1987, 746, 833, 1007, 1065, 1123, 1181, 1988, 747, 834, 1008, 1066, 1124, 1182, 1989, 748, 835, 1009, 1067, 1125, 1183, 1990, 749, 836, 1010, 1068, 1126, 1184, 1991, 750, 837, 1011, 1069, 1127, 1185, 1992, 751, 838, 1012, 1070, 1128, 1186, 1993, 752, 839, 1013, 1071, 1129, 1187, 1994, 753, 840, 1014, 1072, 1130, 1188, 1995, 754, 841, 1015, 1073, 1131, 1189, 1996, 755, 842, 1016, 1074, 1132, 1190, 1997, 756, 843, 1017, 1075, 1133, 1191, 1998, 757, 844, 1018, 1076, 1134, 1192, 1999, 758, 845, 1019, 1077, 1135, 1193, 2000, 759, 846, 1020, 1078, 1136, 1194, 2001, 760, 847, 1021, 1079, 1137, 1195, 2002, 761, 848, 1022, 1080, 1138, 1196, 2003, 762, 849, 1023, 1081, 1139, 1197, 2004, 763, 850, 1024, 1082, 1140, 1198, 2005, 764, 851, 1025, 1083, 1141, 1199, 2006, 765, 852, 1668, 2007, 766, 853, 1670, 2008, 767, 854, 1672, 2009, 768, 855, 1674, 2010, 1026, 1084, 1142, 1200, 1669, 2011, 1027, 1085, 1143, 1201, 1671, 2012, 1028, 1086, 1144, 1202, 1673, 2013, 1029, 1087, 1145, 1203, 1675, 2014, 798, 2015, 799, 2016, 800, 2017, 801, 2018, 802, 2019, 803, 2020, 804, 2021, 805, 2022, 806, 2023, 807, 2024, 808, 2025, 809, 2026, 810, 2027, 811, 2028, 812, 2029, 813, 2030, 814, 2031, 815, 2032, 816, 2033, 817, 2034, 818, 2035, 819, 2036, 820, 2037, 821, 2038, 822, 2039, 823, 2040, 824, 2041, 825, 2042, 826, 2043, 885, 914, 2044, 886, 915, 2045, 887, 916, 2046, 888, 917, 2047, 889, 918, 2048, 890, 919, 2049, 891, 920, 2050, 892, 921, 2051, 893, 922, 2052, 894, 923, 2053, 895, 924, 2054, 896, 925, 2055, 897, 926, 2056, 898, 927, 2057, 899, 928, 2058, 900, 929, 2059, 901, 930, 2060, 902, 931, 2061, 903, 932, 2062, 904, 933, 2063, 905, 934, 2064, 906, 935, 2065, 907, 936, 2066, 908, 937, 2067, 909, 938, 2068, 910, 939, 2069, 911, 940, 2070, 912, 941, 2071, 913, 942, 2072, 856, 943, 2073, 857, 944, 2074, 858, 945, 2075, 859, 946, 2076, 860, 947, 2077, 861, 948, 2078, 862, 949, 2079, 863, 950, 2080, 864, 951, 2081, 865, 952, 2082, 866, 953, 2083, 867, 954, 2084, 868, 955, 2085, 869, 956, 2086, 870, 957, 2087, 871, 958, 2088, 872, 959, 2089, 873, 960, 2090, 874, 961, 2091, 875, 962, 2092, 876, 963, 2093, 877, 964, 2094, 878, 965, 2095, 879, 966, 2096, 880, 967, 2097, 881, 1676, 2098, 882, 1678, 2099, 883, 1680, 2100, 884, 1682, 2101, 968, 1677, 2102, 969, 1679, 2103, 970, 1681, 2104, 971, 1683, 2105, 1030, 1320, 2106, 1031, 1321, 2107, 1032, 1322, 2108, 1033, 1323, 2109, 1034, 1324, 2110, 1035, 1325, 2111, 1036, 1326, 2112, 1037, 1327, 2113, 1038, 1328, 2114, 1039, 1329, 2115, 1040, 1330, 2116, 1041, 1331, 2117, 1042, 1332, 2118, 1043, 1333, 2119, 1044, 1334, 2120, 1045, 1335, 2121, 1046, 1336, 2122, 1047, 1337, 2123, 1048, 1338, 2124, 1049, 1339, 2125, 1050, 1340, 2126, 1051, 1341, 2127, 1052, 1342, 2128, 1053, 1343, 2129, 1054, 1344, 2130, 1055, 1345, 2131, 1056, 1346, 2132, 1057, 1347, 2133, 1058, 1348, 2134, 972, 1233, 2135, 973, 1234, 2136, 974, 1235, 2137, 975, 1236, 2138, 976, 1237, 2139, 977, 1238, 2140, 978, 1239, 2141, 979, 1240, 2142, 980, 1241, 2143, 981, 1242, 2144, 982, 1243, 2145, 983, 1244, 2146, 984, 1245, 2147, 985, 1246, 2148, 986, 1247, 2149, 987, 1248, 2150, 988, 1249, 2151, 989, 1250, 2152, 990, 1251, 2153, 991, 1252, 2154, 992, 1253, 2155, 993, 1254, 2156, 994, 1255, 2157, 995, 1256, 2158, 996, 1257, 2159, 997, 1258, 2160, 998, 1259, 2161, 999, 1260, 2162, 1000, 1261, 2163, 1088, 1407, 2164, 1089, 1408, 2165, 1090, 1409, 2166, 1091, 1410, 2167, 1092, 1411, 2168, 1093, 1412, 2169, 1094, 1413, 2170, 1095, 1414, 2171, 1096, 1415, 2172, 1097, 1416, 2173, 1098, 1417, 2174, 1099, 1418, 2175, 1100, 1419, 2176, 1101, 1420, 2177, 1102, 1421, 2178, 1103, 1422, 2179, 1104, 1423, 2180, 1105, 1424, 2181, 1106, 1425, 2182, 1107, 1426, 2183, 1108, 1427, 2184, 1109, 1428, 2185, 1110, 1429, 2186, 1111, 1430, 2187, 1112, 1431, 2188, 1113, 1432, 2189, 1114, 1433, 2190, 1115, 1434, 2191, 1116, 1435, 2192, 1146, 1494, 2193, 1147, 1495, 2194, 1148, 1496, 2195, 1149, 1497, 2196, 1150, 1498, 2197, 1151, 1499, 2198, 1152, 1500, 2199, 1153, 1501, 2200, 1154, 1502, 2201, 1155, 1503, 2202, 1156, 1504, 2203, 1157, 1505, 2204, 1158, 1506, 2205, 1159, 1507, 2206, 1160, 1508, 2207, 1161, 1509, 2208, 1162, 1510, 2209, 1163, 1511, 2210, 1164, 1512, 2211, 1165, 1513, 2212, 1166, 1514, 2213, 1167, 1515, 2214, 1168, 1516, 2215, 1169, 1517, 2216, 1170, 1518, 2217, 1171, 1519, 2218, 1172, 1520, 2219, 1173, 1521, 2220, 1174, 1522, 2221, 1204, 1581, 2222, 1205, 1582, 2223, 1206, 1583, 2224, 1207, 1584, 2225, 1208, 1585, 2226, 1209, 1586, 2227, 1210, 1587, 2228, 1211, 1588, 2229, 1212, 1589, 2230, 1213, 1590, 2231, 1214, 1591, 2232, 1215, 1592, 2233, 1216, 1593, 2234, 1217, 1594, 2235, 1218, 1595, 2236, 1219, 1596, 2237, 1220, 1597, 2238, 1221, 1598, 2239, 1222, 1599, 2240, 1223, 1600, 2241, 1224, 1601, 2242, 1225, 1602, 2243, 1226, 1603, 2244, 1227, 1604, 2245, 1228, 1605, 2246, 1229, 1606, 2247, 1230, 1607, 2248, 1231, 1608, 2249, 1232, 1609, 2250, 1262, 1291, 2251, 1263, 1292, 2252, 1264, 1293, 2253, 1265, 1294, 2254, 1266, 1295, 2255, 1267, 1296, 2256, 1268, 1297, 2257, 1269, 1298, 2258, 1270, 1299, 2259, 1271, 1300, 2260, 1272, 1301, 2261, 1273, 1302, 2262, 1274, 1303, 2263, 1275, 1304, 2264, 1276, 1305, 2265, 1277, 1306, 2266, 1278, 1307, 2267, 1279, 1308, 2268, 1280, 1309, 2269, 1281, 1310, 2270, 1282, 1311, 2271, 1283, 1312, 2272, 1284, 1313, 2273, 1285, 1314, 2274, 1286, 1315, 2275, 1287, 1316, 2276, 1288, 1317, 2277, 1289, 1318, 2278, 1290, 1319, 2279, 1349, 1378, 2280, 1350, 1379, 2281, 1351, 1380, 2282, 1352, 1381, 2283, 1353, 1382, 2284, 1354, 1383, 2285, 1355, 1384, 2286, 1356, 1385, 2287, 1357, 1386, 2288, 1358, 1387, 2289, 1359, 1388, 2290, 1360, 1389, 2291, 1361, 1390, 2292, 1362, 1391, 2293, 1363, 1392, 2294, 1364, 1393, 2295, 1365, 1394, 2296, 1366, 1395, 2297, 1367, 1396, 2298, 1368, 1397, 2299, 1369, 1398, 2300, 1370, 1399, 2301, 1371, 1400, 2302, 1372, 1401, 2303, 1373, 1402, 2304, 1374, 1403, 2305, 1375, 1404, 2306, 1376, 1405, 2307, 1377, 1406, 2308, 1436, 1465, 2309, 1437, 1466, 2310, 1438, 1467, 2311, 1439, 1468, 2312, 1440, 1469, 2313, 1441, 1470, 2314, 1442, 1471, 2315, 1443, 1472, 2316, 1444, 1473, 2317, 1445, 1474, 2318, 1446, 1475, 2319, 1447, 1476, 2320, 1448, 1477, 2321, 1449, 1478, 2322, 1450, 1479, 2323, 1451, 1480, 2324, 1452, 1481, 2325, 1453, 1482, 2326, 1454, 1483, 2327, 1455, 1484, 2328, 1456, 1485, 2329, 1457, 1486, 2330, 1458, 1487, 2331, 1459, 1488, 2332, 1460, 1489, 2333, 1461, 1490, 2334, 1462, 1491, 2335, 1463, 1492, 2336, 1464, 1493, 2337, 1523, 1552, 2338, 1524, 1553, 2339, 1525, 1554, 2340, 1526, 1555, 2341, 1527, 1556, 2342, 1528, 1557, 2343, 1529, 1558, 2344, 1530, 1559, 2345, 1531, 1560, 2346, 1532, 1561, 2347, 1533, 1562, 2348, 1534, 1563, 2349, 1535, 1564, 2350, 1536, 1565, 2351, 1537, 1566, 2352, 1538, 1567, 2353, 1539, 1568, 2354, 1540, 1569, 2355, 1541, 1570, 2356, 1542, 1571, 2357, 1543, 1572, 2358, 1544, 1573, 2359, 1545, 1574, 2360, 1546, 1575, 2361, 1547, 1576, 2362, 1548, 1577, 2363, 1549, 1578, 2364, 1550, 1579, 2365, 1551, 1580, 2366, 1610, 1639, 2367, 1611, 1640, 2368, 1612, 1641, 2369, 1613, 1642, 2370, 1614, 1643, 2371, 1615, 1644, 2372, 1616, 1645, 2373, 1617, 1646, 2374, 1618, 1647, 2375, 1619, 1648, 2376, 1620, 1649, 2377, 1621, 1650, 2378, 1622, 1651, 2379, 1623, 1652, 2380, 1624, 1653, 2381, 1625, 1654, 2382, 1626, 1655, 2383, 1627, 1656, 2384, 1628, 1657, 2385, 1629, 1658, 2386, 1630, 1659, 2387, 1631, 1660, 2388, 1632, 1661, 2389, 1633, 1662, 2390, 1634, 1663, 2391, 1635, 1664, 2392, 1636, 1665, 2393, 1637, 1666, 2394, 1638, 1667, 2395, 206, 1668, 1670, 1672, 1674, 1676, 1678, 1680, 1682, 207, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 208, 827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 209, 1842, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1727, 1843, 2396, 1730, 2397, 210, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 211, 1839, 212, 1839, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1835, 1840, 2398, 1838, 2399, 213, 1175, 1176, 1177, 1178, 1179, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1198, 1199, 1200, 1201, 1202, 1203, 1819, 1822, 2400, 1820, 1823, 1827, 1835, 2401, 1581, 1582, 1583, 1584, 1585, 1586, 1587, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 1595, 1596, 1597, 1598, 1599, 1600, 1601, 1602, 1603, 1604, 1605, 1606, 1607, 1608, 1609, 2402, 1841, 2403, 1774, 1776, 2404, 214, 1774, 1735, 1775, 2405, 1684, 1686, 1738, 1859, 2406, 215, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 1346, 1347, 1348, 1735, 1858, 216, 1738, 1741, 1744, 2407, 1734, 2408, 217, 1684, 1731, 1782, 1685, 1688, 1692, 2409, 218, 1686, 1693, 1706, 1709, 2410, 219, 1694, 1700, 1706, 1697, 1703, 1709, 1712, 2411, 1713, 1714, 1717, 1860, 2412, 220, 1745, 1756, 1760, 2413, 1720, 1757, 1769, 1850, 2414, 221, 1751, 1754, 1761, 2415, 1726, 1762, 1766, 2416, 222, 1720, 1762, 1723, 1757, 2417, 223, 1766, 1769, 1724, 1772, 2418, 224, 1772, 1773, 2419, 225, 1810, 2420, 226, 1807, 1809, 1808, 1847, 2421, 2422, 1407, 1408, 1409, 1410, 1411, 1412, 1413, 1414, 1415, 1416, 1417, 1418, 1419, 1420, 1421, 1422, 1423, 1424, 1425, 1426, 1427, 1428, 1429, 1430, 1431, 1432, 1433, 1434, 1435, 1494, 1495, 1496, 1497, 1498, 1499, 1500, 1501, 1502, 1503, 1504, 1505, 1506, 1507, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1517, 1518, 1519, 1520, 1521, 1522, 1856, 2423, 1784, 1787, 1789, 1802, 1857, 2424, 227, 1784, 1788, 1793, 2425, 228, 1787, 1862, 1852, 2426, 1802, 1803, 1811, 1831, 1851, 1864, 1872, 2427, 229, 1789, 1790, 1793, 2428, 1796, 1797, 2429, 230, 1791, 231, 1791, 1821, 1792, 1797, 2430, 1800, 1801, 2431, 232, 1803, 1807, 1806, 2432, 233, 1817, 234, 1811, 1814, 1817, 1725, 2433, 1818, 1820, 1847, 2434, 235, 1776, 1777, 1778, 2435, 1781, 2436, 236, 1778, 1782, 1783, 2437, 237, 1714, 1717, 1745, 1848, 2438, 619, 621, 623, 1727, 1731, 1741, 1748, 1751, 1755, 1814, 1849, 1868, 1875, 2439, 238, 1831, 1854, 1834, 2440, 239, 1853, 2441, 240, 1821, 241, 1827, 0, 1830, 2442, 242, 1823, 1844, 1826, 2443, 243, 1842, 1728, 1729, 2444, 1732, 1733, 2445, 1736, 1737, 2446, 1739, 1740, 2447, 1742, 1743, 2448, 1707, 1708, 2449, 1710, 1711, 2450, 1715, 1716, 2451, 1718, 1719, 2452, 1752, 1753, 2453, 1763, 1764, 2454, 1721, 1722, 1765, 2455, 1767, 1768, 2456, 1770, 1771, 2457, 1779, 1780, 2458, 1785, 1786, 2459, 1794, 1795, 2460, 1798, 1799, 2461, 1804, 1805, 2462, 1812, 1813, 2463, 1815, 1816, 2464, 2465, 1824, 1825, 2466, 1828, 1829, 2467, 1832, 1833, 2468, 1836, 1837, 2469, 1746, 1747, 2470, 1749, 1750, 2471, 1758, 1759, 2472, 1687, 1690, 1692, 2473, 1691, 1700, 1703, 2474, 1701, 1702, 2475, 1704, 1705, 2476, 2, 3, 4, 5, 6, 11, 12, 13, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 244, 248, 1844, 1845, 2477, 1689, 1694, 1697, 2478, 1695, 1696, 2479, 1698, 1699, 2480, 1855, 2481, 1860, 1861, 2482, 1862, 1863, 2483, 248, 251, 249, 252, 245, 429, 246, 425, 427, 247, 321, 423, 457, 324, 458, 326, 459, 460, 461, 462, 463, 464, 465
};

void dwdx_rowvals_SPARCED(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED