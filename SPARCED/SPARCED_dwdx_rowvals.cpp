#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 2995> dwdx_rowvals_SPARCED_ = {
    1, 143, 254, 255, 257, 256, 258, 272, 292, 313, 144, 254, 258, 261, 262, 263, 311, 1830, 145, 256, 265, 267, 255, 261, 267, 268, 255, 262, 270, 271, 274, 273, 276, 275, 278, 277, 280, 279, 282, 281, 284, 283, 286, 285, 288, 287, 290, 3, 289, 291, 294, 293, 296, 295, 298, 297, 300, 299, 302, 301, 304, 303, 306, 305, 308, 307, 310, 4, 309, 146, 315, 147, 316, 148, 317, 266, 315, 269, 316, 317, 1831, 1832, 312, 1833, 149, 311, 312, 150, 313, 314, 151, 152, 21, 153, 319, 323, 21, 328, 333, 430, 431, 21, 154, 319, 328, 155, 321, 340, 322, 323, 325, 342, 399, 341, 343, 400, 344, 401, 156, 331, 346, 348, 332, 333, 336, 346, 347, 350, 362, 402, 318, 346, 349, 351, 358, 362, 370, 403, 157, 348, 354, 355, 372, 346, 352, 362, 404, 350, 357, 360, 405, 158, 381, 382, 159, 358, 359, 160, 161, 362, 364, 365, 346, 362, 363, 367, 406, 318, 346, 362, 366, 368, 374, 381, 407, 346, 362, 369, 408, 162, 343, 351, 368, 370, 391, 375, 376, 409, 163, 354, 374, 378, 387, 410, 371, 372, 373, 411, 367, 380, 383, 412, 164, 385, 387, 388, 385, 386, 390, 413, 374, 385, 389, 391, 394, 414, 165, 394, 396, 364, 387, 395, 397, 415, 390, 416, 166, 349, 366, 389, 417, 385, 392, 418, 167, 357, 380, 320, 325, 419, 329, 336, 420, 168, 421, 346, 362, 422, 346, 362, 423, 385, 424, 169, 432, 170, 432, 433, 434, 1845, 435, 437, 1846, 171, 435, 436, 1847, 172, 437, 448, 518, 438, 439, 1848, 440, 442, 457, 1849, 173, 440, 441, 1850, 174, 442, 493, 443, 444, 1851, 445, 451, 454, 1852, 175, 445, 446, 447, 1853, 448, 1854, 449, 450, 1855, 176, 451, 498, 500, 452, 453, 1856, 177, 454, 455, 456, 1857, 1858, 178, 457, 458, 459, 1859, 460, 462, 1860, 179, 460, 502, 504, 506, 511, 516, 461, 1861, 180, 462, 508, 463, 464, 1862, 465, 1863, 466, 467, 469, 1864, 467, 471, 475, 517, 1865, 468, 1866, 470, 471, 473, 1867, 472, 1868, 474, 475, 477, 1869, 476, 1870, 477, 478, 479, 480, 483, 513, 181, 480, 481, 482, 486, 515, 1871, 182, 483, 484, 485, 496, 514, 1872, 487, 488, 1873, 183, 488, 489, 490, 1874, 491, 1875, 184, 491, 492, 493, 498, 1876, 494, 495, 1877, 497, 500, 1878, 499, 1879, 501, 1880, 1881, 185, 502, 1838, 1841, 186, 504, 187, 506, 503, 1882, 505, 1883, 507, 1884, 188, 508, 511, 1834, 509, 510, 1885, 512, 1886, 1835, 1836, 1887, 1837, 1888, 1839, 1840, 1889, 1844, 1890, 1842, 1843, 1891, 189, 519, 525, 549, 571, 641, 643, 645, 649, 190, 531, 539, 543, 553, 561, 565, 575, 581, 647, 651, 653, 655, 657, 659, 191, 591, 595, 663, 192, 599, 193, 605, 610, 667, 194, 669, 671, 195, 673, 675, 196, 519, 523, 547, 569, 585, 623, 625, 627, 629, 586, 1892, 197, 521, 545, 567, 587, 625, 631, 633, 635, 588, 1893, 198, 529, 543, 551, 573, 627, 633, 637, 199, 537, 559, 565, 579, 589, 629, 635, 637, 639, 590, 1894, 200, 535, 557, 577, 201, 591, 593, 661, 202, 599, 203, 605, 612, 665, 204, 613, 205, 615, 624, 641, 1895, 626, 643, 1896, 628, 645, 647, 1897, 630, 649, 651, 1898, 632, 1899, 634, 653, 1900, 636, 655, 1901, 638, 657, 1902, 640, 659, 1903, 662, 663, 1904, 666, 667, 1905, 614, 669, 1906, 616, 673, 1907, 520, 521, 523, 528, 529, 534, 535, 537, 542, 1908, 534, 544, 545, 547, 551, 556, 557, 559, 564, 1909, 542, 564, 566, 567, 569, 573, 577, 579, 584, 1910, 592, 593, 598, 1911, 600, 601, 1912, 606, 607, 612, 1913, 522, 644, 677, 735, 1914, 536, 678, 736, 1915, 524, 525, 642, 679, 737, 1916, 530, 531, 646, 681, 739, 1917, 538, 539, 650, 683, 741, 1918, 546, 654, 685, 743, 1919, 548, 549, 648, 687, 745, 1920, 552, 553, 1921, 558, 686, 744, 1922, 560, 561, 688, 746, 1923, 568, 656, 690, 748, 1924, 578, 691, 749, 1925, 570, 571, 652, 692, 750, 1926, 574, 575, 658, 693, 751, 1927, 580, 581, 660, 694, 752, 1928, 594, 595, 664, 696, 754, 1929, 604, 699, 757, 1930, 610, 611, 668, 701, 759, 1931, 670, 671, 702, 760, 1932, 674, 675, 703, 761, 1933, 526, 527, 680, 738, 1934, 532, 533, 550, 682, 740, 1935, 540, 541, 572, 684, 742, 1936, 554, 555, 1937, 562, 563, 576, 689, 747, 1938, 582, 583, 695, 753, 1939, 596, 597, 697, 755, 1940, 602, 603, 698, 756, 1941, 608, 609, 700, 758, 1942, 672, 704, 762, 1943, 676, 705, 763, 1944, 617, 618, 1945, 619, 620, 1946, 621, 622, 1947, 706, 793, 967, 1025, 1083, 1141, 1948, 707, 794, 968, 1026, 1084, 1142, 1949, 708, 795, 969, 1027, 1085, 1143, 1950, 709, 796, 970, 1028, 1086, 1144, 1951, 710, 797, 971, 1029, 1087, 1145, 1952, 711, 798, 972, 1030, 1088, 1146, 1953, 712, 799, 973, 1031, 1089, 1147, 1954, 713, 800, 974, 1032, 1090, 1148, 1955, 714, 801, 975, 1033, 1091, 1149, 1956, 715, 802, 976, 1034, 1092, 1150, 1957, 716, 803, 977, 1035, 1093, 1151, 1958, 717, 804, 978, 1036, 1094, 1152, 1959, 718, 805, 979, 1037, 1095, 1153, 1960, 719, 806, 980, 1038, 1096, 1154, 1961, 720, 807, 981, 1039, 1097, 1155, 1962, 721, 808, 982, 1040, 1098, 1156, 1963, 722, 809, 983, 1041, 1099, 1157, 1964, 723, 810, 984, 1042, 1100, 1158, 1965, 724, 811, 985, 1043, 1101, 1159, 1966, 725, 812, 986, 1044, 1102, 1160, 1967, 726, 813, 987, 1045, 1103, 1161, 1968, 727, 814, 988, 1046, 1104, 1162, 1969, 728, 815, 989, 1047, 1105, 1163, 1970, 729, 816, 990, 1048, 1106, 1164, 1971, 730, 817, 991, 1049, 1107, 1165, 1972, 731, 818, 1634, 1973, 732, 819, 1636, 1974, 733, 820, 1638, 1975, 734, 821, 1640, 1976, 992, 1050, 1108, 1166, 1635, 1977, 993, 1051, 1109, 1167, 1637, 1978, 994, 1052, 1110, 1168, 1639, 1979, 995, 1053, 1111, 1169, 1641, 1980, 764, 1981, 765, 1982, 766, 1983, 767, 1984, 768, 1985, 769, 1986, 770, 1987, 771, 1988, 772, 1989, 773, 1990, 774, 1991, 775, 1992, 776, 1993, 777, 1994, 778, 1995, 779, 1996, 780, 1997, 781, 1998, 782, 1999, 783, 2000, 784, 2001, 785, 2002, 786, 2003, 787, 2004, 788, 2005, 789, 2006, 790, 2007, 791, 2008, 792, 2009, 851, 880, 2010, 852, 881, 2011, 853, 882, 2012, 854, 883, 2013, 855, 884, 2014, 856, 885, 2015, 857, 886, 2016, 858, 887, 2017, 859, 888, 2018, 860, 889, 2019, 861, 890, 2020, 862, 891, 2021, 863, 892, 2022, 864, 893, 2023, 865, 894, 2024, 866, 895, 2025, 867, 896, 2026, 868, 897, 2027, 869, 898, 2028, 870, 899, 2029, 871, 900, 2030, 872, 901, 2031, 873, 902, 2032, 874, 903, 2033, 875, 904, 2034, 876, 905, 2035, 877, 906, 2036, 878, 907, 2037, 879, 908, 2038, 822, 909, 2039, 823, 910, 2040, 824, 911, 2041, 825, 912, 2042, 826, 913, 2043, 827, 914, 2044, 828, 915, 2045, 829, 916, 2046, 830, 917, 2047, 831, 918, 2048, 832, 919, 2049, 833, 920, 2050, 834, 921, 2051, 835, 922, 2052, 836, 923, 2053, 837, 924, 2054, 838, 925, 2055, 839, 926, 2056, 840, 927, 2057, 841, 928, 2058, 842, 929, 2059, 843, 930, 2060, 844, 931, 2061, 845, 932, 2062, 846, 933, 2063, 847, 1642, 2064, 848, 1644, 2065, 849, 1646, 2066, 850, 1648, 2067, 934, 1643, 2068, 935, 1645, 2069, 936, 1647, 2070, 937, 1649, 2071, 996, 1286, 2072, 997, 1287, 2073, 998, 1288, 2074, 999, 1289, 2075, 1000, 1290, 2076, 1001, 1291, 2077, 1002, 1292, 2078, 1003, 1293, 2079, 1004, 1294, 2080, 1005, 1295, 2081, 1006, 1296, 2082, 1007, 1297, 2083, 1008, 1298, 2084, 1009, 1299, 2085, 1010, 1300, 2086, 1011, 1301, 2087, 1012, 1302, 2088, 1013, 1303, 2089, 1014, 1304, 2090, 1015, 1305, 2091, 1016, 1306, 2092, 1017, 1307, 2093, 1018, 1308, 2094, 1019, 1309, 2095, 1020, 1310, 2096, 1021, 1311, 2097, 1022, 1312, 2098, 1023, 1313, 2099, 1024, 1314, 2100, 938, 1199, 2101, 939, 1200, 2102, 940, 1201, 2103, 941, 1202, 2104, 942, 1203, 2105, 943, 1204, 2106, 944, 1205, 2107, 945, 1206, 2108, 946, 1207, 2109, 947, 1208, 2110, 948, 1209, 2111, 949, 1210, 2112, 950, 1211, 2113, 951, 1212, 2114, 952, 1213, 2115, 953, 1214, 2116, 954, 1215, 2117, 955, 1216, 2118, 956, 1217, 2119, 957, 1218, 2120, 958, 1219, 2121, 959, 1220, 2122, 960, 1221, 2123, 961, 1222, 2124, 962, 1223, 2125, 963, 1224, 2126, 964, 1225, 2127, 965, 1226, 2128, 966, 1227, 2129, 1054, 1373, 2130, 1055, 1374, 2131, 1056, 1375, 2132, 1057, 1376, 2133, 1058, 1377, 2134, 1059, 1378, 2135, 1060, 1379, 2136, 1061, 1380, 2137, 1062, 1381, 2138, 1063, 1382, 2139, 1064, 1383, 2140, 1065, 1384, 2141, 1066, 1385, 2142, 1067, 1386, 2143, 1068, 1387, 2144, 1069, 1388, 2145, 1070, 1389, 2146, 1071, 1390, 2147, 1072, 1391, 2148, 1073, 1392, 2149, 1074, 1393, 2150, 1075, 1394, 2151, 1076, 1395, 2152, 1077, 1396, 2153, 1078, 1397, 2154, 1079, 1398, 2155, 1080, 1399, 2156, 1081, 1400, 2157, 1082, 1401, 2158, 1112, 1460, 2159, 1113, 1461, 2160, 1114, 1462, 2161, 1115, 1463, 2162, 1116, 1464, 2163, 1117, 1465, 2164, 1118, 1466, 2165, 1119, 1467, 2166, 1120, 1468, 2167, 1121, 1469, 2168, 1122, 1470, 2169, 1123, 1471, 2170, 1124, 1472, 2171, 1125, 1473, 2172, 1126, 1474, 2173, 1127, 1475, 2174, 1128, 1476, 2175, 1129, 1477, 2176, 1130, 1478, 2177, 1131, 1479, 2178, 1132, 1480, 2179, 1133, 1481, 2180, 1134, 1482, 2181, 1135, 1483, 2182, 1136, 1484, 2183, 1137, 1485, 2184, 1138, 1486, 2185, 1139, 1487, 2186, 1140, 1488, 2187, 1170, 1547, 2188, 1171, 1548, 2189, 1172, 1549, 2190, 1173, 1550, 2191, 1174, 1551, 2192, 1175, 1552, 2193, 1176, 1553, 2194, 1177, 1554, 2195, 1178, 1555, 2196, 1179, 1556, 2197, 1180, 1557, 2198, 1181, 1558, 2199, 1182, 1559, 2200, 1183, 1560, 2201, 1184, 1561, 2202, 1185, 1562, 2203, 1186, 1563, 2204, 1187, 1564, 2205, 1188, 1565, 2206, 1189, 1566, 2207, 1190, 1567, 2208, 1191, 1568, 2209, 1192, 1569, 2210, 1193, 1570, 2211, 1194, 1571, 2212, 1195, 1572, 2213, 1196, 1573, 2214, 1197, 1574, 2215, 1198, 1575, 2216, 1228, 1257, 2217, 1229, 1258, 2218, 1230, 1259, 2219, 1231, 1260, 2220, 1232, 1261, 2221, 1233, 1262, 2222, 1234, 1263, 2223, 1235, 1264, 2224, 1236, 1265, 2225, 1237, 1266, 2226, 1238, 1267, 2227, 1239, 1268, 2228, 1240, 1269, 2229, 1241, 1270, 2230, 1242, 1271, 2231, 1243, 1272, 2232, 1244, 1273, 2233, 1245, 1274, 2234, 1246, 1275, 2235, 1247, 1276, 2236, 1248, 1277, 2237, 1249, 1278, 2238, 1250, 1279, 2239, 1251, 1280, 2240, 1252, 1281, 2241, 1253, 1282, 2242, 1254, 1283, 2243, 1255, 1284, 2244, 1256, 1285, 2245, 1315, 1344, 2246, 1316, 1345, 2247, 1317, 1346, 2248, 1318, 1347, 2249, 1319, 1348, 2250, 1320, 1349, 2251, 1321, 1350, 2252, 1322, 1351, 2253, 1323, 1352, 2254, 1324, 1353, 2255, 1325, 1354, 2256, 1326, 1355, 2257, 1327, 1356, 2258, 1328, 1357, 2259, 1329, 1358, 2260, 1330, 1359, 2261, 1331, 1360, 2262, 1332, 1361, 2263, 1333, 1362, 2264, 1334, 1363, 2265, 1335, 1364, 2266, 1336, 1365, 2267, 1337, 1366, 2268, 1338, 1367, 2269, 1339, 1368, 2270, 1340, 1369, 2271, 1341, 1370, 2272, 1342, 1371, 2273, 1343, 1372, 2274, 1402, 1431, 2275, 1403, 1432, 2276, 1404, 1433, 2277, 1405, 1434, 2278, 1406, 1435, 2279, 1407, 1436, 2280, 1408, 1437, 2281, 1409, 1438, 2282, 1410, 1439, 2283, 1411, 1440, 2284, 1412, 1441, 2285, 1413, 1442, 2286, 1414, 1443, 2287, 1415, 1444, 2288, 1416, 1445, 2289, 1417, 1446, 2290, 1418, 1447, 2291, 1419, 1448, 2292, 1420, 1449, 2293, 1421, 1450, 2294, 1422, 1451, 2295, 1423, 1452, 2296, 1424, 1453, 2297, 1425, 1454, 2298, 1426, 1455, 2299, 1427, 1456, 2300, 1428, 1457, 2301, 1429, 1458, 2302, 1430, 1459, 2303, 1489, 1518, 2304, 1490, 1519, 2305, 1491, 1520, 2306, 1492, 1521, 2307, 1493, 1522, 2308, 1494, 1523, 2309, 1495, 1524, 2310, 1496, 1525, 2311, 1497, 1526, 2312, 1498, 1527, 2313, 1499, 1528, 2314, 1500, 1529, 2315, 1501, 1530, 2316, 1502, 1531, 2317, 1503, 1532, 2318, 1504, 1533, 2319, 1505, 1534, 2320, 1506, 1535, 2321, 1507, 1536, 2322, 1508, 1537, 2323, 1509, 1538, 2324, 1510, 1539, 2325, 1511, 1540, 2326, 1512, 1541, 2327, 1513, 1542, 2328, 1514, 1543, 2329, 1515, 1544, 2330, 1516, 1545, 2331, 1517, 1546, 2332, 1576, 1605, 2333, 1577, 1606, 2334, 1578, 1607, 2335, 1579, 1608, 2336, 1580, 1609, 2337, 1581, 1610, 2338, 1582, 1611, 2339, 1583, 1612, 2340, 1584, 1613, 2341, 1585, 1614, 2342, 1586, 1615, 2343, 1587, 1616, 2344, 1588, 1617, 2345, 1589, 1618, 2346, 1590, 1619, 2347, 1591, 1620, 2348, 1592, 1621, 2349, 1593, 1622, 2350, 1594, 1623, 2351, 1595, 1624, 2352, 1596, 1625, 2353, 1597, 1626, 2354, 1598, 1627, 2355, 1599, 1628, 2356, 1600, 1629, 2357, 1601, 1630, 2358, 1602, 1631, 2359, 1603, 1632, 2360, 1604, 1633, 2361, 206, 1634, 1636, 1638, 1640, 1642, 1644, 1646, 1648, 207, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 208, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 209, 1808, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 1693, 1809, 2362, 1696, 2363, 210, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 211, 1805, 212, 1805, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1801, 1806, 2364, 1804, 2365, 213, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1785, 1788, 2366, 1786, 1789, 1793, 1801, 2367, 1547, 1548, 1549, 1550, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558, 1559, 1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 1570, 1571, 1572, 1573, 1574, 1575, 2368, 1807, 2369, 1740, 1742, 2370, 214, 1740, 1701, 1741, 2371, 1650, 1652, 1704, 1825, 2372, 215, 1199, 1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 1286, 1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297, 1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308, 1309, 1310, 1311, 1312, 1313, 1314, 1701, 1824, 216, 1704, 1707, 1710, 2373, 1700, 2374, 217, 1650, 1697, 1748, 1651, 1654, 1658, 2375, 218, 1652, 1659, 1672, 1675, 2376, 219, 1660, 1666, 1672, 1663, 1669, 1675, 1678, 2377, 1679, 1680, 1683, 1826, 2378, 220, 1711, 1722, 1726, 2379, 1686, 1723, 1735, 1816, 2380, 221, 1717, 1720, 1727, 2381, 1692, 1728, 1732, 2382, 222, 1686, 1728, 1689, 1723, 2383, 223, 1732, 1735, 1690, 1738, 2384, 224, 1738, 1739, 2385, 225, 1776, 2386, 226, 1773, 1775, 1774, 1813, 2387, 2388, 1373, 1374, 1375, 1376, 1377, 1378, 1379, 1380, 1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1389, 1390, 1391, 1392, 1393, 1394, 1395, 1396, 1397, 1398, 1399, 1400, 1401, 1460, 1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468, 1469, 1470, 1471, 1472, 1473, 1474, 1475, 1476, 1477, 1478, 1479, 1480, 1481, 1482, 1483, 1484, 1485, 1486, 1487, 1488, 1822, 2389, 1750, 1753, 1755, 1768, 1823, 2390, 227, 1750, 1754, 1759, 2391, 228, 1753, 1828, 1818, 2392, 1768, 1769, 1777, 1797, 1817, 1830, 1838, 2393, 229, 1755, 1756, 1759, 2394, 1762, 1763, 2395, 230, 1757, 231, 1757, 1787, 1758, 1763, 2396, 1766, 1767, 2397, 232, 1769, 1773, 1772, 2398, 233, 1783, 234, 1777, 1780, 1783, 1691, 2399, 1784, 1786, 1813, 2400, 235, 1742, 1743, 1744, 2401, 1747, 2402, 236, 1744, 1748, 1749, 2403, 237, 1680, 1683, 1711, 1814, 2404, 585, 587, 589, 1693, 1697, 1707, 1714, 1717, 1721, 1780, 1815, 1834, 1841, 2405, 238, 1797, 1820, 1800, 2406, 239, 1819, 2407, 240, 1787, 241, 1793, 0, 1796, 2408, 242, 1789, 1810, 1792, 2409, 243, 1808, 1694, 1695, 2410, 1698, 1699, 2411, 1702, 1703, 2412, 1705, 1706, 2413, 1708, 1709, 2414, 1673, 1674, 2415, 1676, 1677, 2416, 1681, 1682, 2417, 1684, 1685, 2418, 1718, 1719, 2419, 1729, 1730, 2420, 1687, 1688, 1731, 2421, 1733, 1734, 2422, 1736, 1737, 2423, 1745, 1746, 2424, 1751, 1752, 2425, 1760, 1761, 2426, 1764, 1765, 2427, 1770, 1771, 2428, 1778, 1779, 2429, 1781, 1782, 2430, 2431, 1790, 1791, 2432, 1794, 1795, 2433, 1798, 1799, 2434, 1802, 1803, 2435, 1712, 1713, 2436, 1715, 1716, 2437, 1724, 1725, 2438, 1653, 1656, 1658, 2439, 1657, 1666, 1669, 2440, 1667, 1668, 2441, 1670, 1671, 2442, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 244, 248, 1810, 1811, 2443, 1655, 1660, 1663, 2444, 1661, 1662, 2445, 1664, 1665, 2446, 1821, 2447, 1826, 1827, 2448, 1828, 1829, 2449, 248, 251, 249, 252, 245, 246, 331, 247, 321, 339, 425, 324, 327, 426, 326, 330, 427, 334, 335, 428, 337, 338, 429
};

void dwdx_rowvals_SPARCED(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED