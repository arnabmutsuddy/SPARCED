#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 2921> dwdx_rowvals_SPARCED_ = {
    1, 146, 260, 261, 263, 262, 264, 278, 298, 319, 147, 260, 264, 267, 268, 269, 317, 1779, 148, 262, 271, 273, 261, 267, 273, 274, 261, 268, 276, 277, 280, 279, 282, 281, 284, 283, 286, 285, 288, 287, 290, 289, 292, 291, 294, 293, 296, 3, 295, 297, 300, 299, 302, 301, 304, 303, 306, 305, 308, 307, 310, 309, 312, 311, 314, 313, 316, 4, 315, 149, 321, 150, 322, 151, 323, 272, 321, 275, 322, 323, 1780, 1781, 318, 1782, 152, 317, 318, 153, 319, 320, 154, 155, 156, 325, 329, 334, 339, 373, 374, 157, 325, 334, 158, 327, 328, 329, 331, 346, 347, 159, 337, 324, 338, 339, 342, 348, 160, 349, 350, 161, 162, 163, 164, 324, 351, 352, 165, 353, 166, 354, 355, 356, 167, 357, 168, 358, 359, 169, 360, 361, 170, 326, 331, 362, 335, 342, 363, 171, 364, 365, 366, 367, 172, 381, 173, 381, 382, 383, 1794, 384, 386, 1795, 174, 384, 385, 1796, 175, 386, 397, 467, 387, 388, 1797, 389, 391, 406, 1798, 176, 389, 390, 1799, 177, 391, 442, 392, 393, 1800, 394, 400, 403, 1801, 178, 394, 395, 396, 1802, 397, 1803, 398, 399, 1804, 179, 400, 447, 449, 401, 402, 1805, 180, 403, 404, 405, 1806, 1807, 181, 406, 407, 408, 1808, 409, 411, 1809, 182, 409, 451, 453, 455, 460, 465, 410, 1810, 183, 411, 457, 412, 413, 1811, 414, 1812, 415, 416, 418, 1813, 416, 420, 424, 466, 1814, 417, 1815, 419, 420, 422, 1816, 421, 1817, 423, 424, 426, 1818, 425, 1819, 426, 427, 428, 429, 432, 462, 184, 429, 430, 431, 435, 464, 1820, 185, 432, 433, 434, 445, 463, 1821, 436, 437, 1822, 186, 437, 438, 439, 1823, 440, 1824, 187, 440, 441, 442, 447, 1825, 443, 444, 1826, 446, 449, 1827, 448, 1828, 450, 1829, 1830, 188, 451, 1787, 1790, 189, 453, 190, 455, 452, 1831, 454, 1832, 456, 1833, 191, 457, 460, 1783, 458, 459, 1834, 461, 1835, 1784, 1785, 1836, 1786, 1837, 1788, 1789, 1838, 1793, 1839, 1791, 1792, 1840, 192, 468, 474, 498, 520, 590, 592, 594, 598, 193, 480, 488, 492, 502, 510, 514, 524, 530, 596, 600, 602, 604, 606, 608, 194, 540, 544, 612, 195, 548, 196, 554, 559, 616, 197, 618, 620, 198, 622, 624, 199, 468, 472, 496, 518, 534, 572, 574, 576, 578, 535, 1841, 200, 470, 494, 516, 536, 574, 580, 582, 584, 537, 1842, 201, 478, 492, 500, 522, 576, 582, 586, 202, 486, 508, 514, 528, 538, 578, 584, 586, 588, 539, 1843, 203, 484, 506, 526, 204, 540, 542, 610, 205, 548, 206, 554, 561, 614, 207, 562, 208, 564, 573, 590, 1844, 575, 592, 1845, 577, 594, 596, 1846, 579, 598, 600, 1847, 581, 1848, 583, 602, 1849, 585, 604, 1850, 587, 606, 1851, 589, 608, 1852, 611, 612, 1853, 615, 616, 1854, 563, 618, 1855, 565, 622, 1856, 469, 470, 472, 477, 478, 483, 484, 486, 491, 1857, 483, 493, 494, 496, 500, 505, 506, 508, 513, 1858, 491, 513, 515, 516, 518, 522, 526, 528, 533, 1859, 541, 542, 547, 1860, 549, 550, 1861, 555, 556, 561, 1862, 471, 593, 626, 684, 1863, 485, 627, 685, 1864, 473, 474, 591, 628, 686, 1865, 479, 480, 595, 630, 688, 1866, 487, 488, 599, 632, 690, 1867, 495, 603, 634, 692, 1868, 497, 498, 597, 636, 694, 1869, 501, 502, 1870, 507, 635, 693, 1871, 509, 510, 637, 695, 1872, 517, 605, 639, 697, 1873, 527, 640, 698, 1874, 519, 520, 601, 641, 699, 1875, 523, 524, 607, 642, 700, 1876, 529, 530, 609, 643, 701, 1877, 543, 544, 613, 645, 703, 1878, 553, 648, 706, 1879, 559, 560, 617, 650, 708, 1880, 619, 620, 651, 709, 1881, 623, 624, 652, 710, 1882, 475, 476, 629, 687, 1883, 481, 482, 499, 631, 689, 1884, 489, 490, 521, 633, 691, 1885, 503, 504, 1886, 511, 512, 525, 638, 696, 1887, 531, 532, 644, 702, 1888, 545, 546, 646, 704, 1889, 551, 552, 647, 705, 1890, 557, 558, 649, 707, 1891, 621, 653, 711, 1892, 625, 654, 712, 1893, 566, 567, 1894, 568, 569, 1895, 570, 571, 1896, 655, 742, 916, 974, 1032, 1090, 1897, 656, 743, 917, 975, 1033, 1091, 1898, 657, 744, 918, 976, 1034, 1092, 1899, 658, 745, 919, 977, 1035, 1093, 1900, 659, 746, 920, 978, 1036, 1094, 1901, 660, 747, 921, 979, 1037, 1095, 1902, 661, 748, 922, 980, 1038, 1096, 1903, 662, 749, 923, 981, 1039, 1097, 1904, 663, 750, 924, 982, 1040, 1098, 1905, 664, 751, 925, 983, 1041, 1099, 1906, 665, 752, 926, 984, 1042, 1100, 1907, 666, 753, 927, 985, 1043, 1101, 1908, 667, 754, 928, 986, 1044, 1102, 1909, 668, 755, 929, 987, 1045, 1103, 1910, 669, 756, 930, 988, 1046, 1104, 1911, 670, 757, 931, 989, 1047, 1105, 1912, 671, 758, 932, 990, 1048, 1106, 1913, 672, 759, 933, 991, 1049, 1107, 1914, 673, 760, 934, 992, 1050, 1108, 1915, 674, 761, 935, 993, 1051, 1109, 1916, 675, 762, 936, 994, 1052, 1110, 1917, 676, 763, 937, 995, 1053, 1111, 1918, 677, 764, 938, 996, 1054, 1112, 1919, 678, 765, 939, 997, 1055, 1113, 1920, 679, 766, 940, 998, 1056, 1114, 1921, 680, 767, 1583, 1922, 681, 768, 1585, 1923, 682, 769, 1587, 1924, 683, 770, 1589, 1925, 941, 999, 1057, 1115, 1584, 1926, 942, 1000, 1058, 1116, 1586, 1927, 943, 1001, 1059, 1117, 1588, 1928, 944, 1002, 1060, 1118, 1590, 1929, 713, 1930, 714, 1931, 715, 1932, 716, 1933, 717, 1934, 718, 1935, 719, 1936, 720, 1937, 721, 1938, 722, 1939, 723, 1940, 724, 1941, 725, 1942, 726, 1943, 727, 1944, 728, 1945, 729, 1946, 730, 1947, 731, 1948, 732, 1949, 733, 1950, 734, 1951, 735, 1952, 736, 1953, 737, 1954, 738, 1955, 739, 1956, 740, 1957, 741, 1958, 800, 829, 1959, 801, 830, 1960, 802, 831, 1961, 803, 832, 1962, 804, 833, 1963, 805, 834, 1964, 806, 835, 1965, 807, 836, 1966, 808, 837, 1967, 809, 838, 1968, 810, 839, 1969, 811, 840, 1970, 812, 841, 1971, 813, 842, 1972, 814, 843, 1973, 815, 844, 1974, 816, 845, 1975, 817, 846, 1976, 818, 847, 1977, 819, 848, 1978, 820, 849, 1979, 821, 850, 1980, 822, 851, 1981, 823, 852, 1982, 824, 853, 1983, 825, 854, 1984, 826, 855, 1985, 827, 856, 1986, 828, 857, 1987, 771, 858, 1988, 772, 859, 1989, 773, 860, 1990, 774, 861, 1991, 775, 862, 1992, 776, 863, 1993, 777, 864, 1994, 778, 865, 1995, 779, 866, 1996, 780, 867, 1997, 781, 868, 1998, 782, 869, 1999, 783, 870, 2000, 784, 871, 2001, 785, 872, 2002, 786, 873, 2003, 787, 874, 2004, 788, 875, 2005, 789, 876, 2006, 790, 877, 2007, 791, 878, 2008, 792, 879, 2009, 793, 880, 2010, 794, 881, 2011, 795, 882, 2012, 796, 1591, 2013, 797, 1593, 2014, 798, 1595, 2015, 799, 1597, 2016, 883, 1592, 2017, 884, 1594, 2018, 885, 1596, 2019, 886, 1598, 2020, 945, 1235, 2021, 946, 1236, 2022, 947, 1237, 2023, 948, 1238, 2024, 949, 1239, 2025, 950, 1240, 2026, 951, 1241, 2027, 952, 1242, 2028, 953, 1243, 2029, 954, 1244, 2030, 955, 1245, 2031, 956, 1246, 2032, 957, 1247, 2033, 958, 1248, 2034, 959, 1249, 2035, 960, 1250, 2036, 961, 1251, 2037, 962, 1252, 2038, 963, 1253, 2039, 964, 1254, 2040, 965, 1255, 2041, 966, 1256, 2042, 967, 1257, 2043, 968, 1258, 2044, 969, 1259, 2045, 970, 1260, 2046, 971, 1261, 2047, 972, 1262, 2048, 973, 1263, 2049, 887, 1148, 2050, 888, 1149, 2051, 889, 1150, 2052, 890, 1151, 2053, 891, 1152, 2054, 892, 1153, 2055, 893, 1154, 2056, 894, 1155, 2057, 895, 1156, 2058, 896, 1157, 2059, 897, 1158, 2060, 898, 1159, 2061, 899, 1160, 2062, 900, 1161, 2063, 901, 1162, 2064, 902, 1163, 2065, 903, 1164, 2066, 904, 1165, 2067, 905, 1166, 2068, 906, 1167, 2069, 907, 1168, 2070, 908, 1169, 2071, 909, 1170, 2072, 910, 1171, 2073, 911, 1172, 2074, 912, 1173, 2075, 913, 1174, 2076, 914, 1175, 2077, 915, 1176, 2078, 1003, 1322, 2079, 1004, 1323, 2080, 1005, 1324, 2081, 1006, 1325, 2082, 1007, 1326, 2083, 1008, 1327, 2084, 1009, 1328, 2085, 1010, 1329, 2086, 1011, 1330, 2087, 1012, 1331, 2088, 1013, 1332, 2089, 1014, 1333, 2090, 1015, 1334, 2091, 1016, 1335, 2092, 1017, 1336, 2093, 1018, 1337, 2094, 1019, 1338, 2095, 1020, 1339, 2096, 1021, 1340, 2097, 1022, 1341, 2098, 1023, 1342, 2099, 1024, 1343, 2100, 1025, 1344, 2101, 1026, 1345, 2102, 1027, 1346, 2103, 1028, 1347, 2104, 1029, 1348, 2105, 1030, 1349, 2106, 1031, 1350, 2107, 1061, 1409, 2108, 1062, 1410, 2109, 1063, 1411, 2110, 1064, 1412, 2111, 1065, 1413, 2112, 1066, 1414, 2113, 1067, 1415, 2114, 1068, 1416, 2115, 1069, 1417, 2116, 1070, 1418, 2117, 1071, 1419, 2118, 1072, 1420, 2119, 1073, 1421, 2120, 1074, 1422, 2121, 1075, 1423, 2122, 1076, 1424, 2123, 1077, 1425, 2124, 1078, 1426, 2125, 1079, 1427, 2126, 1080, 1428, 2127, 1081, 1429, 2128, 1082, 1430, 2129, 1083, 1431, 2130, 1084, 1432, 2131, 1085, 1433, 2132, 1086, 1434, 2133, 1087, 1435, 2134, 1088, 1436, 2135, 1089, 1437, 2136, 1119, 1496, 2137, 1120, 1497, 2138, 1121, 1498, 2139, 1122, 1499, 2140, 1123, 1500, 2141, 1124, 1501, 2142, 1125, 1502, 2143, 1126, 1503, 2144, 1127, 1504, 2145, 1128, 1505, 2146, 1129, 1506, 2147, 1130, 1507, 2148, 1131, 1508, 2149, 1132, 1509, 2150, 1133, 1510, 2151, 1134, 1511, 2152, 1135, 1512, 2153, 1136, 1513, 2154, 1137, 1514, 2155, 1138, 1515, 2156, 1139, 1516, 2157, 1140, 1517, 2158, 1141, 1518, 2159, 1142, 1519, 2160, 1143, 1520, 2161, 1144, 1521, 2162, 1145, 1522, 2163, 1146, 1523, 2164, 1147, 1524, 2165, 1177, 1206, 2166, 1178, 1207, 2167, 1179, 1208, 2168, 1180, 1209, 2169, 1181, 1210, 2170, 1182, 1211, 2171, 1183, 1212, 2172, 1184, 1213, 2173, 1185, 1214, 2174, 1186, 1215, 2175, 1187, 1216, 2176, 1188, 1217, 2177, 1189, 1218, 2178, 1190, 1219, 2179, 1191, 1220, 2180, 1192, 1221, 2181, 1193, 1222, 2182, 1194, 1223, 2183, 1195, 1224, 2184, 1196, 1225, 2185, 1197, 1226, 2186, 1198, 1227, 2187, 1199, 1228, 2188, 1200, 1229, 2189, 1201, 1230, 2190, 1202, 1231, 2191, 1203, 1232, 2192, 1204, 1233, 2193, 1205, 1234, 2194, 1264, 1293, 2195, 1265, 1294, 2196, 1266, 1295, 2197, 1267, 1296, 2198, 1268, 1297, 2199, 1269, 1298, 2200, 1270, 1299, 2201, 1271, 1300, 2202, 1272, 1301, 2203, 1273, 1302, 2204, 1274, 1303, 2205, 1275, 1304, 2206, 1276, 1305, 2207, 1277, 1306, 2208, 1278, 1307, 2209, 1279, 1308, 2210, 1280, 1309, 2211, 1281, 1310, 2212, 1282, 1311, 2213, 1283, 1312, 2214, 1284, 1313, 2215, 1285, 1314, 2216, 1286, 1315, 2217, 1287, 1316, 2218, 1288, 1317, 2219, 1289, 1318, 2220, 1290, 1319, 2221, 1291, 1320, 2222, 1292, 1321, 2223, 1351, 1380, 2224, 1352, 1381, 2225, 1353, 1382, 2226, 1354, 1383, 2227, 1355, 1384, 2228, 1356, 1385, 2229, 1357, 1386, 2230, 1358, 1387, 2231, 1359, 1388, 2232, 1360, 1389, 2233, 1361, 1390, 2234, 1362, 1391, 2235, 1363, 1392, 2236, 1364, 1393, 2237, 1365, 1394, 2238, 1366, 1395, 2239, 1367, 1396, 2240, 1368, 1397, 2241, 1369, 1398, 2242, 1370, 1399, 2243, 1371, 1400, 2244, 1372, 1401, 2245, 1373, 1402, 2246, 1374, 1403, 2247, 1375, 1404, 2248, 1376, 1405, 2249, 1377, 1406, 2250, 1378, 1407, 2251, 1379, 1408, 2252, 1438, 1467, 2253, 1439, 1468, 2254, 1440, 1469, 2255, 1441, 1470, 2256, 1442, 1471, 2257, 1443, 1472, 2258, 1444, 1473, 2259, 1445, 1474, 2260, 1446, 1475, 2261, 1447, 1476, 2262, 1448, 1477, 2263, 1449, 1478, 2264, 1450, 1479, 2265, 1451, 1480, 2266, 1452, 1481, 2267, 1453, 1482, 2268, 1454, 1483, 2269, 1455, 1484, 2270, 1456, 1485, 2271, 1457, 1486, 2272, 1458, 1487, 2273, 1459, 1488, 2274, 1460, 1489, 2275, 1461, 1490, 2276, 1462, 1491, 2277, 1463, 1492, 2278, 1464, 1493, 2279, 1465, 1494, 2280, 1466, 1495, 2281, 1525, 1554, 2282, 1526, 1555, 2283, 1527, 1556, 2284, 1528, 1557, 2285, 1529, 1558, 2286, 1530, 1559, 2287, 1531, 1560, 2288, 1532, 1561, 2289, 1533, 1562, 2290, 1534, 1563, 2291, 1535, 1564, 2292, 1536, 1565, 2293, 1537, 1566, 2294, 1538, 1567, 2295, 1539, 1568, 2296, 1540, 1569, 2297, 1541, 1570, 2298, 1542, 1571, 2299, 1543, 1572, 2300, 1544, 1573, 2301, 1545, 1574, 2302, 1546, 1575, 2303, 1547, 1576, 2304, 1548, 1577, 2305, 1549, 1578, 2306, 1550, 1579, 2307, 1551, 1580, 2308, 1552, 1581, 2309, 1553, 1582, 2310, 209, 1583, 1585, 1587, 1589, 1591, 1593, 1595, 1597, 210, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 211, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 212, 1757, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940, 941, 942, 943, 944, 1642, 1758, 2311, 1645, 2312, 213, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1002, 214, 1754, 215, 1754, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1750, 1755, 2313, 1753, 2314, 216, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1734, 1737, 2315, 1735, 1738, 1742, 1750, 2316, 1496, 1497, 1498, 1499, 1500, 1501, 1502, 1503, 1504, 1505, 1506, 1507, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1517, 1518, 1519, 1520, 1521, 1522, 1523, 1524, 2317, 1756, 2318, 1689, 1691, 2319, 217, 1689, 1650, 1690, 2320, 1599, 1601, 1653, 1774, 2321, 218, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1650, 1773, 219, 1653, 1656, 1659, 2322, 1649, 2323, 220, 1599, 1646, 1697, 1600, 1603, 1607, 2324, 221, 1601, 1608, 1621, 1624, 2325, 222, 1609, 1615, 1621, 1612, 1618, 1624, 1627, 2326, 1628, 1629, 1632, 1775, 2327, 223, 1660, 1671, 1675, 2328, 1635, 1672, 1684, 1765, 2329, 224, 1666, 1669, 1676, 2330, 1641, 1677, 1681, 2331, 225, 1635, 1677, 1638, 1672, 2332, 226, 1681, 1684, 1639, 1687, 2333, 227, 1687, 1688, 2334, 228, 1725, 2335, 229, 1722, 1724, 1723, 1762, 2336, 2337, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 1346, 1347, 1348, 1349, 1350, 1409, 1410, 1411, 1412, 1413, 1414, 1415, 1416, 1417, 1418, 1419, 1420, 1421, 1422, 1423, 1424, 1425, 1426, 1427, 1428, 1429, 1430, 1431, 1432, 1433, 1434, 1435, 1436, 1437, 1771, 2338, 1699, 1702, 1704, 1717, 1772, 2339, 230, 1699, 1703, 1708, 2340, 231, 1702, 1777, 1767, 2341, 1717, 1718, 1726, 1746, 1766, 1779, 1787, 2342, 232, 1704, 1705, 1708, 2343, 1711, 1712, 2344, 233, 1706, 234, 1706, 1736, 1707, 1712, 2345, 1715, 1716, 2346, 235, 1718, 1722, 1721, 2347, 236, 1732, 237, 1726, 1729, 1732, 1640, 2348, 1733, 1735, 1762, 2349, 238, 1691, 1692, 1693, 2350, 1696, 2351, 239, 1693, 1697, 1698, 2352, 240, 1629, 1632, 1660, 1763, 2353, 534, 536, 538, 1642, 1646, 1656, 1663, 1666, 1670, 1729, 1764, 1783, 1790, 2354, 241, 1746, 1769, 1749, 2355, 242, 1768, 2356, 243, 1736, 244, 1742, 0, 1745, 2357, 245, 1738, 1759, 1741, 2358, 246, 1757, 1643, 1644, 2359, 1647, 1648, 2360, 1651, 1652, 2361, 1654, 1655, 2362, 1657, 1658, 2363, 1622, 1623, 2364, 1625, 1626, 2365, 1630, 1631, 2366, 1633, 1634, 2367, 1667, 1668, 2368, 1678, 1679, 2369, 1636, 1637, 1680, 2370, 1682, 1683, 2371, 1685, 1686, 2372, 1694, 1695, 2373, 1700, 1701, 2374, 1709, 1710, 2375, 1713, 1714, 2376, 1719, 1720, 2377, 1727, 1728, 2378, 1730, 1731, 2379, 2380, 1739, 1740, 2381, 1743, 1744, 2382, 1747, 1748, 2383, 1751, 1752, 2384, 1661, 1662, 2385, 1664, 1665, 2386, 1673, 1674, 2387, 1602, 1605, 1607, 2388, 1606, 1615, 1618, 2389, 1616, 1617, 2390, 1619, 1620, 2391, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 247, 254, 1759, 1760, 2392, 1604, 1609, 1612, 2393, 1610, 1611, 2394, 1613, 1614, 2395, 1770, 2396, 1775, 1776, 2397, 1777, 1778, 2398, 254, 257, 255, 258, 248, 249, 337, 250, 327, 345, 368, 330, 333, 369, 332, 336, 370, 340, 341, 371, 343, 344, 372, 251, 252, 253, 375, 376, 377, 378, 379, 380
};

void dwdx_rowvals_SPARCED(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED