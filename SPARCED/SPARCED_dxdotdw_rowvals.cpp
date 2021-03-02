#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SPARCED {

static constexpr std::array<sunindextype, 5054> dxdotdw_rowvals_SPARCED_ = {
    0, 0, 1, 3, 4, 37, 38, 39, 42, 42, 42, 43, 43, 43, 46, 46, 48, 52, 51, 53, 54, 57, 59, 63, 65, 68, 71, 74, 769, 770, 771, 771, 79, 80, 80, 83, 85, 85, 88, 90, 90, 93, 97, 99, 102, 105, 105, 105, 107, 120, 123, 127, 130, 137, 138, 139, 143, 151, 152, 158, 160, 162, 163, 165, 166, 153, 167, 167, 154, 639, 640, 641, 644, 644, 645, 645, 645, 645, 646, 646, 646, 646, 649, 655, 655, 658, 658, 658, 659, 662, 664, 666, 666, 669, 672, 672, 672, 672, 675, 677, 679, 681, 683, 688, 690, 690, 693, 696, 697, 700, 702, 703, 706, 706, 706, 706, 709, 711, 711, 714, 716, 718, 719, 719, 721, 723, 34, 35, 168, 168, 155, 155, 757, 638, 638, 156, 156, 169, 28, 27, 29, 170, 157, 1, 3, 4, 27, 28, 29, 34, 35, 37, 38, 39, 42, 43, 46, 48, 51, 52, 53, 54, 57, 59, 63, 65, 68, 71, 74, 79, 80, 83, 85, 88, 90, 93, 97, 99, 102, 105, 107, 120, 123, 127, 130, 137, 138, 139, 143, 151, 152, 153, 154, 155, 156, 157, 158, 160, 162, 163, 165, 166, 167, 168, 169, 170, 638, 639, 640, 641, 644, 645, 646, 649, 655, 658, 659, 662, 664, 666, 669, 672, 675, 677, 679, 681, 683, 688, 690, 693, 696, 697, 700, 702, 703, 706, 709, 711, 714, 716, 718, 719, 721, 723, 757, 769, 770, 771, 757, 767, 768, 757, 767, 768, 767, 767, 768, 1, 1, 1, 2, 1, 2, 1, 2, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 3, 33, 2, 35, 36, 2, 35, 36, 30, 31, 31, 30, 31, 39, 42, 72, 39, 42, 72, 43, 44, 771, 43, 44, 771, 39, 44, 773, 39, 44, 773, 44, 72, 774, 44, 72, 774, 40, 44, 773, 40, 42, 73, 40, 42, 73, 44, 73, 774, 46, 47, 770, 46, 47, 770, 40, 47, 775, 40, 47, 775, 41, 47, 775, 47, 73, 776, 47, 73, 776, 47, 772, 776, 41, 42, 772, 39, 40, 40, 41, 44, 45, 57, 44, 45, 57, 47, 49, 57, 47, 50, 57, 55, 56, 57, 55, 56, 57, 57, 64, 70, 57, 64, 70, 57, 61, 57, 61, 44, 45, 47, 49, 50, 55, 56, 58, 60, 61, 62, 64, 66, 67, 69, 70, 72, 73, 75, 76, 77, 78, 772, 773, 774, 775, 776, 40, 41, 79, 80, 81, 79, 80, 81, 81, 82, 82, 83, 84, 82, 83, 84, 82, 85, 86, 82, 85, 86, 82, 86, 87, 87, 88, 89, 87, 88, 89, 87, 90, 91, 87, 90, 91, 87, 91, 92, 92, 93, 94, 92, 93, 94, 92, 94, 95, 85, 95, 96, 85, 95, 96, 87, 95, 96, 92, 97, 98, 92, 97, 98, 97, 98, 136, 92, 99, 100, 92, 99, 100, 92, 100, 101, 87, 102, 103, 87, 102, 103, 87, 103, 104, 104, 105, 106, 104, 105, 106, 104, 107, 108, 104, 107, 108, 104, 108, 109, 109, 110, 109, 110, 110, 111, 112, 110, 111, 112, 110, 113, 110, 113, 111, 113, 114, 111, 113, 114, 113, 115, 113, 115, 111, 115, 116, 111, 115, 116, 115, 117, 118, 115, 117, 118, 115, 118, 119, 119, 120, 121, 119, 120, 121, 119, 121, 122, 119, 123, 124, 119, 123, 124, 119, 124, 125, 122, 126, 122, 126, 126, 127, 128, 126, 127, 128, 126, 128, 129, 129, 130, 131, 129, 130, 131, 90, 131, 132, 90, 131, 132, 92, 131, 132, 125, 133, 125, 133, 97, 131, 134, 97, 131, 134, 97, 133, 135, 97, 133, 135, 105, 137, 140, 105, 137, 140, 105, 138, 141, 105, 138, 141, 105, 139, 142, 105, 139, 142, 107, 143, 144, 107, 143, 144, 109, 143, 144, 105, 143, 145, 105, 143, 145, 117, 119, 123, 125, 120, 122, 105, 111, 105, 111, 85, 87, 151, 158, 184, 151, 158, 184, 160, 184, 190, 160, 184, 190, 158, 184, 192, 158, 184, 192, 151, 192, 210, 151, 192, 210, 184, 210, 184, 210, 162, 184, 193, 162, 184, 193, 152, 193, 211, 152, 193, 211, 184, 185, 211, 184, 185, 211, 165, 184, 191, 165, 184, 191, 163, 184, 194, 163, 184, 194, 152, 194, 212, 152, 194, 212, 184, 186, 212, 184, 186, 212, 152, 162, 185, 152, 162, 185, 160, 185, 195, 160, 185, 195, 158, 185, 196, 158, 185, 196, 151, 196, 211, 151, 196, 211, 162, 185, 197, 162, 185, 197, 152, 197, 213, 152, 197, 213, 185, 213, 185, 213, 165, 185, 198, 165, 185, 198, 163, 185, 199, 163, 185, 199, 152, 199, 214, 152, 199, 214, 185, 186, 214, 185, 186, 214, 152, 163, 186, 152, 163, 186, 160, 186, 200, 160, 186, 200, 158, 186, 202, 158, 186, 202, 151, 202, 212, 151, 202, 212, 162, 186, 203, 162, 186, 203, 152, 203, 214, 152, 203, 214, 165, 186, 201, 165, 186, 201, 163, 186, 204, 163, 186, 204, 152, 204, 215, 152, 204, 215, 186, 215, 186, 215, 158, 221, 713, 158, 159, 160, 222, 713, 160, 161, 163, 223, 713, 163, 164, 153, 166, 187, 153, 166, 187, 166, 187, 205, 166, 187, 205, 153, 205, 216, 153, 205, 216, 187, 216, 187, 216, 154, 167, 188, 154, 167, 188, 188, 217, 188, 217, 154, 206, 217, 167, 188, 206, 155, 168, 189, 155, 168, 189, 189, 218, 189, 218, 155, 207, 218, 155, 207, 218, 168, 189, 207, 168, 189, 207, 169, 182, 169, 182, 170, 183, 170, 183, 158, 221, 713, 159, 221, 713, 160, 222, 713, 161, 222, 713, 163, 223, 713, 164, 223, 713, 158, 171, 158, 171, 158, 160, 172, 158, 160, 172, 158, 162, 173, 158, 162, 173, 158, 163, 174, 158, 163, 174, 160, 175, 160, 175, 160, 162, 176, 160, 162, 176, 160, 163, 177, 160, 163, 177, 162, 163, 178, 162, 163, 178, 163, 179, 163, 179, 151, 171, 192, 151, 171, 192, 151, 172, 190, 151, 172, 190, 151, 173, 193, 151, 173, 193, 152, 173, 196, 152, 173, 196, 151, 174, 194, 151, 174, 194, 152, 174, 202, 152, 174, 202, 152, 176, 195, 152, 176, 195, 152, 177, 200, 152, 177, 200, 152, 178, 203, 152, 178, 203, 152, 179, 204, 152, 179, 204, 166, 180, 166, 180, 153, 180, 205, 153, 180, 205, 168, 181, 168, 181, 155, 181, 207, 155, 181, 207, 156, 182, 208, 156, 182, 208, 156, 208, 219, 156, 208, 219, 157, 183, 209, 157, 183, 209, 157, 209, 220, 157, 209, 220, 190, 224, 191, 225, 192, 226, 210, 227, 193, 228, 211, 229, 194, 230, 212, 231, 195, 232, 198, 233, 196, 234, 199, 235, 214, 236, 200, 237, 201, 238, 202, 239, 203, 240, 204, 241, 215, 242, 205, 243, 216, 244, 217, 245, 206, 246, 218, 247, 207, 248, 208, 249, 209, 250, 219, 251, 220, 252, 190, 224, 191, 225, 192, 226, 210, 227, 193, 228, 211, 229, 194, 230, 212, 231, 195, 232, 198, 233, 196, 234, 199, 235, 214, 236, 200, 237, 201, 238, 202, 239, 203, 240, 204, 241, 215, 242, 205, 243, 216, 244, 217, 245, 206, 246, 218, 247, 207, 248, 208, 249, 209, 250, 219, 251, 220, 252, 190, 257, 639, 191, 258, 639, 192, 259, 639, 210, 260, 639, 193, 261, 639, 211, 262, 639, 194, 263, 639, 212, 264, 639, 195, 265, 639, 198, 266, 639, 196, 267, 639, 199, 268, 639, 214, 269, 639, 200, 270, 639, 201, 271, 639, 202, 272, 639, 203, 273, 639, 204, 274, 639, 215, 275, 639, 205, 276, 639, 216, 277, 639, 217, 278, 639, 206, 279, 639, 218, 280, 639, 207, 281, 639, 208, 282, 639, 209, 283, 639, 219, 284, 639, 220, 285, 639, 190, 257, 639, 191, 258, 639, 192, 259, 639, 210, 260, 639, 193, 261, 639, 211, 262, 639, 194, 263, 639, 212, 264, 639, 195, 265, 639, 198, 266, 639, 196, 267, 639, 199, 268, 639, 214, 269, 639, 200, 270, 639, 201, 271, 639, 202, 272, 639, 203, 273, 639, 204, 274, 639, 215, 275, 639, 205, 276, 639, 216, 277, 639, 217, 278, 639, 206, 279, 639, 218, 280, 639, 207, 281, 639, 208, 282, 639, 209, 283, 639, 219, 284, 639, 220, 285, 639, 224, 315, 225, 316, 226, 317, 227, 318, 228, 319, 229, 320, 230, 321, 231, 322, 232, 323, 233, 324, 234, 325, 235, 326, 236, 327, 237, 328, 238, 329, 239, 330, 240, 331, 241, 332, 242, 333, 243, 334, 244, 335, 245, 336, 246, 337, 247, 338, 248, 339, 249, 340, 250, 341, 251, 342, 252, 343, 286, 315, 287, 316, 288, 317, 289, 318, 290, 319, 291, 320, 292, 321, 293, 322, 294, 323, 295, 324, 296, 325, 297, 326, 298, 327, 299, 328, 300, 329, 301, 330, 302, 331, 303, 332, 304, 333, 305, 334, 306, 335, 307, 336, 308, 337, 309, 338, 310, 339, 311, 340, 312, 341, 313, 342, 314, 343, 286, 315, 287, 316, 288, 317, 289, 318, 290, 319, 291, 320, 292, 321, 293, 322, 294, 323, 295, 324, 296, 325, 297, 326, 298, 327, 299, 328, 300, 329, 301, 330, 302, 331, 303, 332, 304, 333, 305, 334, 306, 335, 307, 336, 308, 337, 309, 338, 310, 339, 311, 340, 312, 341, 313, 342, 314, 343, 190, 286, 191, 287, 192, 288, 210, 289, 193, 290, 211, 291, 194, 292, 212, 293, 195, 294, 198, 295, 196, 296, 199, 297, 214, 298, 200, 299, 201, 300, 202, 301, 203, 302, 204, 303, 215, 304, 205, 305, 216, 306, 217, 307, 206, 308, 218, 309, 207, 310, 208, 311, 209, 312, 219, 313, 220, 314, 315, 377, 642, 316, 378, 642, 317, 379, 642, 318, 380, 642, 319, 381, 642, 320, 382, 642, 321, 383, 642, 322, 384, 642, 323, 385, 642, 324, 386, 642, 325, 387, 642, 326, 388, 642, 327, 389, 642, 328, 390, 642, 329, 391, 642, 330, 392, 642, 331, 393, 642, 332, 394, 642, 333, 395, 642, 334, 396, 642, 335, 397, 642, 336, 398, 642, 337, 399, 642, 338, 400, 642, 339, 401, 642, 344, 402, 642, 345, 403, 642, 346, 404, 642, 347, 405, 642, 315, 377, 642, 316, 378, 642, 317, 379, 642, 318, 380, 642, 319, 381, 642, 320, 382, 642, 321, 383, 642, 322, 384, 642, 323, 385, 642, 324, 386, 642, 325, 387, 642, 326, 388, 642, 327, 389, 642, 328, 390, 642, 329, 391, 642, 330, 392, 642, 331, 393, 642, 332, 394, 642, 333, 395, 642, 334, 396, 642, 335, 397, 642, 336, 398, 642, 337, 399, 642, 338, 400, 642, 339, 401, 642, 344, 402, 642, 345, 403, 642, 346, 404, 642, 347, 405, 642, 224, 348, 642, 225, 349, 642, 226, 350, 642, 227, 351, 642, 228, 352, 642, 229, 353, 642, 230, 354, 642, 231, 355, 642, 232, 356, 642, 233, 357, 642, 234, 358, 642, 235, 359, 642, 236, 360, 642, 237, 361, 642, 238, 362, 642, 239, 363, 642, 240, 364, 642, 241, 365, 642, 242, 366, 642, 243, 367, 642, 244, 368, 642, 245, 369, 642, 246, 370, 642, 247, 371, 642, 248, 372, 642, 253, 373, 642, 254, 374, 642, 255, 375, 642, 256, 376, 642, 224, 348, 642, 225, 349, 642, 226, 350, 642, 227, 351, 642, 228, 352, 642, 229, 353, 642, 230, 354, 642, 231, 355, 642, 232, 356, 642, 233, 357, 642, 234, 358, 642, 235, 359, 642, 236, 360, 642, 237, 361, 642, 238, 362, 642, 239, 363, 642, 240, 364, 642, 241, 365, 642, 242, 366, 642, 243, 367, 642, 244, 368, 642, 245, 369, 642, 246, 370, 642, 247, 371, 642, 248, 372, 642, 253, 373, 642, 254, 374, 642, 255, 375, 642, 256, 376, 642, 224, 406, 644, 225, 407, 644, 226, 408, 644, 227, 409, 644, 228, 410, 644, 229, 411, 644, 230, 412, 644, 231, 413, 644, 232, 414, 644, 233, 415, 644, 234, 416, 644, 235, 417, 644, 236, 418, 644, 237, 419, 644, 238, 420, 644, 239, 421, 644, 240, 422, 644, 241, 423, 644, 242, 424, 644, 243, 425, 644, 244, 426, 644, 245, 427, 644, 246, 428, 644, 247, 429, 644, 248, 430, 644, 253, 431, 644, 254, 432, 644, 255, 433, 644, 256, 434, 644, 224, 406, 644, 225, 407, 644, 226, 408, 644, 227, 409, 644, 228, 410, 644, 229, 411, 644, 230, 412, 644, 231, 413, 644, 232, 414, 644, 233, 415, 644, 234, 416, 644, 235, 417, 644, 236, 418, 644, 237, 419, 644, 238, 420, 644, 239, 421, 644, 240, 422, 644, 241, 423, 644, 242, 424, 644, 243, 425, 644, 244, 426, 644, 245, 427, 644, 246, 428, 644, 247, 429, 644, 248, 430, 644, 253, 431, 644, 254, 432, 644, 255, 433, 644, 256, 434, 644, 224, 435, 647, 225, 436, 647, 226, 437, 647, 227, 438, 647, 228, 439, 647, 229, 440, 647, 230, 441, 647, 231, 442, 647, 232, 443, 647, 233, 444, 647, 234, 445, 647, 235, 446, 647, 236, 447, 647, 237, 448, 647, 238, 449, 647, 239, 450, 647, 240, 451, 647, 241, 452, 647, 242, 453, 647, 243, 454, 647, 244, 455, 647, 245, 456, 647, 246, 457, 647, 247, 458, 647, 248, 459, 647, 253, 460, 647, 254, 461, 647, 255, 462, 647, 256, 463, 647, 224, 435, 647, 225, 436, 647, 226, 437, 647, 227, 438, 647, 228, 439, 647, 229, 440, 647, 230, 441, 647, 231, 442, 647, 232, 443, 647, 233, 444, 647, 234, 445, 647, 235, 446, 647, 236, 447, 647, 237, 448, 647, 238, 449, 647, 239, 450, 647, 240, 451, 647, 241, 452, 647, 242, 453, 647, 243, 454, 647, 244, 455, 647, 245, 456, 647, 246, 457, 647, 247, 458, 647, 248, 459, 647, 253, 460, 647, 254, 461, 647, 255, 462, 647, 256, 463, 647, 224, 464, 649, 225, 465, 649, 226, 466, 649, 227, 467, 649, 228, 468, 649, 229, 469, 649, 230, 470, 649, 231, 471, 649, 232, 472, 649, 233, 473, 649, 234, 474, 649, 235, 475, 649, 236, 476, 649, 237, 477, 649, 238, 478, 649, 239, 479, 649, 240, 480, 649, 241, 481, 649, 242, 482, 649, 243, 483, 649, 244, 484, 649, 245, 485, 649, 246, 486, 649, 247, 487, 649, 248, 488, 649, 253, 489, 649, 254, 490, 649, 255, 491, 649, 256, 492, 649, 224, 464, 649, 225, 465, 649, 226, 466, 649, 227, 467, 649, 228, 468, 649, 229, 469, 649, 230, 470, 649, 231, 471, 649, 232, 472, 649, 233, 473, 649, 234, 474, 649, 235, 475, 649, 236, 476, 649, 237, 477, 649, 238, 478, 649, 239, 479, 649, 240, 480, 649, 241, 481, 649, 242, 482, 649, 243, 483, 649, 244, 484, 649, 245, 485, 649, 246, 486, 649, 247, 487, 649, 248, 488, 649, 253, 489, 649, 254, 490, 649, 255, 491, 649, 256, 492, 649, 377, 493, 658, 378, 494, 658, 379, 495, 658, 380, 496, 658, 381, 497, 658, 382, 498, 658, 383, 499, 658, 384, 500, 658, 385, 501, 658, 386, 502, 658, 387, 503, 658, 388, 504, 658, 389, 505, 658, 390, 506, 658, 391, 507, 658, 392, 508, 658, 393, 509, 658, 394, 510, 658, 395, 511, 658, 396, 512, 658, 397, 513, 658, 398, 514, 658, 399, 515, 658, 400, 516, 658, 401, 517, 658, 402, 518, 658, 403, 519, 658, 404, 520, 658, 405, 521, 658, 377, 493, 658, 378, 494, 658, 379, 495, 658, 380, 496, 658, 381, 497, 658, 382, 498, 658, 383, 499, 658, 384, 500, 658, 385, 501, 658, 386, 502, 658, 387, 503, 658, 388, 504, 658, 389, 505, 658, 390, 506, 658, 391, 507, 658, 392, 508, 658, 393, 509, 658, 394, 510, 658, 395, 511, 658, 396, 512, 658, 397, 513, 658, 398, 514, 658, 399, 515, 658, 400, 516, 658, 401, 517, 658, 402, 518, 658, 403, 519, 658, 404, 520, 658, 405, 521, 658, 377, 493, 657, 378, 494, 657, 379, 495, 657, 380, 496, 657, 381, 497, 657, 382, 498, 657, 383, 499, 657, 384, 500, 657, 385, 501, 657, 386, 502, 657, 387, 503, 657, 388, 504, 657, 389, 505, 657, 390, 506, 657, 391, 507, 657, 392, 508, 657, 393, 509, 657, 394, 510, 657, 395, 511, 657, 396, 512, 657, 397, 513, 657, 398, 514, 657, 399, 515, 657, 400, 516, 657, 401, 517, 657, 402, 518, 657, 403, 519, 657, 404, 520, 657, 405, 521, 657, 348, 522, 658, 349, 523, 658, 350, 524, 658, 351, 525, 658, 352, 526, 658, 353, 527, 658, 354, 528, 658, 355, 529, 658, 356, 530, 658, 357, 531, 658, 358, 532, 658, 359, 533, 658, 360, 534, 658, 361, 535, 658, 362, 536, 658, 363, 537, 658, 364, 538, 658, 365, 539, 658, 366, 540, 658, 367, 541, 658, 368, 542, 658, 369, 543, 658, 370, 544, 658, 371, 545, 658, 372, 546, 658, 373, 547, 658, 374, 548, 658, 375, 549, 658, 376, 550, 658, 348, 522, 658, 349, 523, 658, 350, 524, 658, 351, 525, 658, 352, 526, 658, 353, 527, 658, 354, 528, 658, 355, 529, 658, 356, 530, 658, 357, 531, 658, 358, 532, 658, 359, 533, 658, 360, 534, 658, 361, 535, 658, 362, 536, 658, 363, 537, 658, 364, 538, 658, 365, 539, 658, 366, 540, 658, 367, 541, 658, 368, 542, 658, 369, 543, 658, 370, 544, 658, 371, 545, 658, 372, 546, 658, 373, 547, 658, 374, 548, 658, 375, 549, 658, 376, 550, 658, 348, 522, 657, 349, 523, 657, 350, 524, 657, 351, 525, 657, 352, 526, 657, 353, 527, 657, 354, 528, 657, 355, 529, 657, 356, 530, 657, 357, 531, 657, 358, 532, 657, 359, 533, 657, 360, 534, 657, 361, 535, 657, 362, 536, 657, 363, 537, 657, 364, 538, 657, 365, 539, 657, 366, 540, 657, 367, 541, 657, 368, 542, 657, 369, 543, 657, 370, 544, 657, 371, 545, 657, 372, 546, 657, 373, 547, 657, 374, 548, 657, 375, 549, 657, 376, 550, 657, 406, 551, 686, 407, 552, 686, 408, 553, 686, 409, 554, 686, 410, 555, 686, 411, 556, 686, 412, 557, 686, 413, 558, 686, 414, 559, 686, 415, 560, 686, 416, 561, 686, 417, 562, 686, 418, 563, 686, 419, 564, 686, 420, 565, 686, 421, 566, 686, 422, 567, 686, 423, 568, 686, 424, 569, 686, 425, 570, 686, 426, 571, 686, 427, 572, 686, 428, 573, 686, 429, 574, 686, 430, 575, 686, 431, 576, 686, 432, 577, 686, 433, 578, 686, 434, 579, 686, 406, 551, 686, 407, 552, 686, 408, 553, 686, 409, 554, 686, 410, 555, 686, 411, 556, 686, 412, 557, 686, 413, 558, 686, 414, 559, 686, 415, 560, 686, 416, 561, 686, 417, 562, 686, 418, 563, 686, 419, 564, 686, 420, 565, 686, 421, 566, 686, 422, 567, 686, 423, 568, 686, 424, 569, 686, 425, 570, 686, 426, 571, 686, 427, 572, 686, 428, 573, 686, 429, 574, 686, 430, 575, 686, 431, 576, 686, 432, 577, 686, 433, 578, 686, 434, 579, 686, 406, 551, 654, 685, 407, 552, 654, 685, 408, 553, 654, 685, 409, 554, 654, 685, 410, 555, 654, 685, 411, 556, 654, 685, 412, 557, 654, 685, 413, 558, 654, 685, 414, 559, 654, 685, 415, 560, 654, 685, 416, 561, 654, 685, 417, 562, 654, 685, 418, 563, 654, 685, 419, 564, 654, 685, 420, 565, 654, 685, 421, 566, 654, 685, 422, 567, 654, 685, 423, 568, 654, 685, 424, 569, 654, 685, 425, 570, 654, 685, 426, 571, 654, 685, 427, 572, 654, 685, 428, 573, 654, 685, 429, 574, 654, 685, 430, 575, 654, 685, 431, 576, 654, 685, 432, 577, 654, 685, 433, 578, 654, 685, 434, 579, 654, 685, 435, 580, 686, 436, 581, 686, 437, 582, 686, 438, 583, 686, 439, 584, 686, 440, 585, 686, 441, 586, 686, 442, 587, 686, 443, 588, 686, 444, 589, 686, 445, 590, 686, 446, 591, 686, 447, 592, 686, 448, 593, 686, 449, 594, 686, 450, 595, 686, 451, 596, 686, 452, 597, 686, 453, 598, 686, 454, 599, 686, 455, 600, 686, 456, 601, 686, 457, 602, 686, 458, 603, 686, 459, 604, 686, 460, 605, 686, 461, 606, 686, 462, 607, 686, 463, 608, 686, 435, 580, 686, 436, 581, 686, 437, 582, 686, 438, 583, 686, 439, 584, 686, 440, 585, 686, 441, 586, 686, 442, 587, 686, 443, 588, 686, 444, 589, 686, 445, 590, 686, 446, 591, 686, 447, 592, 686, 448, 593, 686, 449, 594, 686, 450, 595, 686, 451, 596, 686, 452, 597, 686, 453, 598, 686, 454, 599, 686, 455, 600, 686, 456, 601, 686, 457, 602, 686, 458, 603, 686, 459, 604, 686, 460, 605, 686, 461, 606, 686, 462, 607, 686, 463, 608, 686, 435, 580, 687, 436, 581, 687, 437, 582, 687, 438, 583, 687, 439, 584, 687, 440, 585, 687, 441, 586, 687, 442, 587, 687, 443, 588, 687, 444, 589, 687, 445, 590, 687, 446, 591, 687, 447, 592, 687, 448, 593, 687, 449, 594, 687, 450, 595, 687, 451, 596, 687, 452, 597, 687, 453, 598, 687, 454, 599, 687, 455, 600, 687, 456, 601, 687, 457, 602, 687, 458, 603, 687, 459, 604, 687, 460, 605, 687, 461, 606, 687, 462, 607, 687, 463, 608, 687, 464, 609, 652, 465, 610, 652, 466, 611, 652, 467, 612, 652, 468, 613, 652, 469, 614, 652, 470, 615, 652, 471, 616, 652, 472, 617, 652, 473, 618, 652, 474, 619, 652, 475, 620, 652, 476, 621, 652, 477, 622, 652, 478, 623, 652, 479, 624, 652, 480, 625, 652, 481, 626, 652, 482, 627, 652, 483, 628, 652, 484, 629, 652, 485, 630, 652, 486, 631, 652, 487, 632, 652, 488, 633, 652, 489, 634, 652, 490, 635, 652, 491, 636, 652, 492, 637, 652, 464, 609, 652, 465, 610, 652, 466, 611, 652, 467, 612, 652, 468, 613, 652, 469, 614, 652, 470, 615, 652, 471, 616, 652, 472, 617, 652, 473, 618, 652, 474, 619, 652, 475, 620, 652, 476, 621, 652, 477, 622, 652, 478, 623, 652, 479, 624, 652, 480, 625, 652, 481, 626, 652, 482, 627, 652, 483, 628, 652, 484, 629, 652, 485, 630, 652, 486, 631, 652, 487, 632, 652, 488, 633, 652, 489, 634, 652, 490, 635, 652, 491, 636, 652, 492, 637, 652, 464, 609, 653, 465, 610, 653, 466, 611, 653, 467, 612, 653, 468, 613, 653, 469, 614, 653, 470, 615, 653, 471, 616, 653, 472, 617, 653, 473, 618, 653, 474, 619, 653, 475, 620, 653, 476, 621, 653, 477, 622, 653, 478, 623, 653, 479, 624, 653, 480, 625, 653, 481, 626, 653, 482, 627, 653, 483, 628, 653, 484, 629, 653, 485, 630, 653, 486, 631, 653, 487, 632, 653, 488, 633, 653, 489, 634, 653, 490, 635, 653, 491, 636, 653, 492, 637, 653, 249, 253, 638, 249, 253, 638, 250, 254, 638, 250, 254, 638, 251, 255, 638, 251, 255, 638, 252, 256, 638, 252, 256, 638, 340, 344, 638, 340, 344, 638, 341, 345, 638, 341, 345, 638, 342, 346, 638, 342, 346, 638, 343, 347, 638, 343, 347, 638, 657, 662, 663, 657, 662, 663, 657, 664, 753, 657, 664, 753, 663, 759, 663, 759, 753, 754, 753, 754, 663, 665, 753, 663, 665, 753, 666, 759, 760, 666, 759, 760, 667, 759, 760, 667, 759, 761, 667, 759, 761, 668, 759, 761, 666, 754, 755, 666, 754, 755, 667, 754, 755, 667, 754, 756, 667, 754, 756, 668, 754, 756, 665, 666, 729, 665, 666, 729, 665, 667, 729, 665, 667, 730, 665, 667, 730, 665, 668, 730, 666, 667, 667, 668, 668, 711, 731, 668, 711, 731, 668, 712, 731, 668, 712, 732, 668, 712, 732, 668, 713, 732, 671, 675, 735, 671, 675, 735, 670, 675, 735, 675, 676, 677, 678, 703, 704, 673, 674, 642, 713, 724, 642, 713, 724, 643, 713, 724, 642, 643, 662, 713, 725, 662, 713, 725, 661, 713, 725, 661, 662, 656, 658, 726, 656, 658, 726, 656, 657, 726, 657, 659, 727, 657, 659, 727, 658, 659, 727, 659, 713, 728, 659, 713, 728, 660, 713, 728, 659, 660, 669, 712, 750, 669, 712, 750, 669, 711, 750, 712, 713, 669, 713, 751, 669, 712, 751, 672, 713, 733, 672, 713, 733, 673, 713, 733, 672, 673, 671, 713, 670, 671, 671, 676, 752, 671, 676, 752, 670, 676, 752, 670, 711, 673, 674, 674, 675, 734, 674, 675, 734, 674, 676, 734, 671, 676, 735, 674, 677, 736, 674, 677, 736, 674, 678, 736, 671, 677, 737, 671, 677, 737, 671, 678, 737, 678, 679, 680, 678, 679, 680, 654, 655, 656, 654, 655, 656, 654, 706, 707, 654, 706, 707, 707, 709, 738, 707, 709, 738, 707, 708, 738, 708, 709, 662, 709, 710, 662, 709, 710, 687, 688, 739, 687, 688, 739, 686, 688, 739, 687, 689, 690, 687, 689, 690, 687, 693, 694, 687, 693, 694, 696, 697, 698, 696, 697, 698, 689, 694, 740, 689, 694, 740, 694, 695, 740, 689, 695, 695, 698, 741, 695, 698, 741, 698, 699, 741, 695, 699, 687, 692, 699, 687, 692, 699, 692, 700, 742, 692, 700, 742, 692, 701, 742, 700, 701, 683, 684, 683, 684, 682, 683, 682, 683, 692, 703, 743, 692, 703, 743, 692, 704, 743, 703, 713, 744, 703, 713, 744, 704, 713, 744, 702, 703, 705, 702, 703, 705, 650, 651, 650, 651, 650, 697, 718, 650, 697, 718, 651, 721, 746, 651, 721, 746, 651, 722, 746, 721, 722, 651, 719, 747, 651, 719, 747, 651, 720, 747, 719, 720, 692, 714, 748, 692, 714, 748, 692, 715, 748, 714, 715, 647, 651, 749, 647, 651, 749, 648, 651, 749, 647, 648, 645, 646, 647, 645, 646, 647, 652, 653, 641, 642, 723, 641, 642, 723, 721, 757, 758, 721, 757, 758, 686, 684, 711, 712, 712, 713, 670, 671, 691, 692, 690, 691, 716, 717, 714, 762, 714, 762, 686, 687, 686, 687, 657, 658, 657, 658, 668, 763, 764, 668, 763, 764, 690, 765, 766, 690, 765, 766, 3, 32, 692, 3, 32, 692, 32, 33, 692, 3, 33, 143, 146, 713, 143, 146, 713, 146, 147, 713, 143, 147, 137, 148, 692, 137, 148, 692, 148, 149, 692, 137, 150, 713, 137, 150, 713, 149, 150, 713, 137, 149, 81, 82, 84, 86, 87, 89, 91, 92, 94, 95, 96, 98, 100, 101, 103, 104, 106, 108, 109, 110, 111, 112, 113, 114, 115, 116, 122, 125, 126, 128, 129, 131, 132, 133, 134, 135, 136, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 159, 161, 164, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 642, 643, 647, 648, 650, 651, 652, 653, 654, 656, 657, 660, 661, 663, 665, 667, 668, 670, 671, 673, 674, 676, 678, 680, 682, 684, 685, 686, 687, 689, 691, 692, 694, 695, 698, 699, 701, 704, 705, 707, 708, 710, 712, 713, 715, 717, 720, 722, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 758, 759, 760, 761, 762, 764, 766
};

void dxdotdw_rowvals_SPARCED(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_SPARCED_));
}
} // namespace amici
} // namespace model_SPARCED