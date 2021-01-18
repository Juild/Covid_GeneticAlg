# ifndef __HEADER_UTILS__
# define __HEADER_UTILS__

# include <stdlib.h>
# include <stdio.h>
# include <gsl/gsl_rng.h>

//functions
void exit_error(const char *miss, int errcode);

void init_rng();
void free_rng();
double random_double();
unsigned long random_ulong();
unsigned int random_int(unsigned int n_max);

# define DAYS 101
# define N_PARAMS 5

const double DATA[DAYS][N_PARAMS] = {
	{ 1.00, 1.00, 0.00, 0.00, 0.000 },
	{ 1.84, 1.25, 0.05, 0.34, 0.000 },
	{ 2.28, 1.60, 0.12, 0.73, 0.002 },
	{ 2.57, 2.04, 0.20, 1.16, 0.004 },
	{ 2.81, 2.54, 0.30, 1.67, 0.008 },
	{ 3.05, 3.09, 0.41, 2.24, 0.013 },
	{ 3.33, 3.70, 0.54, 2.91, 0.020 },
	{ 3.65, 4.35, 0.69, 3.68, 0.029 },
	{ 4.01, 5.06, 0.87, 4.55, 0.040 },
	{ 4.42, 5.81, 1.06, 5.54, 0.054 },
	{ 4.88, 6.62, 1.27, 6.67, 0.071 },
	{ 5.40, 7.49, 1.51, 7.94, 0.091 },
	{ 5.97, 8.44, 1.77, 9.36, 0.115 },
	{ 6.60, 9.46, 2.05, 10.95, 0.143 },
	{ 7.29, 10.57, 2.36, 12.73, 0.175 },
	{ 8.06, 11.77, 2.71, 14.72, 0.212 },
	{ 8.91, 13.09, 3.08, 16.93, 0.253 },
	{ 9.84, 14.52, 3.49, 19.38, 0.301 },
	{ 10.87, 16.10, 3.93, 22.10, 0.355 },
	{ 12.00, 17.83, 4.42, 25.11, 0.415 },
	{ 13.25, 19.72, 4.95, 28.45, 0.483 },
	{ 14.63, 21.80, 5.53, 32.14, 0.559 },
	{ 16.15, 24.09, 6.17, 36.23, 0.644 },
	{ 17.82, 26.62, 6.87, 40.74, 0.738 },
	{ 19.66, 29.39, 7.63, 45.73, 0.844 },
	{ 21.70, 32.45, 8.48, 51.24, 0.960 },
	{ 23.94, 35.82, 9.40, 57.33, 1.090 },
	{ 26.41, 39.53, 10.42, 64.05, 1.233 },
	{ 29.14, 43.62, 11.54, 71.46, 1.392 },
	{ 32.15, 48.14, 12.77, 79.65, 1.569 },
	{ 35.46, 53.11, 14.12, 88.69, 1.763 },
	{ 39.12, 58.60, 15.61, 98.66, 1.979 },
	{ 43.16, 64.65, 17.25, 109.67, 2.217 },
	{ 47.61, 71.32, 19.06, 121.81, 2.480 },
	{ 52.52, 78.68, 21.05, 135.21, 2.770 },
	{ 57.94, 86.80, 23.25, 149.99, 3.091 },
	{ 63.91, 95.75, 25.67, 166.29, 3.446 },
	{ 70.50, 105.63, 28.33, 184.28, 3.837 },
	{ 77.77, 116.52, 31.27, 204.13, 4.269 },
	{ 85.79, 128.53, 34.52, 226.03, 4.745 },
	{ 94.63, 141.78, 38.09, 250.18, 5.271 },
	{ 104.39, 156.39, 42.03, 276.82, 5.851 },
	{ 115.14, 172.50, 46.38, 306.21, 6.492 },
	{ 127.01, 190.27, 51.17, 338.63, 7.198 },
	{ 140.09, 209.87, 56.45, 374.39, 7.978 },
	{ 154.52, 231.48, 62.28, 413.83, 8.838 },
	{ 170.44, 255.31, 68.70, 457.34, 9.787 },
	{ 187.99, 281.59, 75.79, 505.32, 10.833 },
	{ 207.34, 310.56, 83.60, 558.25, 11.988 },
	{ 228.68, 342.51, 92.22, 616.62, 13.261 },
	{ 252.21, 377.73, 101.72, 681.00, 14.666 },
	{ 278.16, 416.56, 112.19, 751.99, 16.215 },
	{ 306.76, 459.37, 123.74, 830.29, 17.924 },
	{ 338.29, 506.55, 136.48, 916.63, 19.809 },
	{ 373.05, 558.55, 150.51, 1011.85, 21.887 },
	{ 411.37, 615.86, 165.99, 1116.83, 24.180 },
	{ 453.60, 679.01, 183.05, 1232.60, 26.708 },
	{ 500.14, 748.59, 201.86, 1360.23, 29.496 },
	{ 551.43, 825.25, 222.58, 1500.95, 32.570 },
	{ 607.94, 909.70, 245.43, 1656.08, 35.960 },
	{ 670.19, 1002.70, 270.60, 1827.09, 39.697 },
	{ 738.77, 1105.11, 298.33, 2015.59, 43.818 },
	{ 814.30, 1217.87, 328.89, 2223.34, 48.361 },
	{ 897.47, 1341.98, 362.55, 2452.30, 53.369 },
	{ 989.04, 1478.57, 399.62, 2704.60, 58.889 },
	{ 1089.83, 1628.85, 440.44, 2982.59, 64.974 },
	{ 1200.76, 1794.15, 485.39, 3288.84, 71.680 },
	{ 1322.81, 1975.91, 534.87, 3626.19, 79.069 },
	{ 1457.05, 2175.71, 589.32, 3997.73, 87.212 },
	{ 1604.66, 2395.25, 649.23, 4406.86, 96.183 },
	{ 1766.92, 2636.39, 715.14, 4857.30, 106.065 },
	{ 1945.23, 2901.15, 787.61, 5353.11, 116.949 },
	{ 2141.07, 3191.68, 867.28, 5898.76, 128.936 },
	{ 2356.10, 3510.34, 954.82, 6499.10, 142.133 },
	{ 2592.07, 3859.64, 1050.99, 7159.42, 156.661 },
	{ 2850.89, 4242.29, 1156.59, 7885.53, 172.651 },
	{ 3134.61, 4661.16, 1272.48, 8683.70, 190.246 },
	{ 3445.42, 5119.36, 1399.61, 9560.78, 209.600 },
	{ 3785.68, 5620.13, 1538.98, 10524.19, 230.885 },
	{ 4157.89, 6166.94, 1691.68, 11582.00, 254.286 },
	{ 4564.73, 6763.43, 1858.87, 12742.90, 280.004 },
	{ 5009.00, 7413.38, 2041.80, 14016.31, 308.258 },
	{ 5493.67, 8120.75, 2241.78, 15412.35, 339.286 },
	{ 6021.85, 8889.60, 2460.21, 16941.91, 373.345 },
	{ 6596.77, 9724.09, 2698.56, 18616.67, 410.714 },
	{ 7221.77, 10628.42, 2958.39, 20449.12, 451.691 },
	{ 7900.26, 11606.81, 3241.30, 22452.54, 496.601 },
	{ 8635.72, 12663.38, 3548.97, 24641.08, 545.789 },
	{ 9431.63, 13802.12, 3883.12, 27029.66, 599.628 },
	{ 10291.43, 15026.82, 4245.50, 29634.03, 658.513 },
	{ 11218.49, 16340.93, 4637.90, 32470.69, 722.867 },
	{ 12216.01, 17747.48, 5062.07, 35556.87, 793.138 },
	{ 13286.98, 19248.93, 5519.77, 38910.42, 869.799 },
	{ 14434.06, 20847.08, 6012.67, 42549.78, 953.349 },
	{ 15659.53, 22542.93, 6542.39, 46493.79, 1044.309 },
	{ 16965.18, 24336.51, 7110.40, 50761.62, 1143.224 },
	{ 18352.19, 26226.78, 7718.01, 55372.56, 1250.660 },
	{ 19821.06, 28211.45, 8366.32, 60345.87, 1367.197 },
	{ 21371.50, 30286.92, 9056.19, 65700.54, 1493.434 },
	{ 23002.29, 32448.13, 9788.18, 71455.07, 1629.977 },
	{ 24711.26, 34688.48, 10562.48, 77627.19, 1777.437 }
};

# endif
