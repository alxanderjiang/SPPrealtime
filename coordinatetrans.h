#include<iostream>
#include<math.h>
using namespace std;
//以下常数参量基于WGS-84坐标系
#define a 6378137.0
#define pi 3.1415926535897932384626433832795
#define es2 0.00669437999013//第一偏心率的平方
#define e22 0.006739496742227//第二偏心率的平方
#define clight 299792458.0
//以下参数常量基于CGCS2000坐标系

//以下为直角换大地函数
void xyztoblh(double x, double y, double z, double* str) {
    double tb, h, b, l, n;
    double temp1, temp2;
    double m1, m2, m3, m4;
    l = acos(x / sqrt(pow(x, 2) + pow(y, 2)));
    m1 = 1 / sqrt(pow(x, 2) + pow(y, 2));
    m2 = a * es2;
    m3 = 1 - es2;
    temp1 = z / sqrt(pow(x, 2) + pow(y, 2));
    //以下为迭代
    do {
        temp2 = temp1;
        temp1 = m1 * (z + m2 * temp1 / sqrt(1 + m3 * pow(temp1, 2)));
    } while (temp1 - temp2 > 1e-14);
    tb = temp1; b = atan(tb);
    n = a / sqrt(1 - es2 * pow(sin(b), 2));
    h = sqrt(pow(x, 2) + pow(y, 2)) / cos(b) - a / sqrt(1 - es2 * pow(sin(b), 2));
    if (y < 0) l = -l;
    str[0] = b * 180 / pi, str[1] = l * 180 / pi, str[2] = h;//储存结果
    //printf("大地坐标(%.4lf,%.4lf,%.4lf)\n",str[0],str[1],str[2]);
}

//以下为大地换直角函数
void blhtoxyz(double b, double l, double h, double* str) {
    double x, y, z;
    double n;
    b = b * pi / 180; l = l * pi / 180;
    n = a / sqrt(1 - es2 * pow(sin(b), 2));
    x = (n + h) * cos(b) * cos(l);
    y = (n + h) * cos(b) * sin(l);
    z = (n * (1 - es2) + h) * sin(b);
    str[0] = x, str[1] = y, str[2] = z;//储存结果
    //printf("空间直角坐标(%.10lf,%.10lf,%.10lf)\n",x,y,z);
}

//以下为大地坐标转换为高斯平面坐标
void blhtoGuss(double b, double l, double l0, double* result) {
    //公式常数项
    double X0, A, B, C, D, E, n;
    A = 1 + 3.0 / 4 * es2 + 45.0 / 63 * es2 * es2 + 175.0 / 256 * es2 * es2 * es2 + 11025.0 / 16384 * es2 * es2 * es2 * es2 + 43659.0 / 65536 * es2 * es2 * es2 * es2 * es2;
    B = 3.0 / 4 * es2 + 15.0 / 16 * es2 * es2 + 525.0 / 512 * es2 * es2 * es2 + 2205.0 / 2048 * es2 * es2 * es2 * es2 + 72765.0 / 65536 * es2 * es2 * es2 * es2 * es2;
    C = 15.0 / 64 * es2 * es2 + 105.0 / 256 * es2 * es2 * es2 + 2205.0 / 4096 * es2 * es2 * es2 * es2 + 10395.0 / 16384 * es2 * es2 * es2 * es2 * es2;
    D = 35.0 / 512 * es2 * es2 * es2 + 315.0 / 2048 * es2 * es2 * es2 * es2 + 31185.0 / 131072 * es2 * es2 * es2 * es2 * es2;
    E = 315.0 / 16384 * es2 * es2 * es2 * es2 + 3645.0 / 65536 * es2 * es2 * es2 * es2 * es2;
    X0 = a * (1 - es2) * (A * b - B / 2 * sin(2 * b) + C / 4 * sin(4 * b) - D / 6 * sin(6 * b) + E / 8 * sin(8 * b));
    n = a / sqrt(1 - es2 * pow(sin(b), 2));//卯酉圈曲率半径
    //计算公式中省略简写的参量
    double t, mew, dl;
    dl = l - l0;//中央经差
    t = tan(b);
    mew = e22 * cos(b) * cos(b);
    //主计算公式
    double x, y;
    x = X0 + n / 2 * sin(b) * cos(b) * dl * dl
        + n / 24 * sin(b) * pow(cos(b), 3) * (5 - t * t + 9 * pow(mew, 2) + 4 * pow(mew, 4)) * pow(dl, 4)
        + n / 720 * sin(b) * pow(cos(b), 5) * (61 - 58 * t * t + pow(t, 4)) * pow(dl, 6);
    y = n * cos(b) * dl
        + n / 6 * pow(cos(b), 3) * (1 - t * t + mew * mew) * dl * dl * dl
        + n / 120 * pow(cos(b), 5) * (5 - 18 * t * t + pow(t, 4) + 14 * mew * mew - 58 * mew * mew * t * t) * pow(dl, 5);
    //printf("%.4lf %.4lf\n",x,y);
    result[0] = x;
    result[1] = y;
}

void Gusstoblh(double x, double y, double l0, double* result) {
    double Bf, A, B, C, D, E, n, m;
    //迭代计算Bf底点纬度
    A = 1 + 3.0 / 4 * es2 + 45.0 / 63 * es2 * es2 + 175.0 / 256 * es2 * es2 * es2 + 11025.0 / 16384 * es2 * es2 * es2 * es2 + 43659.0 / 65536 * es2 * es2 * es2 * es2 * es2;
    B = 3.0 / 4 * es2 + 15.0 / 16 * es2 * es2 + 525.0 / 512 * es2 * es2 * es2 + 2205.0 / 2048 * es2 * es2 * es2 * es2 + 72765.0 / 65536 * es2 * es2 * es2 * es2 * es2;
    C = 15.0 / 64 * es2 * es2 + 105.0 / 256 * es2 * es2 * es2 + 2205.0 / 4096 * es2 * es2 * es2 * es2 + 10395.0 / 16384 * es2 * es2 * es2 * es2 * es2;
    D = 35.0 / 512 * es2 * es2 * es2 + 315.0 / 2048 * es2 * es2 * es2 * es2 + 31185.0 / 131072 * es2 * es2 * es2 * es2 * es2;
    E = 315.0 / 16384 * es2 * es2 * es2 * es2 + 3645.0 / 65536 * es2 * es2 * es2 * es2 * es2;
    double temp1, temp2;
    temp1 = x / (a * A * (1 - es2));//迭代初值
    //以下为迭代
    int i = 0;
    do {
        temp2 = temp1;
        temp1 = x / (a * A * (1 - es2)) + (B / 2 * sin(2 * temp1) - C / 4 * sin(4 * temp1) + D / 6 * sin(6 * temp1) - E / 8 * sin(8 * temp1)) / A;
        i++;
    } while (temp1 - temp2 > 0.0001 / 3600.0 / 180 * pi);
    Bf = temp1;//printf("%d %lf\n",i,Bf/pi*180);
    //计算公式参数
    n = a / sqrt(1 - es2 * pow(sin(Bf), 2));
    m = a * (1 - es2) / sqrt(pow(1 - es2 * pow(sin(Bf), 2), 3));
    double t, mew;
    t = tan(Bf);
    mew = e22 * cos(Bf) * cos(Bf);
    //主计算公式
    double b, l;
    b = Bf - y * y * t / 2 / m / n
        * (1 - y * y / 12 / n / n * (5 + mew * mew + 3 * t * t - 9 * mew * mew * t * t)
            + pow(y, 4) * (61 + 90 * t * t + 45 * t * t * t * t) / 360 / n / n / n / n);
    l = l0 + y / n / cos(Bf)
        * (1 - y * y / 6 / n / n * (1 + mew * mew + 2 * t * t)
            + pow(y, 4) / 120 / pow(n, 4) * (5 + 6 * mew * mew + 28 * t * t + 8 * mew * mew * t * t + 24 * t * t * t * t));
    //printf("%lf %lf\n",b,l);

    result[0] = b;
    result[1] = l;
}



//角度制弧度制变换
double degtorad(const char* deg) {
    double x1 = 0, x2 = 0, x3 = 0;
    sscanf(deg, "%lf %lf %lf", &x1, &x2, &x3);
    return (x1 + x2 / 60 + x3 / 60 / 60) / 180 * pi;
}
double degtorad(double deg) {
    return deg / 180.0 * pi;
}

void radtodeg(const double rad, double* result) {
    result[0] = rad / pi * 180;
    result[1] = int(result[0]);
    result[2] = int((result[0] - result[1]) * 60);
    result[3] = (result[0] - result[1] - result[2] / 60.0) * 60 * 60;
    //printf("%.0lf %.0lf %.4lf\n",result[1],result[2],result[3]);
}
