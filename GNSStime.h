//文件名GNSStime.h
#include <stdio.h>
#include <math.h>
#include <iostream>
#define SECOND 315964800
//定义Unix时结构
struct gtime_t {
    time_t time;
    double second;
};
//定义GPS时结构
struct GPSTime {
    unsigned short Week;
    double Second;
};
//定义通用计时结构
struct COMMONTIME {
    unsigned short Year;
    unsigned short Month;
    unsigned short Day;
    unsigned short Hour;
    unsigned short Minute;
    double Second;
};
void printcommtime(struct COMMONTIME p)
{
    printf("%d 年 %d 月 %d 日 %d 时 %d 分 %lf 秒\n", p.Year, p.Month, p.Day, p.Hour, p.Minute, p.Second);
}
//通用时转UnixTime
struct gtime_t str2time(const char* str) {
    struct COMMONTIME t0;
    struct gtime_t t;
    int i = 0;
    //按字符串输入结构读取年月日时分秒至通用时结构体中
    t0.Year = (str[0] - '0') * 1000 + (str[1] - '0') * 100 + (str[2] - '0') * 10 + (str[3] - '0') * 1;
    t0.Month = (str[5] - '0') * 10 + (str[6] - '0') * 1;//printf("%d",t0.Month);
    t0.Day = (str[8] - '0') * 10 + (str[9] - '0') * 1;//printf("%d",t0.Day);
    t0.Hour = (str[11] - '0') * 10 + (str[12] - '0') * 1;//printf("%d",t0.Hour);
    t0.Minute = (str[14] - '0') * 10 + (str[15] - '0') * 1;
    t0.Second = (str[17] - '0') * 10 + (str[18] - '0') * 1;
    //开始计算通用时距Unixtime起点时的Days
    int Days = 0;
    /*
    for(int i=1970;i<=t0.Year;i++){
        if(i!=t0.Year){
            Days+=365;
            if((i%4==0&&i%100!=0)||i%400==0)
                    Days++;//闰年计365天
        }
        else{
            for(int j=1;j<t0.Month;j++){
                if(j==1||j==3||j==5||j==7||j==8||j==10||j==12)
                    Days+=31;
                if(j==4||j==6||j==9||j==11)
                    Days+=30;
                if(j==2&&t0.Year%4==0)
                    Days+=29;
                if(j==2&&t0.Year%4!=0)
                    Days+=28;
            }
            Days+=t0.Day;
            Days--;//减去最后一天
        }
    }*/
    int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
    Days = (t0.Year - 1970) * 365 + (t0.Year - 1969) / 4 + doy[t0.Month - 1] + t0.Day - 2 + (t0.Year % 4 == 0 && t0.Month >= 3 ? 1 : 0);//老师提供算法
    t.time = Days * 86400 + t0.Hour * 3600 + t0.Minute * 60 + floor(t0.Second);
    t.second = t0.Second - floor(t0.Second);
    return  t;
}
//GPS时转Unix time
struct gtime_t gpst2time(int week, double seconds) {
    gtime_t result = {};
    result.second = seconds - floor(seconds);
    result.time = time_t(86400 * 7 * week) + int(seconds) + SECOND;
    return result;
}

double time2gpst(gtime_t t0, int& week) {
    time_t sec = t0.time - SECOND;
    week = int(sec / (86400 * 7));
    double second = (double)(sec - week * 86400 * 7) + t0.second;
    return second;
}
//unix时转通用时
struct COMMONTIME time2epoch(gtime_t t0) {
    int Days = (int)(t0.time) / 86400;
    double secs = (int)(t0.time - (time_t)Days * 86400);

    int mday[] = { 31,28,31,30,31,30,31,31,30,31,30,31,
                31,28,31,30,31,30,31,31,30,31,30,31,
                31,29,31,30,31,30,31,31,30,31,30,31,
                31,28,31,30,31,30,31,31,30,31,30,31 };
    int d, mon;
    for (d = Days % 1461, mon = 0; mon < 48; mon++) {
        if (d >= mday[mon])
            d -= mday[mon];
        else
            break;
    }

    struct COMMONTIME comtime;
    comtime.Year = 1970 + Days / 1461 * 4 + mon / 12;
    comtime.Month = mon % 12 + 1;
    comtime.Day = d + 1;
    comtime.Hour = secs / 3600;
    comtime.Minute = (secs - comtime.Hour * 3600) / 60;
    comtime.Second = secs - comtime.Hour * 3600 - comtime.Minute * 60;
    return comtime;
}
