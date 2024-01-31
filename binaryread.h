#include<iostream>
#include"GNSStime.h"
#define MAXSIZE 1024*1024//最大容忍1M大小的观测文件
#include <fstream>
#include<time.h>
#include<string.h>
#define Frame_IN 0
#define Frame_OUT 1//定义文件指针帧内外
#define EPHYES 1
#define EPHNO 0
#define POLYCRC32 0xEDB88320u//CRC校验用
using namespace std;

unsigned short bit2ushort(unsigned char* p) {
    unsigned short r;
    memcpy(&r, p, 2);
    return r;
}
long bit2long(unsigned char* p) {
    long r;
    memcpy(&r, p, 4);
    return r;
}
double bit2double(unsigned char* p) {
    double r;
    memcpy(&r, p, 8);
    return r;
}
float bit2float(unsigned char* p) {
    float r;
    memcpy(&r, p, 4);
    return r;
}
unsigned long bit2ulong(unsigned char* p) {
    unsigned long r;
    memcpy(&r, p, 4);
    return r;
}

bool bit2bool(unsigned char* p) {
    bool r;
    memcpy(&r, p, 4);
    return r;
}

unsigned int bit2uint(unsigned char* p) {
    unsigned int r;
    memcpy(&r, p, 4);
    return r;
}


/*结构体：消息报告*/
/*包含：消息在数据中的开始、结束索引，消息类型（ID），消息长度，crc校验码*/
typedef struct {
    int start;
    int end;
    unsigned short ID;
    unsigned short Len;
    unsigned int crc;//CRC校验
}breport;

/*结构体：GPS卫星星历*/
/*包括，卫星星历存在指示与星历参数*/
typedef struct {
    int statu = EPHNO;
    double af0, af1, af2;//卫星钟差改正参数
    double tgd;//卫星电路群延时改正
    bool AS;//AS标识
    int svh;//卫星健康状况(svh)
    int week; //GPS周
    gtime_t toe, toc; //卫星星历的参考时刻，卫星钟时间
    double A, e; //轨道长半轴(A)，偏心率(e)
    double i0; //轨道倾角（i0)
    double OMG0; //升交点经度(0MG0)
    double omg; //近地点角距（omg）
    double M0; //平近点角（M0）
    double deln; //平均角速度(deln)
    double OMGD; //升交点速率（OMGd）
    double idot; //轨道倾角变化率（idot）
    double crc, crs; //地心距的摄动改正项(crc,crs)
    double cuc, cus; //升交角距的摄动改正项(cuc,cus)
    double cic, cis; //倾角的摄动改正项（cic，cis）
    double toes; //周内秒(toes)
    double N;//改正的平均角速度
}eph_t;

/*结构体：北斗三代卫星星历*/
/*包括，卫星星历存在指示与星历参数*/
typedef struct {
    int statu = EPHNO;//初始化为该卫星星历不存在
    int prn;//卫星PRN编号（BDS 1-63）
    int health;//卫星健康状态,0=healthey
    unsigned char sattype;//卫星类别(GE0/MEO/IGSO)
    unsigned char sismai;//空间信号监测精度
    unsigned short IODE, IODC;//星历、时钟数据龄期
    int week, Zweek;//以GPS时为准的周数，注意与原始数据的类型区分
    double tow;//子帧1的时间标识
    double toe;//星历参考时刻（基于GPS秒）
    double DeltaA;//参考时刻轨道长半轴相对于参考值的偏差
    double dDeltaA;//长半轴变化率
    double DeltaN;//参考时刻卫星平均角速度与计算值之差
    double dDeltaN;//参考时刻卫星平均角速度与计算值之差的变化率
    double M0;//平近点角
    double e;//偏心率
    double Omg;//近地点幅角
    double cuc, cus;//纬度幅角摄动改正
    double crc, crs;//轨道半径摄动改正
    double cic, cis;//轨道倾角摄动改正
    double I0;//轨道倾角
    double Idot;//轨道倾角变化率
    double OMG0;//升交点赤经
    double OMGD;//升交点赤经变化率
    double toc;//参考时刻的卫星钟时间
    double tgdb1cp;//B1C导频分量时延差
    double tgdB2ap;//B2A导频分量时延差
    double ISCB2ad, ISCb1cd;//B2A/B1C数据分量相对于B2A/B1C导频分量的时延修正项。
    double af0, af1, af2;//卫星钟差参数
    int iTop;//数据预测的周内时刻
    unsigned int FreqType;//信号类型（B1C/B2A）
}eph_bds3;

/*结构体：北斗二代卫星星历*/
/*包括，卫星星历存在指示与星历参数*/
typedef struct {
    int statu = EPHNO;//初始化为该卫星星历不存在
    int prn;//卫星PRN号BDS1-63
    double tow;//子帧 1 的时间标识（基于 GPS 时间），s
    unsigned long health;//卫星健康指示，0=healthey
    unsigned long AODE1, AODE2;//星历数据龄期
    unsigned long week, Zweek;//周与周计数，均以GPS时为准
    double toe;//注意：北斗时的周内秒
    double A;//轨道长半轴
    double deltN;//平均角速度改正值
    double M0;//参考时刻平近点角
    double e;//偏心率
    double omg;//近地点幅角（角距）
    double cuc, cus, crc, crs, cic, cis;//纬度、半径和轨道倾角改正数
    double I0;//参考时刻轨道倾角
    double idot;//轨道倾角变化率
    double OMG0;//升交点赤经
    double OMGD;//升交点赤经变化率
    double AODC;//时钟数据龄期
    double toc;//参考时刻的卫星钟时间，GPS秒为准
    double tgd1, tgd2;//群延迟（B1，B2）
    double af0, af1, af2;//钟差参数
    double N;//改正的平均角速度
    double URA;//用户距离精度
}eph_bds2;

void charreplace(char* str, char c, char b) {
    for (int i = 0; i < strlen(str); i++)
        if (str[i] == c)
            str[i] = b;
}

void eph2file(eph_t eph, FILE* fs, int prn) {
    COMMONTIME comtoc = time2epoch(eph.toc);
    char str[100] = {};
    fprintf(fs, "%2d%3d%3d%3d%3d%3d%5.1lf", prn, comtoc.Year % 100, comtoc.Month, comtoc.Day, comtoc.Hour, comtoc.Minute, comtoc.Second);
    sprintf(str, "%19.12E%19.12E%19.12E\n", eph.af0, eph.af1, eph.af2); charreplace(str, 'E', 'D'); fprintf(fs, str);//第一行
    sprintf(str, "   %19.12E%19.12E%19.12E%19.12E\n", 48.0, eph.crs, eph.deln, eph.M0); charreplace(str, 'E', 'D'); fprintf(fs, str);
    sprintf(str, "   %19.12E%19.12E%19.12E%19.12E\n", eph.cuc, eph.e, eph.cus, sqrt(eph.A)); charreplace(str, 'E', 'D'); fprintf(fs, str);
    sprintf(str, "   %19.12E%19.12E%19.12E%19.12E\n", eph.toes, eph.cic, eph.OMG0, eph.cis); charreplace(str, 'E', 'D'); fprintf(fs, str);
    sprintf(str, "   %19.12E%19.12E%19.12E%19.12E\n", eph.i0, eph.crc, eph.omg, eph.OMGD); charreplace(str, 'E', 'D'); fprintf(fs, str);
    sprintf(str, "   %19.12E%19.12E%19.12E%19.12E\n", eph.idot, 0.0, double(eph.week), 1.0); charreplace(str, 'E', 'D'); fprintf(fs, str);
    sprintf(str, "   %19.12E%19.12E%19.12E%19.12E\n", 2.0, double(eph.svh), eph.tgd, 48.0); charreplace(str, 'E', 'D'); fprintf(fs, str);
    sprintf(str, "   %19.12E%19.12E%19.12E%19.12E\n", 0.0, 0.0, 0.0, 0.0); charreplace(str, 'E', 'D'); fprintf(fs, str);
}


GPSTime breport2GPStime(breport r, unsigned char* data) {
    GPSTime result;
    if (r.ID == 47 || r.ID == 43 || r.ID == 7 || r.ID == 631) {
        unsigned char oweek[2] = { data[r.start + 14],data[r.start + 15] };
        unsigned char osec[4] = { data[r.start + 16],data[r.start + 17],data[r.start + 18],data[r.start + 19] };
        result.Week = bit2ushort(oweek);
        result.Second = double(bit2long(osec)) / 1000.0;
    }
    if (r.ID == 8) {
        unsigned char oweek[2] = { data[r.start + 10],data[r.start + 11] };
        unsigned char osec[4] = { data[r.start + 12],data[r.start + 13],data[r.start + 14],data[r.start + 15] };
        result.Week = bit2ushort(oweek);
        result.Second = double(bit2long(osec)) / 1000.0;
    }
    return result;
}

//CRC校验函数
unsigned int crc32(const unsigned char* buff, int len) {
    int i, j;
    unsigned int crc = 0;
    for (i = 0; i < len; i++) {
        crc ^= buff[i];
        for (j = 0; j < 8; j++) {
            if (crc & 1) crc = (crc >> 1) ^ POLYCRC32; else crc >>= 1;
        }
    }
    return crc;
}

//检查CRC校验结果
int getcrc(unsigned char* fullodata, breport r) {
    return crc32(fullodata + r.start, r.end - r.start - 3) == bit2uint(fullodata + r.end - 3);
}
//从文件中获取消息报告
int binaryfileread(breport* epoch, unsigned char* fullodata, FILE* fp) {
    unsigned char str[100];
    int epochnum = 0;
    int statu = Frame_OUT;//初始化帧状态为帧外
    int fcount = 0, count = 0;
    //中间关键变量存储
    unsigned char headerlen = 0, messagesID[2] = {}, messagesLen[2] = {};
    int flag = 0;//帧头查找标志
    while (!feof(fp)) {
        fread(str, 1, 1, fp);
        fullodata[fcount] = str[0];
        fcount++;//每次读取一位
    }
    //同步帧头
    for (int i = 0; i < fcount; i++)
    {
        if (i < fcount - 2 && fullodata[i] == 0xaa && fullodata[i + 1] == 0x44 && fullodata[i + 2] == 0x12)
            epoch[epochnum].start = i, statu = Frame_IN, headerlen = 28;
        if (i < fcount - 2 && fullodata[i] == 0xaa && fullodata[i + 1] == 0x44 && fullodata[i + 2] == 0xb5)
            epoch[epochnum].start = i, statu = Frame_IN, headerlen = 24;

        if (statu == Frame_IN && headerlen == 28) {
            count++;//记录
            if (count == 5) messagesID[0] = fullodata[i];
            if (count == 6) messagesID[1] = fullodata[i];
            if (count == 9) messagesLen[0] = fullodata[i];
            if (count == 10) messagesLen[1] = fullodata[i];

            if (count > 28 && count == 28 + int(bit2ushort(messagesLen)) + 4) {
                //printf("%d\n",fcount); 
                epoch[epochnum].end = i;
                epoch[epochnum].ID = bit2ushort(messagesID);
                epoch[epochnum].Len = bit2ushort(messagesLen);
                epoch[epochnum].crc = getcrc(fullodata, epoch[epochnum]);
                epochnum++;
                statu = Frame_OUT; count = 0;
            }
        }
        if (statu == Frame_IN && headerlen == 24) {
            count++;//记录
            if (count == 5) messagesID[0] = fullodata[i];
            if (count == 6) messagesID[1] = fullodata[i];
            if (count == 7) messagesLen[0] = fullodata[i];
            if (count == 8) messagesLen[1] = fullodata[i];

            if (count > 24 && count == 24 + int(bit2ushort(messagesLen)) + 4) {
                //printf("%d\n",fcount); 
                epoch[epochnum].end = i;
                epoch[epochnum].ID = bit2ushort(messagesID);
                epoch[epochnum].Len = bit2ushort(messagesLen);
                epoch[epochnum].crc = getcrc(fullodata, epoch[epochnum]);
                epochnum++;
                statu = Frame_OUT; count = 0;
            }
        }
    }
    return epochnum;
}
//从一定长度的二进制数组中获取消息报告
int getbinaryreport(breport* epoch, unsigned char* fullodata, int fulllen) {
    int epochnum = 0;
    int statu = Frame_OUT;//初始化帧状态为帧外
    int fcount = fulllen, count = 0;
    //中间关键变量存储
    unsigned char headerlen = 0, messagesID[2] = {}, messagesLen[2] = {};
    int flag = 0;//帧头查找标志
    //同步帧头
    for (int i = 0; i < fcount; i++)
    {
        if (i < fcount - 2 && fullodata[i] == 0xaa && fullodata[i + 1] == 0x44 && fullodata[i + 2] == 0x12)
            epoch[epochnum].start = i, statu = Frame_IN, headerlen = 28;
        if (i < fcount - 2 && fullodata[i] == 0xaa && fullodata[i + 1] == 0x44 && fullodata[i + 2] == 0xb5)
            epoch[epochnum].start = i, statu = Frame_IN, headerlen = 24;

        if (statu == Frame_IN && headerlen == 28) {
            count++;//记录
            if (count == 5) messagesID[0] = fullodata[i];
            if (count == 6) messagesID[1] = fullodata[i];
            if (count == 9) messagesLen[0] = fullodata[i];
            if (count == 10) messagesLen[1] = fullodata[i];

            if (count > 28 && count == 28 + int(bit2ushort(messagesLen)) + 4) {
                //printf("%d\n",fcount); 
                epoch[epochnum].end = i;
                epoch[epochnum].ID = bit2ushort(messagesID);
                epoch[epochnum].Len = bit2ushort(messagesLen);
                epoch[epochnum].crc = getcrc(fullodata, epoch[epochnum]);
                epochnum++;
                statu = Frame_OUT; count = 0;
            }
        }
        if (statu == Frame_IN && headerlen == 24) {
            count++;//记录
            if (count == 5) messagesID[0] = fullodata[i];
            if (count == 6) messagesID[1] = fullodata[i];
            if (count == 7) messagesLen[0] = fullodata[i];
            if (count == 8) messagesLen[1] = fullodata[i];

            if (count > 24 && count == 24 + int(bit2ushort(messagesLen)) + 4) {
                //printf("%d\n",fcount); 
                epoch[epochnum].end = i;
                epoch[epochnum].ID = bit2ushort(messagesID);
                epoch[epochnum].Len = bit2ushort(messagesLen);
                epoch[epochnum].crc = getcrc(fullodata, epoch[epochnum]);
                epochnum++;
                statu = Frame_OUT; count = 0;
            }
        }
    }
    return epochnum;
}

//二进制星历（GPS）数据读取函数，MessageID 7
unsigned long getsat(eph_t* eph, unsigned char* fullodata, breport epoch) {
    int start = epoch.start, end = epoch.end;
    unsigned long prn = bit2ulong(fullodata + start + 28);
    if (prn > 35) return 0;
    eph[prn].svh = int(bit2ulong(fullodata + start + 28 + 12));
    eph[prn].week = int(bit2ulong(fullodata + start + 28 + 24));
    eph[prn].toes = bit2double(fullodata + start + 28 + 32);
    eph[prn].A = bit2double(fullodata + start + 28 + 40);
    eph[prn].deln = bit2double(fullodata + start + 28 + 48);
    eph[prn].M0 = bit2double(fullodata + start + 28 + 56);
    eph[prn].e = bit2double(fullodata + start + 28 + 64);
    eph[prn].omg = bit2double(fullodata + start + 28 + 72);
    eph[prn].cuc = bit2double(fullodata + start + 28 + 80);
    eph[prn].cus = bit2double(fullodata + start + 28 + 88);
    eph[prn].crc = bit2double(fullodata + start + 28 + 96);
    eph[prn].crs = bit2double(fullodata + start + 28 + 104);
    eph[prn].cic = bit2double(fullodata + start + 28 + 112);
    eph[prn].cis = bit2double(fullodata + start + 28 + 120);
    eph[prn].i0 = bit2double(fullodata + start + 28 + 128);
    eph[prn].idot = bit2double(fullodata + start + 28 + 136);
    eph[prn].OMG0 = bit2double(fullodata + start + 28 + 144);
    eph[prn].OMGD = bit2double(fullodata + start + 28 + 152);
    double etoc = bit2double(fullodata + start + 28 + 164);//卫星钟时间
    eph[prn].tgd = bit2double(fullodata + start + 28 + 172);
    eph[prn].af0 = bit2double(fullodata + start + 28 + 180);
    eph[prn].af1 = bit2double(fullodata + start + 28 + 188);
    eph[prn].af2 = bit2double(fullodata + start + 28 + 196);
    eph[prn].AS = bit2ulong(fullodata + start + 28 + 204) ? true : false;
    eph[prn].N = bit2double(fullodata + start + 28 + 208);
    eph[prn].toe = gpst2time(eph[prn].week, eph[prn].toes);//卫星星历的参考时刻
    eph[prn].toc = gpst2time(eph[prn].week, etoc);
    return prn;//返回卫星prn编号
}

//二进制星历（北斗三代）读取程序，Message ID：3000
int getbds3eph(eph_bds3* eph, unsigned char* data, breport epoch)
{
    int start = epoch.start, end = epoch.end;
    unsigned char prnc = data[start + 28];
    int prn = int(prnc);//卫星PRN码
    eph[prn].prn = prn;
    eph[prn].health = int(data[start + 28 + 1]);
    eph[prn].sattype = data[start + 28 + 2];
    eph[prn].sismai = data[start + 28 + 3];
    eph[prn].IODE = bit2ushort(data + start + 28 + 4);
    eph[prn].IODC = bit2ushort(data + start + 28 + 6);
    eph[prn].week = int(bit2ushort(data + start + 28 + 8));
    eph[prn].Zweek = int(bit2ushort(data + start + 28 + 10));
    eph[prn].tow = bit2double(data + start + 28 + 12);
    eph[prn].toe = bit2double(data + start + 28 + 20);
    eph[prn].DeltaA = bit2double(data + start + 28 + 28);
    eph[prn].dDeltaA = bit2double(data + start + 28 + 36);
    eph[prn].DeltaN = bit2double(data + start + 28 + 44);
    eph[prn].dDeltaN = bit2double(data + start + 28 + 52);
    eph[prn].M0 = bit2double(data + start + 28 + 60);
    eph[prn].e = bit2double(data + start + 28 + 68);
    eph[prn].Omg = bit2double(data + start + 28 + 76);
    eph[prn].cuc = bit2double(data + start + 28 + 84);
    eph[prn].cus = bit2double(data + start + 28 + 92);
    eph[prn].crc = bit2double(data + start + 28 + 100);
    eph[prn].crs = bit2double(data + start + 28 + 108);
    eph[prn].cic = bit2double(data + start + 28 + 116);
    eph[prn].cis = bit2double(data + start + 28 + 124);
    eph[prn].I0 = bit2double(data + start + 28 + 132);
    eph[prn].Idot = bit2double(data + start + 28 + 140);
    eph[prn].OMG0 = bit2double(data + start + 28 + 148);
    eph[prn].OMGD = bit2double(data + start + 28 + 156);
    eph[prn].toc = bit2double(data + start + 28 + 164);
    eph[prn].tgdb1cp = bit2double(data + start + 28 + 172);
    eph[prn].tgdB2ap = bit2double(data + start + 28 + 180);
    eph[prn].ISCB2ad = bit2double(data + start + 28 + 188);
    eph[prn].ISCb1cd = bit2double(data + start + 28 + 196);
    eph[prn].af0 = bit2double(data + start + 28 + 204);
    eph[prn].af1 = bit2double(data + start + 28 + 212);
    eph[prn].af2 = bit2double(data + start + 28 + 220);
    eph[prn].iTop = bit2double(data + start + 28 + 228);
    eph[prn].FreqType = bit2uint(data + start + 28 + 244);
    return prn;
}
//二进制星历读取函数，北斗二代，MessageID：1047
int getbds2eph(eph_bds2* eph, unsigned char* data, breport epoch)
{
    int start = epoch.start, end = epoch.end;
    unsigned long prnul = bit2ulong(data + start + 28);
    int prn = int(prnul);//卫星PRN码
    eph[prn].prn = prn;
    eph[prn].tow = bit2double(data + start + 28 + 4);
    eph[prn].health = bit2ulong(data + start + 28 + 12);
    eph[prn].AODE1 = bit2ulong(data + start + 28 + 16);
    eph[prn].AODE2 = bit2ulong(data + start + 28 + 20);
    eph[prn].week = bit2ulong(data + start + 28 + 24);
    eph[prn].Zweek = bit2ulong(data + start + 28 + 28);
    eph[prn].toe = bit2double(data + start + 28 + 32);
    eph[prn].A = bit2double(data + start + 28 + 40);
    eph[prn].deltN = bit2double(data + start + 28 + 48);
    eph[prn].M0 = bit2double(data + start + 28 + 56);
    eph[prn].e = bit2double(data + start + 28 + 64);
    eph[prn].omg = bit2double(data + start + 28 + 72);
    eph[prn].cuc = bit2double(data + start + 28 + 80);
    eph[prn].cus = bit2double(data + start + 28 + 88);
    eph[prn].crc = bit2double(data + start + 28 + 96);
    eph[prn].crs = bit2double(data + start + 28 + 104);
    eph[prn].cic = bit2double(data + start + 28 + 112);
    eph[prn].cis = bit2double(data + start + 28 + 120);
    eph[prn].I0 = bit2double(data + start + 28 + 128);
    eph[prn].idot = bit2double(data + start + 28 + 136);
    eph[prn].OMG0 = bit2double(data + start + 28 + 144);
    eph[prn].OMGD = bit2double(data + start + 28 + 152);
    eph[prn].AODC = bit2ulong(data + start + 28 + 160);
    eph[prn].toc = bit2double(data + start + 28 + 164);
    eph[prn].tgd1 = bit2double(data + start + 28 + 172);
    eph[prn].tgd2 = bit2double(data + start + 28 + 180);
    eph[prn].af0 = bit2double(data + start + 28 + 188);
    eph[prn].af1 = bit2double(data + start + 28 + 196);
    eph[prn].af2 = bit2double(data + start + 28 + 204);
    eph[prn].N = bit2double(data + start + 28 + 216);
    eph[prn].URA = bit2double(data + start + 28 + 224);
    return prn;
};
//星历归一化函数
void bdseph2gpseph(eph_t* gps, eph_bds2* bds, int prn) {
    
    gps[prn].svh = int(bds[prn].health);
    gps[prn].week = int(bds[prn].week);
    gps[prn].toes = bds[prn].toe;
    gps[prn].A = bds[prn].A;
    gps[prn].deln = bds[prn].deltN;
    gps[prn].M0 = bds[prn].M0;
    gps[prn].e = bds[prn].e;
    gps[prn].omg = bds[prn].omg;
    gps[prn].cuc = bds[prn].cuc;
    gps[prn].cus = bds[prn].cus;
    gps[prn].crc = bds[prn].crc;
    gps[prn].crs = bds[prn].crs;
    gps[prn].cic = bds[prn].cic;
    gps[prn].cis = bds[prn].cis;
    gps[prn].i0 = bds[prn].I0;
    gps[prn].idot = bds[prn].idot;
    gps[prn].OMG0 = bds[prn].OMG0;
    gps[prn].OMGD = bds[prn].OMGD;
    double etoc = bds[prn].toc;//卫星钟时间
    gps[prn].tgd = bds[prn].tgd1;
    gps[prn].af0 = bds[prn].af0;
    gps[prn].af1 = bds[prn].af1;
    gps[prn].af2 = bds[prn].af2;
    gps[prn].N = bds[prn].N;
    gps[prn].toe = gpst2time(gps[prn].week, gps[prn].toes);//卫星星历的参考时刻
    gps[prn].toc = gpst2time(gps[prn].week, etoc);
}
/*结构体：GPS系统的观测值*/
/*包括，观测时间和观测值*/
typedef struct {
    int num = 0;//观测到卫星的数量
    gtime_t rt;//观测时间
    int name[36] = {};//记录卫星编号
    double R0[36] = {};//伪距观测值
    float psrstd[36] = {};//伪距标准差
    double adr[36] = {};//载波相位
    float adrstd[36] = {};//载波相位标准差
    float dopp[36] = {};//瞬时多普勒Hz
    float  Cno[36] = {};//载噪比dB-Hz
    float loctime[36] = {};//连续跟踪时间(无周跳)
    unsigned long trackstatu[36] = {};//连续跟踪状态
}GPSOBS;

int getobs(GPSOBS& R, unsigned char* fullodata, breport epoch) {
    int start = epoch.start, end = epoch.end, snum, name[35];
    double R0[35] = {};
    //观测时间
    COMMONTIME comt = time2epoch(gpst2time(int(breport2GPStime(epoch, fullodata).Week), breport2GPStime(epoch, fullodata).Second));
    snum = bit2long(fullodata + start + 28);
    R.num = snum;
    GPSTime obswas = breport2GPStime(epoch, fullodata);
    R.rt = gpst2time(obswas.Week, obswas.Second);
    for (int j = 0; j < snum; j++)
    {
        //PRN码
        R.name[j] = int(bit2ushort(fullodata + start + 32 + j * 44));
        //伪距观测值
        R.R0[j] = bit2double(fullodata + start + 36 + j * 44);
        R.psrstd[j] = bit2float(fullodata + start + 28 + 16 + j * 44);
        R.adr[j] = bit2double(fullodata + start + 28 + 20 + j * 44);
        R.adrstd[j] = bit2float(fullodata + start + 28 + 28 + j * 44);
        R.dopp[j] = bit2float(fullodata + start + 28 + 32 + j * 44);
        R.Cno[j] = bit2float(fullodata + start + 28 + 36 + j * 44);
        R.loctime[j] = bit2float(fullodata + start + 28 + 40 + j * 44);
    }
    return snum;
}

