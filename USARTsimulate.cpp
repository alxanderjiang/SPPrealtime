#include<Windows.h>//由于链接器问题，Windows.h头文件最好包含在第一行
#include"satpos.h"

typedef struct {
    unsigned char buff[8192] = {};
    int bufflen = 0;
}USARTbuff;


int main() {

    srand((unsigned)time(NULL));
    //创建结果文件
    FILE* fr, * fd;
    const char frname[50] = "BDStest.txt";
    const char fdname[50] = "0114direct.txt";
    //fr = fopen(frname, "w"); 
    //fd = fopen(fdname, "w");
    //fprintf(fr, "timestamp\tB\tL\tH\tDtr\tvX\tvY\tvZ\tvDtr\n");
    //fprintf(fr, "timestamp\tX\tY\tZ\tDtr\tvX\tvY\tvZ\tvDtr\n");
    //fprintf(fd, "timestamp\tB\tL\tH\tX\t\Y\tZ\n");

    //模拟串口数据的文件名
    //char file[100] = "20220927OBS+NAV+PSRPOS+IONUTC.dat";
    char file[100] = "BDStest.txt";
    FILE* fp;
    unsigned char fulldata[MAXSIZE];//最大支持1M的文件模拟
    int maxlen = 0;
    if ((fp = fopen(file, "rb")) == NULL)
    {
        printf("Can not open file\n"); return 0;
    }

    while (!feof(fp))
    {
        char str[100];
        fread(str, 1, 1, fp);
        fulldata[maxlen] = str[0];
        maxlen++;
    }
    printf("读取的文件中一共有%d字节的数据\n", maxlen);

    //模拟串口数据的缓冲区
    USARTbuff buff;//缓冲区
    unsigned char message[2 * 8192];//数据区(上次缓冲区剩余＋本次缓冲区)

    //开始模拟串口传输
    int i = 0; USARTbuff left = {};
    //定义接收机启动后的全局变量
    eph_t eph[36] = {};//GPS星历
    eph_bds2 bds[65] = {};//北斗星历
    eph_t bds_gps[65] = {};//归一化北斗星历
    double ion[8] = {};//电离层延迟改正参数结构体
    int slen = 0; int history = 0;
    int testlen[10] = { 400,512,225,428,654,713,500,481,442,256};
    while (true) {
        //固定模拟的单次长度为8192
        slen = rand()%400;//模拟长度
        memcpy(buff.buff, fulldata + history, slen);
        buff.bufflen = slen;
        i++; history += slen;

        //构建有效消息序列数据区(上次缓冲区剩余＋本次缓冲区)
        int messagex = 0;
        {
            memcpy(message, left.buff, left.bufflen);
            memcpy(message + left.bufflen, buff.buff, buff.bufflen);
            messagex = buff.bufflen + left.bufflen;
        }


        //从有效数据区中获取消息报告
        breport epoch[400] = {};
        int epochnum = 0;
        epochnum = getbinaryreport(epoch, message, messagex);
        memset(left.buff, 0x00, sizeof(left.buff)); left.bufflen = 0;//上一次剩余数据区清零
        left.bufflen = messagex  - epoch[epochnum].start;//重新构建剩余数据区结构(定长)
        memcpy(left.buff, message + epoch[epochnum].start, left.bufflen);//重新构建剩余数据区

        if (epochnum == 0)    continue;

        //定位主循环
        for (int j = 0; j < epochnum; j++) {
            GPSOBS R = {};//观测值结构体

            //北斗二代星历
            if (epoch[j].ID == 1047) {
                int prn;
                prn = getbds2eph(bds, message, epoch[j]);
                bds[prn].statu = EPHYES;
                bdseph2gpseph(bds_gps, bds, prn);
                bds_gps[prn].statu = EPHYES;
            }
            
            
            //GPS星历
            if (epoch[j].ID == 7) {
                int prn;
                prn = getsat(eph, message, epoch[j]);
                eph[prn].statu = EPHYES;
            }
            //电离层，因为老师提供的文件有两种板卡协议，0927数据一定要加上"OEM7"参数调用重载，Com_16不需要
            if (epoch[j].ID == 8) {
                getion(message, ion, epoch[j]);
            }
            //如果读到观测值，则进行定位尝试
            if (epoch[j].ID == 631 || epoch[j].ID == 43) {
                SppResult Result;
                getobs(R, message, epoch[j]);
                SPPpos(R, eph, ion, Result);
                printf("GPS单系统定位结果\n");
                printResult(Result,"BLH");

                SPPpos(R, bds_gps, ion, Result);
                printf("BDS单系统定位结果\n");
                printResult(Result, "BLH");
                //fprintResult(fr, Result);
                //Sleep(900);//模拟串口数据延时
            }
            //输出接收机原始结果
            if (epoch[j].ID == 47) {
                int start = epoch[j].start;
                double lat = bit2double(message + start + 28 + 8);
                double lon = bit2double(message + start + 28 + 16);
                double he = bit2double(message + start + 28 + 24) + bit2float(message + start + 28 + 32);
                //注意接收机直接输出的是海拔高，需要与其内置的水准面模型相加才能得到直接解算的椭球高he
                printf("\n接收机解算结果 WGS-84 %lf %lf %lf\n\n", lat, lon, he);
                //fprintf(fd, "%lld %lf %lf %lf ", gpst2time(breport2GPStime(epoch[j], message).Week, breport2GPStime(epoch[j], message).Second).time, lat, lon, he);
                double xyz[3] = {};
                blhtoxyz(lat, lon, he, xyz);
                //fprintf(fd, "%lf %lf %lf\n", xyz[0], xyz[1], xyz[2]);
                //printcommtime(time2epoch(gpst2time(breport2GPStime(epoch[j], message).Week, breport2GPStime(epoch[j], message).Second)));
            }
        }


        //printf("第%d条缓冲区数据结束于%d位置,剩余%d移入left\n", i,epoch[epochnum-1].end,left.bufflen);

        if (history > maxlen)
            break;
    }
    return 0;
}