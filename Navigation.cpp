#include "satpos.h"

int main()
{
    /*文件读取并生成报表*/
    clock_t tstart = clock(), et;
    FILE* fp;
    breport epoch[3000];//消息报告
    int epochnum = 0;//消息数量/历元索引
    //unsigned char str[100];
    unsigned char fullodata[MAXSIZE];//消息内容存储结构体
    const char filepath[50] = "Com16_1127212705.dat";
    if ((fp = fopen(filepath, "rb")) == NULL)
    {
        printf("Can not open file\n"); return 0;
    }
    epochnum = binaryfileread(epoch, fullodata, fp);//读取文件
    fclose(fp);
    printf("文件%s中共有%d组卫星文件\n", filepath, epochnum);

    /*冷启动初始化星历结构体*/
    eph_t eph[36] = {};//初始化GPS星历结构体
    for (int e = 0; e < epochnum; e++) {
        //星历文件
        if (epoch[e].ID == 7) {
            unsigned long prn = getsat(eph, fullodata, epoch[e]);
            if (prn == -1) continue;
            eph[prn].statu = EPHYES;
        }
    }
    /*初始化电离层延迟改正参数*/
    double ion[8];
    for (int e = 0; e < epochnum; e++) {
        if (epoch[e].ID == 8) {
            int start = epoch[e].start, end = epoch[e].end;
            getion(fullodata, ion, epoch[e]);
            //printf("%d. %.4e %.4e %.4e %.4e\n",e+1,bit2double(fullodata+start+24+32),bit2double(fullodata+start+24+40),bit2double(fullodata+start+24+48),bit2double(fullodata+start+24+56));
        }
    }
    /*逐历元读取伪距观测值并计算位置*/
    //GPSOBS R={};
    for (int e = 0; e < epochnum; e++) {
        //观测值文件
        if (epoch[e].ID == 43 || epoch[e].ID == 631) {
            GPSOBS R = {};
            SppResult Result = {};
            getobs(R, fullodata, epoch[e]);
            //计算单点定位
            SPPpos(R, eph, ion, Result);
            printResult(Result);
        }
        
        if (epoch[e].ID == 47) {
            int start = epoch[e].start, end = epoch[e].end;
            printf("%d sol statue:%d pos type%d %lf %lf %lf %lf satsnum%x\n", e + 1, bit2ulong(fullodata + start + 28), bit2ulong(fullodata + start + 28 + 4), bit2double(fullodata + start + 28 + 8), bit2double(fullodata + start + 28 + 16), bit2double(fullodata + start + 28 + 24), bit2float(fullodata + start + 28 + 32), fullodata[start + 28 + 65]);
        }
        
    }
    
    gtime_t t0 = gpst2time(2290, 1700.32);
    
    et = clock();
    printf("CPUtime %lf sec\n", (et - tstart) / 1000.0);
    return 0;
}

