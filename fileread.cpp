//文件名：fileread.cpp
/*功能：将二进制文件读取出来并转写为RINEX格式的数据文档*/
/*经过这个主函数运行后，同目录下应当生成两个后缀为.21o的RNIEX观测文档*/
#include"binaryread.h"
int main() {
    clock_t tstart = clock(), et;
    FILE* fp;
    breport epoch[3000]; int epochnum = 0;
    //unsigned char str[100];
    unsigned char fullodata[MAXSIZE];
    //int statu=Frame_OUT;//初始化文件指针指示到数据帧外部
    const char filepath[50] = "BDStest.txt";
    if ((fp = fopen(filepath, "rb")) == NULL)
    {
        printf("Can not open file\n"); return 0;
    }
    epochnum = binaryfileread(epoch, fullodata, fp);//读取文件
    fclose(fp);

    printf("文件%s中共有%d组卫星文件\n", filepath, epochnum);
    /*
    for(int i=0;i<epochnum;i++){
        GPSTime temptime;temptime=breport2GPStime(epoch[i],fullodata);
        printf("%d.索引为[%d,%d],Message ID=%d,Message Length=%d\n",i+1,epoch[i].start,epoch[i].end,epoch[i].ID,epoch[i].Len);
        printf("观测时间为GPS时 %d 周 %lf 秒,PRN=G%02d,未跳秒协调世界时: ",temptime.Week,temptime.Second,bit2ulong(fullodata+epoch[i].start+28));
        printcommtime(time2epoch(gpst2time(int(temptime.Week),temptime.Second)));
    }
    */
    //现在开始处理卫星文件
    eph_t eph[36] = {};//GPS系统星历结构体
    FILE* fs1, * fs2;
    fs1 = fopen("220911271555.21o", "w");
    fs2 = fopen("220911271555e.21o", "w");
    for (int i = 0; i < epochnum; i++) {
    int start = epoch[i].start, end = epoch[i].end;
    if (epoch[i].ID == 8)
    {
        char temp[200] = {};
        sprintf(temp, "   %11.4E %11.4E %11.4E %11.4E          ION ALPHA\n", bit2double(fullodata + start + 24), bit2double(fullodata + start + 32), bit2double(fullodata + start + 40), bit2double(fullodata + start + 48));
        charreplace(temp, 'E', 'D'); fprintf(fs2, temp);
        sprintf(temp, "   %11.4E %11.4E %11.4E %11.4E          ION BETA\n", bit2double(fullodata + start + 56), bit2double(fullodata + start + 64), bit2double(fullodata + start + 72), bit2double(fullodata + start + 80));
        charreplace(temp, 'E', 'D'); fprintf(fs2, temp);
        break;
    }
    }
    
    fprintf(fs1, "RINEX                            END OF HEADER\n");
    fprintf(fs2, "RINEX                            END OF HEADER\n");
    
    for (int i = 0; i < epochnum; i++) {
        //首先检查消息是否校验正确
        if (epoch[i].crc == 0) { 
            printf("第%d号消息校验未通过\n", i); 
            continue; 
        }

        if (epoch[i].ID == 7)//卫星星历文件
        {
            unsigned long prn = getsat(eph, fullodata, epoch[i]);
            
            if (!prn) continue;//prn号不属于GPS，跳出
           
            eph[prn].statu = EPHYES;
            eph2file(eph[prn], fs2, prn);
            //printf("%d\n",bit2ulong(fullodata+start+28+20));
        }

        if (epoch[i].ID == 631 || epoch[i].ID == 43)//卫星伪距观测值文件
        {

            int start = epoch[i].start, end = epoch[i].end, snum, name[35];
            double R0[35] = {};
            //观测时间
            COMMONTIME comt = time2epoch(gpst2time(int(breport2GPStime(epoch[i], fullodata).Week), breport2GPStime(epoch[i], fullodata).Second));
            snum = bit2long(fullodata + start + 28);
            for (int j = 0; j < snum; j++)
            {
                //PRN码
                name[j] = bit2ushort(fullodata + start + 32 + j * 44);
                //伪距观测值
                R0[j] = bit2double(fullodata + start + 36 + j * 44);
            }

            fprintf(fs1, "> %04d  %02d  %02d  %02d  %02d  %lf\n", comt.Year, comt.Month, comt.Day, comt.Hour, comt.Minute, comt.Second);
            for (int i = 0; i < snum; i++)
                fprintf(fs1, "G%02d  %lf\n", name[i], R0[i]);
        }
        if (epoch[i].ID == 47) {
            int start = epoch[i].start, end = epoch[i].end;
            printf("%d sol statue:%d pos type%d %lf %lf %lf %lf satsnum%x\n", i + 1, bit2ulong(fullodata + start + 28), bit2ulong(fullodata + start + 28 + 4), bit2double(fullodata + start + 28 + 8), bit2double(fullodata + start + 28 + 16), bit2double(fullodata + start + 28 + 24), bit2float(fullodata + start + 28 + 32), fullodata[start + 28 + 66]);
        }
        if (epoch[i].ID == 8) {
            int start = epoch[i].start, end = epoch[i].end;
            
        }
    }

    
    //星历最后更新值
    for(int i=0;i<36;i++){
        if(eph[i].statu==EPHYES)
            eph2file(eph[i],fs2,i);
    }

    for (int i = 0; i < epochnum-1; i++) {
        printf("%d\n", epoch[i + 1].start - epoch[i].end);
    }
    et = clock();
    printf("\nCPUtime: %f sec", double(et - tstart) / CLOCKS_PER_SEC);
    return 0;
}