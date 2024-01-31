#include<Windows.h>//�������������⣬Windows.hͷ�ļ���ð����ڵ�һ��
#include"satpos.h"

typedef struct {
    unsigned char buff[8192] = {};
    int bufflen = 0;
}USARTbuff;


int main() {

    srand((unsigned)time(NULL));
    //��������ļ�
    FILE* fr, * fd;
    const char frname[50] = "BDStest.txt";
    const char fdname[50] = "0114direct.txt";
    //fr = fopen(frname, "w"); 
    //fd = fopen(fdname, "w");
    //fprintf(fr, "timestamp\tB\tL\tH\tDtr\tvX\tvY\tvZ\tvDtr\n");
    //fprintf(fr, "timestamp\tX\tY\tZ\tDtr\tvX\tvY\tvZ\tvDtr\n");
    //fprintf(fd, "timestamp\tB\tL\tH\tX\t\Y\tZ\n");

    //ģ�⴮�����ݵ��ļ���
    //char file[100] = "20220927OBS+NAV+PSRPOS+IONUTC.dat";
    char file[100] = "BDStest.txt";
    FILE* fp;
    unsigned char fulldata[MAXSIZE];//���֧��1M���ļ�ģ��
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
    printf("��ȡ���ļ���һ����%d�ֽڵ�����\n", maxlen);

    //ģ�⴮�����ݵĻ�����
    USARTbuff buff;//������
    unsigned char message[2 * 8192];//������(�ϴλ�����ʣ�࣫���λ�����)

    //��ʼģ�⴮�ڴ���
    int i = 0; USARTbuff left = {};
    //������ջ��������ȫ�ֱ���
    eph_t eph[36] = {};//GPS����
    eph_bds2 bds[65] = {};//��������
    eph_t bds_gps[65] = {};//��һ����������
    double ion[8] = {};//������ӳٸ��������ṹ��
    int slen = 0; int history = 0;
    int testlen[10] = { 400,512,225,428,654,713,500,481,442,256};
    while (true) {
        //�̶�ģ��ĵ��γ���Ϊ8192
        slen = rand()%400;//ģ�ⳤ��
        memcpy(buff.buff, fulldata + history, slen);
        buff.bufflen = slen;
        i++; history += slen;

        //������Ч��Ϣ����������(�ϴλ�����ʣ�࣫���λ�����)
        int messagex = 0;
        {
            memcpy(message, left.buff, left.bufflen);
            memcpy(message + left.bufflen, buff.buff, buff.bufflen);
            messagex = buff.bufflen + left.bufflen;
        }


        //����Ч�������л�ȡ��Ϣ����
        breport epoch[400] = {};
        int epochnum = 0;
        epochnum = getbinaryreport(epoch, message, messagex);
        memset(left.buff, 0x00, sizeof(left.buff)); left.bufflen = 0;//��һ��ʣ������������
        left.bufflen = messagex  - epoch[epochnum].start;//���¹���ʣ���������ṹ(����)
        memcpy(left.buff, message + epoch[epochnum].start, left.bufflen);//���¹���ʣ��������

        if (epochnum == 0)    continue;

        //��λ��ѭ��
        for (int j = 0; j < epochnum; j++) {
            GPSOBS R = {};//�۲�ֵ�ṹ��

            //������������
            if (epoch[j].ID == 1047) {
                int prn;
                prn = getbds2eph(bds, message, epoch[j]);
                bds[prn].statu = EPHYES;
                bdseph2gpseph(bds_gps, bds, prn);
                bds_gps[prn].statu = EPHYES;
            }
            
            
            //GPS����
            if (epoch[j].ID == 7) {
                int prn;
                prn = getsat(eph, message, epoch[j]);
                eph[prn].statu = EPHYES;
            }
            //����㣬��Ϊ��ʦ�ṩ���ļ������ְ忨Э�飬0927����һ��Ҫ����"OEM7"�����������أ�Com_16����Ҫ
            if (epoch[j].ID == 8) {
                getion(message, ion, epoch[j]);
            }
            //��������۲�ֵ������ж�λ����
            if (epoch[j].ID == 631 || epoch[j].ID == 43) {
                SppResult Result;
                getobs(R, message, epoch[j]);
                SPPpos(R, eph, ion, Result);
                printf("GPS��ϵͳ��λ���\n");
                printResult(Result,"BLH");

                SPPpos(R, bds_gps, ion, Result);
                printf("BDS��ϵͳ��λ���\n");
                printResult(Result, "BLH");
                //fprintResult(fr, Result);
                //Sleep(900);//ģ�⴮��������ʱ
            }
            //������ջ�ԭʼ���
            if (epoch[j].ID == 47) {
                int start = epoch[j].start;
                double lat = bit2double(message + start + 28 + 8);
                double lon = bit2double(message + start + 28 + 16);
                double he = bit2double(message + start + 28 + 24) + bit2float(message + start + 28 + 32);
                //ע����ջ�ֱ��������Ǻ��θߣ���Ҫ�������õ�ˮ׼��ģ����Ӳ��ܵõ�ֱ�ӽ���������he
                printf("\n���ջ������� WGS-84 %lf %lf %lf\n\n", lat, lon, he);
                //fprintf(fd, "%lld %lf %lf %lf ", gpst2time(breport2GPStime(epoch[j], message).Week, breport2GPStime(epoch[j], message).Second).time, lat, lon, he);
                double xyz[3] = {};
                blhtoxyz(lat, lon, he, xyz);
                //fprintf(fd, "%lf %lf %lf\n", xyz[0], xyz[1], xyz[2]);
                //printcommtime(time2epoch(gpst2time(breport2GPStime(epoch[j], message).Week, breport2GPStime(epoch[j], message).Second)));
            }
        }


        //printf("��%d�����������ݽ�����%dλ��,ʣ��%d����left\n", i,epoch[epochnum-1].end,left.bufflen);

        if (history > maxlen)
            break;
    }
    return 0;
}