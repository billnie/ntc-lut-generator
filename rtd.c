
/*
电工委员会标准IEC751的方程式：

    在-78℃到0℃的温度范围内： Rt=100[1+3.90802×10-3×t-0.5802×10-6×t2-4.27350×10-12(t-100)t^3]

    在0℃到+600℃的温度范围内： Rt=100(1+3.90802×10-3×t-0.5802×10-6×t^2)

其中：

    Rt是温度t时的阻值(单位：Ω)

    t是温度(单位：℃) 

PT100厚膜铂电阻温度传感器允许通过的工作电流为≦5mA
2、PT1000热电阻值计算

        R(t)＝R0(1＋At＋Bt2)

    R(t)是温度为t时铂热电阻的电阻值，Ω；
    t为温度，℃ 
    R0是温度为0℃时铂热电阻的电阻值，Ω；
    A、B为分度常数。
分度常数：A＝0.0038623139728     B＝-0.00000065314932626 
PT1000厚膜铂电阻温度传感器的电阻值与温度的关系偏离所给分度表的允差不超过 ±(0.30+0.005∣t∣) 
注：∣t∣为被测温度的绝对值(℃)  
PT1000厚膜铂电阻温度传感器允许通过的工作电流为≦0.5mA  

Pt100温度传感器的使用，Pt100温度传感器是一个模拟信号，它在实际应用中有二种形式：一种是不需要显示的主要采集到plc，
这样的话在使用的时候就是只需要一块pt100的集成电路，要注意的是这个集成电路采集的不是电流信号是电阻值，
pt100的集成电路(需要一个＋－12VDC电源提供工作电压)直接把采集到的电阻变为1-5VDC输入到plc，经过简单的+-
计算就可以得到相应的温度值．(这样的形式可以同时采集多路)，还有一种就是单独的一个pt100温度传感器(工作电源是24VDC)，
产生一个4-20MA的电流，然后再通过一个4-20MA电流电路板把4-20MA的电流变为1-5V电压，这个不一样的就是可以窜连一个电磁指示仪表，
其他的基本一样就不作详细说明了
*/
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
/**********
* Defines *
**********/
/** Absolute Zero. */
#define TABS (-273.15)


#define NTC_NOMINAL_RESISTANCE 50000
#define NTC_NOMINAL_TEMPERATURE 25
#define KELVIN_CONSTANT 273.15
#define NTC_NOMINAL_CONSTANT 3950

#define A 3.9083e-3 
#define B -5.775e-7 
#define C -4.183e-12 
//1、 温度计算电阻 
static float floatCCalcurtd(int fT) {
	float fR;  
	if(fT >= -200 && fT < 0)
	 {
		  fR = 100 * (1 + A*fT + B*fT*fT + C*(fT-100)*fT*fT*fT) ;  
	 } 
	 else if(fT >= 0 && fT <= 850) {
		  fR = 100 * (1 + A*fT + B*fT*fT);  
	} else{
		fR = 0; 
	 } 
	 return fR;
}

static float  CCalcuTfromRtd(float fR) { 
	 float fT  , fT0; 
	 short i ; 
	fT0 = (fR / 100 - 1) / A ;
	  if(fR >= 18.52 && fR < 100) //-200℃- 0℃ 
	  { 
		  for(i = 0 ; i < 50 ; i ++)
		   { 
			   fT = fT0 + (fR - 100*(1 + A*fT0 + B*fT0*fT0 - 100*C*fT0*fT0*fT0 + C*fT0*fT0*fT0*fT0)) / (100 * (A + 2*B*fT0 - 300*C*fT0*fT0 + 4*C*fT0*fT0*fT0)) ;

				if(fabs(fT - fT0) < 0.001) break ; 
				else fT0 = fT ; 
			} 
//			l_strT.Format(_T("%.3f") , fT); 
	} else if(fR >= 100 && fR <= 390.481) //0℃- 850℃ 
	{ 
		for(i = 0 ; i < 50 ; i ++) {
			 fT = fT0 + (fR - 100*(1 + A*fT0 + B*fT0*fT0)) / (100*(A + 2*B*fT0)) ;
			  if(fabs(fT - fT0) < 0.001) break ; 
			  else fT0 = fT ; 
		} 
//		l_strT.Format(_T("%.3f") , fT); 
	} 
	else ; 
	return fT0;
 }
/*公式法计算NTC温度值*/
static float FormulartdRes(int temp)
{

  float result=0.0;
  if(temp < 0){
//Rt=100[1+3.90802×10-3×t-0.5802×10-6×t2-4.27350×10-12(t-100)t^3]
//	result = 1 + 3.90802×10-3*temp-0.5802×10-6*temp *temp-4.27350×10-12*(temp-100)*temp*temp*temp;
	result = 1 + 3.90802*10-3*temp-0.5802*10-6*temp *temp- 4.27350*10-12*(temp-100)*temp*temp*temp; 
  }else{
// Rt=100(1+3.90802×10-3×t-0.5802×10-6×t^2)
	result = 100*(1+3.90802*10-3*temp-0.5802*10-6*temp*temp);
//	result = 100*(1+3.90802×10-3×temp-0.5802×10-6×temp*temp)	  ;
  }
  
  return result;
}

void usage(char *cmd){
	printf("Calculate rtd and t\n");
	printf("using temperature-resstance data points from a file for the thermistor as input.\n");
	printf("Usage: %s [-s -40 -e 150 [-s AD_STEP (default=32)] -v 1 filename\n v = 0, 1,2\n", cmd);	
}
// ttor -s -20 e 150 
/**
 * Main funktion for conversion from temperature to resistance.
 * Calculates resistance for given temperature and prints result to console.
 * If called without parameters, the user is requested for temperature value,
 * otherwise all arguments are used as temperatures, converted to resistance
 * and printed to console.
 * @param argc number of arguments.
 * @param argv argument list.
 * @return 0 indicating no error.
 */
int main(int argc, char *argv[])
{
	int i;
	double t;
	int start , end, step, vb;
	/* getopt */
	int c;
	start = -40; end = 150; step =0; vb =0;
	while((c=getopt(argc, argv, "s:e:o:v:")) != -1) {
		switch(c) {
		case 's':
			start=atoi(optarg);
			break;
		case 'o':
			step=atoi(optarg);
			break;
		case 'e':
			end=atoi(optarg);
			break;
		case 'v':
			vb=atoi(optarg);
			break;	
		case '?':
			usage(argv[0]);
		default:
//			usage(argv[0]);
			return 1;
		}
	}

	printf("rtd library version 1.0\n");
	printf("Copyright (C) 2018 - billnie\n\n");
	usage(argv[0]);
	printf("start=%d end=%d step = %d\n", start, end, step);
	if(step >0 &&end > start){
		float res, vol,vol2;
		int xv,k=0;
		for(i=start; i < end; i+=step){
			res = floatCCalcurtd(i);
			vol = res*12.5 ;
			xv = vol *4096/3300; 
			if(vb==1)
				printf("%d %.1f  %.1f, %d\n", i ,res, vol, xv/* ttor(i)*/);
			else if(vb==2) printf(", %d	//%d	%d\n", 4096 - xv, i, i+20/* ttor(i)*/);
			else printf(",%d	// %.2f,  %d,  %d\n", (int)xv, res, i,k);
			k++;
		}
	}else{
		while(1){
			int t;
			float fr;
			if(vb ==1){
				scanf("%d", &t);
				printf("temp=%d res = %.2f\n", t, floatCCalcurtd(t));
			}else{
				scanf("%f", &fr);
				printf("res=%.3f temp = %.2f\n", fr, CCalcuTfromRtd(fr));
			}
		}

	}
}
