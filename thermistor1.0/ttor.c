/*
 * NTC thermistor library
 * Version 1.0
 * Copyright (C) 2007, 2013 - SoftQuadrat GmbH, Germany
 * Contact: thermistor (at) softquadrat.de
 * Web site: thermistor.sourceforge.net
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 * USA
 */

/** @file ttor.c
 * Program calculating the resistance from temperature value of an NTC thermistor.
 * The Steinhart-Hart polynom allows calculation of absolute temperature
 * from resistance of an NTC thermistor
 *
 * <center><i>
 * 1/t = 1/t<sub>0</sub>
 *       + c<sub>1</sub> &middot; ln(r/r<sub>0</sub>)
 *       + c<sub>2</sub> &middot; ln(r/r<sub>0</sub>)<sup>2</sup>
 *       + c<sub>3</sub> &middot; ln(r/r<sub>0</sub>)<sup>3</sup>
 * </i></center>
 *
 * where (<b>r<sub>0</sub></b>,<b>t<sub>0</sub></b>) is a fixed resistance temperature pair.
 *
 * By substitution
 *
 * <center><i>
 * ln(r/r<sub>0</sub>) = ln(r) - ln(r<sub>0</sub>)
 * </i></center>
 *
 * this leads to a polynom in <b>ln(r)</b>
 *
 * <center><i>
 * 1/t = a<sub>0</sub>
 *       + a<sub>1</sub> &middot; ln r
 *       + a<sub>2</sub> &middot; (ln r)<sup>2</sup>
 *       + a<sub>3</sub> &middot; (ln r)<sup>3</sup>
 * </i></center>
 *
 *
 * With substitutions
 *
 * <center><i>
 * b = a<sub>2</sub>/a<sub>3</sub>
 * </i></center>
 *
 * <center><i>
 * c = a<sub>1</sub>/a<sub>3</sub>
 * </i></center>
 *
 * <center><i>
 * d = (a<sub>0</sub> - 1/t)/a<sub>3</sub>
 * </i></center>
 *
 * <center><i>
 * p = c - 1/3 * b<sup>2</sup>
 * </i></center>
 *
 * <center><i>
 * q = 2/27 &middot; b<sup>3</sup> - 1/3 * b &middot; c + d
 * </i></center>
 *
 * <center><i>
 * u = [ -q/2 + ( q<sup>2</sup>/4 + p<sup>3</sup>/27 ) <sup>1/2</sup> ] <sup>1/3</sup>
 * </i></center>
 *
 * <center><i>
 * v = [ -q/2 - ( q<sup>2</sup>/4 + p<sup>3</sup>/27 ) <sup>1/2</sup> ] <sup>1/3</sup>
 * </i></center>
 *
 * this gives
 *
 * <center><i>
 * r = e<sup>u + v - b/3</sup>
 * </i></center>
 *
 * The program calculates the calculates the temperature from a given resistance value of an NTC
 * according to that formula.
 */

/***********
* Includes *
***********/
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <stdarg.h>

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

/************
* Variables *
************/
/** Coefficients of Steinhart-Hart polynom. */
double a[] = {
	4.524024725919526e-004,
	3.934722516618191e-004,
	-7.642331765196044e-006,
	4.048572707661904e-007,
};
#define NTC_NOMINAL_RESISTANCE 50000
#define NTC_NOMINAL_TEMPERATURE 25
#define KELVIN_CONSTANT 273.15
#define NTC_NOMINAL_CONSTANT 3950
/*公式法计算NTC温度值*/
static float FormulaNTCRes(float temp)
{

  float result=0.0;
  result =-1/(NTC_NOMINAL_TEMPERATURE+KELVIN_CONSTANT) + 1/(KELVIN_CONSTANT+temp);
  result *= NTC_NOMINAL_CONSTANT;
  result = exp(result);
  result *= NTC_NOMINAL_RESISTANCE;
 
  return result;
}
/**************
* Prototyping *
**************/
/* Conversion from temperature to resistance. */
double ttor(double t);

/************
* Functions *
************/
/**
/* Conversion from temperature to resistance.
 * Calculates and returns resistance for given temperature.
 * @param t temperature (in degree Celsius).
 * @return corresponding resistance.
 */
double ttor(double t)
{
	double r;
	double u, v, p, q, b, c, d;

	t = t - TABS;
	d = (a[0] - 1.0 / t) / a[3];
	c = a[1] / a[3];
	b = a[2] / a[3];
	q = 2.0 / 27.0 * b * b * b - 1.0 / 3.0 * b * c + d;
	p = c - 1.0 / 3.0 * b * b;
	v = - pow(q / 2.0 + sqrt(q * q / 4.0 + p * p * p / 27.0), 1.0 / 3.0);
	u =   pow(-q / 2.0 + sqrt(q * q / 4.0 + p * p * p / 27.0), 1.0 / 3.0);
	r  = exp(u + v - b / 3.0);
	return r;
}
void usage(char *cmd){
	printf("Calculate A/D lookup table for NTC thermistor / resistor voltage divider network\n");
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
			// if (optopt == 'r' || optopt == 'm' || optopt == 's'){
			// 	fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			// } else if (isprint (optopt)) {
			// 	fprintf(stderr, "Unknown option '-%c'.\n", optopt);
			// } else {
			// 	fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
			// }
			/* fall through */
		default:
//			usage(argv[0]);
			return 1;
		}
	}

	printf("Thermistor library version 1.0\n");
	printf("Copyright (C) 2007, 2013 - SoftQuadrat GmbH, Germany\n\n");
	usage(argv[0]);
	if(step >0 &&end > start){
		float res, vol;
		int xv;
		for(i=start; i < end; i+=step){
			res = FormulaNTCRes(i);
			vol = res*3300 / (res+10000);
			xv = vol *4096/3300;
			if(vb==1)
				printf("%d %.1f  %.1f, %d\n", i ,res, vol, xv/* ttor(i)*/);
			else if(vb==2) printf(", %d	//%d	%d\n",  4096 - xv, i, i+20/* ttor(i)*/);
			else printf(", %d	//%d	%d\n",  xv, i, i+20/* ttor(i)*/);
		}
	}else{
		if (argc > 1) {
			for (i = 1; i < argc; i++) {
				sscanf(argv[i], "%lf", &t);
				printf("Temperature : %f\tResistance... : %f\n", t, ttor(t));
			}
		}
		else {
			printf("Temperature.. : ");
			scanf("%lf", &t);
			printf("Resistance... : %f, %.1f\n", ttor(t), FormulaNTCRes(t));
		}
		return 0;
	}
	
}
