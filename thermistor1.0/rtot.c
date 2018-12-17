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

/** @file rtot.c
 * Program calculating the temperature from resistance value of an NTC thermistor.
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
 * The program calculates the calculates the temperature from a given resistance value of an NTC
 * according to that formula.
 */

/***********
* Includes *
***********/
#include <stdio.h>
#include <stdlib.h>
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



/**************
* Prototyping *
**************/
/* Evaluates p(x) for a polynom p. */
double poly(double x, int degree, double p[]);
/* Conversion from resistance to temperature. */
double rtot(double r);

/************
* Functions *
************/
/**
 * Evaluates p(x) for a polynom p.
 * Calculates the value of polynom p at x accordings to
 * Horners schema.
 * @param p polynom.
 * @param x value to be inserted into the polynom.
 * @return calculated polynom value.
 */
double poly(double x, int degree, double p[]) {
	double retval = 0.0;
	int i;

	for (i = degree; i >= 0; i--)
		retval = retval * x + p[i];
	return retval;
}

/* Conversion from resistance to temperature.
 * Calculates and returns temperature for given resistance.
 * @param t resistance (in Ohm).
 * @return corresponding temperature.
 */
double rtot(double r)
{
	double ti;

	ti = poly(log(r), 3, a);
	ti = 1.0 / ti + TABS;
	return ti;
}
#define NTC_NOMINAL_RESISTANCE 50000
#define NTC_NOMINAL_TEMPERATURE 25
#define KELVIN_CONSTANT 273.15
#define NTC_NOMINAL_CONSTANT 3950
/*公式法计算NTC温度值*/
static float FormulaNTCTemperature(float resistance)
{
  float temp;
  float result=0.0;
 
  result=resistance/NTC_NOMINAL_RESISTANCE;
  result=(log(result)/NTC_NOMINAL_CONSTANT)+(1/(NTC_NOMINAL_TEMPERATURE+KELVIN_CONSTANT));
  temp=1/result-KELVIN_CONSTANT;
 
  return temp;
}
/**
 * Main funktion for conversion from resistance to temperature.
 * Calculates temperature for given resistance and prints result to console.
 * If called without parameters, the user is requested for resistance value,
 * otherwise all arguments are used as resistances, converted to temperature
 * and printed to console.
 * @param argc number of arguments.
 * @param argv argument list.
 * @return 0 indicating no error.
 */
int main(int argc, char *argv[])
{
	int i;
	double r;

	printf("Thermistor library version 1.0\n");
	printf("Copyright (C) 2007, 2013 - SoftQuadrat GmbH, Germany\n\n");
	if (argc > 1) {
		for (i = 1; i < argc; i++) {
			sscanf(argv[i], "%lf", &r);
			printf("tResistance : %f\tTemperature : %f\n", r, rtot(r));
		}
	}
	else {
		printf("Resistance... : ");
		scanf("%lf", &r);
		printf("Temperature.. : %f, %.1f\n", rtot(r), FormulaNTCTemperature(r));
		
	}
	return 0;
}
