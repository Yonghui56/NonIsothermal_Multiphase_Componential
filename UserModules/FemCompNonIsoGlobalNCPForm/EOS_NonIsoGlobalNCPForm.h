/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file EOS_H2_H2O_ISOT.h
 *
 * Created on 2014-05-12 by Yonghui HUANG
 */

#pragma once


#include "math.h"
#include <algorithm>
#include "AbstractEOS_NonIsoGlobalNCPForm.h"
#include <cmath>

/**
 *  EOS formulation for calculating each secondary variable
 *  PC PG PL N_L N_G Kr_g Kr_l X_m X_M Sg
 */
class EOS_NonIsoGlobalNCPForm : public AbstractEOS_NonIsoGlobalNCPForm
{
public:
    
	/**
	  * Constructor
	  */
	EOS_NonIsoGlobalNCPForm() : AbstractEOS_NonIsoGlobalNCPForm(2)
    {
    };

	/**
	  * DESTRUCTOR
	  */
	virtual ~EOS_NonIsoGlobalNCPForm()
    {

    };

	


	virtual void set_env_condition(ogsChem::LocalVector & env_condition)
	{
		// mean pressure
		P_G = env_condition(0);
		// molar fraction of the light component
		X_L = env_condition(1);
		//SATURATION of the gas phase
		S_G = env_condition(2);

		T_L = env_condition(3);

	}; 
	/**
	* realization of the calc the residual of saturation
	* EOS--Function 1
	* X=Sl*rho_L_h+Sg*rho_G_h
	*/
	/**
	*  return mass fraction of air in gas phase
	*/
	virtual double get_massfraction(double x_G_a)
	{
		double massfractionX(0.0);
		massfractionX = x_G_a / (x_G_a + (1 - x_G_a)*M_L / M_G);
		return massfractionX;
	}

	virtual double Deriv_massfraction(double x_G_a)
	{
		//double dX(0.0);
		//dX = (get_massfraction(x_G_a + eps) - get_massfraction(x_G_a - eps)) / 2 / eps;
		double dX2 = ((x_G_a + (1 - x_G_a)*M_L / M_G) - x_G_a*(1 - M_L / M_G)) / pow((x_G_a + (1 - x_G_a)*M_L / M_G), 2);
		return dX2;
	}
	/**
	* Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
	*/
    virtual double getPcbySg(double Sg)//for waste package case1 n=2.0
    {
        double Pc = 0.0;
        double S_eff = getEffectSat_lbySg(Sg);
        double p0 = 6.324555320336758e+05;
        Pc = p0*0.05878*((1.263*(1 - S_eff) - 2.120)*(1 - S_eff) + 1.417) * (1 - S_eff);
        if (Sg > 1 - S_lr) {
            double y = p0*0.05878*((1.263*1.0 - 2.120)*1.0 + 1.417)*1.0;
            double m = p0*0.05878*((3 * 1.263*1.0 - 2 * 2.120)*1.0 + 1.417);// *(1 / (1 - S_lr));
            return -S_eff*m + y;
        }
        else if (Sg < 0.0) {
            double y = 0.0;
            double m = p0*0.05878*1.417;// *(1 / (1 - S_lr));
            return (1 - S_eff)*m + y;
        }
        return Pc;
    }
    /**
    * Derivation of PC (capillary pressure) in terms of saturation based on van Genuchten model
    */
    virtual double Deriv_dPCdS(double Sg)
    {
        double dPcdSg(0.0);
        double S_eff = getEffectSat_lbySg(Sg);

        double p0 = 6.324555320336758e+05;

        if (Sg > 1.0 - S_lr)
            S_eff = 0.0;
        dPcdSg = p0*0.05878*((3 * 1.263*(1 - S_eff) - 2 * 2.120)*(1 - S_eff) + 1.417)*(1 / (1 - S_lr));
        if (Sg < 0.0) {
            dPcdSg = p0*0.05878*1.417*(1 / (1 - S_lr));
        }

        return dPcdSg;
    }
    /**
    * calculate the effective saturation of Liquid phase by Gas saturation
    */
    virtual double getEffectSat_lbySg(double Sg)
    {
        double EffectSat_l = 0.0;
        // this->res_saturation_g->eval(Sg, res_S_g);

        /*if (Sg < S_gr)
        EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);
        else if (Sg>(1 - S_lr))
        EffectSat_l = 0.0;
        else*/
        EffectSat_l = (1 - Sg - S_lr) / (1 - S_gr - S_lr);

        return EffectSat_l;
    }
	virtual double get_RHO_G_W(double P_gw, double T)
	{
		double rho_G_w(0.0);
		double C_w = M_L / R / T;
		rho_G_w = P_gw*C_w;
		return rho_G_w;
	}


	virtual double get_RHO_G_a(double P_G_a, double X,double T)
	{
		double rho_G_a(0.0);
		double C_v = M_G / R / T;
		rho_G_a = P_G_a*C_v;// std::max(P_G_a*C_v, M_G*X / M_L);
		return rho_G_a;
	}
	virtual double get_P_sat_S(double T)
	{
		// Here unit of T is Celsius;
		double P_sat(0.0);
		double T_0 = 373.15;
		double P_0 = 101325.0;
		double h_wg = 2258000.0;
		P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//if (T>T_0+2)
			//P_sat = P_0;
		return P_sat;
	}

	/*
	*Kelvin equation
	*/
	virtual double get_P_G_w(double /*PG*/, double PC, double T)
	{
		double P_sat(0.0);
		double P_gw(0.0);
		P_sat = get_P_sat_S(T);
		double C_w = M_L / R / T;
		P_gw = P_sat*exp(-(PC)*C_w / rho_l_std);
		//if (get_P_G_w_ORG(PC, T) > PG)
			//P_gw = PG*0.995;
		return P_gw;
	}
	
	virtual double get_P_G_w_ORG(double PC, double T)
	{
		double P_sat(0.0);
		double P_gw(0.0);
		P_sat = get_P_sat_S(T);
		double C_w = M_L / R / T;
		P_gw = P_sat*exp(-PC*C_w / rho_l_std);
		return P_gw;
	}
	

	virtual double Deriv_dPsat_dT(double PG, double T)
	{
		// Here unit of T is Celsius;
		double dPsat_dT(0.0);
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		dPsat_dT = P_0*(M_L*h_wg / R)*(1. / T / T)*exp(((1. / T_0) - (1. / T))*M_L*h_wg / R);
		double P_sat = P_0*exp(((1 / T_0) - (1 / T))*M_L*h_wg / R);
		//if (T > T_0+2)
			//dPsat_dT = 0.0;
		return dPsat_dT;
	}
	virtual double Deriv_dPgw_dT(double PG, double PC, double T)
	{
		double T_0 = 373.15;
		double P_0 = 101325;
		double h_wg = 2258000.0;
		double dPgw_dT(0.0);
		double C_w = M_L / R / T;
		double dPsat_dT(0.0);
		double P_sat(0);
		P_sat = get_P_sat_S(T);
		dPsat_dT = Deriv_dPsat_dT(PG,T);
		dPgw_dT = dPsat_dT*exp(-PC*C_w / rho_l_std) + P_sat*exp(-PC*C_w / rho_l_std)*(PC*M_L / rho_l_std / R / T / T);
		//dPgw_dT = exp((-PC*M_L / (rho_l_std* R*T)) + (((1 / T_0) - (1 / T))*M_L*h_wg / R))*P_0*((h_wg*M_L / (R*T*T)) + (M_L*PC / (R*rho_l_std*T*T)));
		//if (get_P_G_w_ORG(PC, T) > PG)
			//dPgw_dT = 0.0;
		return dPgw_dT;
	}

	virtual double Deriv_dPgw_dPC(double PG, double PC, double T)
	{
		double T_0 = 373.15;
		double P_0 = 101325.0;
		double h_wg = 2258000.0;
		double dPgw_dPC(0.0);
		double C_w = M_L / R / T;
		double P_sat(0.0);
		P_sat = get_P_sat_S(T);
		dPgw_dPC = P_sat*exp(-PC*C_w / rho_l_std)*(-C_w / rho_l_std);
		//if (get_P_G_w_ORG(PC, T) > PG)
			//dPgw_dPC = 0.0;
		return dPgw_dPC;

	}

	/**
	* OVERALL HEAT conductivity
	*/
	virtual double get_overall_Heat_Conductivity(double Sg, double lambda_pm_dry = 0.582, double lambda_pm_wet = 1.14)
	{
		double lambda_pm(0.0);
		lambda_pm = lambda_pm_dry + std::pow(1 - Sg, 0.5)*(lambda_pm_wet - lambda_pm_dry);
		if (Sg < 0)
			lambda_pm = lambda_pm_wet;
		else if (Sg > 1)
			lambda_pm = lambda_pm_dry;
		return lambda_pm;

	}

	/**
	* Calc henry const
	*/
	virtual double Henry_const(double T)//unit [Pa^(-1)]
	{
		return (0.8942 + 1.47*exp(-0.04394*T))*(1e-10);
		//return (0.8942 + 1.47 * exp(-0.04394*T))*1e-10;
	}
	/**
	* Calc derivatives
	*/
	virtual double deriv_henry_T(double T)
	{
		return (-1.47 * exp(-0.04394*T)*(1e-10)*0.04394);
	}
private:

	/**
	  * parameter for calculate the molar density of gas phase
      */
	const double R = 8.314;//J/mol/kg
	
	//const double P_vapor = 0.0;
	const double eps = 1e-12;
	const double rho_l_std = 1000;// 1000.0;
	const double M_G = 0.02896;
	const double M_L = 0.018;
	/**
	* parameters for heatpipe capillary pressure-saturation model
	*/
	const double S_gr = 0.0;//residual saturation of the gas phase
	const double S_lr = 0.15;//residual saturation of the liquid phase

	double P_G;
	
	double X_L;
	double S_G;
	double T_L;
	double _alpha;
	double _beta;
};


