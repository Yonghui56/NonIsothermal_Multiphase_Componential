/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractEOS_TotalDensityForm.h
 *
 * Created on 2015-04-28 by Yonghui Huang
 */
 
#ifndef ABSTRACTEOS_NONISO_PSFORM_H 
#define ABSTRACTEOS_NONISO_PSFORM_H 

#include "ChemLib/chemconst.h"

class AbstractEOS_NonIso_PSForm {
public:
	/**
	  * constructor, input is the number of unknowns
	  */
	AbstractEOS_NonIso_PSForm(std::size_t n_unknowns)
		:N(n_unknowns)
	{};

	/**
	  * destructor will be overriden depending on the real EOS class. 
	  */
	virtual ~AbstractEOS_NonIso_PSForm() {};

	/**
	  * set the environmental condition of P, T, X ... etc.
	  * these values are stored in the env_condition vector. 
	  */
	virtual void set_env_condition(ogsChem::LocalVector & env_conditon) = 0;

	//virtual void eval(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & res) = 0;
	//virtual void calc_Jacobian_fd(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J) = 0;
	/**
	  * calculate the Jacobian matrix based on the values of unknowns given
	  */
	//virtual void calc_Jacobian(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalMatrix & J) = 0;

	/**
	* calculate REST SECONDARY VARIABLES BASED ON THE SATURATION AND p_l
	*/
	//virtual void calc_Rest_SecVar(ogsChem::LocalVector & vec_unknowns, ogsChem::LocalVector & vec_rest_var) = 0;
	/**
	  * get the number of unknowns and number of equations
	  */
	std::size_t get_n() { return N; };

private: 

	/**
	  * number of unknowns and governing equtions
	  */
	const std::size_t N; 


}; 

#endif