/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file Nested_EOS_NRIterationStepInitializer.h
*
* This class is used in Newton-Raphson iterations to include the
* calcualation of local ODE problems
*
* after the file DiscreteNonlinearSolverFactory.h
*
* Created on 2015-03-21 by Yonghui Huang
*/

#ifndef NONISO_CapPressureForm_NESTED_EOS_NR_ITERATION_STEP_INITIALIZER_H
#define NONISO_CapPressureForm_NESTED_EOS_NR_ITERATION_STEP_INITIALIZER_H

#include "NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler.h"
#include "NonLinearCMP_NonIso_CapPressureForm_TimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonIso_CapPressureForm_Nested_EOS_NRIterationStepInitializer
{
public:
	NonIso_CapPressureForm_Nested_EOS_NRIterationStepInitializer(
		NonLinearCMP_NonIso_CapPressureForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, T_NODAL_FUNCTION_SCALAR, T_FUNCTION_DATA>* assemblerR,
		NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler<T_NODAL_FUNCTION_SCALAR, T_FUNCTION_DATA> *assemblerJ)
		: _assemblerR(assemblerR), _assemblerJ(assemblerJ) {};

	template<class T_X, class F_RESIDUALS, class F_DX>
	void pre_process(const T_X &/*dx*/, const T_X &/*x_new*/, F_RESIDUALS & f_residuals, F_DX &/*f_dx*/)
	{
		
	};

	template<class T_X, class F_RESIDUALS, class F_DX>
	void post_process(const T_X &/*dx*/, const T_X & x_new, F_RESIDUALS &/*f_residuals*/, F_DX &/*f_dx*/)
	{
        for (size_t num_p = 0; num_p < x_new.size()-1; num_p+=3)
		{

			_assemblerJ->get_function_data()->set_PG_node_values(num_p/3, x_new[num_p]);
			_assemblerJ->get_function_data()->set_PC_node_values(num_p/3, x_new[num_p+1]);
			_assemblerJ->get_function_data()->set_T_node_values(num_p / 3, x_new[num_p + 2]);
			
		}
		_assemblerJ->get_function_data()->calc_nodal_eos_sys(0.0);
		
	};

private:
	/**
	* pointer to assemblerR
	*/
	NonLinearCMP_NonIso_CapPressureForm_TimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, T_NODAL_FUNCTION_SCALAR, T_FUNCTION_DATA>* _assemblerR;

	/**
	* pointer to assemblerJ
	*/
	NonLinearCMP_NonIso_CapPressureForm_JacobianLocalAssembler<T_NODAL_FUNCTION_SCALAR, T_FUNCTION_DATA>* _assemblerJ;
	
};

#endif  // end of ifndef