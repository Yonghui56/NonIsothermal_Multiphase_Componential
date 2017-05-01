/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file FunctionCMP_2P2C.hpp
*
* Created on 2014-05-12 by Yonghui HUANG
*/

#include "logog.hpp"

#include "MathLib/DataType.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "NumLib/Function/TXFunctionDirect.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "MathLib/ODE/RungeKutta4.h"
#include "EOS_NonIsoGlobalNCPForm.h"

template <class T1, class T2>
bool FunctionCMP_NonIsoGlobalNCPForm<T1, T2>::initialize(const BaseLib::Options &option)
{

	Ogs6FemData* femData = Ogs6FemData::getInstance();
	_msh_id = option.getOptionAsNum<size_t>("MeshID");
	size_t time_id = option.getOptionAsNum<int>("TimeGroupID");
	NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

	//mesh and FE objects
	MeshLib::IMesh* msh = femData->list_mesh[_msh_id];
	// get the number of elements 
	size_t n_ele = msh->getNumberOfElements();
	size_t n_nodes = msh->getNumberOfNodes();
	MyDiscreteSystem* dis = 0;
	dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
	_feObjects = new FemLib::LagrangeFeObjectContainer(msh);

	// set names of the output parameters
	this->setOutputParameterName(0, "GAS_PRESSURE");
	this->setOutputParameterName(1, "MOLAR_FRACTION");
	this->setOutputParameterName(2, "SATURATION");
	this->setOutputParameterName(3, "TEMPERATURE");
	this->setOutputParameterName(4, "Liquid_Pressure");
	this->setOutputParameterName(5, "Capillary_Pressure");
	this->setOutputParameterName(6, "Mass_Density_L_H");
	this->setOutputParameterName(7, "Mass_Density_G_H");
	// also secondary variables
	// TODO: set seconary variable names also as output parameters

	// create the MyCMPPressureForm problem
	_problem = new MyCMPNonIsoGlobalNCPFormProblemType(dis);
	_problem->setTimeSteppingFunction(*tim);  // applying the time stepping function
	
	// creating mean pressure vector
	MyVariableCMPNonIsoGlobalNCPForm* Gas_Pressure = _problem->addVariable("GAS_PRESSURE");
	FemVariableBuilder var_builder;
	var_builder.doit("GAS_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, Gas_Pressure);
	SolutionLib::FemIC* femIC = _problem->getVariable(0)->getIC();
	_P = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_P->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
		femIC->setup(*_P);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_P->initialize(*dis, _problem->getVariable(0)->getCurrentOrder(), 0.0);
	}

	// creating molar fraction

	MyVariableCMPNonIsoGlobalNCPForm* total_mass_density = _problem->addVariable("MOLAR_FRACTION");
	var_builder.doit("MOLAR_FRACTION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, total_mass_density);
	femIC = _problem->getVariable(1)->getIC();
	_X = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_X->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
		femIC->setup(*_X);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_X->initialize(*dis, _problem->getVariable(1)->getCurrentOrder(), 0.0);
	}

	MyVariableCMPNonIsoGlobalNCPForm* saturation = _problem->addVariable("SATURATION");
	var_builder.doit("SATURATION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, saturation);
	femIC = _problem->getVariable(2)->getIC();
	_S = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_S->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
		femIC->setup(*_S);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_S->initialize(*dis, _problem->getVariable(2)->getCurrentOrder(), 0.0);
	}

	MyVariableCMPNonIsoGlobalNCPForm* temperature = _problem->addVariable("TEMPERATURE");
	var_builder.doit("TEMPERATURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, temperature);
	femIC = _problem->getVariable(3)->getIC();
	_T = new MyNodalFunctionScalar();
	if (femIC)
	{
		// FemIC vector is not empty
		_T->initialize(*dis, _problem->getVariable(3)->getCurrentOrder(), 0.0);
		femIC->setup(*_T);
	}
	else
	{
		// FemIC vector is empty
		// initialize the vector with zeros
		_T->initialize(*dis, _problem->getVariable(3)->getCurrentOrder(), 0.0);
	}
	// initialize the local EOS problem
	
	//ogsChem::LocalVector _output1;
	//_output1 = ogsChem::LocalVector::Zero(_LP_EOS->N);

	_PC = new MyNodalFunctionScalar();	
	_PC->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	//_PC->setValue(0, 0.0); //PC_ini
	_PL = new MyNodalFunctionScalar();
	_PL->initialize(*dis, FemLib::PolynomialOrder::Linear, 1e+6);//PL_ini
	_rho_L_h = new MyNodalFunctionScalar();//MASS density of hydrogen in liquid phase
	_rho_L_h->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	_rho_G_h = new MyNodalFunctionScalar();//mass density of hydrogen in gas phase
	_rho_G_h->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	
	_dPcdSg = new MyNodalFunctionScalar();
	_dPcdSg->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
	

	/**
	* initialize the vector of tempVar
	*/
	//_mat_Jacob = LocalMatrix::Zero(4 * n_nodes, 4 * n_nodes);
	_mat_Jacob.resize(
		4 * n_nodes,
		4 * n_nodes);  // Resizing also initializes the matrix with zero
	_vec_Res = LocalVector::Zero(4 * n_nodes);

	//initialize the local M matrix and K matrix and partially Jacobian matrix
	for (size_t index = 0; index < n_ele; index++)
	{
		MathLib::LocalMatrix test1 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_M_matrix.push_back(test1);
		MathLib::LocalMatrix test2 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_K_matrix.push_back(test2);
		MathLib::LocalMatrix test3 = MathLib::LocalMatrix::Zero(9, 9);
		_elem_J_matrix.push_back(test3);
	}
	

	// linear assemblers
	MyNonLinearAssemblerType* nonlin_assembler = new MyNonLinearAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this, msh->getGeometricProperty()->getCoordinateSystem());

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearCMPNonIsoGlobalNCPFormProblemType(dis);
	_non_linear_eqs = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize(nonlin_assembler, non_linear_r_assembler, non_linear_j_assembler);
	//_non_linear_eqs->initialize(nonlin_assembler,non_linear_r_assembler, non_linear_j_assembler);
	_non_linear_problem->setTimeSteppingFunction(*tim);
	_non_linear_problem->addVariable("GAS_PRESSURE");
	_non_linear_problem->addVariable("MOLAR_FRACTION");
	_non_linear_problem->addVariable("SATURATION");
	_non_linear_problem->addVariable("TEMPERATURE");

	SolutionLib::FemIC* P_ic = new SolutionLib::FemIC(msh);
	P_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_P->getDiscreteData()));
	var_builder.doit("GAS_PRESSURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(0));
	_non_linear_problem->getVariable(0)->setIC(P_ic);

	SolutionLib::FemIC* X_ic = new SolutionLib::FemIC(msh);
	X_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_X->getDiscreteData()));
	var_builder.doit("MOLAR_FRACTION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(1));
	_non_linear_problem->getVariable(1)->setIC(X_ic);

	SolutionLib::FemIC* S_ic = new SolutionLib::FemIC(msh);
	S_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_S->getDiscreteData()));
	var_builder.doit("SATURATION", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(2));
	_non_linear_problem->getVariable(2)->setIC(S_ic);

	SolutionLib::FemIC* T_ic = new SolutionLib::FemIC(msh);
	T_ic->addDistribution(femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>(_T->getDiscreteData()));
	var_builder.doit("TEMPERATURE", option, msh, femData->geo, femData->geo_unique_name, _feObjects, _non_linear_problem->getVariable(3));
	_non_linear_problem->getVariable(3)->setIC(T_ic);

	// set up non-linear solution
	myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler);
	myNSolverFactory = new MyDiscreteNonlinearSolverFactory(myNRIterator);
	this->_non_linear_solution = new MyNonLinearSolutionType(dis, this->_non_linear_problem, this,myNSolverFactory);
	this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
	this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
	const BaseLib::Options* optNum = option.getSubGroup("Numerics");

	// linear solver
	MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
	linear_solver->setOption(*optNum);
	// set nonlinear solver options
	this->_non_linear_solution->getNonlinearSolver()->setOption(*optNum);
    // set theta----------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (!optNum->hasOption("TimeTheta"))
    {
        ERR("Time theta setting not found!!!");
        exit(1);
    }
    else
    {
        double tmp_theta(optNum->getOptionAsNum<double>("TimeTheta"));
        non_linear_j_assembler->setTheta(tmp_theta);
        non_linear_r_assembler->setTheta(tmp_theta);
        _non_linear_eqs->getLinearAssembler()->setTheta(tmp_theta);
    }
	


	// get the nonlinear solution dof manager
	this->_nl_sol_dofManager = this->_non_linear_solution->getDofEquationIdTable();
	// set up solution
	_solution = new MyCMPNonIsoGlobalNCPFormSolution(dis, _problem, _non_linear_solution, this);
	/**
	*Calculate the secondary variables on each node
	*/
	this->calc_nodal_eos_sys(0.0);
	

	// set initial output parameter
	OutputVariableInfo var1(this->getOutputParameterName(0), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2(this->getOutputParameterName(1), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3(this->getOutputParameterName(2), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var3.name, var3);
	OutputVariableInfo var4(this->getOutputParameterName(3), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _T);
	femData->outController.setOutput(var4.name, var4);

	OutputVariableInfo var_sec_1(this->getOutputParameterName(4), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_sec_1.name, var_sec_1);

	OutputVariableInfo var_sec_2(this->getOutputParameterName(5), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_sec_2.name, var_sec_2);

	OutputVariableInfo var_sec_3(this->getOutputParameterName(6), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _rho_L_h);
	femData->outController.setOutput(var_sec_3.name, var_sec_3);

	OutputVariableInfo var_sec_4(this->getOutputParameterName(7), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _rho_G_h);
	femData->outController.setOutput(var_sec_4.name, var_sec_4);

	return true;
}

template <class T1, class T2>
void FunctionCMP_NonIsoGlobalNCPForm<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	/*
	size_t i;
	//const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for (i = 0; i < _linear_problems.size(); i++) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}
	*/
}

template <class T1, class T2>
void FunctionCMP_NonIsoGlobalNCPForm<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{

}

template <class T1, class T2>
void FunctionCMP_NonIsoGlobalNCPForm<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
	
	//update data for output
	Ogs6FemData* femData = Ogs6FemData::getInstance();
	this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	// we have 2 primary variables
	OutputVariableInfo var1("GAS_PRESSURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _P);
	femData->outController.setOutput(var1.name, var1);
	OutputVariableInfo var2("MOLAR_FRACTION", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _X);
	femData->outController.setOutput(var2.name, var2);
	OutputVariableInfo var3("SATURATION", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _S);
	femData->outController.setOutput(var3.name, var3);
	OutputVariableInfo var4("TEMPERATURE", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _T);
	femData->outController.setOutput(var4.name, var4);
	// and 8 secondary variables
	// add all seconary variables as well. 
	OutputVariableInfo var_Sec_1("Liquid_Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PL);
	femData->outController.setOutput(var_Sec_1.name, var_Sec_1);

	OutputVariableInfo var_Sec_2("Capillary_Pressure", _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _PC);
	femData->outController.setOutput(var_Sec_2.name, var_Sec_2);



}


template <class T1, class T2>
void FunctionCMP_NonIsoGlobalNCPForm<T1, T2>::calc_nodal_eos_sys(double dt = 0.0)
{
	std::size_t node_id(0);
	std::size_t num_nodes(0);
	num_nodes = _P->getDiscreteData()->getRangeEnd();
	double EPS = 1E-7;

	MathLib::LocalMatrix test = MathLib::LocalMatrix::Zero(1, 1);

	MathLib::LocalVector input = MathLib::LocalVector::Zero(4);
	MathLib::LocalVector output = MathLib::LocalVector::Zero(4);

	double PC_np(0.0);
	double PG_w_np(0.0);

	double dPCdSG_np(0.0);
	double dPG_wdPC_np(0.0);
	double dPG_wdT_np(0.0);
	INFO("Calculating EOS on each node...");

	// loop for each node
	for (node_id = _P->getDiscreteData()->getRangeBegin();
		node_id < _P->getDiscreteData()->getRangeEnd();
		node_id++)
	{
		input(0) = _P->getValue(node_id);
		input(1) = _X->getValue(node_id);
		input(2) = _S->getValue(node_id);//SG gas phase
		input(3) = _T->getValue(node_id);// Temperature
		//This part is for standard newton iteration of local problem
		//solve EOS 
		//This part is for complementary condition
		_EOS->set_env_condition(input);//apply the primary variables as the input variables for Eos
		
		PC_np = _EOS->getPcbySg(input(2));
		PG_w_np = _EOS->get_P_G_w(input(0), PC_np, input(3));
		dPCdSG_np = _EOS->Deriv_dPCdS(input(2));
		dPG_wdPC_np = _EOS->Deriv_dPgw_dPC(input(0), PC_np, input(3));
		dPG_wdT_np = _EOS->Deriv_dPgw_dT(input(0), PC_np, input(3));
		//For minimum function 1
		//_vec_Res(node_id ) = std::min( input(2), C_h*(input(0) + output(2)) - rho_L_std*M_G*input(1) / M_L );//+ 2 * num_nodes
		this->_PL->setValue(node_id, input(0) - PC_np);
		this->_PC->setValue(node_id, PC_np);
		//_vec_Res(3 * num_nodes + node_id) = -std::min(input(2), 1 - (input(1) + (PG_w_np/input(0))));//X_G_air + X_G_vap
        auto test0 = 1 - (input(1) + (PG_w_np / input(0)));
        auto test1= -std::min(input(2), 1 - (input(1) + (PG_w_np / input(0))));
        auto test2= -(2* input(2) - (input(1) + (PG_w_np / input(0))))
            + std::min({ input(2), 1 - (input(1) + (PG_w_np / input(0))), input(2)-1 })
            + std::max({ input(2), 1 - (input(1) + (PG_w_np / input(0))), input(2)-1 });
        auto test3 = -std::max(std::min(input(2), 1 - (input(1) + (PG_w_np / input(0)))),
            std::min(std::max(input(2), 1 - (input(1) + (PG_w_np / input(0)))), input(2)-1));
        _vec_Res(3 * num_nodes + node_id)= -(2 * input(2) - (input(1) + (PG_w_np / input(0))))
            + std::min({ input(2), 1 - (input(1) + (PG_w_np / input(0))), input(2) - 1 })
            + std::max({ input(2), 1 - (input(1) + (PG_w_np / input(0))), input(2) - 1 });
        if (input(2) < input(2) - 1) {
        if(test0<input(2) - 1 && test0>input(2)){
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id) = PG_w_np / pow(input(0), 2);
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 1) = -1.0;
            _mat_Jacob.coeffRef(3 * num_nodes + node_id, 4 * node_id + 2) = -dPG_wdPC_np*dPCdSG_np / input(0);
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 3) = -dPG_wdT_np / input(0);
        }
        else if (test0 > input(2) - 1) {
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id) = 0.0;
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 1) = 0.0;
            _mat_Jacob.coeffRef(3 * num_nodes + node_id, 4 * node_id + 2) = 1.0;
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 3) = 0.0;
        }
        else {
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id) = 0.0;
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 1) = 0.0;
            _mat_Jacob.coeffRef(3 * num_nodes + node_id, 4 * node_id + 2) = 1.0;
            _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 3) = 0.0;
        }
        }
        else {
            if (test0<input(2)&& test0>input(2) - 1) {
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id) = PG_w_np / pow(input(0), 2);
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 1) = -1.0;
                _mat_Jacob.coeffRef(3 * num_nodes + node_id, 4 * node_id + 2) = -dPG_wdPC_np*dPCdSG_np / input(0);
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 3) = -dPG_wdT_np / input(0);
            }
            else if (test0 > input(2)) {
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id) = 0.0;
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 1) = 0.0;
                _mat_Jacob.coeffRef(3 * num_nodes + node_id, 4 * node_id + 2) = 1.0;
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 3) = 0.0;
            }
            else {
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id) = 0.0;
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 1) = 0.0;
                _mat_Jacob.coeffRef(3 * num_nodes + node_id, 4 * node_id + 2) = 1.0;
                _mat_Jacob.coeffRef(node_id + 3 * num_nodes, 4 * node_id + 3) = 0.0;
            }
        }
		
	}
	

		
	

}


