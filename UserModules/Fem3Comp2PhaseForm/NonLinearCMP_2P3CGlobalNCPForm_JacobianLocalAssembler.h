/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NonLinearCMP_2P2C_JacobianLocalAssembler.h
 *
 * 2015-03-21 by Yonghui Huang
 */

/**
  * This file is same as the MassTransportTimeODELocalAssembler.h
  * The difference is, the compound molecular diffusion coefficient is disabled, 
  */

#ifndef NON_LINEAR_CMP_2P3CGLOBALNCPFORM_JACOBIAN_LOCAL_ASSEMBLER_H
#define NON_LINEAR_CMP_2P3CGLOBALNCPFORM_JACOBIAN_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"
#include "Ogs6FemData.h"
#include "NonLinearCMP_2P3CGlobalNCPForm_TimeODELocalAssembler.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA >//,
class NonLinearCMP_2P3CGlobalNCPForm_JacobianLocalAssembler : public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	typedef MathLib::LocalVector LocalVectorType;
	typedef MathLib::LocalMatrix LocalMatrixType;
	NonLinearCMP_2P3CGlobalNCPForm_JacobianLocalAssembler(FemLib::LagrangeFeObjectContainer* feObjects, T_FUNCTION_DATA* func_data, const MeshLib::CoordinateSystem &problem_coordinates)
		: _feObjects(*feObjects), _function_data(func_data), _problem_coordinates(problem_coordinates)
    {
		_EOS = new EOS_2P3CGlobalNCPForm();
	};

	virtual ~NonLinearCMP_2P3CGlobalNCPForm_JacobianLocalAssembler()
    {
		_function_data = NULL; 
		BaseLib::releaseObject(_EOS);
    };

	
	T_FUNCTION_DATA* get_function_data(void)
	{
		return _function_data;
	}
	
	void setTheta(double v)
	{
		assert(v >= .0 && v <= 1.0);
		_Theta = v;
	}
	/*
	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
    {
		size_t el(0);
		const size_t n_nodes = e.getNumberOfNodes();
		const size_t mat_id = e.getGroupID();
		const std::size_t elem_id = e.getID();
		//clear the local Jacobian Matrix

		MathLib::LocalMatrix TMP_M(3*n_nodes,3* n_nodes);
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		double dt = time.getTimeStepSize();
		LocalMatrixType _M = LocalMatrixType::Zero(3 * n_nodes, 3 * n_nodes);
		LocalMatrixType _K = LocalMatrixType::Zero(3 * n_nodes, 3 * n_nodes);
		LocalMatrixType _J = LocalMatrixType::Zero(3 * n_nodes, 3 * n_nodes);
		_M = _function_data->get_elem_M_matrix().at(elem_id);
		_K = _function_data->get_elem_K_matrix().at(elem_id);
		//_J = _function_data->get_elem_J_matrix().at(elem_id);
		
		double eps = 1e-9;// can also be defined as numerical way
		TMP_M = (1.0 / dt)*_M + _Theta*_K;
		//std::cout << TMP_M;
		//LocalMatrixType Local_LHS=(1/dt)*M+
		//LocalVectorType e_vec = LocalVectorType::Zero(2 * n_nodes);
		//_local_assembler=assembleODE(time, e, u1, u0, M, K);
		for (size_t u_idx = 0; u_idx < u1.size(); u_idx++)
		{
			//clear the epsilon vectoe 
			LocalVectorType e_vec = LocalVectorType::Zero(3 * n_nodes);
			e_vec(u_idx) = eps*(1 + std::abs(u1(u_idx)));
			localJ.block(0, u_idx, u1.size(), 1) = (TMP_M*(u1 - e_vec) - TMP_M*(u1 + e_vec)) / 2 / eps / (1 + std::abs(u1(u_idx)));
			//debugging--------------------------
			//std::cout << "localJ: \n";
			//std::cout << localJ << std::endl;
			//end of debugging-------------------
			//localJ=
		}

		//localJ.noalias() -= _J;
    }  // end of function assembly
	*/
	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
	{
		size_t el(0);
		const size_t n_nodes = e.getNumberOfNodes();
		const size_t mat_id = e.getGroupID();
		const std::size_t elem_id = e.getID();
		const size_t n_dof = u1.size();
		double dt = time.getTimeStepSize();
		TMP_M_1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
		TMP_M_2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
		//FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		//MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

		double eps = 1e-7;

		for (size_t u_idx = 0; u_idx < n_dof; u_idx++)
		{
			//clear M1,K1,F1,M2,K2,F2
			M1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			K1 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			F1 = MathLib::LocalVector::Zero(n_dof);
			M2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			K2 = MathLib::LocalMatrix::Zero(n_dof, n_dof);
			F2 = MathLib::LocalVector::Zero(n_dof);
			//clear local_r
			local_r_1 = MathLib::LocalVector::Zero(n_dof);
			local_r_2 = MathLib::LocalVector::Zero(n_dof);
			//clear the epsilon vectoe 
			MathLib::LocalVector e_vec = MathLib::LocalVector::Zero(n_dof);
			e_vec(u_idx) = eps*(1 + std::abs(u1(u_idx)));
			assemble_jac_ODE(time, e, u1 - e_vec, u0, M1, K1, F1);
			TMP_M_1 = (1.0 / dt)*M1 + _Theta*K1;
			local_r_1 = TMP_M_1*(u1 - e_vec);
			TMP_M_1 = (1.0 / dt) * M1 - (1. - _Theta) * K1;
			local_r_1.noalias() -= TMP_M_1 * u0;
			local_r_1.noalias() -= F1;
			assemble_jac_ODE(time, e, u1 + e_vec, u0, M2, K2, F2);
			TMP_M_2 = (1.0 / dt)*M2 + _Theta*K2;
			local_r_2 = TMP_M_2*(u1 + e_vec);
			TMP_M_2 = (1.0 / dt) * M2 - (1. - _Theta) * K2;
			local_r_2.noalias() -= TMP_M_2 * u0;
			local_r_2.noalias() -= F2;
			//localJ.block(0, u_idx, u1.size(), 1) += (local_r_1 - local_r_2) / 2 / eps / (1 + std::abs(u1(u_idx)));
			localJ.col(u_idx) += (local_r_1 - local_r_2) / 2 / eps / (1 + std::abs(u1(u_idx)));
			//debugging--------------------------
			//std::cout << "localJ: \n";
			//std::cout << localJ << std::endl;
			//end of debugging-------------------
			//localJ=
		}
	}
protected:
	virtual void assemble_jac_ODE(const NumLib::TimeStep & /*time*/, const MeshLib::IElement &e, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix  & localM, MathLib::LocalMatrix & localK, MathLib::LocalMatrix & localF)
	{
		// -------------------------------------------------
		// current element
		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
		// integration method 
		FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
		// now get the sizes
		const std::size_t n_dim = e.getDimension();
		const std::size_t mat_id = e.getGroupID();
		const std::size_t ele_id = e.getID();
		const std::size_t n_nodes = e.getNumberOfNodes(); // number of connecting nodes
		const std::size_t n_gsp = q->getNumberOfSamplingPoints(); // number of Gauss points
		std::size_t node_id(0);  // index of the node

		// get the instance of porous media class
		MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
		// get the fluid property for gas phase
		MaterialLib::Fluid* fluid1 = Ogs6FemData::getInstance()->list_fluid[0];
		// get the fluid property for liquid phase
		MaterialLib::Fluid* fluid2 = Ogs6FemData::getInstance()->list_fluid[1];
		// get the first (light) component - hydrogen
		MaterialLib::Compound* component1 = Ogs6FemData::getInstance()->list_compound[0];
		// get the second (heavy) component - water
		MaterialLib::Compound* component2 = Ogs6FemData::getInstance()->list_compound[1];

		const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());
		double geo_area = 1.0;
		pm->geo_area->eval(e_pos, geo_area);
		const bool hasGravityEffect = _problem_coordinates.hasY();//detect the gravity term based on the mesh structure

		if (hasGravityEffect) {
			vec_g = MathLib::LocalVector::Zero(_problem_coordinates.getDimension());
			vec_g[_problem_coordinates.getIndexOfY()] = 9.81;
		}
		/*
		*SECONDARY VARIABLES initialize
		*/

		
		PC_gp = MathLib::LocalMatrix::Zero(1, 1);
		PL_gp = MathLib::LocalVector::Zero(1, 1);
		dPC_dSg_gp = MathLib::LocalMatrix::Zero(1, 1);
		
		X_L_h_gp = MathLib::LocalMatrix::Zero(1, 1);
		X_L_w_gp = MathLib::LocalMatrix::Zero(1, 1);
		X_L_c_gp = MathLib::LocalMatrix::Zero(1, 1);

		rho_mol_L_gp = MathLib::LocalVector::Zero(1, 1);
		rho_mol_G_gp = MathLib::LocalVector::Zero(1, 1);
		rho_mass_G_gp = MathLib::LocalVector::Zero(1, 1);
		rho_mass_L_gp = MathLib::LocalVector::Zero(1, 1);

		P_sat_gp = MathLib::LocalVector::Zero(1, 1);
		// Loop over each node on this element e
		// calculate the secondary variables' values on each node 
		// calculate the derivative of secondary variables based on P and X on each node
		// store the value in the vector respectively

		M = MathLib::LocalMatrix::Zero(3, 5);//method 2
		D = MathLib::LocalMatrix::Zero(3, 5);//method 2
		H = MathLib::LocalVector::Zero(3);//for gravity term
		tmp = MathLib::LocalMatrix::Zero(1, 1);
		Input = MathLib::LocalVector::Zero(3);
		localMass_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); // for method2
		localDispersion_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); //for method2
		localDispersion = MathLib::LocalMatrix::Zero(n_nodes, n_nodes);
		LocalJacb = MathLib::LocalMatrix::Zero(3 * n_nodes, 3 * n_nodes);
		localGravity_tmp = MathLib::LocalVector::Zero(n_nodes);//tmp vect for gravity 
		localRes = MathLib::LocalVector::Zero(n_nodes);//tmp vect for gravity 
		/*
		* Loop for each Gauss Point
		*/
		Kr_L_gp = 0.0;
		Kr_G_gp = 0.0; lambda_L = 0.0; lambda_G = 0.0;
		Lambda_h = 0.0;
		//-------------------Begin assembly on each gauss point---------------------------
		for (j = 0; j < n_gsp; j++)
		{
			// calc on each gauss point
			// get the sampling point
			q->getSamplingPoint(j, gp_x);
			// compute basis functions
			fe->computeBasisFunctions(gp_x);
			// get the coordinates of current location 
			fe->getRealCoordinates(real_x);
			// transfer to Gauss point position
			NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);
			double fac = geo_area*fe->getDetJ()*q->getWeight(j);
			// get the porosity
			pm->porosity->eval(gp_pos, poro);
			// get the intrinsic permeability
			pm->permeability->eval(gp_pos, K_perm);
			// evaluate viscosity of each fluid phase
			fluid1->dynamic_viscosity->eval(gp_pos, mu_G);
			fluid2->dynamic_viscosity->eval(gp_pos, mu_L);

			//fluid1->molar_mass->eval(gp_pos, M_G);
			//fluid2->molar_mass->eval(gp_pos, M_L);

			fluid1->density->eval(gp_pos, rho_G_std);
			fluid2->density->eval(gp_pos, rho_L_std);

			// evaluate componential diffusion coefficients
			component1->molecular_diffusion->eval(gp_pos, D_G);
			component2->molecular_diffusion->eval(gp_pos, D_L);

			// evaluation of the shape function
			MathLib::LocalMatrix &Np = *fe->getBasisFunction();
			MathLib::LocalMatrix &dNp = *fe->getGradBasisFunction();


			P_gp = Np*u1.head(n_nodes);
			X1_gp = Np*u1.block(n_nodes, 0, n_nodes, 1);//H2
			X2_gp = Np*u1.block(2 * n_nodes, 0, n_nodes, 1);//CO2
			X3_gp = Np*u1.block(3 * n_nodes, 0, n_nodes, 1);//WATER VAPOR
			S_gp = Np*u1.tail(n_nodes);

			P_sat_gp(0, 0) = _EOS->get_P_sat(T_0);
			X_L_w_gp(0, 0) = P_gp(0, 0)*X3_gp(0, 0) / P_sat_gp(0, 0);
			X_L_h_gp(0, 0) = P_gp(0, 0)*X1_gp(0, 0) / Hen_L_h;
			X_L_c_gp(0, 0) = P_gp(0, 0)*X2_gp(0, 0) / Hen_L_c;
			PC_gp(0, 0) = pm->getPc_bySat(S_gp(0, 0));
			PL_gp(0, 0) = P_gp(0, 0) - PC_gp(0, 0);
			dPC_dSg_gp(0, 0) = pm->Deriv_dPCdS(S_gp(0, 0));
			rho_mol_G_gp(0, 0) = P_gp(0, 0) / R / T_0;
			rho_mol_L_gp(0, 0) = rho_L_std / M_L / X_L_w_gp(0, 0);

			rho_mass_G_gp(0, 0) = rho_mol_G_gp(0, 0)*(X1_gp(0, 0)*M_H + X2_gp(0, 0)*M_C + X3_gp(0, 0)*M_L);
			rho_mass_L_gp(0, 0) = rho_mol_L_gp(0, 0)*(X_L_h_gp(0, 0)*M_H + X_L_c_gp(0, 0)*M_C + X_L_w_gp(0, 0)*M_L);

			isinf = std::isfinite(dPC_dSg_gp(0, 0));
			if (isinf == 0)
			{
				dPC_dSg_gp(0, 0) = 0.0;
			}
			else
			{
				dPC_dSg_gp(0, 0) = dPC_dSg_gp(0, 0);
			}



			M(0, 0) = poro*(S_gp(0, 0)*X1_gp(0, 0) / R / T_0);//dPG
			M(0, 1) = poro*(rho_mol_G_gp(0, 0)*S_gp(0, 0) + (1 - S_gp(0, 0))*rho_l_std*P_sat_gp(0, 0) / M_L / Hen_L_h / X3_gp(0, 0));//dX1 h2
			M(0, 2) = 0.0;//Dx2 co2
			M(0, 3) = -poro*((1 - S_gp(0, 0))*rho_l_std*P_sat_gp(0, 0)*X1_gp(0, 0) / M_L / Hen_L_h / pow(X3_gp(0, 0), 2));//dX2
			M(0, 4) = poro*(rho_mol_G_gp(0, 0)*X1_gp(0, 0) - rho_mol_L_gp(0, 0)*X_L_h_gp(0, 0));

			M(1, 0) = poro*(S_gp(0, 0)*X2_gp(0, 0) / R / T_0);//dPG
			M(1, 1) = 0.0;
			M(1, 2) = poro*(rho_mol_G_gp(0, 0)*S_gp(0, 0) + (1 - S_gp(0, 0))*rho_l_std*P_sat_gp(0, 0) / M_L / Hen_L_c / X3_gp(0, 0));//dX1 h2
			M(1, 3) = -poro*((1 - S_gp(0, 0))*rho_l_std*P_sat_gp(0, 0)*X2_gp(0, 0) / M_L / Hen_L_c / pow(X3_gp(0, 0), 2));//dX2
			M(1, 4) = poro*(rho_mol_G_gp(0, 0)*X2_gp(0, 0) - rho_mol_L_gp(0, 0)*X_L_c_gp(0, 0));

			M(2, 0) = poro*(S_gp(0, 0)*X3_gp(0, 0) / R / T_0);//dPG
			M(2, 1) = 0.0;
			M(2, 2) = 0.0;
			M(2, 3) = poro*rho_mol_G_gp(0, 0)*S_gp(0, 0);
			M(2, 4) = poro*(rho_mol_G_gp(0, 0)*X3_gp(0, 0) - rho_l_std / M_L);
			//-------------debugging------------------------
			//std::cout << "M=" << std::endl;
			//std::cout << M << std::endl;
			//--------------end debugging-------------------

			//assembly the mass matrix
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 5; jj++){
					tmp(0, 0) = M(ii, jj);
					localMass_tmp.setZero();
					fe->integrateWxN(j, tmp, localMass_tmp);
					localM.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localMass_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localM=" << std::endl;
			//std::cout << localM << std::endl;
			//-------------End debugging------------------------
			//calculate the relative permeability 
			Kr_L_gp = pm->getKr_L_bySg(S_gp(0, 0)); // change the Relative permeability calculation into porous media file			
			Kr_G_gp = pm->getKr_g_bySg(S_gp(0, 0));
			lambda_L = K_perm*Kr_L_gp / mu_L;
			lambda_G = K_perm*Kr_G_gp / mu_G;


			//Calc each entry of the Laplace Matrix
			D(0, 0) = lambda_G*rho_mol_G_gp(0, 0)*X1_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X_L_h_gp(0, 0)
				+ poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*X1_gp(0, 0) / Hen_L_h;
			D(0, 1) = poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*P_gp(0, 0) / Hen_L_h
				+ poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);
			D(0, 2) = 0.0;
			D(0, 3) = 0.0;
			D(0, 4) = -lambda_L*rho_mol_L_gp(0, 0)*X_L_h_gp(0, 0)*dPC_dSg_gp(0, 0);

			D(1, 0) = lambda_G*rho_mol_G_gp(0, 0)*X2_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X_L_c_gp(0, 0)
				+ poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0) / Hen_L_c;
			D(1, 1) = 0.0;
			D(1, 2) = poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*P_gp(0, 0) / Hen_L_c
				+ poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);
			D(1, 3) = 0.0;
			D(1, 4) = -lambda_L*rho_mol_L_gp(0, 0)*X_L_c_gp(0, 0)*dPC_dSg_gp(0, 0);

			D(2, 0) = lambda_G*rho_mol_G_gp(0, 0)*X2_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X_L_w_gp(0, 0)
				+ poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*X2_gp(0, 0) / P_sat_gp(0, 0);
			D(2, 1) = 0.0;
			D(2, 2) = 0.0;
			D(2, 3) = poro*(1 - S_gp(0, 0))*D_L*rho_mol_L_gp(0, 0)*P_gp(0, 0) / P_sat_gp(0, 0)
				+ poro*D_G*S_gp(0, 0)*rho_mol_G_gp(0, 0);
			D(2, 4) = -lambda_L*rho_mol_L_gp(0, 0)*X_L_w_gp(0, 0)*dPC_dSg_gp(0, 0);

			//-------------debugging------------------------
			//std::cout << "D=" << std::endl;
			//std::cout << D << std::endl;
			//--------------end debugging-------------------
			//
			for (ii = 0; ii < 3; ii++){
				for (jj = 0; jj < 5; jj++){
					tmp(0, 0) = D(ii, jj);
					localDispersion_tmp.setZero();
					fe->integrateDWxDN(j, tmp, localDispersion_tmp);
					localK.block(n_nodes*ii, n_nodes*jj, n_nodes, n_nodes) += localDispersion_tmp;
				}
			}
			//-------------debugging------------------------
			//std::cout << "localK=" << std::endl;
			//std::cout << localK << std::endl;
			//--------------end debugging-------------------
			//----------------assembly the gravity term--------------------------------
			H(0) = lambda_G*rho_mol_G_gp(0, 0)*X1_gp(0, 0)*rho_mass_G_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X_L_h_gp(0, 0)*rho_mass_L_gp(0, 0);
			H(1) = lambda_G*rho_mol_G_gp(0, 0)*X2_gp(0, 0)*rho_mass_G_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X_L_c_gp(0, 0)*rho_mass_L_gp(0, 0);
			H(2) = lambda_G*rho_mol_G_gp(0, 0)*X3_gp(0, 0)*rho_mass_G_gp(0, 0) + lambda_L*rho_mol_L_gp(0, 0)*X_L_w_gp(0, 0)*rho_mass_L_gp(0, 0);
			//-------------debugging------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << H << std::endl;
			//--------------end debugging-------------------
			if (hasGravityEffect) {
				// F += dNp^T * H* gz

				for (int idx = 0; idx < 3; idx++){
					tmp(0, 0) = H(idx);
					localGravity_tmp.setZero();
					//fe->integrateDWxvec_g(j, tmp, localGravity_tmp, vec_g);
					localGravity_tmp = fac * dNp.transpose()*tmp(0, 0)*vec_g;
					localF.block(n_nodes*idx, 0, n_nodes, 1) += localGravity_tmp;

				}
			}
			//----------------end assembly--------------------------------
			//std::cout << "H=" << std::endl;
			//std::cout << localF << std::endl;

		}
		//----------------end assemby gaus---------------------------------
		//----------------calc jacobian for complementary contrains--------------
		/*
		for (i = 0; i < n_nodes; i++)
		{
		node_id = e.getNodeID(i);
		localF[i + 2 * n_nodes] = -std::min(u1[i+2*n_nodes],C_h*(u1[i] + PC[i]) - u1[i + n_nodes] * M_G*rho_L_std / M_L);
		if (u1[i + 2 * n_nodes] <=C_h*(u1[i] + PC[i]) - u1[i + n_nodes] * M_G*rho_l_std / M_L){//then rho_L^h=C_h*PG
		//Calc each entry of the mass matrix
		LocalJacb(2 * n_nodes+i, 2 * n_nodes + i) = 1;
		}
		else if (u1[i + 2 * n_nodes] > C_h*(u1[i] + PC[i]) - u1[i + n_nodes] * M_G*rho_l_std / M_L){//then rho_L^h=M^h*rho_L^std*X_L^h/M^w
		LocalJacb(i + 2 * n_nodes, i) = C_h;
		LocalJacb(i + 2 * n_nodes, n_nodes + i) = -M_G*rho_L_std / M_L;
		LocalJacb(i + 2 * n_nodes, 2 * n_nodes + i) = C_h*dPC_dSg[i];
		}
		}
		*/
		//std::cout << "H=" << std::endl;
		//std::cout << LocalJacb<< std::endl;
		//_function_data->get_elem_M_matrix()[ele_id] = localM;
		//_function_data->get_elem_K_matrix()[ele_id] = localK;

		//_function_data->get_elem_J_matrix()[ele_id] = LocalJacb;
	}
private:
    /**
      * FEM object
      */ 
    FemLib::LagrangeFeObjectContainer _feObjects;
	MeshLib::CoordinateSystem _problem_coordinates;
	/**
	  * pointer to the function data class
	  */
	T_FUNCTION_DATA* _function_data; 
	//LocalMatrixType M = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
	//LocalMatrixType K = LocalMatrixType::Zero(2 * n_nodes, 2 * n_nodes);
	double _Theta;
	std::size_t i, j, ii, jj;

	double real_x[3], gp_x[3];
	//LocalMatrixType* _M;
	//LocalMatrixType* _K;
	MathLib::LocalVector local_r_1;
	MathLib::LocalVector local_r_2;
	MathLib::LocalMatrix M1, M2, K1, K2, F1, F2, TMP_M_1, TMP_M_2;

	MathLib::LocalVector vec_g; // gravity term 
	MathLib::LocalVector H;//this is for local gravity vector 
	MathLib::LocalVector localGravity_tmp;
	MathLib::LocalVector localRes;

	MathLib::LocalVector PC;
	MathLib::LocalVector rho_L_h;
	MathLib::LocalVector rho_G_h;
	MathLib::LocalVector dPC_dSg;

	MathLib::LocalVector Input;


	MathLib::LocalMatrix rho_mol_L_gp;
	MathLib::LocalMatrix rho_mol_G_gp;
	MathLib::LocalMatrix PC_gp;
	MathLib::LocalMatrix dPC_dSg_gp;
	//LocalMatrixType PG_gp;
	MathLib::LocalVector PL_gp;
	MathLib::LocalVector P_sat_gp;

	//primary variable on gauss point
	MathLib::LocalMatrix P_gp;
	MathLib::LocalMatrix X1_gp;
	MathLib::LocalMatrix X2_gp;
	MathLib::LocalMatrix X3_gp;
	MathLib::LocalMatrix S_gp;


	MathLib::LocalMatrix M;//Local Mass Matrix
	MathLib::LocalMatrix D;//Local Laplace Matrix
	MathLib::LocalMatrix tmp;//method 2
	MathLib::LocalMatrix X_L_w_gp;
	MathLib::LocalMatrix X_L_h_gp;
	MathLib::LocalMatrix X_L_c_gp;

	MathLib::LocalMatrix localMass_tmp; //for method2
	MathLib::LocalMatrix localDispersion_tmp; //for method2
	MathLib::LocalMatrix localDispersion;
	MathLib::LocalMatrix matN;
	MathLib::LocalMatrix LocalJacb;

	MathLib::LocalVector rho_mass_L_gp;
	MathLib::LocalVector rho_mass_G_gp;

	double Kr_L_gp, Kr_G_gp;
	double lambda_L, lambda_G;
	double poro, K_perm, mu_L, mu_G, D_G, D_L;
	double rho_G_std;
	double rho_L_std;
	double Lambda_h;
	double Var_a;
	const double M_H = 0.002;
	const double M_L = 0.01;
	const double M_C = 0.016;

	double RHO_L;// RHO_L=rho_L_std+rho_L_h
	double RHO_G;// RHO_G=rho_G_h+rho_G_w

	const double R = 8.314;
	const double Hen_L_h = 7.26e+9; //Henry constant in [Pa]
	const double Hen_L_c =4.13e+9; //Henry constant in [Pa]
	const double T_0 = 303.15; //in [K]
	const double rho_l_std = 1000.0;
	EOS_2P3CGlobalNCPForm * _EOS;
	std::size_t isinf;
};

#endif  // end of ifndef
