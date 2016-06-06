/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LocalProblem__EOS->cpp
 *
 * Created on 2014-05-14 by Yonghui HUANG & Haibing Shao
 */


#include "logog.hpp"
#include <algorithm>
#include "math.h"
#include "LocalProblem_EOS_NonIso_LocalNCP.h"
#include "EOS_NonIso_HeatPipeProb.h"

/**
  * constructor of the class
  */
LocalProblem_EOS_NonIso_LocalNCP::LocalProblem_EOS_NonIso_LocalNCP()
{
	// initialize the corresponding EOS
	// _EOS = new EOS_H2_H2O_ISOT();
	_EOS = new EOS_NonIso_HeatPipeProb();

	// get the degree of freedom
	N = _EOS->get_n(); 
	Res = ogsChem::LocalVector::Ones(N); 
	Matrix_Jacobian = ogsChem::LocalMatrix::Zero(N, N);
	/*test */
	
}

/**
  * destructor of the class
  */
LocalProblem_EOS_NonIso_LocalNCP::~LocalProblem_EOS_NonIso_LocalNCP(void)
{
	BaseLib::releaseObject(_EOS);
}

/**
 * solve the EOS
 * In our case, we have 2 inputs, P and X
 * and 8 outputs, which are the secondary variables.
 * the initial guess will be embedded into the Output vector.
 */
void LocalProblem_EOS_NonIso_LocalNCP::solve(ogsChem::LocalVector & Input, ogsChem::LocalVector & Output)
{
	if (Input.size() == 3 && Output.size() == 4)
	{		
		std::size_t m_flag = 1; 
		//This part is for semi-smooth newton iteration of local problem with complementary condition
		_EOS->set_env_condition(Input);
		U_ini = Output.head(4);
		
		this->solve_LocalProblem_Newton_LineSearch(m_flag);
		if (m_flag == 0)
		{
			Output.head(4) = U_cur;
			//std::cout << U_cur << std::endl;
		}
		else
		{
			WARN("Solving local EOS problem does not converge! \n Using old values as seoncdary varibales. \n"); 
			U_ini << 0.0, 0.0, 0.0,0.0;
			_EOS->set_env_condition(Input);
			this->solve_LocalProblem_Newton_LineSearch(m_flag);
			if (m_flag == 0)
			{
				Output.head(4) = U_cur;
			}
			//std::cout << Input << std::endl;
			//std::cout << U_ini << std::endl;
			//std::cout << U_cur << std::endl;
			//exit(1);
		}
	}
	else
	{
		ERR("When solving local EOS problem, the size of input and out vector is not correct! ");
		exit(1); 
	}

}
//void LocalProblem_EOS::der
void LocalProblem_EOS_NonIso_LocalNCP::calc_Deriv_xx(ogsChem::LocalVector & INPUT, ogsChem::LocalVector & OUTPUT, MathLib::LocalMatrix & matSecDer)
{
	//CALCULATE THE DERIVATIVE OF THE SMALL VALUE ON P

	MathLib::LocalVector output_p, output_x, output_T;
	MathLib::LocalVector output_p1, output_x1, output_T1;
	MathLib::LocalVector input_p, input_x, input_T;
	MathLib::LocalVector input_p1, input_x1, input_T1;

	output_p = OUTPUT;
	output_p1 = OUTPUT;
	input_p = INPUT;
	input_p1 = INPUT;
	input_p(0) = input_p(0) + eps;
	// directly use the output value as the initial guess of the next calculation
	solve(input_p, output_p);
	//store the output value
	input_p1(0) = input_p1(0) - eps;
	solve(input_p1, output_p1);

	////CALCULATE THE DERIVATIVE OF THE SMALL VALUE ON X
	output_x = OUTPUT;
	output_x1 = OUTPUT;
	input_x = INPUT;
	input_x1 = INPUT;
	input_x(1) = input_x(1) + eps;
	solve(input_x, output_x);
	input_x1(1) = input_x1(1) - eps;
	solve(input_x1, output_x1);

	////CALCULATE THE DERIVATIVE OF THE SMALL VALUE ON T
	output_T = OUTPUT;
	output_T1 = OUTPUT;
	input_T = INPUT;
	input_T1 = INPUT;
	input_T(2) = input_T(2) + eps;
	solve(input_T, output_T);
	input_T1(2) = input_T1(2) - eps;
	solve(input_T1, output_T1);
	// here is the derivative operations. 
	matSecDer.col(0) = (output_p - output_p1) /2/ eps;
	matSecDer.col(1) = (output_x - output_x1) / 2/eps;
	matSecDer.col(2) = (output_T - output_T1) /2/ eps;
	//std::cout << matSecDer << std::endl;
};
void LocalProblem_EOS_NonIso_LocalNCP::calc_Deriv_aa(ogsChem::LocalVector & INPUT, ogsChem::LocalVector & OUTPUT, MathLib::LocalMatrix & matSecDer)
{
	//CALCULATE THE DERIVATIVE OF THE SMALL VALUE ON P
	MathLib::LocalMatrix Jac_sec = MathLib::LocalMatrix::Zero(N, N);
	MathLib::LocalMatrix Jac_prior = MathLib::LocalMatrix::Zero(N, N-1);
	_EOS->calc_Jacobian(OUTPUT, Jac_sec);
	_EOS->calc_Jacobian_loc_Prior(OUTPUT, Jac_prior);
	
	matSecDer = Jac_sec.inverse()*Jac_prior;
    //std::cout << matSecDer << std::endl;
};
/**
  * TODO: describe this function
  */
void LocalProblem_EOS_NonIso_LocalNCP::solve_minimization(ogsChem::LocalMatrix & J,
	                                      ogsChem::LocalVector & b,
	                                      ogsChem::LocalVector & dx)
{
	// step 0: variable definition and initialization
	int n_r, n_rows_M;
	ogsChem::LocalMatrix Q, R, P, B, RB, V, M, tmp;
	ogsChem::LocalVector z, y, y1, y2;
	Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr_decomp;
	y = ogsChem::LocalVector::Zero(dx.rows());

	// step 1: perform Housholder QR decomposition on J, 
	// so that Q^T * J * P = { R  B }
	//                       { 0  0 }
	qr_decomp.compute(J);
	Q = qr_decomp.matrixQ();
	P = qr_decomp.colsPermutation();
	n_r = qr_decomp.rank();
	n_rows_M = J.cols() - n_r;
	RB = Q.transpose() * J * P;

	if (n_r == J.cols())
	{
		// if n_rank == n_cols, directly solve
		dx = qr_decomp.solve(-b);
	}
	else
	{
		// step 2: split R and B
		R = RB.topLeftCorner(n_r, n_r);
		B = RB.topRightCorner(n_r, RB.cols() - n_r);

		// step 3: if n_rank < n_cols, calculate V, z and y based on R and B. 
		// solve R*V = B
		qr_decomp.compute(R);
		V = qr_decomp.solve(B);
		// Rz = (Q^T *(-b))
		// (I + V^TV)*y2 = V^T * z
		M = ogsChem::LocalMatrix::Identity(n_rows_M, n_rows_M) + V.transpose() * V;
		tmp = (Q.transpose() * (-1.0 * b)).topRows(n_r);
		z = qr_decomp.solve(tmp);
		y2 = M.fullPivHouseholderQr().solve(V.transpose() * z);
		// y1 = z - V*y2
		y1 = z - V * y2;
		// formulate y
		y.head(n_r) = y1;
		y.tail(J.rows() - n_r) = y2;
		// apply permuation
		dx = P * y;
	}
	return;
}

/**
  * Newton iteration with line search
  */
void LocalProblem_EOS_NonIso_LocalNCP::solve_LocalProblem_Newton_LineSearch(std::size_t & flag)
{
	// flag = 0: converged within max_iter; 
	//      = 1: hit max_iter without converg;
	flag = 1; 

	ogsChem::LocalVector U_pre;// vec_residual;
	ogsChem::LocalVector deltaU;

	// number of iterations
	size_t j(0), Iter(0);//j for line search
	const double neta(1.0);// damping factor for ..
	const double alpha(0.5);// damping factor for line search
	double d_norm(0.0), d1_norm(0.0);

	/**
	*  Initialize
	*/
	U_cur = ogsChem::LocalVector::Zero(N);
	U_pre = ogsChem::LocalVector::Zero(N);
	deltaU = ogsChem::LocalVector::Ones(N);
	//Matrix_Jacobian = ogsChem::LocalMatrix::Zero(N, N);
	//Res = ogsChem::LocalVector::Zero(N);

	// start solving the system


#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "Residual Vector: \n";
	// std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

	// save the previous values
	U_cur = U_ini;
	U_pre = U_ini;
	_EOS->eval(U_pre, Res);
	d_norm = Res.norm();
	while (true)
	{

		if (Iter > Iter_tol)
		{
			flag = 1;
			//INFO("residuals : %1.3e", d_norm);
			break; 
		}
		if (d_norm < tot)
		{
			flag = 0; 
			
			break; 
		}

		// construct Jacobian matrix
		_EOS->calc_Jacobian(U_pre, Matrix_Jacobian);
		//_EOS->calc_Jacobian_fd(U_pre, Matrix_Jacobian);
        #ifdef _DEBUG
		// debugging--------------------------
		//std::cout << "Matrix_Jacobian: \n";
		//std::cout << Matrix_Jacobian << std::endl;
		//std::cout << "Residual_Vector: \n";
		//std::cout << Res << std::endl;
		//end of debugging-------------------
        #endif

		//deltaU = Res / Matrix_Jacobian;
		this->solve_minimization(Matrix_Jacobian, Res, deltaU);

        #ifdef _DEBUG
		// debugging--------------------------
		//std::cout << "dx Vector: \n";
		//std::cout << deltaU << std::endl;
		// end of debugging-------------------
        #endif
		
		U_cur = U_pre + neta*deltaU;//  set the damping factor as 0.5

		// evaluate residual with x_new
		_EOS->eval(U_cur, Res);
		// std::cout << "Norm of Residual is "<< std::endl;
		// std::cout << Res.norm() << std::endl;
		j = 0;
		// line search begins
		while (j < Iter_tol)
		{
			// d1_norm = norm(res,inf);
			d1_norm = Res.norm();

			if (d1_norm < d_norm)
				break;

			// updating dx
			deltaU = deltaU * alpha;//

			// increment of unknowns
			//U_cur = U_pre  + deltaU;
			U_cur = U_cur  - deltaU;

			// evaluate residual with x_new
			_EOS->eval(U_cur, Res);

            #ifdef _DEBUG
			// std::cout << "vec_residual: \n";
			// std::cout << Res << std::endl;
			// display the residual
			// std::cout << "Line Search Iteration #" << Iter << "||res|| = " << d1_norm << "||delta_x|| = " << Res.norm() << std::endl;
            #endif

			j++;
		}  // end of while
		d_norm = d1_norm;
		U_pre = U_cur;
		
		// increase the iteration count
		Iter++;
	}
}


