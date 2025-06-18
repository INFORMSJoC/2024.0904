#pragma once
#include "gurobi_c++.h"
#include "ReadIns.h"
#include<chrono>

// use gurobi to compute the complete model


// ----------------------------- the information of the solution
class SolutionInfo
{
public:
	int num_car;
	bool OP_status;	// if get the optimal solution
	double F1;
	double F2;
	double F3;
	double opt_obj;

	double total_time;
	double gap;
	std::vector<int> Body;
	std::vector<int> Color;
	std::vector<int> Config;
	std::vector<int> E;
	std::vector<int> A;

	//SolutionInfo() : num_car(0), OP_status(0), F1(-1), F2(-1), F3(-1), opt_obj(-1), total_time(0.0), gap(-1), Body(num_car),
	//	Color(num_car), Config(num_car), E(num_car), A(num_car) {}

	SolutionInfo() : num_car(0), OP_status(0), F1(-1), F2(-1), F3(-1), opt_obj(-1), total_time(0.0), gap(-1) {}
};



// ----------------------------- define the function to compute the complete model
void CompleteModel(Instance ins, std::string const& ins_name, std::string const& csvfile, std::string const& txtfile,
	double time_limit, SolutionInfo* sol)
{
	// open files for writting
	std::ofstream out_txt_f;
	out_txt_f.open(txtfile, std::ios::out);

	std::ofstream out_csv;
	out_csv.open(csvfile, std::ios::app);


	GRBEnv env = GRBEnv(true);
	env.start();

	auto t1 = std::chrono::high_resolution_clock::now();

	GRBModel model = GRBModel(env);
	// model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_DoubleParam_TimeLimit, time_limit);

	// add variables
	std::vector< std::vector<GRBVar> > U(ins.d_T);     //U[p][w]: for linearization
	for (int p = 0; p < ins.V.size(); p++)
	{
		std::vector<GRBVar> row(ins.W.size());
		for (int w = 0; w < ins.W.size(); w++)
		{
			std::string var_name = "u[" + std::to_string(ins.V[p]) + "][" + std::to_string((int)ins.W[w]) + "]";
			GRBVar u = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
			row[w] = u;
		}
		U[p] = row;
	}



	std::map<std::tuple<int, int, int, int>, GRBVar> X;
	for (int i = 0; i < ins.d_T; i++)
		for (int b = 0; b < ins.B.size(); b++)
			for (int m = 0; m < ins.M.size(); m++)
				for (int j = 0; j < ins.J.size(); j++)
				{
					std::string var_name = "x[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.B[b]) + "][" +
						std::to_string(ins.M[m]) + "][" + std::to_string(ins.J[j]) + "]";
					GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
					X[std::make_tuple(ins.V[i], ins.B[b], ins.M[m], ins.J[j])] = x;
				}




	std::map<std::tuple<int, int>, GRBVar> F;
	for (int i = 0; i < ins.d_T - 1; i++)
	{
		std::string var_name = "f[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.V[i] + 1) + "]";
		GRBVar f = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
		F[std::make_tuple(ins.V[i], ins.V[i] + 1)] = f;
	}




	std::map<std::tuple<int, int>, GRBVar> G;
	for (int i = 0; i < ins.d_T - 1; i++)
	{
		std::string var_name = "g[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.V[i] + 1) + "]";
		GRBVar g = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
		G[std::make_tuple(ins.V[i], ins.V[i] + 1)] = g;
	}





	std::map<std::tuple<int, int>, GRBVar> E;	// y
	for (int i = 0; i < ins.d_T; i++)
		for (int p = 0; p < ins.d_T; p++)
		{
			std::string var_name = "e[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.V[p]) + "]";
			GRBVar e = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
			E[std::make_tuple(ins.V[i], ins.V[p])] = e;
		}





	std::map<std::tuple<int, int>, GRBVar> A;	// h
	for (int i = 0; i < ins.d_T; i++)
		for (int l = 0; l < ins.L.size(); l++)
		{
			std::string var_name = "a[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.L[l]) + "]";
			GRBVar a = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
			A[std::make_tuple(ins.V[i], ins.L[l])] = a;
		}



	std::map<std::tuple<int, int, int>, GRBVar> Z;
	for (int i = 1; i <= ins.d_T; i++)
		for (int v = 1; v <= i; v++)
			for (int l = 1; l <= ins.L.size(); l++)
			{
				std::string var_name = "z[" + std::to_string(i) + "][" + std::to_string(v) + "][" + std::to_string(l) + "]";
				GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
				Z[std::make_tuple(i, v, l)] = z;
			}



	std::map<std::tuple<int, int, int, int, int>, GRBVar> D;
	for (int i = 1; i <= ins.d_T; i++)
		for (int b = 1; b <= ins.B.size(); b++)
			for (int m = 1; m <= ins.M.size(); m++)
				for (int j = 1; j <= ins.J.size(); j++)
					for (int p = 1; p <= ins.d_T; p++)
					{
						std::string var_name = "d[" + std::to_string(i) + "][" + std::to_string(b) + "][" +
							std::to_string(m) + "][" + std::to_string(j) + "][" + std::to_string(p) + "]";
						GRBVar d = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
						D[std::make_tuple(i, b, m, j, p)] = d;
					}



	// add constriants
	for (int p = 1; p <= ins.d_T; p++)
		for (int w = 1; w <= ins.W.size(); w++)
		{
			GRBLinExpr con1 = 0;
			for (int j = 1; j <= ins.J.size(); j++)
				for (int p_ = 1; p_ <= p; p_++)
					for (int i = 1; i <= ins.d_T; i++)
						for (int b = 1; b <= ins.B.size(); b++)
							for (int m = 1; m <= ins.M.size(); m++)
								con1 += ins.r_j_w[j][w] * D[std::make_tuple(i, b, m, j, p_)];

			model.addConstr(U[p - 1][w - 1] - con1 + p * ins.sigma_w[w], GRB_GREATER_EQUAL, 0.0);
			model.addConstr(U[p - 1][w - 1] + con1 - p * ins.sigma_w[w], GRB_GREATER_EQUAL, 0.0);
		}



	for (int i = 1; i <= ins.d_T; i++)
		for (int p = 1; p <= ins.d_T; p++)
			for (int b = 1; b <= ins.B.size(); b++)
				for (int m = 1; m <= ins.M.size(); m++)
					for (int j = 1; j <= ins.J.size(); j++)
					{
						GRBLinExpr con_1 = D[std::make_tuple(i, b, m, j, p)] - E[std::make_tuple(i, p)];
						GRBLinExpr con_2 = D[std::make_tuple(i, b, m, j, p)] - X[std::make_tuple(i, b, m, j)];
						GRBLinExpr con_3 = D[std::make_tuple(i, b, m, j, p)] - E[std::make_tuple(i, p)] -
							X[std::make_tuple(i, b, m, j)] + 1;

						model.addConstr(con_1, GRB_LESS_EQUAL, 0.0);
						model.addConstr(con_2, GRB_LESS_EQUAL, 0.0);
						model.addConstr(con_3, GRB_GREATER_EQUAL, 0.0);
					}



	for (int i = 0; i < ins.d_T; i++)
	{
		GRBLinExpr con2 = 0;
		for (int b = 0; b < ins.B.size(); b++)
			for (int m = 0; m < ins.M.size(); m++)
				for (int j = 0; j < ins.J.size(); j++)
					con2 += X[std::make_tuple(ins.V[i], ins.B[b], ins.M[m], ins.J[j])];

		model.addConstr(con2, GRB_EQUAL, 1.0);
	}





	// demand
	for (int i = 0; i < ins.d_bmj.size(); i++)
	{
		GRBLinExpr con3 = 0;
		int b = ins.d_bmj[i].body;
		int c = ins.d_bmj[i].color;
		int j = ins.d_bmj[i].config;
		for (int p = 0; p < ins.d_T; p++)
			con3 += X[std::make_tuple(ins.V[p], b, c, j)];
		model.addConstr(con3, GRB_EQUAL, ins.d_bmj[i].demand);
	}




	// color batch
	for (int p = 1; p <= ins.d_T - ins.c0; p++)
		for (int m = 0; m < ins.M.size(); m++)
		{
			GRBLinExpr con4 = 0;
			for (int i = p; i <= p + ins.c0; i++)
				for (int b = 0; b < ins.B.size(); b++)
					for (int j = 0; j < ins.J.size(); j++)
						con4 += X[std::make_tuple(i, ins.B[b], ins.M[m], ins.J[j])];
			model.addConstr(con4, GRB_LESS_EQUAL, ins.c0);
		}




	for (int i = 0; i < ins.d_T - 1; i++)
		for (int b = 0; b < ins.B.size(); b++)
		{
			GRBLinExpr con5 = 0;

			for (int m = 0; m < ins.M.size(); m++)
				for (int j = 0; j < ins.J.size(); j++)
					con5 += X[std::make_tuple(ins.V[i], ins.B[b], ins.M[m], ins.J[j])];
			for (int m = 0; m < ins.M.size(); m++)
				for (int j = 0; j < ins.J.size(); j++)
					con5 -= X[std::make_tuple(ins.V[i] + 1, ins.B[b], ins.M[m], ins.J[j])];

			model.addConstr(F[std::make_tuple(i + 1, i + 2)] - con5, GRB_GREATER_EQUAL, 0.0);
		}



	for (int i = 0; i < ins.d_T - 1; i++)
		for (int m = 0; m < ins.M.size(); m++)
		{
			GRBLinExpr con6 = 0;


			for (int b = 0; b < ins.B.size(); b++)
				for (int j = 0; j < ins.J.size(); j++)
					con6 += X[std::make_tuple(ins.V[i], ins.B[b], ins.M[m], ins.J[j])];


			for (int b = 0; b < ins.B.size(); b++)
				for (int j = 0; j < ins.J.size(); j++)
					con6 -= X[std::make_tuple(ins.V[i] + 1, ins.B[b], ins.M[m], ins.J[j])];


			model.addConstr(G[std::make_tuple(i + 1, i + 2)] - con6, GRB_GREATER_EQUAL, 0.0);
		}



	for (int i = 0; i < ins.d_T; i++)
	{
		GRBLinExpr con7 = 0;
		for (int p = 0; p < ins.d_T; p++)
			con7 += E[std::make_tuple(ins.V[i], ins.V[p])];
		model.addConstr(con7, GRB_EQUAL, 1.0);
	}




	for (int p = 0; p < ins.d_T; p++)
	{
		GRBLinExpr con8 = 0;
		for (int i = 0; i < ins.d_T; i++)
			con8 += E[std::make_tuple(ins.V[i], ins.V[p])];
		model.addConstr(con8, GRB_EQUAL, 1.0);
	}



	for (int i = 0; i < ins.d_T; i++)
	{
		GRBLinExpr con9 = 0;
		for (int l = 0; l < ins.L.size(); l++)
			con9 += A[std::make_tuple(ins.V[i], ins.L[l])];
		model.addConstr(con9, GRB_EQUAL, 1.0);
	}


	for (int i1 = 1; i1 < ins.d_T; i1++)
		for (int p1 = 1; p1 <= ins.d_T; p1++)
			for (int i2 = i1 + 1; i2 <= ins.d_T; i2++)
				for (int l = 1; l <= ins.L.size(); l++)
				{
					GRBLinExpr con10 = A[std::make_tuple(i1, l)] + A[std::make_tuple(i2, l)];
					for (int p = p1; p <= ins.d_T; p++)
						con10 += E[std::make_tuple(i1, p)];
					for (int p = 1; p <= p1; p++)
						con10 += E[std::make_tuple(i2, p)];
					model.addConstr(con10, GRB_LESS_EQUAL, 3.0);
				}





	for (int v = 2; v <= ins.d_T; v++)
		for (int l = 1; l <= ins.L.size(); l++)
			for (int k = 1; k < v; k++)
				for (int i = v; i <= ins.d_T; i++)
				{
					GRBLinExpr con11 = A[std::make_tuple(k, l)] + A[std::make_tuple(v, l)];
					for (int l_ = 1; l_ <= ins.L.size(); l_++)
						con11 += Z[std::make_tuple(i, v, l_)];
					for (int i_ = k; i_ < i; i_++)
						for (int l_ = 1; l_ <= ins.L.size(); l_++)
							con11 -= Z[std::make_tuple(i_, k, l_)];
					model.addConstr(con11, GRB_LESS_EQUAL, 2.0);
				}



	for (int i = 1; i <= ins.d_T; i++)
	{
		GRBLinExpr con12 = 0;
		for (int v = 1; v <= i; v++)
			for (int l = 1; l <= ins.L.size(); l++)
				con12 += Z[std::make_tuple(i, v, l)];
		model.addConstr(con12, GRB_LESS_EQUAL, 1.0);
	}



	for (int v = 1; v <= ins.d_T; v++)
	{
		GRBLinExpr con13 = 0;
		for (int i = v; i <= ins.d_T; i++)
			for (int l = 1; l <= ins.L.size(); l++)
				con13 += Z[std::make_tuple(i, v, l)];
		model.addConstr(con13, GRB_LESS_EQUAL, 1.0);
	}




	for (int i = 1; i <= ins.d_T; i++)
	{
		GRBLinExpr con14 = 0;
		for (int p = 1; p <= ins.d_T; p++)
			con14 += p * E[std::make_tuple(i, p)];
		con14 -= i;
		for (int l = 1; l <= ins.L.size(); l++)
			con14 += i * Z[std::make_tuple(i, i, l)];
		model.addConstr(con14, GRB_GREATER_EQUAL, 1.0);
	}



	for (int i = 2; i <= ins.d_T; i++)
	{
		GRBLinExpr con15 = ins.d_T - 1;
		for (int p = 1; p <= ins.d_T; p++)
			con15 += p * E[std::make_tuple(i, p)];

		for (int i_ = 1; i_ < i; i_++)
			for (int v = 1; v <= i_; v++)
				for (int l = 1; l <= ins.L.size(); l++)
					con15 -= Z[std::make_tuple(i_, v, l)];

		for (int l = 1; l <= ins.L.size(); l++)
			con15 -= ins.d_T * Z[std::make_tuple(i, i, l)];

		model.addConstr(con15, GRB_GREATER_EQUAL, 0.0);
	}




	for (int i = 2; i <= ins.d_T; i++)
	{
		GRBLinExpr con16 = 0;
		for (int p = 1; p <= ins.d_T; p++)
			con16 += p * E[std::make_tuple(i, p)];

		for (int i_ = 1; i_ < i; i_++)
			for (int v = 1; v <= i_; v++)
				for (int l = 1; l <= ins.L.size(); l++)
					con16 -= Z[std::make_tuple(i_, v, l)];

		for (int l = 1; l <= ins.L.size(); l++)
			con16 += ins.d_T * Z[std::make_tuple(i, i, l)];

		model.addConstr(con16, GRB_LESS_EQUAL, ins.d_T + 1);
	}




	for (int i = 2; i <= ins.d_T; i++)
		for (int v = 1; v < i; v++)
		{
			GRBLinExpr con17 = 0;

			for (int p = 1; p <= ins.d_T; p++)
				con17 += p * E[std::make_tuple(i, p)];

			for (int p = 1; p <= ins.d_T; p++)
				con17 -= p * E[std::make_tuple(v, p)];

			for (int l = 1; l <= ins.L.size(); l++)
				con17 -= ins.d_T * Z[std::make_tuple(i, v, l)];
			model.addConstr(con17, GRB_GREATER_EQUAL, 1 - ins.d_T);
		}




	for (int i = 1; i <= ins.d_T; i++)
		for (int l = 1; l <= ins.L.size(); l++)
		{
			GRBLinExpr con18 = 0;
			for (int i_ = 1; i_ <= i; i_++)
				con18 += A[std::make_tuple(i_, l)];
			for (int i_ = 1; i_ <= i; i_++)
				for (int v = 1; v <= i_; v++)
					con18 -= Z[std::make_tuple(i_, v, l)];
			model.addConstr(con18, GRB_LESS_EQUAL, ins.q0);
		}



	for (int v = 1; v <= ins.d_T; v++)
		for (int l = 1; l <= ins.L.size(); l++)
		{
			GRBLinExpr con19 = 0;
			for (int i = v; i <= ins.d_T; i++)
				con19 += Z[std::make_tuple(i, v, l)];
			con19 -= A[std::make_tuple(v, l)];
			model.addConstr(con19, GRB_LESS_EQUAL, 0.0);
		}



	// set objective
	GRBLinExpr obj = 0;
	for (int i = 1; i < ins.d_T; i++)
		obj += ins.alpha * F[std::make_tuple(i, i + 1)];
	for (int i = 1; i < ins.d_T; i++)
		obj += ins.beta * G[std::make_tuple(i, i + 1)];
	for (int p = 1; p <= ins.d_T; p++)
		for (int w = 1; w <= ins.W.size(); w++)
			obj += ins.gamma * U[p - 1][w - 1];

	model.setObjective(obj, GRB_MINIMIZE);

	std::cout << "Finish building MILP" << std::endl;

	// optimize
	model.optimize();
	auto t2 = std::chrono::high_resolution_clock::now();

	auto t_dur = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

	sol->total_time = model.get(GRB_DoubleAttr_Runtime);

	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		// get the optimal solution
		sol->gap = 0.0;
		sol->OP_status = 1;
		//sol->total_time = t_dur.count();
		//sol->total_time = model.get(GRB_DoubleAttr_Runtime);

		double f1 = 0.0;
		double f2 = 0.0;
		double f3 = 0.0;

		for (int i = 1; i <= ins.d_T - 1; i++)
		{
			f1 += ins.alpha * F[std::make_tuple(i, i + 1)].get(GRB_DoubleAttr_X);
			f2 += ins.beta * G[std::make_tuple(i, i + 1)].get(GRB_DoubleAttr_X);
		}

		for (int p = 1; p <= ins.d_T; p++)
			for (int w = 1; w <= ins.W.size(); w++)
				f3 += ins.gamma * U[p - 1][w - 1].get(GRB_DoubleAttr_X);

		sol->F1 = f1;
		sol->F2 = f2;
		sol->F3 = f3;
		sol->opt_obj = model.get(GRB_DoubleAttr_ObjVal);


		for (int i = 1; i <= ins.d_T; i++)
			for (int b = 1; b <= ins.B.size(); b++)
				for (int m = 1; m <= ins.M.size(); m++)
					for (int j = 1; j <= ins.J.size(); j++)
						if (X[std::make_tuple(i, b, m, j)].get(GRB_DoubleAttr_X) > 0.05)
						{
							//sol->Body[i - 1] = b;
							//sol->Color[i - 1] = m;
							//sol->Config[i - 1] = j;

							sol->Body.push_back(b);
							sol->Color.push_back(m);
							sol->Config.push_back(j);

						}

		for (int i = 1; i <= ins.d_T; i++)
			for (int l = 1; l <= ins.L.size(); l++)
				if (A[std::make_tuple(i, l)].get(GRB_DoubleAttr_X) > 0.05)
					sol->A.push_back(l);
					//sol->A[i - 1] = l;

		for (int p = 1; p <= ins.d_T; p++)
			for (int i = 1; i <= ins.d_T; i++)
				if (E[std::make_tuple(i, p)].get(GRB_DoubleAttr_X) > 0.05)
					sol->E.push_back(i);
					//sol->E[p - 1] = i;

		std::cout << "Get the optimal solution!" << std::endl;
		std::cout << "F1: " << f1 << ", F2: " << f2 << ", F3: " << f3 << std::endl;

		out_txt_f << "Get the optimal solution!" << std::endl;

		out_txt_f << "The body type of each car: " << std::endl;
		for (int bb = 0; bb < sol->Body.size(); bb++) out_txt_f << std::to_string(sol->Body[bb]) << " - ";
		out_txt_f << std::endl;

		out_txt_f << "The color of each car: " << std::endl;
		for (int cc = 0; cc < sol->Color.size(); cc++) out_txt_f << std::to_string(sol->Color[cc]) << " - ";
		out_txt_f << std::endl;

		out_txt_f << "The configuration of each car: " << std::endl;
		for (int oo = 0; oo < sol->Config.size(); oo++) out_txt_f << std::to_string(sol->Config[oo]) << " - ";
		out_txt_f << std::endl;

		out_txt_f << "The lane to which each car is allcoated: " << std::endl;
		for (int ll = 0; ll < sol->A.size(); ll++) out_txt_f << std::to_string(sol->A[ll]) << " - ";
		out_txt_f << std::endl;

		out_txt_f << "The upstream car_id in each downstream position: " << std::endl;
		for (int ee = 0; ee < sol->E.size(); ee++) out_txt_f << std::to_string(sol->E[ee]) << " - ";
		out_txt_f << std::endl;

		out_csv << ins_name << ',' << std::to_string(model.get(GRB_IntAttr_Status)) << ',' 
			<< std::to_string(sol->OP_status) << ',' << std::to_string(sol->gap) << ','
			<< std::to_string(sol->total_time) << ','
			<< std::to_string(sol->opt_obj) << ',' << std::to_string(sol->F1) << ','
			<< std::to_string(sol->F2) << ',' << std::to_string(sol->F3) << std::endl;

	}
	else
	{
		// do not get the optimal solution
		out_txt_f << "Do not Get the optimal solution!" << std::endl;
		std::cout << "Do not Get the optimal solution!" << std::endl;

		sol->total_time = t_dur.count();
		sol->gap = model.get(GRB_DoubleAttr_MIPGap);

		if (sol->gap > 10)
			out_csv << ins_name << ',' << std::to_string(model.get(GRB_IntAttr_Status)) << ',' 
			<< std::to_string(sol->OP_status) << ',' << std::to_string(-1) << ','
			<< std::to_string(sol->total_time * 0.001) << ','
			<< std::to_string(sol->opt_obj) << ',' << std::to_string(sol->F1) << ','
			<< std::to_string(sol->F2) << ',' << std::to_string(sol->F3) << std::endl;
		else 
		{
			out_csv << ins_name << ',' << std::to_string(model.get(GRB_IntAttr_Status)) << ',' 
				<< std::to_string(sol->OP_status) << ',' << std::to_string(sol->gap) << ','
				<< std::to_string(sol->total_time * 0.001) << ','
				<< std::to_string(sol->opt_obj) << ',' << std::to_string(sol->F1) << ','
				<< std::to_string(sol->F2) << ',' << std::to_string(sol->F3) << std::endl;
		}

	}

	out_txt_f.close();
	out_csv.close();

	std::cout << "Finish solving " << ins_name << std::endl;

	return;

}

