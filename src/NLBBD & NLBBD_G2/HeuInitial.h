#pragma once
#include<map>
#include "SolveSP.h"
#include "Elements.h"



void StatConBodyColor(Instance ins, std::map<int, std::vector<std::tuple<int, int>>>& Config_BodyColorMap)
{
	for (int i = 0; i < ins.d_bmj.size(); i++)
		Config_BodyColorMap[ins.d_bmj[i].config].push_back(std::make_tuple(ins.d_bmj[i].body, ins.d_bmj[i].color));
}



bool HeuDone(GRBEnv& env, Instance ins, std::map<int, std::vector<std::tuple<int, int>>>& Con_BodyColor,
	std::vector<int>& BodyList, std::vector<int>& ColorList, LBBD_Info& lbbd, int sp2_method)
{
	// sp2_method: use which method to solve SP2
	bool is_done = false;
	lbbd.Heu_sp = 0;
	lbbd.Kahn_Heu = 0;
	lbbd.Kahn_SP2 = 0;

	std::cout << "H1" << std::endl;
	auto ini_sp_begin = std::chrono::high_resolution_clock::now();

	std::map<std::tuple<int, int>, std::vector<int>> BodyColor_Car;
	std::map<int, std::vector<int>> Con_car;

	for (int i = 0; i < ins.d_T; i++)
		BodyColor_Car[std::make_tuple(BodyList[i], ColorList[i])].push_back(i + 1);		// car_id


	for (auto it = Con_BodyColor.begin(); it != Con_BodyColor.end(); it++)		// j: (b1, m1), (b2, m3)
	{
		int con_j = it->first;	// config
		for (int i = 0; i < it->second.size(); i++)
		{
			// (body,color) = it-second[i]
			for (int car_j = 0; car_j < BodyColor_Car[it->second[i]].size(); car_j++)	// Con_car[j1]: car1, car5, car9
				Con_car[con_j].push_back(BodyColor_Car[it->second[i]][car_j]);
		}
	}

	std::vector<int> better_initial_seq;	// E = []
	std::map<int, int> better_car_pos;		// key = pos, value = car_id

	for (int i = 0; i < lbbd.LB_Config_Seq.size(); i++)
	{
		int con_j = lbbd.LB_Config_Seq[i];		// configuration is con_j
		auto smallest = std::min_element(std::begin(Con_car[con_j]), std::end(Con_car[con_j]));
		int small = *smallest;
		better_initial_seq.push_back(small);		// upstream car id
		better_car_pos[small] = i + 1;				// car small's downstream position is i + 1
		// remove the smallest car_id from other config_car
		for (auto it = Con_car.begin(); it != Con_car.end(); it++)
		{
			std::vector<int>::iterator itr;
			for (itr = it->second.begin(); itr != it->second.end(); itr++)
			{
				if (*itr == small)
				{
					it->second.erase(itr);
					break;
				}
			}
		}
	}

	// get the initial downstream sequence constructed by Heuristic
	// judge if the constructed solution is feasible


	std::vector<std::tuple<int, int>> Ini_upDownSeq;		// upstream carid - downstream position
	for (int n = 0; n < better_initial_seq.size(); n++)
		Ini_upDownSeq.push_back(std::make_tuple(n + 1, better_car_pos[n + 1]));

	bool ini_E_feasible = is_E_feasible(Ini_upDownSeq, ins);
	if (ini_E_feasible == true)
	{
		// do not valid three valid inequalities
		// allocate each cars to a certain lane
		std::map<int, int> ICar_Seq1;		// car_id - pos
		for (int n1 = 0; n1 < ins.d_T; n1++)
			ICar_Seq1[better_initial_seq[n1]] = n1 + 1;

		std::vector<std::vector<int>> IA_sol1(ins.num_l);
		std::vector<int> ILane_allocated1(ins.d_T);
		IA_sol1[0].push_back(1);
		ILane_allocated1[0] = 1;		// put the first upstream car to the first lane

		for (int n1 = 2; n1 <= ins.d_T; n1++)	// car n1
		{
			int Ic_pos = better_car_pos[n1];
			for (int l1 = 0; l1 < ins.num_l; l1++)
			{
				bool Ican_enter = true;
				for (int n2 = 0; n2 < IA_sol1[l1].size(); n2++)
				{
					int Icar_pos = better_car_pos[IA_sol1[l1][n2]];
					if ((n1 < IA_sol1[l1][n2] && Ic_pos > Icar_pos) || (n1 > IA_sol1[l1][n2] && Ic_pos < Icar_pos))
					{
						Ican_enter = false;
						break;
					}
				}

				if (Ican_enter == true)
				{
					IA_sol1[l1].push_back(n1);
					ILane_allocated1[n1 - 1] = l1 + 1;	// car n1 go to lane l1+1
					break;
				}
			}
		}   // finish assigning all cars into PBS


		// check whether the solution A is feasible
		std::map<std::tuple<int, int>, int> IA1;
		for (int Ic = 1; Ic <= ins.d_T; Ic++)
			for (int Il = 1; Il <= ins.num_l; Il++)
				IA1[std::make_tuple(Ic, Il)] = 0;

		for (int Il = 0; Il < IA_sol1.size(); Il++)
			for (int Ic = 0; Ic < IA_sol1[Il].size(); Ic++)
				IA1[std::make_tuple(IA_sol1[Il][Ic], Il + 1)] = 1;

		// get parameters z
		std::map<std::tuple<int, int>, int> Iz1;
		for (int ii = 1; ii <= ins.d_T; ii++)
			for (int pp = 1; pp <= ins.d_T; pp++)
				Iz1[std::make_tuple(ii, pp)] = 0;
		int Isb1 = 1;
		int Icn1 = 0;
		for (int sp = 1; sp <= ins.d_T; sp++)
		{
			int sum_e1 = 0;
			for (int ii = 1; ii <= sp; ii++)
				if (Isb1 == ICar_Seq1[ii]) sum_e1 += 1;

			if (sum_e1 < 1e-4)
			{
				Isb1 = Isb1;
				Icn1 = Icn1;

				if (Icn1 < 1e-4)
					for (int pp = 1; pp < sp + 1; pp++) Iz1[std::make_tuple(sp, pp)] = 0;
				else
				{
					for (int pp = 1; pp < Icn1 + 1; pp++) Iz1[std::make_tuple(sp, pp)] = 1;
					for (int pp = Icn1 + 1; pp < sp + 1; pp++) Iz1[std::make_tuple(sp, pp)] = 0;
				}
			}
			else
			{
				Isb1++;
				Icn1++;
				for (int pp = 1; pp < Icn1 + 1; pp++) Iz1[std::make_tuple(sp, pp)] = 1;
				for (int pp = Icn1 + 1; pp < sp + 1; pp++) Iz1[std::make_tuple(sp, pp)] = 0;
			}
		}

		bool IA_is_feasible1 = true;
		for (int ii = 1; ii <= ins.d_T; ii++)
		{
			for (int ll = 1; ll <= ins.num_l; ll++)
			{
				int Q_il = 0;
				for (int i_ = 1; i_ < ii + 1; i_++) Q_il += IA1[std::make_tuple(i_, ll)];
				for (int pp = 1; pp <= ii; pp++)
					for (int k = 1; k <= ii; k++)
						if (pp == ICar_Seq1[k]) Q_il -= Iz1[std::make_tuple(ii, pp)] * IA1[std::make_tuple(k, ll)];
				if (Q_il > ins.q0)
				{
					IA_is_feasible1 = false;
					break;
				}
			}
			if (IA_is_feasible1 == false) break;
		}

		if (IA_is_feasible1 == true)
		{
			is_done = true;
			lbbd.Kahn_Heu = true;


			lbbd.lA.assign(ILane_allocated1.begin(), ILane_allocated1.end());
			for (int ss = 0; ss < ins.d_T; ss++)
			{
				lbbd.lE.push_back(better_initial_seq[ss]);
			}


			for (int ss = 0; ss < ins.d_T; ss++)
			{
				lbbd.lO.push_back(lbbd.LB_Config_Seq[ss]);
			}

			return is_done;

		}
		else
		{
			// assign A by heuristic is infeasible
			// solve SP-Stage2 model to find a feasible solution if exists
			std::map<std::tuple<int, int>, double> stage1_sol;
			for (int i = 1; i <= ins.d_T; i++)
				for (int p = 1; p <= ins.d_T; p++)
					stage1_sol[std::make_tuple(i, p)] = 0.0;

			for (int i = 0; i < Ini_upDownSeq.size(); i++)
				stage1_sol[std::make_tuple(std::get<0>(Ini_upDownSeq[i]), std::get<1>(Ini_upDownSeq[i]))] = 1.0;

			// solve SP-Stage 2 model
			SP2_Sol stage2_sol = Solve_SP2(env, stage1_sol, ins.d_T, ins.num_l, ins.q0, sp2_method);      // the last parameter means use which method to solve SP-S2

			auto ini_sp_endl1 = std::chrono::high_resolution_clock::now();

			lbbd.total_num_SP2_solved += 1;
			lbbd.total_time_SP2 += stage2_sol.sp2_solve_time;
			lbbd.total_num_SP2_NogoodCut += stage2_sol.num_nogood_cut;

			lbbd.lastSP->num_nogood_cut = stage2_sol.num_nogood_cut;
			lbbd.lastSP->num_sp2_solved = 1;
			lbbd.lastSP->time_solve_sp2 = stage2_sol.sp2_solve_time;

			if (stage2_sol.sp2_status == 1)
			{
				// update the information of lbbd
				is_done = true;
				lbbd.Kahn_SP2 = true;

				for (int ss = 0; ss < ins.d_T; ss++)
				{
					lbbd.lE.push_back(better_initial_seq[ss]);
					lbbd.lA.push_back(stage2_sol.A[ss]);
				}
				for (int ss = 0; ss < ins.d_T; ss++)
				{
					lbbd.lO.push_back(lbbd.LB_Config_Seq[ss]);
				}

				return is_done;

			}
		}

	}


	return is_done;


}