#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>

// newLiftCon2
void LiftCon2(std::vector<std::tuple<int, int>>& liftSet, int num_cars)
{
	for (int i = 1; i <= num_cars; i++)
		for (int p = 1; p <= num_cars; p++)
		{
			auto it = std::find(liftSet.begin(), liftSet.end(), std::make_tuple(i, p));

			if (!(it != liftSet.end()))
			{
				bool in_conflict = true;
				for (int s = 0; s < liftSet.size(); s++)
				{
					int k1 = std::get<0>(liftSet[s]);
					int p1 = std::get<1>(liftSet[s]);

					if (i == k1 || p == p1) continue;
					else if ((i > k1 && p > p1) || (i < k1 && p < p1))
					{
						in_conflict = false;
						break;
					}
				}

				if (in_conflict == true)
					liftSet.push_back(std::make_tuple(i, p));
			}
		}
}



void liftCon3(std::vector<std::tuple<int, int>>& liftSet, std::vector<std::tuple<int, int>> invalidClique, std::vector<std::tuple<int, int>> all_matches)
{
	for (int i = 0; i < invalidClique.size(); i++)
		liftSet.push_back(invalidClique[i]);

	//std::cout << "Original lazy3 tuple: " << std::endl;
	//for (int i = 0; i < liftSet.size(); i++)
	//	std::cout << "(" << std::get<0>(liftSet[i]) << ',' << std::get<1>(liftSet[i]) << ") - ";
	//std::cout << std::endl;

	for (int i = 0; i < all_matches.size(); i++)
	{
		auto it = std::find(liftSet.begin(), liftSet.end(), all_matches[i]);
		if (!(it != liftSet.end()))
		{
			bool in_conflict = true;
			int k1 = std::get<0>(all_matches[i]);
			int p1 = std::get<1>(all_matches[i]);

			for (int j = 0; j < liftSet.size(); j++)
			{
				int k2 = std::get<0>(liftSet[j]);
				int p2 = std::get<1>(liftSet[j]);

				if (k1 == k2 || p1 == p2) continue;
				else if ((k1 > k2 && p1 > p2) || (k1 < k2 && p1 < p2))
				{
					in_conflict = false;
					break;
				}
			}

			if (in_conflict) liftSet.push_back(all_matches[i]);
		}
	}

	/*std::cout << "After lifting: " << std::endl;
	for (int i = 0; i < liftSet.size(); i++)
		std::cout << "(" << std::get<0>(liftSet[i]) << ',' << std::get<1>(liftSet[i]) << ") - ";
	std::cout << std::endl;*/

}