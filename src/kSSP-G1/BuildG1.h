#pragma once
#pragma warning(disable:4267)
#pragma warning(disable:4244)
#include "solution2.h"
#include <cmath>
#include <algorithm>



int DEMAND_BITS = 6;


//Generate the key from a set of demands where one car is added
uint64_t GetKey3(int car_to_add, int color_count, std::vector<uint8_t>& demands)
{
	uint64_t key = car_to_add;
	key = (key << CAR_BITS) + color_count;
	for (size_t i = 0; i < demands.size(); i++)
		key = (key << DEMAND_BITS) + demands[i] + (i == car_to_add ? 1 : 0);
	return key;
}


class Node
{
public:
	int id;
	uint8_t car;			// car type (body, color) id is "car"
	uint8_t color_count;
	std::vector<uint8_t> demands;

	Node(int nid, int ca, int cc, int car_count) : id(nid), car(ca), color_count(cc), demands(car_count) {}
	Node(Node* n) : id(n->id), car(n->car), color_count(n->color_count), demands(n->demands) {}
	void Show()
	{
		printf("Id:%d pos:%d car:%d color_count:%d demands:", id, GetSumDemands(), car, color_count);
		for (size_t i = 0; i < demands.size(); i++)
			printf("%d ", (int)demands[i]);
		printf("\n");
	}
	int GetSumDemands()
	{
		int dmd = 0;
		for (size_t i = 0; i < demands.size(); i++)
			dmd += demands[i];
		return dmd;
	}
};

class Arc1
{
public:
	int from;
	int to;
	int car;		// which (body, color) type is added - car
	int cost;
	uint8_t body;
	uint8_t color;
	Arc1() :from(-1), to(-1), cost(0), car(-1), body(0), color(0) {}

	void Show() { printf("arc cost:%3d car:%d from:%d to:%d\n", cost, car, from, to); }
};

class Graph1
{
public:
	std::vector<Node*> nodes;
	std::vector<Arc1*> arcs;
	std::vector<Arc1*> ShortestPathArc;
	std::vector<int> ShortestPathNode;
	double tt1;			// time to build G1
	double tt2;			// time to find the shortest path
	double ShortestPathCost;
	int source_id;
	int terminal_id;	// the id of the terminal node
	std::map<std::tuple<int, int>, Arc1*> from_to_arc;   // key: from node XX to node XX; value: arc
	bool is_built;

	Graph1() : tt1(0.0), tt2(0.0), ShortestPathCost(0.0), source_id(-1), terminal_id(-1), is_built(1) {}

	~Graph1() { Clear(); }
	void Clear()
	{
		for (size_t i = 0; i < nodes.size(); i++)
			delete nodes[i];
		for (size_t i = 0; i < arcs.size(); i++)
			delete arcs[i];

		nodes.clear();
		arcs.clear();
		ShortestPathArc.clear();
		from_to_arc.clear();
	}


};


void BuildGraph1(Graph1& g, Instance ins, const std::string& g_file_name, double re_time)
{
	// g_file_name = 10_1 for example
	std::string g_w_file = "../GBS/" + g_file_name + ".txt";
	std::ofstream out_g_f;
	out_g_f.open(g_w_file, std::ios::out);


	int d_T = 0; int max_body = 0; int max_color = 0;
	DEMAND_BITS = 0;
	for (int i = 0; i < ins.cars.size(); i++)
	{
		ins.cars[i].id = i;
		ins.cars[i].Show();
		max_body = std::max(ins.cars[i].body, max_body);
		max_color = std::max(ins.cars[i].color, max_color);
		d_T += ins.cars[i].demand;
		DEMAND_BITS = std::max(DEMAND_BITS, (IsPowerOf2(ins.cars[i].demand) ? 1 : 0) + (int)std::ceil(log2(ins.cars[i].demand)));
	}
	COLOR_COUNT_BITS = (int)std::ceil(log2(ins.c0)) + (IsPowerOf2(ins.c0) ? 1 : 0);
	CAR_BITS = (int)std::ceil(log2(ins.cars.size())) + (IsPowerOf2(ins.cars.size()) ? 1 : 0);
	int total_bits = COLOR_COUNT_BITS + CAR_BITS + DEMAND_BITS * ins.cars.size();
	if (total_bits > sizeof(uint64_t) * 8)
	{
		printf("The code can't manage that many bits\n");
		printf("Exiting\n");
		exit(1);
	}

	// printf("Demand bits:%d Car bits:%d Color Count bits:%d total:%d\n", DEMAND_BITS, CAR_BITS, COLOR_COUNT_BITS, total_bits);
	// rintf("Color Max:%d Alpha:%d Beta:%d\n", ins.c0, ins.alpha, ins.beta);


	// --------------------------- create the source node
	std::vector< std::vector<Node*> > node_created(ins.d_T + 1);			// for loop
	Node* root = new Node(0, 0, 0, ins.cars.size());
	node_created[0].push_back(root);

	std::vector< std::map<uint64_t, int> > node_map_list(d_T + 1);


	std::vector<Arc1*> arcs;							// the set of arcs in the basic graph
	std::vector<Node*> nodes;							// the set of nodes, each element is a Node object
	arcs.reserve(10000000);
	nodes.reserve(10000000);
	nodes.push_back(root);  // the set of nodes
	int node_count = 1;
	g.source_id = 0;		// the id of the source node is 1
	//std::cout << "cost: " << ins.alpha << " " << ins.beta << std::endl;

	// -------- build Graph1
	auto t1 = std::chrono::steady_clock::now();
	for (int pos = 0; pos < d_T; pos++)
		for (int j = 0; j < node_created[pos].size(); j++)
		{
			Node* prev_node = node_created[pos][j];
			int prev_color = ins.cars[prev_node->car].color;
			int prev_body = ins.cars[prev_node->car].body;


			for (size_t i = 0; i < ins.cars.size(); i++)
			{
				Car* car = &ins.cars[i];
				if (prev_node->demands[car->id] + 1 > car->demand) continue;

				int new_color_count = prev_color == car->color ? prev_node->color_count + 1 : 1;
				if (new_color_count > ins.c0) continue;

				uint64_t key = pos + 1 == d_T ? -1 : GetKey3(car->id, new_color_count, prev_node->demands);

				Node* n;
				int node_pos = node_map_list[pos + 1][key];
				if (node_pos != 0)
					n = nodes[node_pos];
				else
				{
					node_map_list[pos + 1][key] = node_count;
					n = new Node(prev_node);
					n->id = node_count++;	// first give the value of node_count to n-id, then node_cout ++
					n->car = car->id;		// add which (body, color) type
					n->color_count = new_color_count;
					n->demands[car->id]++;
					nodes.push_back(n);
					node_created[pos + 1].push_back(n);
				}


				auto atime = std::chrono::steady_clock::now();
				auto a_dur = std::chrono::duration_cast<std::chrono::milliseconds>(atime - t1);
				double a_dur1 = a_dur.count();
				if (a_dur1 * 0.001 > re_time)
				{
					g.is_built = false;
					g.tt1 = a_dur1;
					break;
				}



				Arc1* a = new Arc1();
				a->car = car->id;
				a->cost = 0;
				a->from = prev_node->id;
				a->to = n->id;
				a->body = car->body;
				a->color = car->color;
				if (pos != 0)
				{
					if (prev_body != car->body) a->cost += ins.alpha;
					if (prev_color != car->color) a->cost += ins.beta;
				}

				arcs.push_back(a);
				//std::cout << "arcs: " << a->from << " to: " << a->to << " cost: " << a->cost << std::endl;
				g.from_to_arc[std::make_tuple(a->from, a->to)] = a;



			}	// for loop i on the cars

		}	// for loop j on the nodes

	auto t2 = std::chrono::steady_clock::now();

	// -------------------------- create the arcs from the last nodes to the terminal node
	Node* terminal = node_created.back().back();
	g.terminal_id = terminal->id;

	std::cout << "********* Source node id: " << g.source_id << std::endl;
	std::cout << "********* Terminal node id: " << g.terminal_id << std::endl;

	std::cout << "********* Finish building original G1" << std::endl;
	// ------------------------- output the graph into the g_file_name
	for (int i = 0; i < arcs.size(); i++)
	{
		std::string write_s = "a " + std::to_string(arcs[i]->from) + " " + std::to_string(arcs[i]->to) + " " + std::to_string(arcs[i]->cost);
		out_g_f << write_s << std::endl;
	}

	out_g_f.close();
	std::cout << "********* Finish writing original G1" << std::endl;


	// ------------------------- find the shortest path by Kahn algorithm
	auto t3 = std::chrono::steady_clock::now();
	std::vector<int> degree(nodes.size(), 0);
	std::vector<int> list_degree0;
	degree[terminal->id] = 100000;

	for (int a = 0; a < arcs.size(); a++)
		degree[arcs[a]->to]++;

	std::vector< std::vector<Arc1*> > arcs_outgoing(nodes.size());
	for (int a = 0; a < arcs.size(); a++)
		arcs_outgoing[arcs[a]->from].push_back(arcs[a]);


	// list of nodes with degree 0
	for (size_t i = 0; i + 1 < nodes.size(); i++)
		if (degree[i] == 0)
			list_degree0.push_back(i);

	std::vector<int> distances(nodes.size(), 99999999);
	std::vector<int> parent(nodes.size(), -1);
	std::vector<Arc1*> arc_from(nodes.size());
	distances[0] = 0;

	// iterate through list_degree_0 to get the nodes without a link to the terminal mode
	for (int l = 0; l < list_degree0.size(); l++)
	{
		int i = list_degree0[l];
		for (int a = 0; a < arcs_outgoing[i].size(); a++)
		{
			Arc1* arc = arcs_outgoing[i][a];
			if (distances[arc->from] + arc->cost < distances[arc->to])
			{
				distances[arc->to] = distances[arc->from] + arc->cost;
				parent[arc->to] = arc->from;
				arc_from[arc->to] = arc;
			}
			degree[arc->to]--;
			if (degree[arc->to] == 0)
				list_degree0.push_back(arc->to);
		}
	}
	// printf("Shortest path:%d\n", distances[terminal->id]);
	// std::cout << "Nodes in Shortest path: " << terminal->id << " ";

	Arc1* a = arc_from[terminal->id];
	std::cout << "*----- Shortest path by Kahn: " << a->to << " ";
	g.ShortestPathNode.push_back(a->to);

	while (a != NULL)
	{
		g.ShortestPathArc.push_back(a);
		g.ShortestPathCost += a->cost;
		g.ShortestPathNode.push_back(a->from);

		std::cout << a->from << " ";
		int c = ins.cars[a->car].color;
		int b = ins.cars[a->car].body;
		printf("arc cost:%3d body:%d color:%d from:%d to:%d\n",a->cost,b,c, a->from, a->to);
			//printf("\t"); nodes[ a->to ]->Show();
			//printf("\t"); nodes[ a->from ]->Show();
		a = arc_from[a->from];
	}

	auto t4 = std::chrono::steady_clock::now();
	std::cout << std::endl;


	auto d1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	auto d2 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);

	g.nodes = nodes;
	g.arcs = arcs;
	g.tt1 = d1.count();
	g.tt2 = d2.count();
	std::reverse(g.ShortestPathNode.begin(), g.ShortestPathNode.end());
	std::cout << "********* Finish finding Shorest path by Kahn algorithm" << std::endl;
	std::cout << "----------------------------------------" << std::endl;

}


