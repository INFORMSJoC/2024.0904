#pragma once
#pragma warning(disable:4267)
#pragma warning(disable:4244)
#include "ReadIns.h"
#include <chrono>
#include <cmath>

int COLOR_COUNT_BITS = 4;
int CAR_BITS = 4;
int DEMAND_BITS = 6;
int POSITION_BITS = 6;

bool IsPowerOf2(int x)
{
    return x > 0 && !(x & (x - 1));
}

// Generate the key for each node
uint64_t GetKey(int position, int car_to_add, int color_count)
{
    uint64_t key = position;
    key = (key << POSITION_BITS) + car_to_add;
    key = (key << CAR_BITS) + color_count;
    return key;
}

class Node
{
public:
    int id;
    uint8_t position;
    uint8_t car;
    uint8_t color_count;

    Node(int nid, int pos, int ca, int cc) : id(nid), position(pos), car(ca), color_count(cc) {}
    Node(Node* n) : id(n->id), position(n->position), car(n->car), color_count(n->color_count) {}

    void Show()
    {
        printf("Id:%d pos:%d car:%d color_count:%d", id, position, car, color_count);
        printf("\n");
    }
};

class Arc_
{
public:
    int from;
    int to;
    uint8_t car;
    int cost;
    int body;
    int color;
    Arc_() : from(-1), to(-1), cost(0), car(-1), body(-1), color(-1) {}
    void Show()
    {
        printf("arc cost:%3d car:%d from:%d to:%d\n", cost, car, from, to);
    }
};

class Graph2
{
public:
    std::vector<Node*> nodes;
    std::vector<Arc_*> arcs;
    double tt1;     // time for building G2

    std::vector<Arc_*> arcs_0;
    std::vector<Node*> nodes_terminal;
    std::vector<Node*> nodes_inter;

    std::map< int, std::vector<Arc_*> > in_A;
    std::map< int, std::vector<Arc_*> > out_A;

    std::map<std::tuple<int, int>, Arc_*> from_to_arc;  // key: from node XX to node XX; value: arc

    Graph2() : tt1(0.0) {}

    ~Graph2() { Clear(); }
    void Clear()
    {
        for (size_t i = 0; i < nodes.size(); i++)
            delete nodes[i];
        for (size_t i = 0; i < arcs.size(); i++)
            delete arcs[i];

        nodes.clear();
        arcs.clear();
        arcs_0.clear();
        nodes_terminal.clear();
        nodes_inter.clear();

        in_A.clear();
        out_A.clear();
        from_to_arc.clear();
    }
};

void BuildGraph2(Graph2& g, Instance ins)
{
    int max_body = 0; int max_color = 0;
    POSITION_BITS = 0;
    DEMAND_BITS = 0;
    for (int i = 0; i < ins.cars.size(); i++)
    {
        ins.cars[i].id = i;
        ins.cars[i].Show();
        max_body = std::max(ins.cars[i].body, max_body);
        max_color = std::max(ins.cars[i].color, max_color);
        DEMAND_BITS = std::max(DEMAND_BITS, (IsPowerOf2(ins.cars[i].demand) ? 1 : 0) + (int)std::ceil(log2(ins.cars[i].demand)));
    }
    POSITION_BITS = (int)std::ceil(log2(ins.d_T)) + (IsPowerOf2(ins.d_T) ? 1 : 0);
    CAR_BITS = (int)std::ceil(log2(ins.cars.size())) + (IsPowerOf2(ins.cars.size()) ? 1 : 0);
    COLOR_COUNT_BITS = (int)std::ceil(log2(ins.c0)) + (IsPowerOf2(ins.c0) ? 1 : 0);
    int total_bits = POSITION_BITS + CAR_BITS + COLOR_COUNT_BITS + DEMAND_BITS;
    if (total_bits > sizeof(uint64_t) * 8)
    {
        printf("The code can't manage that many bits.\n");
        printf("Exiting\n");
        exit(1);
    }

    printf("Position bits:%d  Car bits:%d  Color_count bits:%d\n", POSITION_BITS, CAR_BITS, COLOR_COUNT_BITS);
    printf("Color Max:%d Alpha:%d Beta:%d\n", ins.c0, ins.alpha, ins.beta);


    // ----------------------------------- create the source node
    std::vector< std::vector<Node*> > node_created(ins.d_T + 1);      // for loop
    Node* root = new Node(0, 0, 0, 0);  // the source node id is 0
    node_created[0].push_back(root);

    std::vector< std::map<uint64_t, int> > node_map_list(ins.d_T + 1);

    std::vector<Arc_*> arcs;                 // the set of arcs
    std::vector<Node*> nodes;                // the set of nodes
    arcs.reserve(10000000);
    nodes.reserve(10000000);
    std::vector<Arc_*> arcs_0;      // the set of arcs whose source node is the source node
    std::vector<Arc_*> arcs_t;      // the set of arcs whose terminal node points to the last position
    std::map< int, std::vector<Arc_*> > incoming_arcs;  // key: node, value: the set of arcs whose end node is this node
    std::map< int, std::vector<Arc_*> > outgoing_arcs;  // key: node, value: the set of arcs whose begin node is this node

    std::vector<Node*> nodes_terminal;  // the set of nodes which point to the last position
    std::vector<Node*> nodes_inter;     // the set of nodes between the source node and the ternimal nodes


    nodes.push_back(root);
    int node_count = 1;


    // build graph2
    auto t1 = std::chrono::steady_clock::now();

    for (int pos = 0; pos < ins.d_T; pos++)
    {
        for (int j = 0; j < node_created[pos].size(); j++)
        {
            Node* prev_node = node_created[pos][j];
            int prev_body = ins.cars[prev_node->car].body;
            int prev_color = ins.cars[prev_node->car].color;

            for (size_t i = 0; i < ins.cars.size(); i++)
            {
                Car* car = &ins.cars[i];
                int new_position = prev_node->position + 1;
                if (new_position > ins.d_T) break;

                int new_color_count = prev_color == car->color ? prev_node->color_count + 1 : 1;
                if (new_color_count > ins.c0) continue;

                uint64_t key = GetKey(new_position, car->id, new_color_count);

                Node* n;
                int node_pos = node_map_list[pos + 1][key];
                if (node_pos != 0)
                    n = nodes[node_pos];
                else
                {
                    // create a new node
                    node_map_list[pos + 1][key] = node_count;
                    n = new Node(prev_node);
                    n->id = node_count++;
                    n->position = new_position;
                    n->car = car->id;
                    n->color_count = new_color_count;

                    nodes.push_back(n);
                    node_created[pos + 1].push_back(n);

                    if (n->position == ins.d_T)
                        nodes_terminal.push_back(n);
                    else
                        nodes_inter.push_back(n);

                }

                // add an arc
                Arc_* a = new Arc_();
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
                if (a->from == 0) arcs_0.push_back(a);
                incoming_arcs[a->to].push_back(a);
                outgoing_arcs[a->from].push_back(a);

                g.from_to_arc[std::make_tuple(a->from, a->to)] = a;

            } // for loop i on the cars

        } // for loop j on the nodes
    }

    //t2 = get_cpu_time();
    auto t2 = std::chrono::steady_clock::now();
    auto d1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    std::cout << "Num of nodes: " << nodes.size() << std::endl;
    std::cout << "Num of arcs: " << arcs.size() << std::endl;
    std::cout << "Num of ternimal nodes: " << nodes_terminal.size() << std::endl;
    std::cout << "Time for building G2: " << d1.count() * 0.001 << std::endl;

    g.nodes = nodes;
    g.arcs = arcs;
    g.arcs_0 = arcs_0;
    g.tt1 = d1.count();


    g.nodes_inter = nodes_inter;
    g.nodes_terminal = nodes_terminal;

    g.in_A = incoming_arcs;
    g.out_A = outgoing_arcs;

}


