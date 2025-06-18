#pragma once
#include <iostream>
#include <cassert>
#include <vector>
#include <time.h>
#include <algorithm>
#include <map>
#include <tuple>
#include <cstdint>
#include "ReadIns.h"
#include <chrono>
#define NOMINMAX
#include <Windows.h>

using namespace std;
typedef uint64_t data_t;

int CAR_COUNT = 1;//how many different types of cars
int COLOR_COUNT_BITS = 4;
int CAR_BITS = 4;
data_t COLOR_MASK = 1;

int DEMAND_BITS_ARRAY[30];
data_t DEMAND_MASK_ARRAY[30];
int MAX_TIME = 3600;    // time limit

//void print_pages()
//{
//    static size_t pagesize = sysconf(_SC_PAGESIZE);
//    int64_t bytes = getCurrentRSS();
//    assert((bytes % pagesize) == 0);
//    size_t pages = bytes / pagesize;
//    //std::cout << "page size: " << pagesize << "\t";
//    std::cout << "Usage bytes: " << bytes << std::endl;
//    //std::cout << "pages: " << pages << std::endl;
//}

//  Windows
// #ifdef _WIN32
// #include <Windows.h>
double get_wall_time() {
    LARGE_INTEGER time, freq;
    if (!QueryPerformanceFrequency(&freq)) {
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)) {
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time() {
    FILETIME a, b, c, d;
    if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0) {
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
                ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }
    else {
        //  Handle error
        return 0;
    }
}

////  Posix/Linux
//#else
//#include <time.h>
//#include <sys/time.h>
//double get_wall_time() {
//    struct timeval time;
//    if (gettimeofday(&time, NULL)) {
//        //  Handle error
//        return 0;
//    }
//    return (double)time.tv_sec + (double)time.tv_usec * .000001;
//}
//double get_cpu_time() {
//    return (double)clock() / CLOCKS_PER_SEC;
//}
//#endif

bool IsPowerOf2(int x)
{
    return x > 0 && !(x & (x - 1));
}

data_t GetMask(int bits)
{
    data_t mask = 0;
    for (int i = 0; i < bits; i++)
        mask = (mask << 1) + 1;
    return mask;
}

//Generate the key from a set of demands where one car is added
data_t GetKey(int car_to_add, int color_count, std::vector<uint8_t>& demands)
{
    data_t key = car_to_add;
    key = (key << COLOR_COUNT_BITS) + color_count;
    for (size_t i = 0; i < demands.size(); i++)
        key = (key << DEMAND_BITS_ARRAY[i]) + demands[i] + (i == car_to_add ? 1 : 0);
    return key;
}
data_t GetKey2(int car_id, int color_count, std::vector<uint8_t>& demands)
{
    data_t key = car_id;
    key = (key << COLOR_COUNT_BITS) + color_count;
    for (size_t i = 0; i < demands.size(); i++)
        key = (key << DEMAND_BITS_ARRAY[i]) + demands[i];
    return key;
}
void GetData(data_t key, int& car_id, int& color_count, std::vector<uint8_t>& demands)
{
    for (int i = 0; i < CAR_COUNT; i++)
    {
        demands[CAR_COUNT - i - 1] = key & DEMAND_MASK_ARRAY[CAR_COUNT - i - 1];
        key = key >> DEMAND_BITS_ARRAY[CAR_COUNT - i - 1];
    }
    color_count = key & COLOR_MASK;
    key = key >> COLOR_COUNT_BITS;
    car_id = key;
}


class Node2
{
public:
    uint32_t id;
    uint32_t parent;
    data_t data;
    uint64_t distance;
    //uint16_t distance;
    // int distance;   // depends on the value of alpha and beta

    Node2() :id(-1), data(0), distance(9999999), parent(-1) {}
    void Show()
    {
        int car_id; int color_count;
        std::vector<uint8_t> demands(CAR_COUNT, 0);
        GetData(data, car_id, color_count, demands);
        int dmd = 0;
        for (size_t i = 0; i < demands.size(); i++)
            dmd += demands[i];
        printf("Id:%d pos:%d car:%d color_count:%d demands:", id, dmd, car_id, color_count);
        for (size_t i = 0; i < demands.size(); i++)
            printf("%d ", (int)demands[i]);
        printf("\n");
    }
};


struct solution_newG1
{
    double time_find_solution;      // milliseconds
    std::vector<int> solution;      // the sequence of car_type id
    double obj_shortest_path;       // the cost of the shortest path
    bool finish_build = true;
};

void Solve(std::vector<Car>& cars, int c0, int alpha, int beta, solution_newG1 & solution_g1)
// void Solve(std::vector<Car>& cars, int c0, int alpha, int beta, std::vector<int>& solution)
{
    int d_T = 0;
    int min_required_bits = 0;
    for (int i = 0; i < cars.size(); i++)
    {
        cars[i].id = i;
        d_T += cars[i].demand;
        int nb_bits = (int)std::ceil(log2(cars[i].demand)) + (IsPowerOf2(cars[i].demand) ? 1 : 0);
        DEMAND_BITS_ARRAY[i] = nb_bits;
        DEMAND_MASK_ARRAY[i] = GetMask(nb_bits);
        min_required_bits += nb_bits;

        printf("Car:%d body:%d color:%d demand:%d ", cars[i].id, cars[i].body, cars[i].color, cars[i].demand);
        printf("bits:%d mask:%d\n", DEMAND_BITS_ARRAY[i], (int)DEMAND_MASK_ARRAY[i]);
    }
    COLOR_COUNT_BITS = (int)std::ceil(log2(c0)) + (IsPowerOf2(c0) ? 1 : 0);
    COLOR_MASK = GetMask(COLOR_COUNT_BITS);
    CAR_BITS = (int)std::ceil(log2(cars.size())) + (IsPowerOf2(cars.size()) ? 1 : 0);
    min_required_bits += COLOR_COUNT_BITS + CAR_BITS;
    if (min_required_bits > sizeof(data_t) * 8)
    {
        printf("The code can't manage that many bits\n");
        printf("Required:%d available:%d\n", min_required_bits, (int)sizeof(data_t) * 8);
        printf("Exiting\n");
        exit(1);
    }

    printf("Car bits:%d\n", CAR_BITS);
    printf("Color bits:%d mask:%d\n", COLOR_COUNT_BITS, (int)COLOR_MASK);
    printf("Color Max:%d Alpha:%d Beta:%d\n", c0, alpha, beta);
    //printf("Size data_t:%lu min_required_bits:%d\n", sizeof(data_t) * 8, min_required_bits);
    //printf("NodeSize:%lu\n", sizeof(Node2));

    vector< vector<Node2*> > node_created(d_T + 1);			// for loop
    Node2* root = new Node2();
    root->id = 0;
    root->data = 0;
    root->distance = 0;
    //root->Show();
    node_created[0].push_back(root);

    std::vector< map<data_t, int> > node_map_list(d_T + 1);

    vector<Node2*> nodes;
    nodes.reserve(100000000);
    nodes.push_back(root);  // the set of nodes
    int node_count = 1;

    bool has_timed_out = false;

    // build the basic graph
    // auto t1 = std::chrono::steady_clock::now();
    double t1, t2;
    t1 = get_cpu_time();
    for (int pos = 0; pos < d_T && !has_timed_out; pos++)
    {
        if (pos >= 1)
        {
            node_created[pos - 1].resize(0);
            node_map_list[pos - 1].clear();
        }

        for (int j = 0; j < node_created[pos].size() && !has_timed_out; j++)
        {
            Node2* prev_node = node_created[pos][j];
            int prev_car; int prev_color_count;
            std::vector<uint8_t> prev_demands(CAR_COUNT, 0);
            GetData(prev_node->data, prev_car, prev_color_count, prev_demands);
            int prev_color = cars[prev_car].color;
            int prev_body = cars[prev_car].body;

            for (size_t i = 0; i < cars.size(); i++)
            {
                Car* car = &cars[i];
                if (prev_demands[car->id] + 1 > car->demand) continue;

                int new_color_count = prev_color == car->color ? prev_color_count + 1 : 1;
                if (new_color_count > c0) continue;

                data_t key = pos + 1 == d_T ? -1 : GetKey(car->id, new_color_count, prev_demands);

                Node2* n;
                int node_pos = node_map_list[pos + 1][key];
                if (node_pos != 0)
                    n = nodes[node_pos];
                else
                {
                    node_map_list[pos + 1][key] = node_count;
                    n = new Node2();
                    n->data = key;
                    n->id = node_count++;

                    //prev_node->Show();
                    //n->Show();
                    //getchar();
                    nodes.push_back(n);
                    node_created[pos + 1].push_back(n);

                    // auto t2 = std::chrono::steady_clock::now();
                    // auto dur1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
                    // double dur11 = dur1.count() * 0.001;

                    // if (nodes.size() % 10000000 == 0 && dur11 > MAX_TIME)
                    if (nodes.size() % 10000000 == 0 && (get_cpu_time() - t1) > MAX_TIME)
                    {
                        has_timed_out = true;
                        i = cars.size();
                        j = node_created[pos].size();
                        pos = d_T;
                        solution_g1.finish_build = false;
                        break;
                    }
                }

                int cost = 0;
                if (pos != 0)
                {
                    if (prev_body != car->body) cost += alpha;
                    if (prev_color != car->color) cost += beta;
                }
                if (n->distance > prev_node->distance + cost)
                {
                    n->distance = prev_node->distance + cost;
                    n->parent = prev_node->id;
                }

            }// for loop i on the cars
        }// for loop j on the nodes
    }

    // auto t3 = std::chrono::steady_clock::now();
    // auto dur2 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t1);
    // double dur22 = dur2.count();

    t2 = get_cpu_time();			// time for basic Graph 1
    cout << "Nodes: " << nodes.size() << endl;
    cout << "Time Building: " << t2 - t1 << endl;
    solution_g1.time_find_solution = (t2 - t1) * 1000;
    // cout << "Time Building: " << std::to_string(dur22 * 0.001) << " sec" << endl;
    // solution_g1.time_find_solution = dur22;

    if (has_timed_out)
    {
        printf("Has timed out\nReturns\n");
        return;
    }
    //print_pages();

      //create the arcs from the last nodes to the terminal node
    Node2* terminal = node_created.back().back();
    printf("Optimal:%d\n", terminal->distance);
    terminal->Show();
    solution_g1.obj_shortest_path = terminal->distance;

    printf("Solution:");
    int cur = terminal->id;
    int prev = -1;
    while (cur != -1)
    {
        Node2* n = nodes[cur];
        int car_id; int color_count; int terminal_car_id = -1;
        std::vector<uint8_t> demands(CAR_COUNT, 0);
        GetData(n->data, car_id, color_count, demands);
        if (prev == terminal->id)
        {
            for (int i = 0; i < CAR_COUNT; i++)
                if (cars[i].demand != demands[i])
                {
                    terminal_car_id = i;
                    break;
                }
        }

        // if (prev == terminal->id) { printf("%d ", terminal_car_id); solution.push_back(terminal_car_id); }
        // if (cur != terminal->id) { printf("%d ", car_id); solution.push_back(car_id); }

        if (prev == terminal->id) { printf("%d ", terminal_car_id); solution_g1.solution.push_back(terminal_car_id); }
        if (cur != terminal->id) { printf("%d ", car_id); solution_g1.solution.push_back(car_id); }

        prev = cur;
        cur = n->parent;
    }
    printf("\n");


    //print_pages();
}