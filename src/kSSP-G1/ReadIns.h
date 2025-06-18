#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <iostream>
#include <vector>
#include <time.h>
#include <algorithm>
#include <map>
#include <tuple>
#include <cstdint>
#include <string>


class CarT      // d_bmj (body b, color m, config j, demand dmd)
{
public:
    int body;
    int color;
    int config;
    int demand;
    //CarT() {}
    CarT(int b, int m, int j, int d) : body(b), color(m), config(j), demand(d) {}
    void Show()
    {
        printf("Car_type: body:%d color:%d config:%d demand:%d\n", body, color, config, demand);
    }
};

class Car       // d_bm (body b, color m, demand dmd)
{
public:
    int id;
    int body;
    int color;
    int demand;
    //Car() {}
    Car(int b, int c, int d) :id(-1), body(b), color(c), demand(d) {}
    void Show()
    {
        printf("Car:%d body:%d color:%d demand:%d\n", id, body, color, demand);
    }
};


class Instance
{
public:
    int alpha;
    int beta;
    int gamma;
    int q0;         // the capacity of the lane
    int c0;         // color batch limit
    int d_T;
    std::vector<CarT> d_bmj;
    std::vector<Car> cars;          // d_bm
    std::map <int, std::map<int, int> > r_j_w;      // j, w, demand

    int num_b;
    int num_m;
    int num_j;
    int num_l;
    int num_w;
    std::map<int, float> sigma_w;

    Instance() :alpha(1), beta(1), gamma(1), q0(0), c0(0), d_T(0), num_b(0), num_m(0), num_j(0), num_l(0), num_w(0) {}

    void ShowCT()
    {
        printf("Car type:\n");
        for (int i = 0; i < d_bmj.size(); i++)
            d_bmj[i].Show();
        std::cout << std::endl;
    }
    void ShowCar()
    {
        printf("Cars:\n");
        for (int i = 0; i < cars.size(); i++)
            cars[i].Show();
        std::cout << std::endl;
    }
    void ShowR()
    {
        printf("r_j_w:\n");
        for (auto it1 = r_j_w.begin(); it1 != r_j_w.end(); it1++)
        {
            for (auto it2 = r_j_w[it1->first].begin(); it2 != r_j_w[it1->first].end(); it2++)
                printf("j:%d w:%d dmd:%d\n", it1->first, it2->first, it2->second);
        }

        std::cout << std::endl;
    }
    void ShowS()
    {
        std::cout << "sigma:" << std::endl;
        for (auto it = sigma_w.begin(); it != sigma_w.end(); it++)
        {
            std::cout << it->first << ": " << it->second << std::endl;
        }
        std::cout << std::endl;
    }
    void Show()
    {
        std::cout << "alpha: " << alpha << std::endl;
        std::cout << "beta: " << beta << std::endl;
        std::cout << "gamma: " << gamma << std::endl;
        std::cout << "q_0: " << q0 << std::endl;
        std::cout << "c_0: " << c0 << std::endl;
        std::cout << "d_T: " << d_T << std::endl;
        ShowCT();
        ShowCar();
        ShowR();
        ShowS();
    }
};

void readIns(const char* file_name, Instance& ins)
{
    FILE* ff;
    fopen_s(&ff, file_name, "r");
    if (ff == NULL)
    {
        printf("Error in the input filename:%s\n", file_name);
        exit(1);
    }

    fscanf_s(ff, "%d\n", &ins.alpha);
    fscanf_s(ff, "%d\n", &ins.beta);
    fscanf_s(ff, "%d\n", &ins.gamma);
    fscanf_s(ff, "%d\n", &ins.num_b);
    fscanf_s(ff, "%d\n", &ins.num_m);
    fscanf_s(ff, "%d\n", &ins.num_j);
    fscanf_s(ff, "%d\n", &ins.num_w);
    fscanf_s(ff, "%d\n", &ins.num_l);
    fscanf_s(ff, "%d\n", &ins.q0);
    fscanf_s(ff, "%d\n", &ins.c0);

    //ins.num_l = 2;
    //ins.q0 = 1;

    if (ins.d_T < 35)
    {
        ins.num_l = 3;
        ins.q0 = 3;
    }
    else if (ins.d_T < 65 && ins.d_T > 35)
    {
        ins.num_l = 5;
        ins.q0 = 5;
    }
    else
    {
        ins.num_l = 6;
        ins.q0 = 6;
    }

    // ins.num_l = 2;
    // ins.q0 = 1;

    ins.alpha = 1;
    ins.beta = 10000;
    ins.gamma = 100;

    //ins.num_l = insl;
    //ins.q0 = insq;
    //ins.alpha = insal;
    //ins.beta = insbe;
    //ins.gamma = insga;


    int b1, m1, j1, d1, n_v = 0;
    for (int i = 0; i < ins.num_b * ins.num_m * ins.num_j; i++)
    {
        int r = fscanf_s(ff, "%d %d %d %d\n", &b1, &m1, &j1, &d1);
        if (r == 0) break;
        if (d1 == 0) continue;
        n_v += d1;
        ins.d_bmj.push_back(CarT(b1, m1, j1, d1));
    }

    ins.d_T = n_v;
    // for (int v = 1; v <= n_v; v++) ins.V.push_back(v);


    int b2, m2, d2;
    for (int i = 0; i < ins.num_b * ins.num_m; i++)
    {
        int r = fscanf_s(ff, "%d %d %d \n", &b2, &m2, &d2);
        if (r == 0) break;
        if (d2 == 0) continue;
        ins.cars.push_back(Car(b2, m2, d2));
    }


    int j2, w1, d3;
    for (int i = 0; i < ins.num_j * ins.num_w; i++)
    {
        int r = fscanf_s(ff, "%d %d %d\n", &j2, &w1, &d3);
        ins.r_j_w[j2][w1] = d3;         // j2 is configuration, w1 is option, d3 is the demand
    }


    int w2; float sig;
    for (int i = 0; i < ins.num_w; i++)
    {
        int r = fscanf_s(ff, "%d %f \n", &w2, &sig);
        ins.sigma_w[w2] = sig;
    }

    std::cout << "********* Finish Reading Ins" << std::endl;



}