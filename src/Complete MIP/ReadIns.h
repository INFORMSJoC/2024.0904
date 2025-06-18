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
    CarT() {}
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
    Car() {}
    Car(int b, int c, int d) :body(b), color(c), demand(d) {}
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
    std::vector<uint8_t> B;
    std::vector<uint8_t> M;
    std::vector<uint8_t> J;
    std::vector<uint8_t> L;
    std::vector<uint8_t> W;
    std::vector<int> V;
    std::vector<CarT> d_bmj;
    std::vector<Car> cars;          // d_bm
    std::map <int, std::map<int, int> > r_j_w;


    std::map<int, float> sigma_w;

    void ShowB()
    {
        printf("B:\n");
        for (int i = 0; i < B.size(); i++) printf("%d ", B[i]);
        std::cout << std::endl;
    }
    void ShowM()
    {
        printf("M:\n");
        for (int i = 0; i < M.size(); i++) printf("%d ", M[i]);
        std::cout << std::endl;
    }
    void ShowJ()
    {
        printf("J:\n");
        for (int i = 0; i < J.size(); i++) printf("%d ", J[i]);
        std::cout << std::endl;
    }
    void ShowL()
    {
        printf("L:\n");
        for (int i = 0; i < L.size(); i++) printf("%d ", L[i]);
        std::cout << std::endl;
    }
    void ShowW()
    {
        printf("W:\n");
        for (int i = 0; i < W.size(); i++) printf("%d ", W[i]);
        std::cout << std::endl;
    }
    void ShowV()
    {
        printf("V:\n");
        for (int i = 0; i < V.size(); i++) printf("%d ", V[i]);
        std::cout << std::endl;
    }
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
        ShowB();
        ShowM();
        ShowJ();
        ShowL();
        ShowW();
        ShowV();
        ShowCT();
        ShowCar();
        ShowR();
        ShowS();
    }
};

void readIns(const char* file_name, Instance& ins)
{
    FILE* ff = fopen(file_name, "r");
    if (!ff)
    {
        printf("Error in the input filename:%s\n", file_name);
        exit(1);
    }


    std::fscanf(ff, "%d\n", &ins.alpha);
    std::fscanf(ff, "%d\n", &ins.beta);
    std::fscanf(ff, "%d\n", &ins.gamma);
    int bb; int mm; int jj; int ww; int ll;
    std::fscanf(ff, "%d\n", &bb);
    std::fscanf(ff, "%d\n", &mm);
    std::fscanf(ff, "%d\n", &jj);
    std::fscanf(ff, "%d\n", &ww);
    std::fscanf(ff, "%d\n", &ll);
    std::fscanf(ff, "%d\n", &ins.q0);
    std::fscanf(ff, "%d\n", &ins.c0);

    //std::cout << ins.alpha << std::endl;
    //std::cout << ins.beta << std::endl;
    //std::cout << ins.gamma << std::endl;
    //std::cout << bb << std::endl;
    //std::cout << mm << std::endl;
    //std::cout << jj << std::endl;
    //std::cout << ww << std::endl;
    //std::cout << ll << std::endl;
    //std::cout << ins.q0 << std::endl;
    //std::cout << ins.q0 << std::endl;

    for (int b = 1; b <= bb; b++) ins.B.push_back(b);
    for (int m = 1; m <= mm; m++) ins.M.push_back(m);
    for (int j = 1; j <= jj; j++) ins.J.push_back(j);
    for (int w = 1; w <= ww; w++) ins.W.push_back(w);
    for (int l = 1; l <= ll; l++) ins.L.push_back(l);
    // for (int l = 1; l <= 2; l++) ins.L.push_back(l);
    // ins.q0 = 2;


    int b1, m1, j1, d1, n_v = 0;
    for (int i = 0; i < bb * mm * jj; i++)
    {
        int r = std::fscanf(ff, "%d %d %d %d\n", &b1, &m1, &j1, &d1);
        if (r == 0) break;
        if (d1 == 0) continue;
        n_v += d1;
        ins.d_bmj.push_back(CarT(b1, m1, j1, d1));
    }

    ins.d_T = n_v;
    for (int v = 1; v <= n_v; v++) ins.V.push_back(v);


    int b2, m2, d2;
    for (int i = 0; i < bb * mm; i++)
    {
        int r = std::fscanf(ff, "%d %d %d \n", &b2, &m2, &d2);
        if (r == 0) break;
        if (d2 == 0) continue;
        ins.cars.push_back(Car(b2, m2, d2));
    }


    int j2, w1, d3;
    for (int i = 0; i < jj * ww; i++)
    {
        int r = std::fscanf(ff, "%d %d %d\n", &j2, &w1, &d3);
        if (r == 0) break;
        ins.r_j_w[j2][w1] = d3;         // j2 is configuration, w1 is option, d3 is the demand
    }


    int w2; float sig;
    for (int i = 0; i < ww; i++)
    {
        int r = std::fscanf(ff, "%d %f \n", &w2, &sig);
        if (r == 0) break;
        ins.sigma_w[w2] = sig;
    }

    std::cout << "********* Finish Reading Ins" << std::endl;
    // ins.ShowCT();
    // ins.ShowCar();
    // ins.ShowV();
}