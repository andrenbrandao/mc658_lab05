/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#ifndef PRIZECOLLECTINGPATH_H
#define PRIZECOLLECTINGPATH_H
#include<string>
#include<algorithm>
#include<iostream>
#include<vector>
#include<lemon/list_graph.h>
#include<lemon/graph_to_eps.h>
#include<gurobi_c++.h>
#include"mygraphlib.h"
using namespace lemon;
int read_pcpath(std::string);
int show_usage();
int show_input();
int make_eps_graph(ListBpGraph::ArcMap<int>& color, std::string name);
bool is_feasible_solution(int& cost, int verbose);
void show_graph_mygraphlib(std::string text);
int prize_collecting_st_path_pli(ListDigraph& g, ListDigraph::NodeMap<double> &prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax);
int prize_collecting_st_path_heuristic(ListDigraph& g, ListDigraph::NodeMap<double> &prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax);
#endif
