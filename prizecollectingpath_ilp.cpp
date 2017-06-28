/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "prizecollectingpath.h"
#include <lemon/dijkstra.h>

///Preencher aqui para facilitar a correcao.
// Nome1: André Nogueira Brandão
// RA1: 116130
// Nome2: Felipe Matias Camargo
// RA2: 180878

///
// PLI function
///
int prize_collecting_st_path_pli(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double>& cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax){
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  GRBLinExpr expr;
  model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
  model.set(GRB_DoubleParam_TimeLimit, tMax);
  ListDigraph::ArcMap<GRBVar> x(g); // variavel x indica se a aresta esta ou nao na solucao
  ListDigraph::Node last_node, new_node; // nos utilizados para calcular o path
  double prize_value; // guarda o valor do premio
  double cutoff=0.0;

  // Adiciona uma variavel x para cada aresta
  // Nao utiliza os valores dos premios de s e de t
  for (ArcIt e(g); e!=INVALID; ++e) {
    prize_value = prize[g.source(e)];
    x[e] = model.addVar(0.0, 1.0, prize_value - 1.0 * cost[e],GRB_BINARY,"");
  }
  model.addVar(1.0, 1.0, prize[t],GRB_BINARY,"");
  model.update();

  // para cada no diferente de s e t
  // verifique que o numero que entra eh igual ao que sai
  for ( DNodeIt v(g); v!=INVALID; ++v ) {
    if(v != s && v != t) {
      GRBLinExpr in_edges = 0;
      GRBLinExpr out_edges = 0;

      for ( InArcIt e(g, v); e!=INVALID; ++e ){
          in_edges += x[e];
      }
       for ( OutArcIt e(g, v); e!=INVALID; ++e ){
          out_edges += x[e];
      }
      model.addConstr(in_edges - out_edges == 0);
    }
  }

  // Para o vertice s: nenhuma aresta incide no no e exatamente uma aresta sai do no
  GRBLinExpr in_edges = 0;
  GRBLinExpr out_edges = 0;
  for ( InArcIt e(g, s); e!=INVALID; ++e ){
    in_edges += x[e];
  }
  for ( OutArcIt e(g, s); e!=INVALID; ++e ){
    out_edges += x[e];
  }
  model.addConstr(in_edges == 0);
  model.addConstr(out_edges == 1);

  // Para o vertice t: uma aresta incide no no e nenhuma aresta sai dele
  in_edges = 0;
  out_edges = 0;
  for ( InArcIt e(g, t); e!=INVALID; ++e ){
    in_edges += x[e];
  }
  for ( OutArcIt e(g, t); e!=INVALID; ++e ){
    out_edges += x[e];
  }
  model.addConstr(in_edges == 1);
  model.addConstr(out_edges == 0);

  // para cada no diferente de s e t
  // ate uma aresta sai e ate uma aresta entra
  for ( DNodeIt v(g); v!=INVALID; ++v ) {
    if(v != s && v != t) {
      GRBLinExpr in_edges = 0;
      GRBLinExpr out_edges = 0;

      for ( InArcIt e(g, v); e!=INVALID; ++e ){
          in_edges += x[e];
      }
       for ( OutArcIt e(g, v); e!=INVALID; ++e ){
          out_edges += x[e];
      }
      model.addConstr(in_edges <= 1);
      model.addConstr(out_edges <= 1);
    }
  }

  try {
    cout << "Before DIJKSTRA" << endl;
    cutoff = prize_collecting_st_path_heuristic(g, prize, cost, s, t, path, LB, UB, tMax);
    if (cutoff > 0) {
      cout << "Setting cutoff " << cutoff << endl;
      model.set(GRB_DoubleParam_Cutoff, cutoff );
    }

    cout<< "TIME LIMIT: " <<model.get(GRB_DoubleParam_TimeLimit) << endl;

    model.update();
    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_OPTIMAL) {
      if(cutoff < model.get(GRB_DoubleAttr_ObjVal)) {
        // cria o caminho de s a t
        path.clear();
        path.push_back(s);
        last_node = s;
        while(last_node != t) {
          for ( ArcIt e(g); e!=INVALID; ++e ){
            if(g.source(e) == last_node) {
              if(BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
                new_node = g.target(e);
                path.push_back(new_node);
                last_node = new_node;
              }
            }
          }
        }
      }
      return 1;
    }
    else {
      cout << "Heuristic completed" << endl;
      return 2;
    }

  }
  catch(GRBException e) {
    std::cerr << "Codigo de erro = " << e.getErrorCode() << std::endl;
    std::cerr << e.getMessage();
    return 0;
  }

  return true;
}


///
//
// Heuristic function
///
int prize_collecting_st_path_heuristic(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax){
    ListDigraph::ArcMap<double> distance(g);

	int minDistance = 0;

	for(ArcIt e(g); e!=INVALID; ++e){
      if(cost[e] - prize[g.target(e)] < minDistance) minDistance = cost[e] - prize[g.target(e)];
    }

    for(ArcIt e(g); e!=INVALID; ++e){
      if(g.target(e) != t) distance[e] = cost[e] - prize[g.target(e)] - minDistance;
      else distance[e] = cost[e];
    }

    Dijkstra<ListDigraph, ListDigraph::ArcMap<double>> dijkstra(g, distance);
    dijkstra.run(s);

    for (ListDigraph::Node v=t;v != s; v=dijkstra.predNode(v)) {
      path.insert(path.begin(), v);
    }
    
    cout << "The distance of node t from node s: "
              << -1 * dijkstra.dist(t) + minDistance*(path.size()-1) << endl;

    path.insert(path.begin(), s);

    return  -1 * dijkstra.dist(t)+ minDistance*(path.size()-1);

	// return 0;
}
