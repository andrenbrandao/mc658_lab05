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

class ConnectivityCuts: public GRBCallback
{
    ListDigraph& g;
    ListDigraph::NodeMap<double>& prize;
    ListDigraph::ArcMap<double> &cost;
    ListDigraph::Node s;
    ListDigraph::Node t;
    std::vector<ListDigraph::Node> &path;
    ListDigraph::ArcMap<GRBVar>& x;
    double (GRBCallback::*solution_value)(GRBVar);
  public:
    ConnectivityCuts(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double>& cost, ListDigraph::Node s,
      ListDigraph::Node t, std::vector<ListDigraph::Node> &path,
      ListDigraph::ArcMap<GRBVar>& x) : g(g), prize(prize), cost(cost), s(s), t(t), path(path), x(x)
    {    }
  protected:
    void callback()
    {
      if (where==GRB_CB_MIPSOL){ solution_value = &ConnectivityCuts::getSolution;}
      else if (where==GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS)==GRB_OPTIMAL) {
        solution_value = &ConnectivityCuts::getNodeRel;
      } else return;
      try {
        ArcValueMap capacity(g);
        DCutMap cut(g);
        double vcut;
        for (ArcIt a(g); a!=INVALID; ++a)
          capacity[a] = (this->*solution_value)(x[a]);  // or getSolution(x[a]);

        for (DNodeIt v(g); v!=INVALID; ++v) {
          if(v != s && v!= t) {
            GRBLinExpr expr;
            // find a mincut between root s and other terminal
            vcut = 1.0;
            for (ArcIt a(g); a!=INVALID; ++a) {
              if(g.target(a) == v && getSolution(x[a]) >= 1.0-MY_EPS ) {
                vcut = DiMinCut(g,capacity, s , v, cut);
                break;
              }
            }
            if (vcut >= 1.0-MY_EPS) continue;

            // found violated cut
            for (ArcIt a(g); a!=INVALID; ++a)
             if ((cut[g.source(a)]==cut[s]) && (cut[g.target(a)]!=cut[s]))
               expr += x[a];
            // for(ArcIt a(g); a!=INVALID; ++a) {
              // if(x[a] >= 1.0-MY_EPS )
                addLazy( expr >= 1.0 ); // or addLazy(expr,GRB_GREATER_EQUAL,1.0);
            // }
          }
        }
  } catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during callback**" << endl;
  }
  }
};

void find_other_path(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double>& cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax, ListDigraph::ArcMap<GRBVar>& x) {
  ListDigraph::Node node, last_node, new_node;

  for (ArcIt e(g); e!=INVALID; ++e) {
    if(BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
      if(find(path.begin(), path.end(), g.source(e)) != path.end()) {

      }
      else {
        node = g.source(e);
        break;
      }
    }
  }

  cout << g.id(node) << endl;
  last_node = node;
        while(last_node != INVALID) {
          for ( OutArcIt e(g, last_node); e!=INVALID; ++e ){
            if(BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
              new_node = g.target(e);
              cout << g.id(new_node) << endl;
              last_node = new_node;
              // opt += prize[new_node] - cost[e];
            }
          }
        }
}

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
  double opt = 0.0;

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
    // prize_collecting_st_path_heuristic(g, prize, cost, s, t, path, LB, UB, tMax);
    // if (LB > 0) {
    //   cout << "Setting cutoff " << LB << endl;
    //   model.set(GRB_DoubleParam_Cutoff, LB );
    // }

    model.set(GRB_IntParam_LazyConstraints, 1);
    ConnectivityCuts cb = ConnectivityCuts(g, prize, cost, s, t, path, x);
    model.setCallback(&cb);

    model.update();
    model.optimize();

    cout << "LB: " << LB << endl;
    cout << "OPTIMAL PLI: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    int optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_OPTIMAL) {
      if(LB < model.get(GRB_DoubleAttr_ObjVal)) {
        // cria o caminho de s a t
        path.clear();
        path.push_back(s);
        opt += prize[s];
        last_node = s;
        while(last_node != t) {
          for ( OutArcIt e(g, last_node); e!=INVALID; ++e ){
            if(BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
              new_node = g.target(e);
              path.push_back(new_node);
              last_node = new_node;
              opt += prize[new_node] - cost[e];
            }
          }
        }
        cout << "SIZE S->T: " << path.size() << endl;

          int count = 0;

          for ( ArcIt e(g); e!=INVALID; ++e ){
            if(BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
              cout << g.id(g.source(e)) << " ";
              count+=1;
            }
          }
          cout << "NUMBER VERTICES: " << count << endl;

      }

      // find_other_path(g, prize, cost, s, t, path, LB, UB, tMax, x);
      cout << "OPT FOUND IN PATH S->T: " << opt << endl;
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
  double min_dist = 0.0;
  double opt = 0.0;
  int n_arcs = 0;

  for(ArcIt e(g); e!=INVALID; ++e){
    distance[e] = cost[e] - prize[g.target(e)];
    // else distance[e] = cost[e];

    if(distance[e] < min_dist)
      min_dist = distance[e];
  }

  // for(ArcIt e(g); e!=INVALID; ++e){
  //   if(min_dist < 0)
  //     distance[e] -= min_dist;
  // }

  Dijkstra<ListDigraph, ListDigraph::ArcMap<double>> dijkstra(g, distance);
  dijkstra.run(s);

  cout << "Dijkstra | Distancia de s a t: "
            << dijkstra.dist(t) << endl;

  for (ListDigraph::Node v=t;v != s; v=dijkstra.predNode(v)) {
    n_arcs++;
    path.insert(path.begin(), v);
  }

  path.insert(path.begin(), s);

  // if(min_dist < 0)
  //   LB =  dijkstra.dist(t) - n_arcs*min_dist + prize[s];
  // else
    LB =  -1 * (dijkstra.dist(t) - prize[s]);

  return 0;
}
