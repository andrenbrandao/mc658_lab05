/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "prizecollectingpath.h"
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>

///Preencher aqui para facilitar a correcao.
// Nome1: André Nogueira Brandão
// RA1: 116130
// Nome2: Felipe Matias Camargo
// RA2: 180878

class ConnectivityCuts: public GRBCallback
{
    ListDigraph& g;
    ListDigraph::Node s;
    ListDigraph::Node t;
    ListDigraph::ArcMap<GRBVar>& xa;
    ListDigraph::NodeMap<GRBVar>& xv;
    double (GRBCallback::*solution_value)(GRBVar);
  public:
    ConnectivityCuts(ListDigraph& g, ListDigraph::Node s, ListDigraph::Node t,
      ListDigraph::ArcMap<GRBVar>& xa,ListDigraph::NodeMap<GRBVar>& xv) : g(g), s(s), t(t), xa(xa), xv(xv)
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
          capacity[a] = (this->*solution_value)(xa[a]);  // or getSolution(x[a]);

        for (DNodeIt v(g); v!=INVALID; ++v) {
          if(v != s && v!= t) {
            GRBLinExpr expr;
            // find a mincut between root s and other terminal
            vcut = 1.0;
            if(getSolution(xv[v]) >= 1.0-MY_EPS)
              vcut = DiMinCut(g,capacity, s , v, cut);
            if (vcut >= 1.0-MY_EPS) continue;

            int count=0;
            // found violated cut
            for (ArcIt a(g); a!=INVALID; ++a)
             if ((cut[g.source(a)]==cut[s]) && (cut[g.target(a)]!=cut[s])) {
                if(getSolution(xv[g.target(a)]) >= 1.0-MY_EPS) {
                  count++;
                  expr += xa[a];
                }
             }
             if(count > 0)
              addLazy( expr >= 1.0  ); // or addLazy(expr,GRB_GREATER_EQUAL,1.0);
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

double greedy_opt;

///
// PLI function
///
int prize_collecting_st_path_pli(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double>& cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax){
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  GRBLinExpr expr;
  model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
  model.set(GRB_DoubleParam_TimeLimit, tMax);
  ListDigraph::ArcMap<GRBVar> xa(g); // variavel x indica se a aresta esta ou nao na solucao
  ListDigraph::NodeMap<GRBVar> xv(g); // variavel x indica se a aresta esta ou nao na solucao
  ListDigraph::Node last_node, new_node; // nos utilizados para calcular o path

  // Adiciona uma variavel x para cada aresta
  // Nao utiliza os valores dos premios de s e de t
  for (ArcIt e(g); e!=INVALID; ++e) {
    xa[e] = model.addVar(0.0, 1.0, -1.0 * cost[e],GRB_BINARY,"");
  }
  for(DNodeIt v(g); v!=INVALID; ++v) {
    xv[v] = model.addVar(0.0, 1.0, prize[v], GRB_BINARY, "");
  }
  model.update();

  // para cada no diferente de s e t
  // verifique que o numero que entra eh igual ao que sai
  for ( DNodeIt v(g); v!=INVALID; ++v ) {
    if(v != s && v != t) {
      GRBLinExpr in_edges = 0;
      GRBLinExpr out_edges = 0;

      for ( InArcIt e(g, v); e!=INVALID; ++e ){
          in_edges += xa[e];
      }
       for ( OutArcIt e(g, v); e!=INVALID; ++e ){
          out_edges += xa[e];
      }
      model.addConstr(in_edges - out_edges == 0);
    }
  }

  // Para o vertice s: nenhuma aresta incide no no e exatamente uma aresta sai do no
  GRBLinExpr in_edges = 0;
  GRBLinExpr out_edges = 0;
  for ( InArcIt e(g, s); e!=INVALID; ++e ){
    in_edges += xa[e];
  }
  for ( OutArcIt e(g, s); e!=INVALID; ++e ){
    out_edges += xa[e];
  }
  model.addConstr(in_edges == 0);
  model.addConstr(out_edges == 1);
  model.addConstr(xv[s] <= in_edges + out_edges);

  // Para o vertice t: uma aresta incide no no e nenhuma aresta sai dele
  in_edges = 0;
  out_edges = 0;
  for ( InArcIt e(g, t); e!=INVALID; ++e ){
    in_edges += xa[e];
  }
  for ( OutArcIt e(g, t); e!=INVALID; ++e ){
    out_edges += xa[e];
  }
  model.addConstr(in_edges == 1);
  model.addConstr(out_edges == 0);
  model.addConstr(xv[t] <= in_edges + out_edges);

  // para cada no diferente de s e t
  // ate uma aresta sai e ate uma aresta entra
  for ( DNodeIt v(g); v!=INVALID; ++v ) {
    if(v != s && v != t) {
      GRBLinExpr in_edges = 0;
      GRBLinExpr out_edges = 0;

      for ( InArcIt e(g, v); e!=INVALID; ++e ){
          in_edges += xa[e];
      }
       for ( OutArcIt e(g, v); e!=INVALID; ++e ){
          out_edges += xa[e];
      }
      model.addConstr(in_edges <= 1);
      model.addConstr(out_edges <= 1);
      model.addConstr(xv[v] <= in_edges + out_edges);
    }
  }

  try {
    prize_collecting_st_path_heuristic(g, prize, cost, s, t, path, LB, UB, tMax);
    if (LB > 0) {
      cout << "Setting cutoff " << LB << endl;
      model.set(GRB_DoubleParam_Cutoff, LB );
    }

    model.set(GRB_IntParam_LazyConstraints, 1);
    ConnectivityCuts cb = ConnectivityCuts(g, s, t, xa, xv);
    model.setCallback(&cb);

    model.update();
    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_OPTIMAL) {
      if(LB < model.get(GRB_DoubleAttr_ObjVal)) {
        // cria o caminho de s a t
        path.clear();
        path.push_back(s);
        last_node = s;
        while(last_node != t) {
          for ( OutArcIt e(g, last_node); e!=INVALID; ++e ){
            if(BinaryIsOne(xa[e].get(GRB_DoubleAttr_X))) {
              new_node = g.target(e);
              path.push_back(new_node);
              last_node = new_node;
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

  return 1;
}

//Greedy Alg
int greedy_dijkstra(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, std::vector<ListDigraph::Arc> blacklist){
    ListDigraph::ArcMap<double> distance(g);
    int n_arcs = 0;
	double minDistance = 0.0;

	for(ArcIt e(g); e!=INVALID; ++e){
      if(cost[e] - prize[g.target(e)] < minDistance) minDistance = cost[e] - prize[g.target(e)];
    }

    for(ArcIt e(g); e!=INVALID; ++e){
      if(g.target(e) != t)
      {
		  distance[e] = cost[e] - prize[g.target(e)] - minDistance;
	  }
      else 
      {
		  distance[e] = cost[e];
	  }
    }
    
    //Gives blacklisted edges an very big weight
    for(unsigned int i = 0; i<blacklist.size(); i++)
    {
		distance[blacklist[i]] = 5000000;
	}

    Dijkstra<ListDigraph, ListDigraph::ArcMap<double>> dijkstra(g, distance);
    dijkstra.run(s);
    

    for (ListDigraph::Node v=t;v != s; v=dijkstra.predNode(v)) {
	  n_arcs++;
      path.insert(path.begin(), v);
    }
    
    path.insert(path.begin(), s);

    return  -1 *dijkstra.dist(t) - n_arcs*minDistance + prize[s];
}

int prize_collecting_st_path_heuristic(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax, std::vector<ListDigraph::Arc> blacklist){
  ListDigraph::ArcMap<double> distance(g);
  std::vector<ListDigraph::Node> temp_path;
  double min_dist = 0.0;
  int n_arcs = 0;

  for(ArcIt e(g); e!=INVALID; ++e)
  {
    if(g.target(e) != t)
      distance[e] = cost[e] - prize[g.target(e)];
    else distance[e] = cost[e];

    if(distance[e] < min_dist)
      min_dist = distance[e];
  }

  for(ArcIt e(g); e!=INVALID; ++e)
  {
    if(min_dist < 0)
      distance[e] -= min_dist;
  }
  for(unsigned int i = 0; i<blacklist.size(); i++)
  {
    distance[blacklist[i]] = 50000000;
  }


  Dijkstra<ListDigraph, ListDigraph::ArcMap<double>> dijkstra(g, distance);
  dijkstra.run(s);
 
  std::vector<ListDigraph::Arc> arc_list;
  for (ListDigraph::Node v=t;v != s; v=dijkstra.predNode(v))
  {
    n_arcs++;
    path.insert(path.begin(), v);
    ListDigraph::Arc prev = dijkstra.predArc(v);
    arc_list.push_back(prev);
  } 
  path.insert(path.begin(), s);
    //cout<< n_arcs;

  if(min_dist < 0)
    LB =  dijkstra.dist(t) - n_arcs*min_dist + prize[s];
  else
    LB =  -1 * (dijkstra.dist(t) - prize[s]);
    
  //Explores neightboors solutions
  for(unsigned int i = 0; i < arc_list.size(); i++)
  {
	  //Blacklists a path of the current solution and tries to found better neightboors solutions
	  blacklist.push_back(arc_list[i]);
	  
	  int local_opt = greedy_dijkstra(g, prize, cost, s, t, temp_path, blacklist);
		//cout << local_opt<<endl;
		if(local_opt > LB)
		{
			LB = local_opt;
			//Tries to explore neighboors from best found node
			double exploring_opt = prize_collecting_st_path_heuristic(g, prize, cost, s, t, temp_path, LB, UB, tMax, blacklist);
			cout << exploring_opt<<endl;
			if(exploring_opt > LB)
			{
				LB = exploring_opt;
				path =  temp_path;
			}
		}
  }
  
  return 2;
}

///
//
// Heuristic function
///
int prize_collecting_st_path_heuristic(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax){
  std::vector<ListDigraph::Arc> blacklist;
  prize_collecting_st_path_heuristic(g, prize, cost, s, t, path, LB, UB, tMax, blacklist);
  return 2;
}



