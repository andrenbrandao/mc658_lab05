/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "prizecollectingpath.h"
#include <typeinfo>
#include <list>
#include <set>
#include <vector>
#include <cstdio>

///Preencher aqui para facilitar a correcao.
// Nome1: André Nogueira Brandão
// RA1: 116130
// Nome2: Felipe Matias Camargo
// RA2: 180878

using namespace std;
# define INF 0x3f3f3f3f

// A classe Graph representa um grafo orientado
class Graph
{
	int V;  // N de vertices do grafo

	list< pair<int, double> > *adj; // lista de adjacencia do grafo

	list<pair<int, double>> *profit;

	double minWeight;

public:
	Graph(int V); //Construtor

	void addEdge(int u, int v, double w, double p); //Adicionar aresta ao grafo

	int shortestPath(int s, int t); // Calcular caminho mínimo de s a t
};

Graph::Graph(int V)
{
	this->V = V;
	adj = new list< pair<int, double> >[V];
	profit = new list<pair<int, double>>[V];
	minWeight = 0;
}

void Graph::addEdge(int u, int v, double w, double p)
{
	adj[u].push_back(make_pair(v, w-p));

	profit[u].push_back(make_pair(v, p - w));

	if (w - p < minWeight)
	{
		minWeight = w - p;
	}
}

/**Time Complexity : Set in C++ are typically implemented using Self-balancing binary search trees. 
Therefore, time complexity of set operations like insert, delete is logarithmic and time complexity of above solution is O(ELogV)).
*/
int Graph::shortestPath(int s, int t)
{
	if (minWeight<0)
	{
		for (int u = 0; u < V; ++u)
		{
			list< pair<int, double> >::iterator i;
			for (i = adj[u].begin(); i != adj[u].end(); ++i)
			{
				(*i).second += -minWeight;
			}
		}
	}

	// O set sera usado como fila de prioridade do dijkstra
	set<pair<double,int>> queue;

	// Cria um vetor de distancias para cada um dos vertices, iniciado elas com INF
	vector<double> dist(V, INF);

	// Cria um vetor para guardar o pai de cada um dos vertices
	vector<int> parent(V, -1);

	// Insere s ao set, e a distancia de s até s é igual a 0
	queue.insert(make_pair(0, s));
	dist[s] = 0;

	/* Itera até todas as menores distancias sejam achadas, no momento em que a fila de prioridade estaja vazia */
	while (!queue.empty())
	{
		// O primeiro vertice do set será o com a menor distacia até s que ainda nao foi finalizado
		pair<double, int> tmp = *(queue.begin());
		queue.erase(queue.begin());

		// O primeiro valor do par deve ser usado para a distancia e a segunda para o vertice, para a ordenação correta da fila
		int u = tmp.second;

		// percorre i, onde i é um vertice adjacente ao vértice u
		list< pair<int, double> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i)
		{
			// pega o vertice adjacente e sua distância até u
			int v = (*i).first;
			double weight = (*i).second;

			//  Checa se tem um caminho mais proximo até v passando por u
			if (dist[v] > dist[u] + weight)
			{
				if (dist[v] != INF)
					queue.erase(queue.find(make_pair(dist[v], v)));

				// Atualiza a distância até v
				dist[v] = dist[u] + weight;
				parent[v] = u;
				queue.insert(make_pair(dist[v], v));
			}
		}
	}

	int u = t;
	int result = 0;
	while (u != s)
	{
		// percorre i, onde i é um vertice adjacente ao vértice u
		list< pair<int, double> >::iterator i;
		for (i = profit[parent[u]].begin(); i != profit[parent[u]].end(); ++i)
		{
			if ((*i).first == u)
			{
				result += (*i).second;
				break;
			}
		}
		u = parent[u];
	}

	return result;
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
  double opt = 0.0; // utilizado para calcular o otimo

  // Adiciona uma variavel x para cada aresta
  // Nao utiliza os valores dos premios de s e de t
  for (ArcIt e(g); e!=INVALID; ++e) {
    prize_value = prize[g.source(e)];
    if(g.source(e) == s) {
      prize_value = 0.0;
    }
    x[e] = model.addVar(0.0, 1.0, prize_value - 1.0 * cost[e],GRB_BINARY,"");
  }
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
    model.update();
    model.optimize();

    opt += prize[s];

    // cria o caminho de s a t
    path.push_back(s);
    last_node = s;
    while(last_node != t) {
      for ( ArcIt e(g); e!=INVALID; ++e ){
        if(g.source(e) == last_node) {
          new_node = g.target(e);
          if(BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
            opt += prize[new_node];
            path.push_back(new_node);
            last_node = new_node;
          }
        }
      }
    }

    cout << "OPT: " << opt << endl;

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
  
  	int V = 9; //Mudar para o tamanho do grafo
	Graph graph(V);
  
	//for (ArcIt e(g); e!=INVALID; ++e) 
	//{
		//prize_value = prize[g.source(e)];
		//if(g.source(e) == s) 
		//{
		  //prize_value = 0.0;
		//}
		//g.addEdge(g.id(g.source(e)),g.id(g.target(e)),1.0 * cost[e],prize_value);
	 //}
	 
	 //int lowerBound = g.shortestPath(g.id(s),g.id(t));
  
	return 0;
}
