/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "prizecollectingpath.h"
#include<fstream>
#include<iostream>
#include<unistd.h>
#include<cstring>
#include<map>
using namespace lemon;
int n, m;
int tMax;     //Maximum execution time, in seconds
ListDigraph g;      //Graph
ListDigraph::NodeMap<double> prizes(g);
ListDigraph::ArcMap<double> costs(g);
ListDigraph::NodeMap<std::string> node_names(g);
ListDigraph::NodeMap<double> px(g), py(g);//For drawing the graph
ListDigraph::Node s;
ListDigraph::Node t;
std::vector<ListDigraph::Node> path;
double UB, LB;      //UpperBound and LowerBound
int main(int argc, char* argv[]){
  tMax = 600;
  int verbose = 0;
  std::string input_file;
  int option, i_f = 0;
  if ( argc == 1 ){
    show_usage();
    return 0;
  }
  while ((option = getopt(argc, argv, "t:i:v"))!=-1)
    switch(option){
      case 't':
        tMax=atoi(optarg);
        break;
      case 'i':
        i_f = 1;
        input_file.assign(optarg);
        break;
      case 'v':
        verbose=1;
        break;
      default:
        break;
    }

  if ( i_f == 0 ){
    std::cout << "-i  mandatory argument" << std::endl;
    return 1;
  }

  if ( !read_pcpath(input_file) ) return 1;

  for ( ListDigraph::NodeIt u(g); u!=INVALID; ++u ){
    if ( node_names[u].compare("s")==0 )
      s=u;
    if ( node_names[u].compare("t")==0 )
      t=u;
  }

  //if ( verbose )
  //  show_input();

  prize_collecting_st_path_pli(g, prizes, costs, s, t, path, UB, LB, tMax);
  if ( verbose ){
    //make_eps_graph(path, "sol");
    //set_pdfreader("okular");
    //set_pdfreader("open");
    //set_pdfreader("xpdf");
    show_graph_mygraphlib(input_file);
  }

  for ( int i=0; i<(int)path.size(); i++ )
    std::cout << g.id(path[i]) << " ";
  std::cout << std::endl;
  return 0;
}

int read_pcpath(string input_file){
  std::ifstream kinput;
  kinput.open(input_file.c_str()); if (!kinput) return 0;
  kinput >> n >> m;
  std::map<std::string, ListDigraph::Node> ref;

  for ( int i=0; i<n; i++ ){
    string tmp;
    double r;
    kinput >> tmp >> r;
    ListDigraph::Node n = g.addNode();
    node_names[n] = tmp;
    prizes[n] = r;
    ref[tmp]=n;
  }

  for ( int i=0; i<m; i++){
    string v1, v2;
    double c_tmp;
    kinput >> v1 >> v2 >> c_tmp;
    ListDigraph::Arc a = g.addArc(ref[v1], ref[v2]);  //source, target
    costs[a] = c_tmp;
  }

  return 1;
}

int show_usage(){
  std::cout << "Uso: ./prize_collecting_path.e -i <input_file> -t <max_time>"  << std::endl
       << "input_time eh o arquivo de entrada e max_time eh o tempo maximo de execucao." << std::endl
       << "voce pode usar o argumento -v para imprimir a entrada e debug. " << std::endl;
  return 0;
}

int show_input(){
  ListDigraph::ArcMap<int> e_color(g,0);
  std::string fname("prize_collecting_path");
  //make_eps_graph(e_color, fname);
  return 0;
}

void show_graph_mygraphlib(std::string text){
  ListDigraph::NodeMap<std::string> v_name(g);
  ListDigraph::ArcMap<std::string> e_name(g);
  ListDigraph::NodeMap<int> v_color(g);
  ListDigraph::ArcMap<int> e_color(g);
  for ( ListDigraph::NodeIt nb(g); nb!=INVALID; ++nb ) {
    if ( node_names[nb].compare("s") == 0 ){
      v_color[nb] = 3;
    } else if ( node_names[nb].compare("t") == 0 ) {
      v_color[nb] = 3;
    } else if ( std::find(path.begin(), path.end(), nb) != path.end() )
      v_color[nb] = 4;
    else
      v_color[nb]=5;

    v_name[nb] = node_names[nb] + "-" + std::to_string(prizes[nb]);
  }
  for ( ListDigraph::ArcIt e(g); e!=INVALID; ++e ){
    e_name[e] = to_string(costs[e]);
    e_color[e] = 2;
  }

  ViewListDigraph(g, node_names, e_name, v_color, e_color, text);
}

