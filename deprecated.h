// =============================================================
#ifndef MYGRAPHLIBDEPRECATED_DEFINE
#define MYGRAPHLIBDEPRECATED_DEFINE


//---------------------------------------------------------------------
// Avoid using these types

typedef ListGraph::EdgeMap<double> EdgeWeight; // Use EdgeValueMap  instead
typedef ListGraph::EdgeMap<string> EdgeName;   // Use EdgeStringMap instead
typedef ListGraph::NodeMap<double> NodeWeight; // Use NodeValueMap  instead
typedef ListGraph::NodeMap<string> NodeName;   // Use NodeStringMap instead
typedef Digraph::ArcMap<double> ArcWeight;     // Use ArcValueMap   instead
typedef Digraph::ArcMap<string> ArcName;       // Use ArcStringMap
typedef Digraph::NodeMap<string> DiNodeName;   // Use DiNodeStringMap

typedef Digraph::Node DiNode;
typedef Digraph::NodeIt DiNodeIt;


#endif
