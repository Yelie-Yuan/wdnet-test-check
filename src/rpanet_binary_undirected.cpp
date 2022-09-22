#include<iostream>
#include<queue>
#include<R.h>
#include "funcPtrUnd.h"
#include<Rcpp.h>
using namespace std;
funcPtrUnd prefFuncCpp;

/**
 * Node structure in undirected networks.
 * id: node id
 * strength: node strength
 * p: preference of being chosen from the existing nodes
 * totalp: sum of preference of current node and its children
 * *left, *right, *parent: pointers to its left, right and parent
 */
struct node_und {
  int id;
  double strength;
  double p, totalp;
  node_und *left, *right, *parent;
};

/**
 * Default preference function.
 *
 * @param strength Node strength.
 * @param params Parameters passed to the preference function.
 * 
 * @return Preference of a node.
 */
double prefFuncDefault(double strength, Rcpp::NumericVector params) {
  return pow(strength, params[0]) + params[1];
}

/**
 * Update total preference from current node to root.
 *
 * @param current_node The current node.
 */
void updateTotalp(node_und *current_node) {
  if (current_node->left == NULL) {
    current_node->totalp = current_node->p;
  }
  else if (current_node->right == NULL) {
    current_node->totalp = current_node->p + current_node->left->totalp;
  }
  else {
    current_node->totalp = current_node->p + current_node->left->totalp + current_node->right->totalp;
  }
  while(current_node->id > 0) {
    return updateTotalp(current_node->parent);
  }
}

/**
 * Update node preference and total preference from the sampled node to root.
 *
 * @param temp_node The sampled/new node.
 * @param func_type Default or customized preference function.
 * @param params Parameters passed to the default preference function.
 * @param prefFuncCpp Pointer of the customized preference function.

 */
void updatePrefUnd(node_und *temp_node, int func_type, 
                         Rcpp::NumericVector params,
                         funcPtrUnd prefFuncCpp) {
  if (func_type == 1) {
    temp_node->p = prefFuncDefault(temp_node->strength, params);
  }
  else {
    temp_node->p = prefFuncCpp(temp_node->strength);
  }
  updateTotalp(temp_node);
}

/**
 * Create a new node.
 *
 * @param id Node ID.
 * 
 * @return The new node.
 */
node_und *createNodeUnd(int id) {
  node_und *new_node = new node_und();
  new_node->id = id;
  new_node->strength = 0;
  new_node->p = new_node->totalp = 0;
  new_node->left = new_node->right = new_node->parent = NULL;
  return new_node;
}

/**
 * Insert a new node to the tree.
 *
 * @param q Sequence of nodes that have less than 2 children.
 * @param id New node ID.
 * 
 * @return The new node.
 */
node_und *insertNodeUnd(queue<node_und*> &q, int id) {
  node_und *new_node = createNodeUnd(id);
  node_und *temp_node = q.front();
  // check left
  if(temp_node->left == NULL) {
    temp_node->left = new_node;
  }
  // check right
  else if (temp_node->right == NULL) {
    temp_node->right = new_node;
    q.pop();
  }
  new_node->parent = temp_node;
  q.push(new_node);
  return new_node;
}

/**
 * Find a node with a given cutoff point w.
 *
 * @param root Root node of the tree.
 * @param w A cutoff point.
 * 
 * @return Sampled node.
 */
node_und *findNode(node_und *root, double w) {
  if (w > root->totalp) {
    // numerical error
    // Rprintf("Numerical error. Diff %f.\n", (w - root->totalp) * pow(10, 10));
    w = root->totalp;
  }
  w -= root->p;
  if (w <= 0) {
    return root;
  }
  else {
    if (w > root->left->totalp) {
      return findNode(root->right, w - root->left->totalp);
    }
    else {
      return findNode(root->left, w);
    }
  }
}

/**
 * Sample a node from the tree.
 *
 * @param root Root node of the tree.
 * 
 * @return Sampled node.
 */
node_und *sampleNodeUnd(node_und *root) {
  double w;
  w = 1;
  while (w == 1) {
    w = unif_rand();
  }
  w *= root->totalp;
  return findNode(root, w);
}

//' Preferential attachment algorithm.
//'
//' @param nstep Number of steps.
//' @param m Number of new edges in each step.
//' @param new_node_id New node ID.
//' @param new_edge_id New edge ID.
//' @param node_vec1 Sequence of nodes in the first column of edgelist.
//' @param node_vec2 Sequence of nodes in the second column of edgelist.
//' @param strength Sequence of node strength.
//' @param edgeweight Weight of existing and new edges.
//' @param scenario Scenario of existing and new edges.
//' @param pref Sequence of node preference.
//' @param control List of controlling arguments.
//' @return Sampled network.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List rpanet_binary_undirected_cpp(
    int nstep, 
    Rcpp::IntegerVector m,
    int new_node_id, 
    int new_edge_id, 
    Rcpp::IntegerVector node_vec1, 
    Rcpp::IntegerVector node_vec2, 
    Rcpp::NumericVector strength, 
    Rcpp::NumericVector edgeweight, 
    Rcpp::IntegerVector scenario,
    Rcpp::NumericVector pref, 
    Rcpp::List control) {
  Rcpp::List scenario_ctl = control["scenario"];
  double alpha = scenario_ctl["alpha"];
  double beta = scenario_ctl["beta"];
  double gamma = scenario_ctl["gamma"];
  double xi = scenario_ctl["xi"];
  bool beta_loop = scenario_ctl["beta.loop"];
  Rcpp::List newedge_ctl = control["newedge"];
  bool node_unique = ! newedge_ctl["node.replace"];
  Rcpp::List preference_ctl = control["preference"];
  Rcpp::NumericVector params(2);
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type) {
  case 1: 
    params = preference_ctl["params"];
    break;
  case 2: {
      SEXP pref_func_ptr = preference_ctl["pref.pointer"];
      prefFuncCpp = *Rcpp::XPtr<funcPtrUnd>(pref_func_ptr);
      break;
    }
  }
  
  double u, temp_p;
  bool m_error;
  int i, j, n_existing, current_scenario;
  node_und *node1, *node2;
  // initialize a tree from seed graph
  node_und *root = createNodeUnd(0);
  root->strength = strength[0];
  updatePrefUnd(root, func_type, params, prefFuncCpp);
  queue<node_und*> q, q1;
  q.push(root);
  for (i = 1; i < new_node_id; i++) {
    node1 = insertNodeUnd(q, i);
    node1->strength = strength[i];
    updatePrefUnd(node1, func_type, params, prefFuncCpp);
  }
  // sample edges
  GetRNGstate();
  for (i = 0; i < nstep; i++) {
    m_error = false;
    n_existing = new_node_id;
    for (j = 0; j < m[i]; j++) {
      u = unif_rand();
      if (u <= alpha) {
        current_scenario = 1;
      }
      else if (u <= alpha + beta) {
        current_scenario = 2;
      }
      else if (u <= alpha + beta + gamma) {
        current_scenario = 3;
      }
      else if (u <= alpha + beta + gamma + xi) {
        current_scenario = 4;
      }
      else {
        current_scenario = 5;
      }
      switch (current_scenario) {
      case 1:
        if (root->totalp == 0) {
          m_error = true;
          break;
        }
        node1 = insertNodeUnd(q, new_node_id);
        new_node_id++;
        node2 = sampleNodeUnd(root);
        break;
      case 2:
        if (root->totalp == 0) {
          m_error = true;
          break;
        }
        node1 = sampleNodeUnd(root);
        if (! beta_loop) {
          if (node1->p == root->totalp) {
            m_error = true;
            break;
          }
          else {
            temp_p = node1->p;
            node1->p = 0;
            updateTotalp(node1);
            node2 = sampleNodeUnd(root);
            node1->p = temp_p;
            updateTotalp(node1);
          }
        }
        else {
          node2 = sampleNodeUnd(root);
        }
        break;
      case 3:
        if (root->totalp == 0) {
          m_error = true;
          break;
        }
        node1 = sampleNodeUnd(root);
        node2 = insertNodeUnd(q, new_node_id);
        new_node_id++;
        break;
      case 4:
        node1 = insertNodeUnd(q, new_node_id);
        new_node_id++;
        node2 = insertNodeUnd(q, new_node_id);
        new_node_id++;
        break;
      case 5:
        node1 = node2 = insertNodeUnd(q, new_node_id);
        new_node_id++;
        break;
      }
      if (m_error) {
        break;
      }
      // sample without replacement
      if (node_unique) {
        if (node1->id < n_existing) {
          node1->p = 0;
          updateTotalp(node1);
        }
        if ((node2->id < n_existing) && (node1 != node2)) {
          node2->p = 0;
          updateTotalp(node2);
        }
      }
      node1->strength += edgeweight[new_edge_id];
      node2->strength += edgeweight[new_edge_id];
      node_vec1[new_edge_id] = node1->id;
      node_vec2[new_edge_id] = node2->id;
      scenario[new_edge_id] = current_scenario;
      q1.push(node1);
      q1.push(node2);
      new_edge_id++;
    }
    if (m_error) {
      m[i] = j;
      // need to print this info
      Rprintf("No enough unique nodes for a scenario %d edge at step %d. Added %d edge(s) at current step.\n", current_scenario, i + 1, j);
    }
    while (! q1.empty()) {
      updatePrefUnd(q1.front(), func_type, params, prefFuncCpp);
      q1.pop();
    }
  }
  PutRNGstate();
  // free memory (queue)
  queue<node_und*>().swap(q);
  queue<node_und*>().swap(q1);
  // save strength and preference
  q.push(root);
  node_und *temp_node;
  j = 0;
  while (! q.empty()) {
    temp_node = q.front();
    q.pop();
    if (temp_node->right != NULL) {
      q.push(temp_node->left);
      q.push(temp_node->right);
    }
    else if (temp_node->left != NULL) {
      q.push(temp_node->left);
    }
    strength[j] = temp_node->strength;
    pref[j] = temp_node->p;
    // free memory (node and tree)
    delete temp_node;
    j++;
  }
  // free memory (queue)
  queue<node_und*>().swap(q);

  Rcpp::List ret;
  ret["m"] = m;
  ret["nnode"] = new_node_id;
  ret["nedge"] = new_edge_id;
  ret["node_vec1"] = node_vec1;
  ret["node_vec2"] = node_vec2;
  ret["pref"] = pref;
  ret["strength"] = strength;
  ret["scenario"] = scenario;
  return ret;
}