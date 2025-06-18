/*
 * Postponed node classification algorithm, Based on the Node Classification algorithm proposed by Feng in [1]
 *
 * [1] Gang Feng, "Finding k shortest simple paths in directed graphs: A node
 *     classification algorithm". Networks 64(1):6-17, 2014.
 *     doi:10.1002/net.21552
 */


#ifndef POSTPONED_NODE_CLASSIFICATION_H
#define POSTPONED_NODE_CLASSIFICATION_H

#include <cstdint>
#include <iostream>
#include <limits>       /* used to get maximum possible value of type TV */
#include "digraph.h"
#include "tools.h"
#include "dijkstra.h"
#include "pairing_heap.h"

namespace kssp {

template<
typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
typename TV   // type of edge weights (int, unsigned int, float, etc.)
>
class PostponedNodeClassification
{
public:
    PostponedNodeClassification() = default;

    // Constructor
    explicit PostponedNodeClassification(directed_graph::DirectedGraph<TI,TV> *g,
                                         const TI ssource,
                                         const TI ttarget);

    // Destructor
    virtual ~PostponedNodeClassification();

    // Check whether there is a candidate path
    bool empty();

    // Get path from/to u
    std::pair<std::vector<TI>,TV> next_path();

    size_t cpt_yielded_paths;
    size_t cpt_used_trees;

private:
    directed_graph::DirectedGraph<TI,TV> *graph;
    TI n;
    TI source;
    TI target;
    TV MAX_WEIGHT;

    dijkstra::Dijkstra<TI,TV> *DTo;

    heap::PairingHeap<tools_for_kssp::CandidatePath<TI,TV> *,TV> *heap_sorted_paths;
    std::vector<tools_for_kssp::CandidatePath<TI,TV> *> yielded_paths;

    std::set<TI> forbidden_vertices;
    std::vector<std::set<std::pair<TI,TI> > > forbidden_edges;
    std::set<std::pair<TI,TI> >  repair_forbidden_edges;

    // data structures for colors
    std::vector<std::set<TI> > succ_g_inv;
    size_t *color;

    std::vector<std::vector<std::pair<TI,TV> > > out_neighbors_residual;

    // data structures for shortest path computation
    TI *pred;
    TV *weight;

    // Reset some data structures
    void reset();

    // Initialize data structures for the manipulation of colors and residual weights
    void init_data_structures();

    // Compute and store deviations of this path
    void compute_deviations(tools_for_kssp::CandidatePath<TI,TV> *candidate_path);

    void repair_path(tools_for_kssp::CandidatePath<TI,TV> *candidate_path);

    // Special dijsktra from spur to green node or target
    void special_dijkstra(const TI spur_node);
};


// Constructor
template<typename TI, typename TV>
PostponedNodeClassification<TI, TV>::PostponedNodeClassification(directed_graph::DirectedGraph<TI,TV> *g,
                                                                 const TI ssource, const TI ttarget):
graph(g), n(g->n), source(ssource), target(ttarget)
{
    heap_sorted_paths = new heap::PairingHeap<tools_for_kssp::CandidatePath<TI,TV> *, TV>();
    yielded_paths.clear();
    forbidden_vertices.clear();
    forbidden_edges.clear();
    succ_g_inv.clear();
    out_neighbors_residual.clear();
    color = nullptr;
    DTo = nullptr;
    MAX_WEIGHT = std::numeric_limits<TV>::max();
    cpt_yielded_paths = 0;
    cpt_used_trees = 0;

    // Compute shortest path from source to target and store it
    if (source == target)
    {
        // trivial case
        std::vector<TI> path;
        path.clear();
        path.push_back(source);

        //the shortest path is always simple
        tools_for_kssp::CandidatePath<TI,TV> *C = new tools_for_kssp::CandidatePath<TI,TV>(path , 0, 0, true);
        heap_sorted_paths->push(C, 0);
    }
    else
    {
        // We call Dijkstra to target
        DTo = new dijkstra::Dijkstra<TI,TV>(graph, target, true);
        DTo->run(n);
        cpt_used_trees++;

        // and store computed path, if any
        if (DTo->successor(source) != source)
        {
            std::vector<TI> path = DTo->get_path(source);
            tools_for_kssp::CandidatePath<TI,TV> *C = new tools_for_kssp::CandidatePath<TI,TV>(path , DTo->weight(source), 0, true);
            heap_sorted_paths->push(C, 0);

            init_data_structures();
        }
    }
}


// Destructor
template<typename TI, typename TV>
PostponedNodeClassification<TI, TV>::~PostponedNodeClassification()
{
    if (color != nullptr)
    {
        free(color);
        free(weight);
        free(pred);
    }
    if (DTo != nullptr)
        delete DTo;
    while (not heap_sorted_paths->empty())
    {
        delete heap_sorted_paths->top_item();
        heap_sorted_paths->pop();
    }
    delete heap_sorted_paths;
    while (not yielded_paths.empty())
    {
        delete yielded_paths.back();
        yielded_paths.pop_back();
    }
    forbidden_vertices.clear();
    forbidden_edges.clear();
    succ_g_inv.clear();
    out_neighbors_residual.clear();
}


// Initialize data structures for the manipulation of colors and residual weights
template<typename TI, typename TV>
void PostponedNodeClassification<TI, TV>::init_data_structures()
{
    // Initialize data structure for assigning colors
    succ_g_inv.resize(n);
    for (TI u = 0; u < n; u++)
    {
        if (u == target)
            continue;
        succ_g_inv[DTo->successor(u)].insert(u);
    }

    // Compute residual edge weights
    out_neighbors_residual.clear();
    out_neighbors_residual.resize(n);
    for (TI u = 0; u < n; u++)
    {
        for (auto const& value: graph->out_neighbors[u])
        {
            TI v = value.first;
            TV w = value.second;
            out_neighbors_residual[u].push_back(std::make_pair(v, DTo->weight(v) + w - DTo->weight(u)));
        }
    }

    color = (size_t*)malloc(n * sizeof(size_t));
    pred = (TI*)malloc(n * sizeof(TI));
    weight = (TV*)malloc(n * sizeof(TV));
}


// Check whether there is a candidate path
template<typename TI, typename TV>
inline bool PostponedNodeClassification<TI, TV>::empty()
{
    return heap_sorted_paths->empty();
}


// Get next path
template<typename TI, typename TV>
std::pair<std::vector<TI>,TV> PostponedNodeClassification<TI, TV>::next_path()
{
    std::vector<TI> path;
    path.clear();
    if (heap_sorted_paths->empty())
        // no more paths
        return std::make_pair(path, 0);

    // Extract the next best path from the heap
    tools_for_kssp::CandidatePath<TI,TV> *prev_path = heap_sorted_paths->top_item();
    heap_sorted_paths->pop();

    while(not prev_path->simple){
        // if the extracted path is not simple, we repair it and re insert it to the heap
        repair_path(prev_path);
        // bug of last path fixed
        if (not heap_sorted_paths->empty()){
            prev_path = heap_sorted_paths->top_item();
            heap_sorted_paths->pop();
        }
        else
            return std::make_pair(path, 0);
    }

    if ((prev_path->weight == 0) || (prev_path->weight == MAX_WEIGHT) )
    {
        return std::make_pair(prev_path->path, 0);
    }
    // Compute deviations
    compute_deviations(prev_path);
    // Store and return the extracted path
    yielded_paths.push_back(prev_path);
    cpt_yielded_paths++;
    return std::make_pair(prev_path->path, prev_path->weight);
}


// Compute and store deviations of this path
template<typename TI, typename TV>
void PostponedNodeClassification<TI, TV>::compute_deviations(tools_for_kssp::CandidatePath<TI,TV> *candidate_path)
{
    /* Compute the sets of forbidden edges for each prefix.
     * We care only of edges incident to prev_path[i], as the edges incident
     * to other vertices in the prefix are forbidden by forbidden vertices
     */
    std::vector<TI> prev_path = candidate_path->path;
    size_t prev_path_size = prev_path.size();
    size_t dev_idx = candidate_path->dev_idx;

    forbidden_edges.clear();
    forbidden_edges.resize(prev_path_size);

    for (size_t i = dev_idx; i < prev_path_size - 1; i++)
        forbidden_edges[i].insert(std::make_pair(prev_path[i], prev_path[i + 1]));
    for (auto const& c_path: yielded_paths)
    {
        size_t j = 1;
        size_t jmax = (prev_path_size<c_path->path.size()?prev_path_size:c_path->path.size()) - 1;
        while ((j < jmax) && (prev_path[j] == c_path->path[j]))
            j++;
        if ((dev_idx < j) && (j < c_path->path.size()))
        {
            // j is the first index at which the paths differ
            //std::cout<<"adding to forbidden_edges the arc "<<c_path->path[j - 1]<<" "<<c_path->path[j]<<std::endl;
            forbidden_edges[j - 1].insert(std::make_pair(c_path->path[j - 1], c_path->path[j]));
        }
    }

    /* Set initial colors:
     * - red    = 0 (fobidden)
     * - yellow = 1 (default)
     * - green  = 2 (we can reach target)
     */
    for (size_t i = 0; i < n; i++)
        color[i] = 1;
    for (size_t i = 0; i <= dev_idx; i++)
        color[prev_path[i]] = 0;

    std::vector<TI> stack;
    stack.clear();
    stack.push_back(target);
    while (not stack.empty())
    {
        TI u = stack.back();
        stack.pop_back();
        color[u] = 2;
        for (TI const& v: succ_g_inv[u])
        {
            if (color[v] == 0)
                continue;
            stack.push_back(v);
        }
    }

    // Compute the cost of prefixes
    std::vector<TV> prefix_cost;
    prefix_cost.clear();
    prefix_cost.push_back(0);
    for (size_t i = 1; i < prev_path_size; i++)
        prefix_cost.push_back(prefix_cost[i - 1] + graph->edges[std::make_pair(prev_path[i - 1], prev_path[i])]);

    /*
     * Deviate from the previous path to find the candidate paths
     */
    for (size_t i = dev_idx; i < prev_path_size - 1; i++)
    {
        TI spur_node = prev_path[i];

        // update colors - set many vertices to yellow (1)
        stack.push_back(spur_node);
        while (not stack.empty())
        {
            TI u = stack.back();
            stack.pop_back();
            color[u] = 1;
            for (TI const& v: succ_g_inv[u])
            {
                if ((color[v] == 0) || (forbidden_edges[i].count(std::make_pair(v, u)) > 0))
                    continue;
                stack.push_back(v);
            }
        }
        color[spur_node] = 0;  // set red

        // Check if we can still find a path from this spur_node
        if (forbidden_edges[i].size() == out_neighbors_residual[spur_node].size())
            continue;

        // We check if the min residual weight of outgoing arcs from
        // spur_node leads to a green node. If so, we save a call to Dijkstra
        TI vv = prev_path[0];
        TV ww = MAX_WEIGHT;
        for (auto const& value: out_neighbors_residual[spur_node])
        {
            TI v = value.first;
            if ((color[v] == 0) || (forbidden_edges[i].count(std::make_pair(spur_node, v)) > 0))
                continue;
            TV w = value.second;
            if (w < ww)
            {
                vv = v;
                ww = w;
            }
        }

        if (ww == MAX_WEIGHT)
        { // We cannot follow any arc from spur_node
            continue;
        }
        else
        {
            // We add the path from spur to t, maybe we can skip this step/

            std::vector<TI> new_path;
            new_path.clear();
            for (size_t j = 0; j <= i; j++)
                new_path.push_back(prev_path[j]);

            new_path.push_back(vv);

            std::vector<TI> tmp = DTo->get_path(DTo->successor(vv));
            for (TI const& v: tmp)
                new_path.push_back(v);

            ww = graph->edges[std::make_pair(prev_path[i], vv)] + DTo->weight(vv);

            bool is_simple = false;
            if (color[vv] == 2)
                is_simple = true; // The best path is via a green node !!!

            // So, its extraction will be postponed if it is not simple!
            tools_for_kssp::CandidatePath<TI,TV> *C = new tools_for_kssp::CandidatePath<TI,TV>(new_path, prefix_cost[i] + ww, i, is_simple);
            heap_sorted_paths->push(C, C->weight);
        }
    }  // Go to next spur node
}


// Compute and store deviations of this path
template<typename TI, typename TV>
void PostponedNodeClassification<TI, TV>::repair_path(tools_for_kssp::CandidatePath<TI,TV> *candidate_path)
{
    // first we re-initialise the colors: to be tested w / o this step
    for (size_t u = 0; u < n; u++)
        color[u] = 1;

    for (size_t i = 0; i <= candidate_path->dev_idx; i++)
        color[candidate_path->path[i]] = 0; //forbidding the prefix

    // preparing the graph for a special dijkstra call
    repair_forbidden_edges.clear();

    for (auto const& c_path: yielded_paths)
    {
        size_t j = 1;
        size_t jmax = (candidate_path->path.size() < c_path->path.size() ? candidate_path->path.size() : c_path->path.size() ) - 1; // jmax = min (prev_path_size,c_path->path.size()) -1
        while ((j < jmax) && (candidate_path->path[j] == c_path->path[j]))
            j++;
        if ((candidate_path->dev_idx < j) && (j < c_path->path.size())){
            // j is the first index at which the paths differ
            repair_forbidden_edges.insert(std::make_pair(c_path->path[j - 1], c_path->path[j])); // forbidding arcs leading to old deviations
        }
    }

    TI spur_node = candidate_path->path[candidate_path->dev_idx];

    // Computing a suffix path using dijkstra
    special_dijkstra(spur_node);

    // tracing the total path

    TV prefix_cost = 0;
    for (size_t i = 1; i <= candidate_path->dev_idx; i++)
        prefix_cost += graph->edges[std::make_pair(candidate_path->path[i - 1], candidate_path->path[i])];

    if (pred[target] != target)
    { // We have a new path !!!

        // The new path starts with the prefix
        std::vector<TI> new_path;
        new_path.clear();
        TV ww;

        for (size_t j = 0; j <= candidate_path->dev_idx; j++)
            new_path.push_back(candidate_path->path[j]); // prefix is in new_path

        // We add the path from spur_node to a target
        std::vector<TI> tmp_path;
        ww = 0;
        TI u = target;
        while (u != spur_node)
        {
            tmp_path.push_back(u);
            ww += graph->edges[std::make_pair(pred[u], u)];
            u = pred[u];
        }
        while (not tmp_path.empty())
        {
            new_path.push_back(tmp_path.back());
            tmp_path.pop_back();
        }

        // We can now store the new path
        tools_for_kssp::CandidatePath<TI,TV> *C = new tools_for_kssp::CandidatePath<TI,TV>(new_path, prefix_cost + ww, candidate_path->dev_idx, true);
        heap_sorted_paths->push(C, C->weight);
    }
    delete candidate_path;
}


// Special dijsktra from spur to green node or target
template<typename TI, typename TV>
void PostponedNodeClassification<TI, TV>::special_dijkstra(const TI spur_node)
{
    for (size_t u = 0; u < n; u++)
    {
        weight[u] = MAX_WEIGHT;
        pred[u] = u;
    }

    heap::PairingHeap<TI,TV> *pq = new heap::PairingHeap<TI,TV>();
    pq->push(spur_node, 0);
    weight[spur_node] = 0;

    while (not pq->empty())
    {
        TI u = pq->top_item();
        if (u == target){
            break;
        }

        pq->pop();

        for (auto const& value: out_neighbors_residual[u])
        {
            TI v = value.first;
            TV w = value.second;
            if ((color[v] == 0) || (repair_forbidden_edges.count(std::make_pair(u, v)) > 0)){
                continue;
            }
            if (weight[u] + w < weight[v])
            {
                weight[v] = weight[u] + w;
                if (pred[v] == v)
                    pq->push(v, weight[v]);
                else
                    pq->decrease(v, weight[v]);
                pred[v] = u;
            }
        }
    }
    delete pq;
    cpt_used_trees++;
}

}
#endif
