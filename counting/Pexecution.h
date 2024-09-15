#pragma once

#include <tbb/tbb.h>
#include <tbb/task_group.h>
#include "graph.h"
#include "decompose.h"
#include "equation.h"
#include "compute_set_intersection.h"
#include "forest.h"
#include "triangle.h"


class WorkerBFS_Node {

public:
    int _num_threads;
    tbb::task_group& _task_group;
    const DataGraph& _din;
    const DataGraph& _dout;
    const DataGraph& _dun;
    const std::vector<VertexID>& _child;
    int _orbitType;
    HashTable* _H;
    HashTable _h;
    bool _isRoot;
    const Node& _tau;
    const std::vector<int>& _aggrePos;
    const std::vector<bool>& _nodeInterPos;
    const std::vector<std::vector<int>>& _nodeInPos;
    const std::vector<std::vector<int>>& _nodeOutPos;
    const std::vector<std::vector<int>>& _nodeUnPos;
    const std::vector<std::vector<int>>& _greaterPos;
    const std::vector<std::vector<int>>& _lessPos;
    const std::vector<std::vector<int>>& _childKeyPos;
    const std::vector<VertexID>& _aggreV;
    const std::vector<int>& _aggreWeight;
    EdgeID* _inOffset;
    VertexID* _inNbors;
    EdgeID* _outOffset;
    VertexID* _outNbors;
    EdgeID* _unOffset;
    VertexID* _unNbors;

    tbb::enumerable_thread_specific<int> _thread_id_ets;
    std::atomic<int> _next_thread_id;

    HashTable* _total_hash_table;
    ui** _total_candidates;
    bool** _total_visited_vertices;
    std::list<VertexID*>* _total_extented_matches;
    VertexID** _tmp;

    WorkerBFS_Node(int num_threads, tbb::task_group& task_group_, VertexID nID, const Tree &t, const std::vector<VertexID> &child, 
            VertexID **candidate, ui *candCount, HashTable *H, const DataGraph &din, const DataGraph &dout, 
            const DataGraph &dun, const Pattern &p, bool isRoot, EdgeID *outID, EdgeID *unID, EdgeID *reverseID, 
            EdgeID *startOffset, VertexID *patternV, VertexID *dataV, int mappingSize, bool *visited, ui *pos, ui *keyPos, 
            ui &keyPosSize, ui sizeBound, VertexID *&tmp, VertexID *allV);
    
    Worker_BFS_Node(int num_threads, );

    set()

    ~WorkerBFS_Node();

    void operator()(std::list<ui*>& partial_matches, ui mappingSize);

    void PgenerateCandidate(VertexID *dataV, int mappingSize, VertexID *&candidate, ui& candCount, VertexID *&tmp); 
};

std::vector<std::list<ui*>> evenly_splice(std::list<ui*>& partial_matches, int n);