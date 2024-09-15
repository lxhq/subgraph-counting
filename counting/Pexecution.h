#pragma once

#include <tbb/tbb.h>
#include <tbb/task_group.h>
#include "graph.h"
#include "decompose.h"
#include "equation.h"
#include "compute_set_intersection.h"
#include "forest.h"
#include "triangle.h"

std::vector<std::list<ui*>> evenly_splice(std::list<ui*>& partial_matches, int n);
class ParallelProcessingMeta  {
public:
    tbb::enumerable_thread_specific<int> _thread_id_ets;
    std::atomic<int> _next_thread_id;

    int _num_threads;
    tbb::task_group& _task_group;
    const DataGraph& _din;
    const DataGraph& _dout;
    const DataGraph& _dun;
    HashTable* _H;
    EdgeID* _inOffset;
    VertexID* _inNbors;
    EdgeID* _outOffset;
    VertexID* _outNbors;
    EdgeID* _unOffset;
    VertexID* _unNbors;

    // allocate memory for each thread
    HashTable* _total_hash_table;
    ui** _total_candidates;
    bool** _total_visited_vertices;
    std::list<VertexID*>* _total_extented_matches;
    VertexID** _tmp;

    ParallelProcessingMeta(int num_threads, 
                    tbb::task_group& task_group, 
                    const DataGraph& din,
                    const DataGraph& dout, 
                    const DataGraph& dun, 
                    HashTable* H);
    ~ParallelProcessingMeta();
};

void PexecuteNode(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        bool isRoot,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV,
        ParallelProcessingMeta &pMeta
);

void PexecuteNodeEdgeKey(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        bool isRoot,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV
);

// class ParallelExecuteNode {

// public:
//     std::vector<VertexID>& _child;
//     int _orbitType;
//     bool _isRoot;
//     Node& _tau;
//     std::vector<int>& _aggrePos;
//     std::vector<bool>& _nodeInterPos;
//     std::vector<std::vector<int>>& _nodeInPos;
//     std::vector<std::vector<int>>& _nodeOutPos;
//     std::vector<std::vector<int>>& _nodeUnPos;
//     std::vector<std::vector<int>>& _greaterPos;
//     std::vector<std::vector<int>>& _lessPos;
//     std::vector<std::vector<int>>& _childKeyPos;
//     std::vector<VertexID>& _aggreV;
//     std::vector<int>& _aggreWeight;


//     ParallelExecuteNode(VertexID nID, const Tree &t, const std::vector<VertexID> &child, 
//             VertexID **candidate, ui *candCount, HashTable *H, const DataGraph &din, const DataGraph &dout, 
//             const DataGraph &dun, const Pattern &p, bool isRoot, EdgeID *outID, EdgeID *unID, EdgeID *reverseID, 
//             EdgeID *startOffset, VertexID *patternV, VertexID *dataV, int mappingSize, bool *visited, ui *pos, ui *keyPos, 
//             ui &keyPosSize, ui sizeBound, VertexID *&tmp, VertexID *allV);

//     ~ParallelExecuteNode();

//     void operator()(std::list<ui*>& partial_matches, ui mappingSize);

//     void generateCandidate(VertexID *dataV, int mappingSize, VertexID *&candidate, ui& candCount, VertexID *&tmp); 
// };

