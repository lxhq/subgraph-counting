#pragma once

#include <tbb/tbb.h>
#include <tbb/task_group.h>
#include "graph.h"
#include "decompose.h"
#include "equation.h"
#include "compute_set_intersection.h"
#include "forest.h"
#include "triangle.h"

class ParallelProcessingMeta {
public:
    int _num_threads;
    tbb::enumerable_thread_specific<int> _thread_id_ets;
    std::atomic<int> _next_thread_id;

    // allocate memory for each thread
    HashTable** _total_hash_table; 
    ui**** _total_candidates;
    ui*** _total_candidates_cnt;
    EdgeID** _total_start_offset;
    ui*** _total_partition_candidates;
    ui** _total_partition_candidates_cnt;
    ui** _total_partition_candidates_pos;
    ui*** _total_key_pos;
    ui** _total_key_pos_size;
    bool** _total_visited_vertices;
    VertexID** _total_tmp;
    ParallelProcessingMeta(int num_threads, 
                    const DataGraph& din,
                    const DataGraph& dout, 
                    const DataGraph& dun);
    void setCandidates(const Tree& t, const DataGraph& dout);
    void clearCandidates(const Tree& t);
    void setPartitionCandidates(const Tree& t, const std::vector<VertexID>& partitionOrder, 
                                const std::vector<bool> &partitionCandPos, const DataGraph& dout, 
                                const std::vector<VertexID> &postOrder, int startPos, int endPos, ui n, ui m);
    void clearPartitionCandidates(const std::vector<VertexID>& partitionOrder, const std::vector<bool> &partitionCandPos, 
                                  const std::vector<VertexID> &postOrder, int startPos, int endPos);
    ~ParallelProcessingMeta();
};

// void PexecuteNode(
//         VertexID nID,
//         const Tree &t,
//         const std::vector<VertexID> &child,
//         VertexID **candidate,
//         ui *candCount,
//         HashTable *H,
//         const DataGraph &din,
//         const DataGraph &dout,
//         const DataGraph &dun,
//         const Pattern &p,
//         bool isRoot,
//         EdgeID *outID,
//         EdgeID *unID,
//         EdgeID *reverseID,
//         EdgeID *startOffset,
//         VertexID *patternV,
//         VertexID *dataV,
//         int mappingSize,
//         bool *visited,
//         ui *pos,
//         ui *keyPos,
//         ui &keyPosSize,
//         ui sizeBound,
//         VertexID *&tmp,
//         VertexID *allV,
//         ParallelProcessingMeta &pMeta
// );

// void PexecuteNodeEdgeKey(
//         VertexID nID,
//         const Tree &t,
//         const std::vector<VertexID> &child,
//         VertexID **candidate,
//         ui *candCount,
//         HashTable *H,
//         const DataGraph &din,
//         const DataGraph &dout,
//         const DataGraph &dun,
//         const Pattern &p,
//         bool isRoot,
//         EdgeID *outID,
//         EdgeID *unID,
//         EdgeID *reverseID,
//         EdgeID *startOffset,
//         VertexID *patternV,
//         VertexID *dataV,
//         int mappingSize,
//         bool *visited,
//         ui *pos,
//         ui *keyPos,
//         ui &keyPosSize,
//         ui sizeBound,
//         VertexID *&tmp,
//         VertexID *allV
// );

void PexecuteParation(ParallelProcessingMeta& pMeta);