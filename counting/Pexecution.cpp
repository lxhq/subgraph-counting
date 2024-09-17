#include "Pexecution.h"

ParallelProcessingMeta::ParallelProcessingMeta(
    int num_threads, 
    const DataGraph& din,
    const DataGraph& dout, 
    const DataGraph& dun) : 
    _num_threads(num_threads),
    _thread_id_ets(-1),
    _next_thread_id(0) {
    
    ui n = dun.getNumVertices(), m = dun.getNumEdges();

    _total_hash_table = new HashTable*[num_threads];
    _total_partition_candidates_pos = new ui*[num_threads];
    _total_start_offset = new EdgeID*[num_threads];
    _total_visited_vertices = new bool*[num_threads];
    _total_tmp = new VertexID*[num_threads];

    // allocate memory for each thread
    for (ui i = 0; i < num_threads; i++) {
        _total_hash_table[i] = new HashTable[MAX_NUM_NODE];
        _total_start_offset[i] = new EdgeID[MAX_PATTERN_SIZE];
        std::fill(_total_start_offset[i], _total_start_offset[i] + MAX_PATTERN_SIZE, 0);

        for (ui j = 0; j < MAX_NUM_NODE; j++) {
            _total_hash_table[i][j] = new Count[m];
            std::fill(_total_hash_table[i][j], _total_hash_table[i][j] + m, 0);
        }
        _total_partition_candidates_pos[i] = new ui[MAX_PATTERN_SIZE];
        std::fill(_total_partition_candidates_pos[i], _total_partition_candidates_pos[i] + MAX_PATTERN_SIZE, 0);
        _total_visited_vertices[i] = new bool[dun.getNumVertices()];
        std::fill(_total_visited_vertices[i], _total_visited_vertices[i] + dun.getNumVertices(), false);
        _total_tmp[i] = new VertexID[dun.getNumVertices()];
    }
}

ParallelProcessingMeta::~ParallelProcessingMeta() {
    // deallocate the memory for each thread
    for (ui i = 0; i < _num_threads; i++) {
        for (ui j = 0; j < MAX_NUM_NODE; j++) {
            delete [] _total_hash_table[i][j];
        }
        delete [] _total_hash_table[i];
        delete [] _total_start_offset[i];
        delete [] _total_partition_candidates_pos[i];
        delete [] _total_visited_vertices[i];
        delete [] _total_tmp[i];
    }
    delete [] _total_hash_table;
    delete [] _total_start_offset;
    delete [] _total_partition_candidates_pos;
    delete [] _total_visited_vertices;
    delete [] _total_tmp;
}

void ParallelProcessingMeta::setCandidates(const Tree& t, const DataGraph& dout) {
    int numNodes = (int)t.getNumNodes();
    _total_candidates = new VertexID ***[_num_threads];
    _total_candidates_cnt = new ui **[_num_threads];
    for (ui thread = 0; thread < _num_threads; thread++) {
        _total_candidates[thread] = new VertexID **[numNodes];
        _total_candidates_cnt[thread] = new ui *[numNodes];
        for (VertexID nID = 0; nID < numNodes; ++nID) {
            const Node &tau = t.getNode(nID);
            int partOrderLength = int(tau.nodeOrder.size() - tau.localOrder.size());
            _total_candidates_cnt[thread][nID] = new ui[tau.nodeOrder.size()];
            _total_candidates[thread][nID] = new VertexID *[tau.nodeOrder.size()];
            const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
            for (int i = partOrderLength; i < nodeCandPos.size(); ++i) {
                if (nodeCandPos[i])
                    _total_candidates[thread][nID][i] = new VertexID[dout.getNumVertices()];
            }
        }
    }
}

void ParallelProcessingMeta::clearCandidates(const Tree& t) {
    int numNodes = (int)t.getNumNodes();
    for (ui thread = 0; thread < _num_threads; thread++) {
        for (VertexID nID = 0; nID < numNodes; ++nID) {
            const Node &tau = t.getNode(nID);
            int partOrderLength = int(tau.nodeOrder.size() - tau.localOrder.size());
            const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
            for (int i = partOrderLength; i < nodeCandPos.size(); ++i) {
                if (nodeCandPos[i])
                    delete[] _total_candidates[thread][nID][i];
            }
            delete[] _total_candidates[thread][nID];
            delete[] _total_candidates_cnt[thread][nID];
        }
        delete[] _total_candidates_cnt[thread];
        delete[] _total_candidates[thread];
    }
    delete[] _total_candidates;
    delete[] _total_candidates_cnt;
}

void ParallelProcessingMeta::setPartitionCandidates(const Tree& t, const std::vector<VertexID>& partitionOrder, 
                                                    const std::vector<bool> &partitionCandPos, const DataGraph& dout, 
                                                    const std::vector<VertexID> &postOrder, int startPos, int endPos, ui n, ui m) {
    _total_partition_candidates = new ui **[_num_threads];
    _total_partition_candidates_cnt = new ui *[_num_threads];
    _total_key_pos = new ui **[_num_threads];
    _total_key_pos_size = new ui*[_num_threads];

    for (ui thread = 0; thread < _num_threads; thread++) {
        _total_partition_candidates[thread] = new VertexID *[partitionOrder.size()];
        for (int i = 0; i < partitionOrder.size(); ++i) {
            if (partitionCandPos[i])
                _total_partition_candidates[thread][i] = new VertexID[dout.getNumVertices()];
        }
        _total_partition_candidates_cnt[thread] = new ui[partitionOrder.size()];
        // for each node, store the keys whose values are not 0
        // those keys will be cleared to 0 before the next visit of this node begins
        // If the used entries is larger m / 8. we just clear the whole hashtable
        _total_key_pos[thread] = new ui *[t.getNumNodes()];
        for (int i = startPos; i < endPos; ++i) {
            VertexID nID = postOrder[i];
            if (t.getNode(nID).keySize == 0)
                _total_key_pos[thread][nID] = new ui[1];
            else if (t.getNode(nID).keySize == 1)
                _total_key_pos[thread][nID] = new ui[n / 8 + 1];
            else if (t.getNode(nID).keySize == 2)
                _total_key_pos[thread][nID] = new ui[m / 8 + 1];
        }
        _total_key_pos_size[thread] = new ui[t.getNumNodes()];
        memset(_total_key_pos_size[thread], 0, sizeof(ui) * (t.getNumNodes()));
    }
}

void ParallelProcessingMeta::clearPartitionCandidates(const std::vector<VertexID>& partitionOrder, const std::vector<bool> &partitionCandPos, 
                                                      const std::vector<VertexID> &postOrder, int startPos, int endPos) {
    for (ui thread = 0; thread < _num_threads; thread++) {
        for (int j = 0; j < partitionOrder.size(); ++j) {
            if (partitionCandPos[j])
                delete[] _total_partition_candidates[thread][j];
        }
        delete[] _total_partition_candidates[thread];
        delete[] _total_partition_candidates_cnt[thread];
        for (int j = startPos; j < endPos; ++j) {
            VertexID nID = postOrder[j];
            delete[] _total_key_pos[thread][postOrder[j]];
        }
        delete[] _total_key_pos[thread];
        delete[] _total_key_pos_size[thread];
    }
    delete [] _total_partition_candidates;
    delete [] _total_partition_candidates_cnt;
    delete [] _total_key_pos;
    delete [] _total_key_pos_size;
}

// std::vector<std::list<ui*>> evenly_splice(std::list<ui*>& partial_matches, int n);
// void generateCandidate(
//         const DataGraph &din,
//         const DataGraph &dout,
//         const DataGraph &dun,
//         VertexID *dataV,
//         const std::vector<int> &inPos,
//         const std::vector<int> &outPos,
//         const std::vector<int> &unPos,
//         VertexID *&candidate,
//         ui &candCount,
//         VertexID *&tmp
// );
// void worker(
//         ParallelProcessingMeta& pMeta, 
//         std::list<ui*>& partial_matches, 
//         const DataGraph& din,
//         const DataGraph& dout,
//         const DataGraph& dun,
//         const std::vector<bool> &nodeInterPos,
//         const std::vector<std::vector<int>> &nodeInPos,
//         const std::vector<std::vector<int>> &nodeOutPos,
//         const std::vector<std::vector<int>> &nodeUnPos,
//         const std::vector<std::vector<int>> &greaterPos,
//         const std::vector<std::vector<int>> &lessPos,
//         EdgeID *inOffset,
//         VertexID *inNbors,
//         EdgeID *outOffset,
//         VertexID *outNbors,
//         EdgeID *unOffset,
//         VertexID *unNbors,
//         const std::vector<VertexID> &child,
//         const std::vector<std::vector<int>> &childKeyPos,
//         const std::vector<VertexID> &aggreV,
//         const std::vector<int> &aggrePos,
//         const std::vector<int> &aggreWeight,
//         const Node &tau,
//         HashTable* H,
//         bool isRoot,
//         int orbitType,
//         ui mappingSize
// );

// std::vector<std::list<ui*>> evenly_splice(std::list<ui*>& partial_matches, int n) {
//     std::vector<std::list<ui*>> sublists(n);  // Create n empty sub-lists
//     if (partial_matches.empty()) {
//         return sublists;
//     }

//     int total_size = partial_matches.size();
//     int base_size = total_size / n;   // Base number of elements per sub-list
//     int remainder = total_size % n;   // Extra elements to distribute

//     auto it = partial_matches.begin();
    
//     for (int i = 0; i < n; ++i) {
//         int current_size = base_size + (i < remainder ? 1 : 0);  // Extra element for first 'remainder' lists

//         auto next_it = it;
//         std::advance(next_it, current_size);  // Move iterator to the end of the current sub-list slice

//         sublists[i].splice(sublists[i].end(), partial_matches, it, next_it);  // Move elements to the sub-list

//         it = next_it;  // Update iterator to the next starting point
//     }

//     return sublists;
// }

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
// ) {
//     int orbitType = t.getOrbitType();
//     HashTable h = H[nID];
//     const Node &tau = t.getNode(nID);
//     const std::vector<int> &aggrePos = t.getAggrePos(nID);
//     const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
//     const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
//     const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
//     const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
//     const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
//     const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
//     const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
//     const std::vector<VertexID> &aggreV = t.getAggreV();
//     const std::vector<int> &aggreWeight = t.getAggreWeight();
//     EdgeID *inOffset = din.getOffsets();
//     VertexID *inNbors = din.getNbors();
//     EdgeID *outOffset = dout.getOffsets();
//     VertexID *outNbors = dout.getNbors();
//     EdgeID *unOffset = dun.getOffsets();
//     VertexID *unNbors = dun.getNbors();

//     ui n = din.getNumVertices();
//     if (mappingSize == 0) {
//         candidate[0] = allV;
//         candCount[0] = n;
//         pos[mappingSize] = 0;
//     }
//     else {
//         if (!nodeInterPos[mappingSize]) {
//             if (!nodeOutPos[mappingSize].empty()) {
//                 VertexID v = dataV[nodeOutPos[mappingSize][0]];
//                 startOffset[mappingSize] = outOffset[v];
//                 candidate[mappingSize] = outNbors + outOffset[v];
//                 candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
//             }
//             else if (!nodeInPos[mappingSize].empty()){
//                 VertexID v = dataV[nodeInPos[mappingSize][0]];
//                 startOffset[mappingSize] = inOffset[v];
//                 candidate[mappingSize] = inNbors + inOffset[v];
//                 candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
//             }
//             else {
//                 VertexID v = dataV[nodeUnPos[mappingSize][0]];
//                 startOffset[mappingSize] = unOffset[v];
//                 candidate[mappingSize] = unNbors + unOffset[v];
//                 candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
//             }
//         }
//         else {
//             generateCandidate(din, dout, dun, dataV, nodeInPos[mappingSize],
//                               nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate[mappingSize], candCount[mappingSize], tmp);
//         }
//         pos[mappingSize] = 0;
//         // apply the symmetry breaking rules
//         if (!greaterPos[mappingSize].empty()) {
//             ui maxTarget = dataV[greaterPos[mappingSize][0]];
//             for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
//                 if (dataV[greaterPos[mappingSize][i]] > maxTarget)
//                     maxTarget = dataV[greaterPos[mappingSize][i]];
//             }
//             pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
//         }
//         if (!lessPos[mappingSize].empty()) {
//             ui minTarget = dataV[lessPos[mappingSize][0]];
//             for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
//                 if (dataV[lessPos[mappingSize][i]] < minTarget)
//                     minTarget = dataV[lessPos[mappingSize][i]];
//             }
//             candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
//         }
//     }
//     if (candCount[mappingSize] == 0) return;
//     std::list<ui*> partial_matches;
//     for (ui i = 0; i < candCount[mappingSize]; i++) {
//         ui* partial_match = new ui[tau.nodeOrder.size()];
//         for (ui j = 0; j < mappingSize; j++) {
//             partial_match[j] = dataV[j];
//         }
//         partial_match[mappingSize] = candidate[mappingSize][i];
//         partial_matches.push_back(partial_match);
//     }
//     tbb::task_group task_group;
//     while (mappingSize < tau.nodeOrder.size() - 1) {
//         mappingSize += 1;
//         patternV[mappingSize] = tau.nodeOrder[mappingSize];
//         std::vector<std::list<ui*>> sublists = evenly_splice(partial_matches, pMeta._num_threads);
//         for (std::list<ui*>& sublist : sublists) {
//             if (sublist.empty()) continue;
//             task_group.run([
//                 &pMeta, 
//                 &partial_matches, 
//                 &din,
//                 &dout,
//                 &dun,
//                 &nodeInterPos,
//                 &nodeInPos,
//                 &nodeOutPos,
//                 &nodeUnPos,
//                 &greaterPos,
//                 &lessPos,
//                 inOffset,
//                 inNbors,
//                 outOffset,
//                 outNbors,
//                 unOffset,
//                 unNbors,
//                 &child,
//                 &childKeyPos,
//                 &aggreV,
//                 &aggrePos,
//                 &aggreWeight,
//                 &tau,
//                 H,
//                 isRoot,
//                 orbitType,
//                 mappingSize
//             ]() {
//                 worker(pMeta, 
//                     partial_matches, 
//                     din, 
//                     dout, 
//                     dun, 
//                     nodeInterPos, 
//                     nodeInPos, 
//                     nodeOutPos, 
//                     nodeUnPos, 
//                     greaterPos, 
//                     lessPos, 
//                     inOffset, 
//                     inNbors, 
//                     outOffset, 
//                     outNbors, 
//                     unOffset, 
//                     unNbors, 
//                     child, 
//                     childKeyPos, 
//                     aggreV, 
//                     aggrePos, 
//                     aggreWeight, 
//                     tau, 
//                     H, 
//                     isRoot, 
//                     orbitType, 
//                     mappingSize
//                 );
//             });
//         }   
//         task_group.wait();

//         for (ui i = 0; i < pMeta._num_threads; i++) {
//             partial_matches.splice(partial_matches.end(), pMeta._total_extented_matches[i]);
//         }
//     }

//     // update the hash table and keyPos from the hashtable of each threads
//     // for (ui i = 0; i < pMeta._num_threads; i++) {
//     //     for (ui j = 0; j < dun.getNumEdges(); j++) {
//     //         h[j] += pMeta._total_hash_table[i][j];
//     //         pMeta._total_hash_table[i][j] = 0;
//     //     }
//     // }
//     // for (ui i = 0; i < pMeta._num_threads; i++) {
//     //     if (h[i] != 0 && keyPosSize < sizeBound) {
//     //         keyPos[keyPosSize] = i;
//     //         keyPosSize += 1;
//     //     }
//     // }
//     keyPosSize = 0;

//     for (ui* partial_match : partial_matches) {
//         delete [] partial_match;
//     }
// }

// void worker(
//         ParallelProcessingMeta& pMeta, 
//         std::list<ui*>& partial_matches, 
//         const DataGraph& din,
//         const DataGraph& dout,
//         const DataGraph& dun,
//         const std::vector<bool> &nodeInterPos,
//         const std::vector<std::vector<int>> &nodeInPos,
//         const std::vector<std::vector<int>> &nodeOutPos,
//         const std::vector<std::vector<int>> &nodeUnPos,
//         const std::vector<std::vector<int>> &greaterPos,
//         const std::vector<std::vector<int>> &lessPos,
//         EdgeID *inOffset,
//         VertexID *inNbors,
//         EdgeID *outOffset,
//         VertexID *outNbors,
//         EdgeID *unOffset,
//         VertexID *unNbors,
//         const std::vector<VertexID> &child,
//         const std::vector<std::vector<int>> &childKeyPos,
//         const std::vector<VertexID> &aggreV,
//         const std::vector<int> &aggrePos,
//         const std::vector<int> &aggreWeight,
//         const Node &tau,
//         HashTable* H,
//         bool isRoot,
//         int orbitType,
//         ui mappingSize) {
//     if (pMeta._thread_id_ets.local() == -1) {
//         pMeta._thread_id_ets.local() = pMeta._next_thread_id.fetch_add(1);
//     }
//     ui* candidates = pMeta._total_candidates[pMeta._thread_id_ets.local()];
//     ui candCount = 0;
//     bool* visited_vertices = pMeta._total_visited_vertices[pMeta._thread_id_ets.local()];
//     std::list<ui*>& extented_matches = pMeta._total_extented_matches[pMeta._thread_id_ets.local()];
//     HashTable h = pMeta._total_hash_table[pMeta._thread_id_ets.local()];
//     VertexID* tmp = pMeta._tmp[pMeta._thread_id_ets.local()];
//     while (!partial_matches.empty()) {
//         ui* dataV = partial_matches.front();
//         partial_matches.pop_front();

//         if (!nodeInterPos[mappingSize]) {
//             if (!nodeOutPos[mappingSize].empty()) {
//                 VertexID w = dataV[nodeOutPos[mappingSize][0]];
//                 candidates = outNbors + outOffset[w];
//                 candCount = outOffset[w + 1] - outOffset[w];
//             }
//             else if (!nodeInPos[mappingSize].empty()){
//                 VertexID w = dataV[nodeInPos[mappingSize][0]];
//                 candidates = inNbors + inOffset[w];
//                 candCount = inOffset[w + 1] - inOffset[w];
//             }
//             else {
//                 VertexID w = dataV[nodeUnPos[mappingSize][0]];
//                 candidates = unNbors + unOffset[w];
//                 candCount = unOffset[w + 1] - unOffset[w];
//             }
//         }
//         else {
//             generateCandidate(din, dout, dun, dataV, nodeInPos[mappingSize],
//                               nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidates, candCount, tmp);
//         }
//         // apply the symmetry breaking rules
//         ui start_pos = 0;
//         if (!greaterPos[mappingSize].empty()) {
//             ui maxTarget = dataV[greaterPos[mappingSize][0]];
//             for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
//                 if (dataV[greaterPos[mappingSize][i]] > maxTarget)
//                     maxTarget = dataV[greaterPos[mappingSize][i]];
//             }
//             start_pos = firstPosGreaterThan(candidates, 0, candCount, maxTarget);
//         }
//         if (!lessPos[mappingSize].empty()) {
//             ui minTarget = dataV[lessPos[mappingSize][0]];
//             for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
//                 if (dataV[lessPos[mappingSize][i]] < minTarget)
//                     minTarget = dataV[lessPos[mappingSize][i]];
//             }
//             candCount = firstPosGreaterThan(candidates, 0, candCount, minTarget);
//         }
//         for (ui i = 0; i < mappingSize; ++i) {
//             visited_vertices[dataV[i]] = true;
//         }
//         for (ui i = start_pos; i < candCount; i++) {
//             VertexID v = candidates[i];
//             if (visited_vertices[v]) {
//                 continue;
//             }
//             dataV[mappingSize] = v;
//             if (mappingSize == tau.nodeOrder.size() - 1) {
//                 Count cnt = 1;
//                 for (int j = 0; j < child.size(); j++) {
//                     VertexID cID = child[j];
//                     cnt *= H[cID][dataV[childKeyPos[j][0]]];
//                 }
//                 if (isRoot) {
//                     if (orbitType == 0) h[0] += cnt * aggreWeight[0];
//                     else if (orbitType == 1) {
//                         for (int j = 0; j < aggreV.size(); j++) {
//                             VertexID key = dataV[aggrePos[j]];
//                             h[key] += cnt * aggreWeight[j];
//                         }
//                     }
//                 } else {
//                     h[dataV[aggrePos[0]]] += cnt;
//                 }
//             } else {
//                 ui* new_embedding = new ui[tau.nodeOrder.size()];
//                 std::copy(dataV, dataV + tau.nodeOrder.size(), new_embedding);
//                 extented_matches.push_back(new_embedding);
//             }
//         }
//         for (ui i = 0; i < mappingSize; i++) {
//             visited_vertices[dataV[i]] = false;
//         }
//         delete [] dataV;
//     }
// }

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
// ) {

// }

// void generateCandidate(
//         const DataGraph &din,
//         const DataGraph &dout,
//         const DataGraph &dun,
//         VertexID *dataV,
//         const std::vector<int> &inPos,
//         const std::vector<int> &outPos,
//         const std::vector<int> &unPos,
//         VertexID *&candidate,
//         ui &candCount,
//         VertexID *&tmp
// ) {
//     int num = inPos.size() + outPos.size() + unPos.size();
//     ui tmpCount;
//     if (num > 2) {
//         VertexID **arrays = new VertexID *[num];
//         ui *counts = new ui[num];
//         int p = 0;
//         for (int position: inPos) {
//             arrays[p] = din.getNeighbors(dataV[position], counts[p]);
//             ++p;
//         }
//         for (int position: outPos) {
//             arrays[p] = dout.getNeighbors(dataV[position], counts[p]);
//             ++p;
//         }
//         for (int position: unPos) {
//             arrays[p] = dun.getNeighbors(dataV[position], counts[p]);
//             ++p;
//         }
//         ComputeSetIntersection::LeapfrogJoin(arrays, counts, num, tmp, tmpCount);
//         candCount = tmpCount;
//         std::swap(candidate, tmp);
//         delete[] counts;
//         delete[] arrays;
//         return;
//     }
//     VertexID *firstNeighbor;
//     ui firstCnt;
//     bool firstFlag = true, secondFlag = true;
//     for (int position: inPos) {
//         ui nbrCnt;
//         VertexID *neighbor = din.getNeighbors(dataV[position], nbrCnt);
//         if (firstFlag) {
//             firstFlag = false;
//             firstCnt = nbrCnt;
//             firstNeighbor = neighbor;
//         }
//         else {
//             if (secondFlag) {
//                 secondFlag = false;
// #ifdef COLLECT_STATISTICS
//                 ++gNumIntersect;
// #endif
//                 ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
//                 candCount = tmpCount;
//                 std::swap(candidate, tmp);
//             }
//             else {
// #ifdef COLLECT_STATISTICS
//                 ++gNumIntersect;
// #endif
//                 ComputeSetIntersection::ComputeCandidates(candidate, candCount, neighbor, nbrCnt, tmp, tmpCount);
//                 candCount = tmpCount;
//                 std::swap(candidate, tmp);
//             }
//         }
//     }
//     for (int position: outPos) {
//         ui nbrCnt;
//         VertexID *neighbor = dout.getNeighbors(dataV[position], nbrCnt);
//         if (firstFlag) {
//             firstFlag = false;
//             firstCnt = nbrCnt;
//             firstNeighbor = neighbor;
//         }
//         else {
//             if (secondFlag) {
//                 secondFlag = false;
// #ifdef COLLECT_STATISTICS
//                 ++gNumIntersect;
// #endif
//                 ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
//                 candCount = tmpCount;
//                 std::swap(candidate, tmp);
//             }
//             else {
// #ifdef COLLECT_STATISTICS
//                 ++gNumIntersect;
// #endif
//                 ComputeSetIntersection::ComputeCandidates(candidate, candCount, neighbor, nbrCnt, tmp, tmpCount);
//                 candCount = tmpCount;
//                 std::swap(candidate, tmp);
//             }
//         }
//     }
//     for (int position: unPos) {
//         ui nbrCnt;
//         VertexID *neighbor = dun.getNeighbors(dataV[position], nbrCnt);
//         if (firstFlag) {
//             firstFlag = false;
//             firstCnt = nbrCnt;
//             firstNeighbor = neighbor;
//         }
//         else {
//             if (secondFlag) {
//                 secondFlag = false;
// #ifdef COLLECT_STATISTICS
//                 ++gNumIntersect;
// #endif
//                 ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
//                 candCount = tmpCount;
//                 std::swap(candidate, tmp);
//             }
//             else {
// #ifdef COLLECT_STATISTICS
//                 ++gNumIntersect;
// #endif
//                 ComputeSetIntersection::ComputeCandidates(candidate, candCount, neighbor, nbrCnt, tmp, tmpCount);
//                 candCount = tmpCount;
//                 std::swap(candidate, tmp);
//             }
//         }
//     }
// }