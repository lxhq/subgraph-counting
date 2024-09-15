#include "Pexecution.h"

std::vector<std::list<ui*>> evenly_splice(std::list<ui*>& partial_matches, int n);
void generateCandidate(
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        VertexID *dataV,
        int pos,
        const std::vector<int> &inPos,
        const std::vector<int> &outPos,
        const std::vector<int> &unPos,
        VertexID **candidate,
        ui *&candCount,
        VertexID *&tmp
);
void worker(
        std::list<ui*>& partial_matches, 
        ui mappingSize
);

ParallelProcessingMeta::ParallelProcessingMeta(
    int num_threads, 
    tbb::task_group& task_group, 
    const DataGraph& din,
    const DataGraph& dout, 
    const DataGraph& dun, 
    HashTable* H) : _thread_id_ets(-1),
    _next_thread_id(0),
    _num_threads(num_threads),
    _task_group(task_group),
    _din(din),
    _dout(dout),
    _dun(dun),
    _H(H) {

    _inOffset = din.getOffsets();
    _inNbors = din.getNbors();
    _outOffset = dout.getOffsets();
    _outNbors = dout.getNbors();
    _unOffset = dun.getOffsets();
    _unNbors = dun.getNbors();
    
    // each thread must have their own memory location for:
    // 1. candidates
    // 2. visited nodes
    _total_candidates = new ui*[num_threads];
    _total_visited_vertices = new bool*[num_threads];

    // 3. hash table
    _total_hash_table = new HashTable[num_threads];

    // 4. extended partial matches
    _total_extented_matches = new std::list<ui*>[num_threads];

    // 5. tmp used for Intersection
    _tmp = new VertexID*[num_threads];

    // allocate memory for each thread
    for (ui i = 0; i < num_threads; i++) {
        _total_hash_table[i] = new Count[dun.getNumEdges()];
        _total_candidates[i] = new ui[dun.getNumVertices()];
        _total_visited_vertices[i] = new bool[dun.getNumVertices()];
        std::fill(_total_visited_vertices[i], _total_visited_vertices[i] + dun.getNumVertices(), false);
        std::fill(_total_hash_table[i], _total_hash_table[i] + dun.getNumEdges(), 0);
        _tmp[i] = new VertexID[dun.getNumVertices()];
    }
}

ParallelProcessingMeta::~ParallelProcessingMeta() {
    // deallocate the memory for each thread
    for (ui i = 0; i < _num_threads; i++) {
        delete [] _total_candidates[i];
        delete [] _total_visited_vertices[i];
        delete [] _total_hash_table[i];
        delete [] _tmp[i];
    }
    delete [] _total_candidates;
    delete [] _total_visited_vertices;
    delete [] _total_hash_table;
    delete [] _total_extented_matches;
    delete [] _tmp;
}

std::vector<std::list<ui*>> evenly_splice(std::list<ui*>& partial_matches, int n) {
    std::vector<std::list<ui*>> sublists(n);  // Create n empty sub-lists
    if (partial_matches.empty()) {
        return sublists;
    }

    int total_size = partial_matches.size();
    int base_size = total_size / n;   // Base number of elements per sub-list
    int remainder = total_size % n;   // Extra elements to distribute

    auto it = partial_matches.begin();
    
    for (int i = 0; i < n; ++i) {
        int current_size = base_size + (i < remainder ? 1 : 0);  // Extra element for first 'remainder' lists

        auto next_it = it;
        std::advance(next_it, current_size);  // Move iterator to the end of the current sub-list slice

        sublists[i].splice(sublists[i].end(), partial_matches, it, next_it);  // Move elements to the sub-list

        it = next_it;  // Update iterator to the next starting point
    }

    return sublists;
}

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
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();

    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                              nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    if (candCount[mappingSize] == 0) return;

    std::list<ui*> partial_matches;
    for (ui i = 0; i < candCount[mappingSize]; i++) {
        ui* partial_match = new ui[tau.nodeOrder.size()];
        for (ui j = 0; j < mappingSize; j++) {
            partial_match[j] = dataV[j];
        }
        partial_match[mappingSize] = candidate[mappingSize][i];
        partial_matches.push_back(partial_match);
    }

    while (mappingSize < tau.nodeOrder.size() - 1) {
        mappingSize += 1;
        patternV[mappingSize] = tau.nodeOrder[mappingSize];
        std::vector<std::list<ui*>> sublists = evenly_splice(partial_matches, pMeta._num_threads);
        for (std::list<ui*>& sublist : sublists) {
            pMeta._task_group.run([&pMeta, &sublist, mappingSize]() {
                worker(pMeta, sublist, mappingSize);
            });
        }   
        pMeta._task_group.wait();

        for (ui i = 0; i < pMeta._num_threads; i++) {
            partial_matches.splice(partial_matches.end(), pMeta._total_extented_matches[i]);
        }
    }

    // update the hash table and keyPos from the hashtable of each threads
    for (ui i = 0; i < pMeta._num_threads; i++) {
        for (ui j = 0; j < dun.getNumEdges(); j++) {
            h[j] += pMeta._total_hash_table[i][j];
            pMeta._total_hash_table[i][j] = 0;
        }
    }
    for (ui i = 0; i < pMeta._num_threads; i++) {
        if (h[i] != 0 && keyPosSize < sizeBound) {
            keyPos[keyPosSize] = i;
            keyPosSize += 1;
        }
    }

    for (ui* partial_match : partial_matches) {
        delete [] partial_match;
    }
}

void worker(
        ParallelProcessingMeta& pMeta, 
        std::list<ui*>& partial_matches, 
        ui mappingSize) {
    if (pMeta._thread_id_ets.local() == -1) {
        pMeta._thread_id_ets.local() = pMeta._next_thread_id.fetch_add(1);
    }
    ui* candidates = _total_candidates[_thread_id_ets.local()];
    ui candCount = 0;
    bool* visited_vertices = _total_visited_vertices[_thread_id_ets.local()];
    std::list<ui*>& extented_matches = _total_extented_matches[_thread_id_ets.local()];
    HashTable h = _total_hash_table[_thread_id_ets.local()];
    VertexID* tmp = _tmp[_thread_id_ets.local()];
    while (!partial_matches.empty()) {
        ui* dataV = partial_matches.front();
        partial_matches.pop_front();

        if (!_nodeInterPos[mappingSize]) {
            if (!_nodeOutPos[mappingSize].empty()) {
                VertexID w = dataV[_nodeOutPos[mappingSize][0]];
                candidates = _outNbors + _outOffset[w];
                candCount = _outOffset[w + 1] - _outOffset[w];
            }
            else if (!_nodeInPos[mappingSize].empty()){
                VertexID w = dataV[_nodeInPos[mappingSize][0]];
                candidates = _inNbors + _inOffset[w];
                candCount = _inOffset[w + 1] - _inOffset[w];
            }
            else {
                VertexID w = dataV[_nodeUnPos[mappingSize][0]];
                candidates = _unNbors + _unOffset[w];
                candCount = _unOffset[w + 1] - _unOffset[w];
            }
        }
        else {
            generateCandidate(dataV, mappingSize, candidates, candCount, tmp);
        }
        // apply the symmetry breaking rules
        ui start_pos = 0;
        if (!_greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[_greaterPos[mappingSize][0]];
            for (int i = 1; i < _greaterPos[mappingSize].size(); ++i) {
                if (dataV[_greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[_greaterPos[mappingSize][i]];
            }
            start_pos = firstPosGreaterThan(candidates, 0, candCount, maxTarget);
        }
        if (!_lessPos[mappingSize].empty()) {
            ui minTarget = dataV[_lessPos[mappingSize][0]];
            for (int i = 1; i < _lessPos[mappingSize].size(); ++i) {
                if (dataV[_lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[_lessPos[mappingSize][i]];
            }
            candCount = firstPosGreaterThan(candidates, 0, candCount, minTarget);
        }
        for (ui i = 0; i < mappingSize; ++i) {
            visited_vertices[dataV[i]] = true;
        }
        for (ui i = start_pos; i < candCount; i++) {
            VertexID v = candidates[i];
            if (visited_vertices[v]) {
                continue;
            }
            dataV[mappingSize] = v;
            if (mappingSize == _tau.nodeOrder.size() - 1) {
                Count cnt = 1;
                for (int j = 0; j < _child.size(); j++) {
                    VertexID cID = _child[j];
                    cnt *= _H[cID][dataV[_childKeyPos[j][0]]];
                }
                if (_isRoot) {
                    if (_orbitType == 0) h[0] += cnt * _aggreWeight[0];
                    else if (_orbitType == 1) {
                        for (int j = 0; j < _aggreV.size(); j++) {
                            VertexID key = dataV[_aggrePos[j]];
                            h[key] += cnt * _aggreWeight[j];
                        }
                    }
                } else {
                    h[dataV[_aggrePos[0]]] += cnt;
                }
            } else {
                ui* new_embedding = new ui[_tau.nodeOrder.size()];
                std::copy(dataV, dataV + _tau.nodeOrder.size(), new_embedding);
                extented_matches.push_back(new_embedding);
            }
        }
        for (ui i = 0; i < mappingSize; i++) {
            visited_vertices[dataV[i]] = false;
        }
        delete [] dataV;
    }
}

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
) {

}

void generateCandidate(
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        VertexID *dataV,
        int pos,
        const std::vector<int> &inPos,
        const std::vector<int> &outPos,
        const std::vector<int> &unPos,
        VertexID **candidate,
        ui *&candCount,
        VertexID *&tmp
) {
    int num = inPos.size() + outPos.size() + unPos.size();
    ui tmpCount;
    if (num > 2) {
        VertexID **arrays = new VertexID *[num];
        ui *counts = new ui[num];
        int p = 0;
        for (int position: inPos) {
            arrays[p] = din.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        for (int position: outPos) {
            arrays[p] = dout.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        for (int position: unPos) {
            arrays[p] = dun.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        ComputeSetIntersection::LeapfrogJoin(arrays, counts, num, tmp, tmpCount);
        candCount[pos] = tmpCount;
        std::swap(candidate[pos], tmp);
        delete[] counts;
        delete[] arrays;
        return;
    }
    VertexID *firstNeighbor;
    ui firstCnt;
    bool firstFlag = true, secondFlag = true;
    for (int position: inPos) {
        ui nbrCnt;
        VertexID *neighbor = din.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstNeighbor = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
    for (int position: outPos) {
        ui nbrCnt;
        VertexID *neighbor = dout.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstNeighbor = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
    for (int position: unPos) {
        ui nbrCnt;
        VertexID *neighbor = dun.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstNeighbor = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
}