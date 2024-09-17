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