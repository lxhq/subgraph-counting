data_graph="soc-gowalla"

moniter="/usr/bin/time"
prefix="/home/ubuntu/Documents/workspace/subgraph-counting"
binary="$prefix/build/executable/scope.out"
data_graph_path="$prefix/exp/data_graph/"$data_graph".txt"
pattern_graph_dir="$prefix/exp/pattern_graph/executed_tree_patterns/"

# single thread
mode="single"
echo "$moniter -v $binary -d $data_graph_path -q $pattern_graph_dir -m $mode -b > ""$data_graph"_single.txt" 2>&1 "

# multi-thread
mode="parallel"
thread_nums=("1" "5" "10" "15" "19")
node_partition="100"
prefix_partition="50"
patterns_parallel_size="10"

for n in "${thread_nums[@]}"; do
    echo "$moniter -v $binary -d $data_graph_path -q $pattern_graph_dir -m $mode -n $n -np $node_partition -pp $prefix_partition -patp $patterns_parallel_size -b > ""$data_graph"_multi_"$n".txt" 2>&1"
done