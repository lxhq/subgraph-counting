data_graph="soc-gowalla"
vertex_set="5voc"

moniter="/usr/bin/time"
prefix="/home/ubuntu/Documents/workspace/subgraph-counting"
binary="$prefix/build/executable/scope.out"
data_graph="$prefix/exp/data_graph/"$data_graph".txt"
pattern_graph="$prefix/exp/pattern_graph/"$vertex_set""
result_dir="$prefix/result/"$data_graph"-parallel/"$vertex_set""
mode="parallel"

# thread_nums=("19" "15" "10" "5" "1")
# node_partition=("100" "500" "1000" "5000")
# prefix_partition=("100" "500" "1000" "5000")
thread_nums=("1")
node_partition=("1000")
prefix_partition=("500")

for n in "${thread_nums[@]}"; do
    for np in "${node_partition[@]}"; do
        for pp in "${prefix_partition[@]}"; do
            echo "Thread: $n, Node: $np, Prefix: $pp"
            $moniter $binary -d $data_graph -q $pattern_graph -r $result_dir -m $mode -n $n -np $np -pp $pp -b > ""$n"_"$np"_"$pp"_"$vertex_set".txt"
        done
    done
done

