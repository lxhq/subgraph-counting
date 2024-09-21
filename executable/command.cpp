//
// Created by Qiyan LI on 2022/8/30.
//

#include "command.h"

Command::Command(int argc, char **argv) : CommandParser(argc, argv){
    optionsKey[OptionKeyword::QueryGraphPath] = "-q";
    optionsKey[OptionKeyword::DataGraphPath] = "-d";
    optionsKey[OptionKeyword::BatchQuery] = "-b";
    optionsKey[OptionKeyword::ResultPath] = "-r";
    optionsKey[OptionKeyword::ShareNode] = "-share";
    optionsKey[OptionKeyword::TrianglePath] = "-t";
    booleanOptionValue[OptionKeyword::BatchQuery] = false;
    booleanOptionValue[OptionKeyword::ShareNode] = false;
    optionsKey[OptionKeyword::ExecutionMode] = "-m";
    optionsKey[OptionKeyword::NumThreads] = "-n";
    optionsKey[OptionKeyword::NodePartitionSize] = "-np";
    optionsKey[OptionKeyword::PrefixPartitionSize] = "-pp";
    optionsKey[OptionKeyword::PatternsParallelSize] = "-patp";
    processOptions();
}

void Command::processOptions() {
    optionsValue[OptionKeyword::QueryGraphPath] = getCommandOption(optionsKey[OptionKeyword::QueryGraphPath]);
    optionsValue[OptionKeyword::DataGraphPath] = getCommandOption(optionsKey[OptionKeyword::DataGraphPath]);
    optionsValue[OptionKeyword::TrianglePath] = getCommandOption(optionsKey[OptionKeyword::TrianglePath]);
    optionsValue[OptionKeyword::ResultPath] = getCommandOption(optionsKey[OptionKeyword::ResultPath]);
    booleanOptionValue[OptionKeyword::BatchQuery] = commandOptionExists(optionsKey[OptionKeyword::BatchQuery]);
    booleanOptionValue[OptionKeyword::ShareNode] = commandOptionExists(optionsKey[OptionKeyword::ShareNode]);
    optionsValue[OptionKeyword::ExecutionMode] = getCommandOption(optionsKey[OptionKeyword::ExecutionMode]);
    std::string num_thread = getCommandOption(optionsKey[OptionKeyword::NumThreads]);
    if (num_thread.empty()) {
        intOptionValue[OptionKeyword::NumThreads] = 10;
    }
    else {
        intOptionValue[OptionKeyword::NumThreads] = std::stoi(num_thread);
    }
    std::string node_partition_size = getCommandOption(optionsKey[OptionKeyword::NodePartitionSize]);
    if (node_partition_size.empty()) {
        intOptionValue[OptionKeyword::NodePartitionSize] = 500;
    }
    else {
        intOptionValue[OptionKeyword::NodePartitionSize] = std::stoi(node_partition_size);
    }
    std::string prefix_partition_size = getCommandOption(optionsKey[OptionKeyword::PrefixPartitionSize]);
    if (prefix_partition_size.empty()) {
        intOptionValue[OptionKeyword::PrefixPartitionSize] = 100;
    }
    else {
        intOptionValue[OptionKeyword::PrefixPartitionSize] = std::stoi(prefix_partition_size);
    }
    std::string patterns_parallel_size = getCommandOption(optionsKey[OptionKeyword::PatternsParallelSize]);
    if (patterns_parallel_size.empty()) {
        intOptionValue[OptionKeyword::PatternsParallelSize] = 10;
    } else {
        intOptionValue[OptionKeyword::PatternsParallelSize] = std::stoi(patterns_parallel_size);
    }
}