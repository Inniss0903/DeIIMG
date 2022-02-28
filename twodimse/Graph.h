//
// Created by GemWang on 2022/1/7.
//

#ifndef DEIIMG2_GRAPH_H
#define DEIIMG2_GRAPH_H

#include <unordered_map>
#include <set>
#include <vector>
#include "PairNode.h"

class Graph {
public:
    int numNodes;
    double sumDegrees;

    std::unordered_map<PairNode, double> weights;
    std::unordered_map<int, std::set<int>> connection;

    //此处使用数组而非HashMap的原因是：图的节点的编号连续，如从1-55555，共55555个节点
    std::vector<double> nodeDegree;

    Graph(int numNodes) : weights(3 * numNodes / 4 + 1), connection(3 * numNodes / 4 + 1), nodeDegree(numNodes + 1) {
        this->numNodes = numNodes;
    }

    void write2File(std::string file) {

    }
};


#endif //DEIIMG2_GRAPH_H
