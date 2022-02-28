//
// Created by GemWang on 2022/1/8.
//

#ifndef DEIIMG2_IO_H
#define DEIIMG2_IO_H

#include <Graph.h>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <set>

namespace io {
    Graph getUndirGraphFromFile(const std::string &filePath) {
        std::ifstream fs(filePath);
        int numNodes;
        fs >> numNodes;
        Graph g(numNodes);
        std::unordered_map<PairNode, double> &weights = g.weights;
        std::unordered_map<int, std::set<int>> &connection = g.connection;
        auto &nodeDegree = g.nodeDegree;
        double sumDegrees = 0.;
        while (fs) {
            int start, end;
            double weight;
            fs >> start >> end >> weight;
            PairNode pair(start, end);
            if (weights.find(pair) == weights.end()) {
                weights.insert({pair, weight});
                connection[start].insert(end);
                connection[end].insert(start);
                nodeDegree[start] += weight;
                nodeDegree[end] += weight;
                sumDegrees += 2 * weight;
            }
        }
        g.sumDegrees = sumDegrees;
        printf("graph sumofdeg %f \t weights size %lu \t con size %lu\n", g.sumDegrees, g.weights.size(), g.connection.size());
        return g;
    }
}

#endif //DEIIMG2_IO_H
