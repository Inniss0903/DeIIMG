//
// Created by GemWang on 2022/1/7.
//

#ifndef DEIIMG2_EDGE_H
#define DEIIMG2_EDGE_H

#include<cmath>
#include <functional>

#define DIFF 1e-40

struct Edge {
    int start;
    int end;
    double weight;
    int seqID;

    Edge(int start, int end, double weight) {
        this->start = start;
        this->end = end;
        this->weight = weight;
    }

    Edge(int start, int end, double weight, int seqID) {
        this->start = start;
        this->end = end;
        this->weight = weight;
        this->seqID = seqID;
    }

    bool operator==(const Edge &e) const {
        //自反性
        if (this == &e) {
            return true;
        }

        if (start != e.start) {
            return false;
        } else if (end != e.end) {
            return false;
        }
        return fabs(weight - e.weight) < DIFF;
    }

    bool operator<(const Edge &edge) const {
        if (std::fabs(this->weight - edge.weight) > DIFF)
            return this->weight < edge.weight;
        if (this->start == edge.start) {
            return this->end < edge.end;
        }
        return this->start < edge.start;
    }
};

namespace std {

    template<>
    struct hash<Edge> {
        std::size_t operator()(const Edge &k) const {
            return k.end;
        }
    };

}

#endif //DEIIMG2_EDGE_H
