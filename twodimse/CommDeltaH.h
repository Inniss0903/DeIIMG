//
// Created by GemWang on 2022/1/7.
//

#ifndef DEIIMG2_COMMDELTAH_H
#define DEIIMG2_COMMDELTAH_H

#include "PairNode.h"
#include <algorithm>
#include <functional>

#define DIFF 1e-40

struct CommDeltaH {
    PairNode pairComms;
    double deltaH;

    CommDeltaH(PairNode pairComms, double deltaH) : pairComms(pairComms) {
        this->deltaH = deltaH;
    }

    bool operator==(const CommDeltaH &c) const {
        if (this == &c) {
            return true;
        }

        return pairComms == c.pairComms && fabs(deltaH - c.deltaH) < DIFF;
    }

    bool operator<(const CommDeltaH &c) const {
        if (std::fabs(deltaH - c.deltaH) >= DIFF)
            return deltaH < c.deltaH;
        return pairComms < c.pairComms;
    }
};

namespace std {

    template<>
    struct hash<CommDeltaH> {
        std::size_t operator()(const CommDeltaH &k) const {
            return std::hash<PairNode>()(k.pairComms);
        }
    };

}


#endif //DEIIMG2_COMMDELTAH_H
