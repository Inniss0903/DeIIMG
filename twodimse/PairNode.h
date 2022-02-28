//
// Created by GemWang on 2022/1/7.
//

#ifndef DEIIMG2_PAIRNODE_H
#define DEIIMG2_PAIRNODE_H

#include <algorithm>

class PairNode {
public:
    int p1;
    int p2;

    bool isValid() const {
        return p1 != p2;
    }

    PairNode(int p1, int p2) {
        if (p1 > p2) {
            std::swap(p1, p2);
        }
        this->p1 = p1;
        this->p2 = p2;
    }

    bool operator==(const PairNode & o) const {
        if (this == &o) {
            return true;
        }
        return p1 == o.p1 && p2 == o.p2;
    }

    bool operator<(const PairNode &o) const {
        return p1 == o.p1 ? p2 < o.p2 : p1 < o.p1;
    }
};

namespace std {

    template<>
    struct hash<PairNode> {
        std::size_t operator()(const PairNode &k) const {
            return k.p2 * 2 + k.p1;
        }
    };

}


#endif //DEIIMG2_PAIRNODE_H
