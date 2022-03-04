//
// Created by GemWang on 2022/1/7.
//

#ifndef DEIIMG2_SE_H
#define DEIIMG2_SE_H

#include <cstdio>
#include <unordered_map>
#include <set>
#include <string>
#include <cmath>
#include <vector>
#include <gmpxx.h>

#include "PairNode.h"
#include "Graph.h"
#include "CommDeltaH.h"

using namespace std;

struct TwoDimSE {

    const double sumDegrees;

    const std::size_t initializeCap;

    double oneDimSE;

    double twoDimSE;
    //社区及其包含的节点

    std::unordered_map<int, std::set<int>> communities;
    //节点的度和割边数
    std::vector<double> &volumes;
    std::vector<double> gs;
    //节点之间的割

    std::unordered_map<PairNode, double> &cuts;
    //社区与其相联系的社区
    //connection是双向的，不同于pairNode，更新时应注意双向更新

    std::unordered_map<int, std::set<int>> &connections;
    // 节点之间的△H。使用了额外的存储，待优化（TreeMap）。
    // Treestd::set可以为类按照一定规则排序，排序效率为O(logn)
    //  TreeMap<PairNode, Double> commDeltaHTreeMap;

    std::unordered_map<PairNode, CommDeltaH> commDeltaHMap;

    std::set<CommDeltaH> commDeltaHSet;

    /**
     * 初始化二维结构熵所需的编码树信息
     * 非构造树方法，仅适用于二维结构熵，
     * 因此将graph中保存的原始节点信息直接拿来使用，不再额外使用存储空间
     *
     * @param graph
     */


    TwoDimSE(Graph &graph) : oneDimSE(0.0), twoDimSE(0.0), sumDegrees(graph.sumDegrees),
                             initializeCap(3 * graph.numNodes / 4 + 1),
                             volumes(graph.nodeDegree), gs(graph.nodeDegree), communities(initializeCap),
                             cuts(graph.weights),
                             connections(graph.connection) {
    }


    /**
     * 二维结构熵极小化
     */

#define ASSERT() do {\
    if(commDeltaHSet.size() != commDeltaHMap.size()) { \
        printf("set size: %lu, map size: %lu\n", commDeltaHSet.size(), commDeltaHMap.size()); \
        assert(false); \
    } \
} while (0)


    void min2dSE(string saveFilePath, bool doPrintNDI, bool doSave) {
        initEncodingTree();
        twoDimSE = oneDimSE;

        auto maxCommDeltaH = --commDeltaHSet.end();

        //找到最大的△H，merge这两个节点，直到不满足merge的条件
        while (maxCommDeltaH->deltaH > 0 && !commDeltaHSet.empty()) {
            const PairNode &comms = maxCommDeltaH->pairComms;
            double deltaH = maxCommDeltaH->deltaH;
            twoDimSE -= deltaH;
            updateCommunities(*maxCommDeltaH);
//            ASSERT();
            maxCommDeltaH = --commDeltaHSet.end();
        }

        printf("done!");
        //完成划分后的其他操作
        if (doSave)
            saveResult(saveFilePath);

        //输出解码信息
        if (doPrintNDI)
            ndiInfo();
    }


    void min2dSE(bool doPrintNDI) {
        min2dSE(" ", doPrintNDI, false);
    }


    void min2dSE(string saveFilePath, bool doPrintNDI) {
        min2dSE(saveFilePath, doPrintNDI, true);
    }


    /**
     * 节点merge以后更新社区信息
     *
     * @param commDeltaH
     */

    void updateCommunities(const CommDeltaH &commDeltaH) {
        const PairNode &comms = commDeltaH.pairComms;
        double deltaH = commDeltaH.deltaH;
        int commLeft = comms.p1;
        int commRight = comms.p2;

        //更新新社区的度和割
        double vi = volumes[commLeft];
        double gi = gs[commLeft];
        double vj = volumes[commRight];
        double gj = gs[commRight];

        volumes[commLeft] = vi + vj;
        gs[commLeft] = gi + gj - 2 * cuts[comms];
        volumes[commRight] = 0.0;
        gs[commRight] = 0.0;

        //两社区融合为一个新的社区，同时切断两社区之间的联系
        communities[commLeft].insert(communities[commRight].begin(), communities[commRight].end());
        communities.erase(commRight);
        //printf("communities size: %lu\n", communities.size());
        commDeltaHMap.erase(comms);
        commDeltaHSet.erase(commDeltaH);
//        ASSERT();
        connections[commLeft].erase(commRight);
        connections[commRight].erase(commLeft);
        cuts.erase(comms);

        //更新与comLeft和comRight相关社区的△H和cut
        updateCutAndDeltaH(commLeft, commRight);

    }


    /**
     * 更新与comLeft和comRight相关社区的△H以及cut
     * 同时也会更新connection的信息
     *
     * @param commLeft
     * @param commRight
     */

    void updateCutAndDeltaH(int commLeft, int commRight) {
        const std::set<int> &connLeft = connections[commLeft];
        std::set<int> &connRight = connections[commRight];

        double Vi = volumes[commLeft];
        double Gi = gs[commLeft];
        double Gk;
        double Vk;
        double Gx;
        double newDelta;
        //遍历与社区commLeft相关联的社区
        for (int k : connLeft) {
            double cutIk;
            PairNode *pairLeftAndK = new PairNode(commLeft, k);
            if (connRight.find(k) != connRight.end()) {    //若社区k与commLeft和commRight均有关联
                auto *pairRightAndK = new PairNode(commRight, k);
                cutIk = cuts[*pairLeftAndK] + cuts[*pairRightAndK];
                connRight.erase(k);
//                commDeltaHSet.erase(new CommDeltaH(pairRightAndK, commDeltaHMap[pairRightAndK))];
                auto found = commDeltaHMap.find(*pairRightAndK);
                if (found != commDeltaHMap.end()) {
                    commDeltaHSet.erase(found->second);
                    commDeltaHMap.erase(found);
//                    ASSERT();
                }
                cuts.erase(*pairRightAndK);
                connections[k].erase(commRight);
            } else {
                cutIk = cuts[*pairLeftAndK];
            }
            Gk = gs[k];
            Vk = volumes[k];
            Gx = Gi + Gk - 2 * cutIk;
            newDelta = computeDeltaH(Vi, Vk, Gi, Gk, Gx, sumDegrees);

            //更新与cuts和△H相关的存储
            cuts[*pairLeftAndK] = cutIk;
            auto found = commDeltaHMap.find(*pairLeftAndK);
            if (found != commDeltaHMap.end()) {
                commDeltaHSet.erase(found->second);
                commDeltaHMap.erase(found);
//                ASSERT();
            }
            CommDeltaH *newDeltaH = new CommDeltaH(*pairLeftAndK, newDelta);
            commDeltaHSet.insert(*newDeltaH);
            commDeltaHMap.insert({*pairLeftAndK, *newDeltaH});
//            ASSERT();
            //此处connection不用更新

        }
        //遍历融合前与社区commRight相关联但与commLeft不关联的社区
        for (int k : connRight) {
            PairNode *pairRightAndK = new PairNode(commRight, k);
            double cutJk = cuts[*pairRightAndK];
            Vk = volumes[k];
            Gk = gs[k];
            Gx = Gi + Gk - 2 * cutJk;
            newDelta = computeDeltaH(Vi, Vk, Gi, Gk, Gx, sumDegrees);
            //更新与cuts和△H相关的存储
            PairNode *pairLeftAndK = new PairNode(commLeft, k);
            cuts[*pairLeftAndK] = cutJk;
            cuts.erase(*pairRightAndK);
            CommDeltaH *commDeltaH = new CommDeltaH(*pairLeftAndK, newDelta);
            auto found = commDeltaHMap.find(*pairRightAndK);
            if (found != commDeltaHMap.end()) {
                commDeltaHSet.erase(found->second);
                commDeltaHMap.erase(found);
//                ASSERT();
            }
            commDeltaHSet.insert(*commDeltaH);
            commDeltaHMap.insert({*pairLeftAndK, *commDeltaH});
//            ASSERT();
            connections[commLeft].insert(k);
            connections[k].insert(commLeft);
            connections[k].erase(commRight);
        }

        connRight.clear();
    }

    /**
     * 初始化编码树
     * 根节点下有n个社区，每个社区只包含一个自身节点
     */

    void initEncodingTree() {
        //计算社区节点之间的deltaH
        for (auto &p : cuts) {
            double vi = volumes[p.first.p1];
            double vj = volumes[p.first.p2];
            double gi = vi;
            double gj = vj;
            double gx = vi + vj - 2 * cuts[p.first];
            double deltaH = computeDeltaH(vi, vj, gi, gj, gx, sumDegrees);
            CommDeltaH *commDeltaH = new CommDeltaH(p.first, deltaH);
            commDeltaHMap.insert({p.first, *commDeltaH});
            commDeltaHSet.insert(*commDeltaH);
//            ASSERT();
        }

        //计算一维结构熵，并初始化社区
        for (int i = 1; i < volumes.size(); i++) {
            if (volumes[i] > 0.0) {
                //一维编码树的的每个社区就是图中的每个顶点
                int constI = i;
                communities[i] = std::set<int>{constI};
                //printf("communities size: %lu\n", communities.size());
                oneDimSE -= (volumes[i] / sumDegrees) * log2(volumes[i] / sumDegrees);
            }
        }

        printf("eliminate uncertainty.....\n");

    }

    /**
     * 输出ndi的相关信息
     */


    void ndiInfo() const {
        printf("*****************************************\n");
        printf("The One and Two dimension SE: %f, %f\nDecoding Information : %f\n",
               oneDimSE, twoDimSE, oneDimSE - twoDimSE);
        double ndi = (oneDimSE - twoDimSE) / oneDimSE;
        printf("The Normalized Decoding Information is %f\n", ndi);
        printf("The comms size is %lu\n", communities.size());
    }

    /**
     * 将二维结构熵极小化划分的社区保存起来
     *
     * @param fileName
     * @throws IOException
     */


    void saveResult(string fileName) {
        FILE *f = fopen(fileName.c_str(), "w");
        for (auto &res: communities) {
            for (int i : res.second) {
                fprintf(f, "%d\t", i);
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }

    /**
     * 计算两个节点之间的△H
     * 论文《Structural information and dynamical complexity of networks》37页
     *
     * @param vi
     * @param vj
     * @param gi
     * @param gj
     * @param gx
     * @param sum_d
     * @return
     */
    double computeDeltaH(double vi, double vj, double gi, double gj, double gx, double sum_d) {
        mpf_class a1(vi * log2(vi));
        mpf_class a2(vj * log2(vj));
        mpf_class a3((vi + vj) * log2(vi + vj));
        mpf_class a4(gi * log2(vi / sum_d));
        mpf_class a5(gj * log2(vj / sum_d));
        mpf_class a6(gx * log2((vi + vj) / sum_d));

        mpf_class b1 = a1 + a2;
        mpf_class b2 = b1 - a3;
        mpf_class b3 = b2 - a4;
        mpf_class b4 = b3 - a5;
        mpf_class b5 = b4 + a6;

        return b5.get_d() / sum_d;
    }

    void clear() {
        communities.clear();
        volumes.clear();
        gs.clear();
        cuts.clear();
        connections.clear();
        commDeltaHMap.clear();
        commDeltaHSet.clear();
    }
};

#endif