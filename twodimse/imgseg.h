//
// Created by GemWang on 2022/1/11.
//

#ifndef DEIIMG2_IMGSEG_H
#define DEIIMG2_IMGSEG_H

#include <opencv2/opencv.hpp>
#include <unordered_map>
#include <set>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <gmpxx.h>
#include <gmp.h>
#include <math.h>

#include <Graph.h>
#include <Edge.h>
#include <dirent.h>
#include "se.h"

#define N_CHANNELS 3
#define DIFF 1e-40
#define RANDOM(x) rand()%(x)
#define CONTAINS(boundaries, p) find(boundaries.begin(), boundaries.end(), p) != boundaries.end()
#define COMMUNITY_CONTAINS(p) community.find(p + 1) != community.end()

using namespace std;
using namespace cv;


extern const vector<int> COLOR_RED;
extern int k;
extern double t1;
extern double t2;
extern double cRatio;


struct imgseg {

    struct edgeDescComparator {
        bool operator()(const Edge &e1, const Edge &e2) {
            if (std::fabs(e1.weight - e2.weight) > DIFF)
                return e2.weight < e1.weight;
            return e1.seqID < e2.seqID;
        }
    };

    vector<vector<int>> commAverage;
    vector<set<Edge, edgeDescComparator>> edgePQs;


    /**
     * 把图像转化为图结构的数据
     *
     * @param img
     * @return
     */
    Graph constructGraph(const Mat &img) {
        int height = img.rows;
        int width = img.cols;
        int boundary[4] = {0, 0, width, height};

        Graph graph(height * width);
        //图像中的坐标和对应到数学中的坐标轴是不一样的
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
//                System.out.printf("Processing at %d , %d%n", y, x);
                pixelINAV(graph, img, boundary, x, y);
            }
        }

//        printf("graph weights size %d \t con size %d\n", graph.weights.size(), graph.connection.size());
        return graph;
    }

    /**
     * 图像构建为图结构的数据时，每个像素点的交互
     * 一般我们以图像左上角为原点建立坐标系，水平为x，垂直为y
     *
     * @param graph
     * @param img
     * @param boundary
     * @param x
     * @param y
     */
    void pixelINAV(Graph &graph, const Mat &img, const int boundary[], int x, int y) {
        int leftBDY = boundary[0];
        int topBDY = boundary[1];
        int rightBDY = boundary[2];
        int bottomBDY = boundary[3];

        int width = img.cols;
        int fromNode = y * width + x;
        int seqID = 0;

        set<Edge, edgeDescComparator> edges;

        for (int i = 1; i <= k; i++) {
            int toNode;
            double weight;
            // x + i, y
            if (x + i < rightBDY) {
                toNode = y * width + (x + i);
                weight = metric(img, x, y, x + i, y);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

            // x, y + i
            if (y + i < bottomBDY) {
                toNode = (y + i) * width + x;
                weight = metric(img, x, y, x, y + i);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

            // x + i, y + i
            if ((x + i < rightBDY) && (y + i < bottomBDY)) {
                toNode = (y + i) * width + (x + i);
                weight = metric(img, x, y, x + i, y + i);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

            // x + i, y - i
            if ((x + i < rightBDY) && (y - i > topBDY - 1)) {
                toNode = (y - i) * width + (x + i);
                weight = metric(img, x, y, x + i, y - i);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

            // x , y - i
            if (y - i > topBDY - 1) {
                toNode = (y - i) * width + x;
                weight = metric(img, x, y, x, y - i);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

            // x - i, y - i
            if ((x - i > leftBDY - 1) && (y - i > topBDY - 1)) {
                toNode = (y - i) * width + (x - i);
                weight = metric(img, x, y, x - i, y - i);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

            // x - i, y
            if (x - i > leftBDY - 1) {
                toNode = y * width + (x - i);
                weight = metric(img, x, y, x - i, y);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

            // x - i, y + i
            if ((x - i > leftBDY - 1) && (y + i < bottomBDY)) {
                toNode = (y + i) * width + (x - i);
                weight = metric(img, x, y, x - i, y + i);
                edges.insert(Edge(fromNode, toNode, weight, seqID++));
            }

        }

        selectNMostSimilar(graph, edges, edges.size() / 2);
    }

    /**
     * 选取权重较大的边，数量为n
     *
     * @param graph 信息系统G
     * @param edges 某个节点在图中的边集合
     * @param n     选取边的数量
     */
    void selectNMostSimilar(Graph &graph, set<Edge, edgeDescComparator> edges, int n) {
        double &sumDegrees = graph.sumDegrees;
        unordered_map<PairNode, double> &weights = graph.weights;
        unordered_map<int, set<int>> &connection = graph.connection;
        auto &nodeDegree = graph.nodeDegree;

        auto iter = edges.begin();
        for (int i = 0; i < n && i < edges.size(); ++iter, ++i) {
            Edge edge = *iter;
            //加一为了使其从1开始编码
            int start = edge.start + 1;
            int end = edge.end + 1;
            double weight = edge.weight;

//            printf("%d -> %d : %f\n", start, end, weight);
            PairNode pair(start, end);
            if (weights.find(pair) == weights.end()) {
                //边及其对应的权重
                weights.insert({pair, edge.weight});
                connection[start].insert(end);
                connection[end].insert(start);
                //每一个节点的度数
                nodeDegree[start] += weight;
                nodeDegree[end] += weight;
                //无向图，所以权重翻倍
                sumDegrees += 2 * weight;
            }
        }

        graph.sumDegrees = sumDegrees;
    }

    /**
     * 构建G**
     *
     * @param communities2D G*的划分结果
     * @param n             构建图时每个节点有其他节点建立联系的数量
     * @return 返回G**
     */
    Graph constructGraphBy2D(const unordered_map<int, set<int>> &communities2D, int n) {
        cout << "construct graph G**" << endl;

        int size = communities2D.size();
        Graph graph(size);
        if (edgePQs.empty()) {
            edgePQs.resize(size);
        }

        for (int fromNode = 0; fromNode < size; fromNode++) {
            if (!edgePQs[fromNode].empty()) {
                //经过首次计算后过无需再次计算，大大缩短运行时间。以空间换时间
                selectNMostSimilar(graph, edgePQs[fromNode], n);
                continue;
            }

            set<Edge, edgeDescComparator> edges;
            for (int toNode = 0; toNode < size; toNode++) {
                if (fromNode != toNode) {
                    double weight = metric(fromNode, toNode);
                    edges.insert(Edge(fromNode, toNode, weight));
                }
            }
            //todo deleted?
            edgePQs[fromNode] = edges;
            selectNMostSimilar(graph, edges, n);
        }

        cout << "done!" << endl;
        return graph;
    }


    double metric(const Mat &img, int x1, int y1, int x2, int y2) {
        auto r1 = img.at<Vec3b>(y1, x1)[0];
        auto r2 = img.at<Vec3b>(y2, x2)[0];
        auto g1 = img.at<Vec3b>(y1, x1)[1];
        auto g2 = img.at<Vec3b>(y2, x2)[1];
        auto b1 = img.at<Vec3b>(y1, x1)[2];
        auto b2 = img.at<Vec3b>(y2, x2)[2];
        double dst = sqrt(
                square(r1 - r2)
                + square(g1 - g2)
                + square(b1 - b2)
        );

        return pow(2, -dst / t1);
    }

    double metric(int node1, int node2) {
        vector<int> c1 = commAverage[node1];
        vector<int> c2 = commAverage[node2];

        double dst = sqrt(
                square(c1[0] - c2[0])
                + square(c1[1] - c2[1])
                + square(c1[2] - c2[2])
        );

        return pow(2, -square(dst / t2));
    }

    int square(int value) {
        return value * value;
    }

    double square(double value) {
        return value * value;
    }


    Mat deIIMG_2D(const unordered_map<int, set<int>> &communities, const Mat &img) {
        return deIIMG_2D(communities, img, false, true, "");
    }

    Mat deIIMG_2D(const unordered_map<int, set<int>> &communities, const Mat &img, bool randomColor, bool showBDY,
                  string segFile) {
        Mat outputImage(img.rows, img.cols, img.type());
        commAverage.resize(communities.size());
        vector<vector<int>> segMap(img.rows, vector<int>(img.cols, 0));
        int index = -1;
        for (auto &comm: communities) {
            vector<int> boundaries = findBoundary(comm.second, img.cols);
            //分割区域的颜色表示
            commAverage[++index] = computeCommAve(img, comm.second);
            vector<int> regionColor = randomColor ? randomRegionColor() : commAverage[index];

            for (auto &p: comm.second) {
                int loc = p - 1;
                int y = loc / img.cols;
                int x = loc % img.cols;
                Vec3b pix;
                //是否需要展示区域的边界条件
                if (showBDY && CONTAINS(boundaries, p)) {
                    for (int c = 0; c < N_CHANNELS; ++c)
                        pix[c] = COLOR_RED[c];
                } else {
                    for (int c = 0; c < N_CHANNELS; ++c)
                        pix[c] = regionColor[c];
                }
                segMap[y][x] = index + 1;
                outputImage.at<Vec3b>(y, x) = pix;
            }
        }

        if (segFile != "") {
            char spliteChar = ' ';
            write2file(segMap, segFile, spliteChar);
        }

        return outputImage;
    }


    Mat deIIMG_3D(const unordered_map<int, set<int>> &communities2D, const unordered_map<int, set<int>> &communities3D,
                  const Mat &img) {
        return deIIMG_3D(communities2D, communities3D, img, false, false, false, "");
    }

    Mat deIIMG_3D(const unordered_map<int, set<int>> &communities2D, const unordered_map<int, set<int>> &communities3D,
                  const Mat &img, bool showBDY, bool is2S, bool randomColor, string segFile) {
        Mat outputImage(img.rows, img.cols, img.type());
        vector<vector<int>> segMap(img.rows, vector<int>(img.cols, 0));
        vector<set<int>> comms2D;
        int label = 0;

        for (auto &comm: communities2D) {
            comms2D.push_back(comm.second);
        }

        for (auto &community: communities3D) {
            //对应到图像像素点的社区
            set<int> realCommunity;
            //每个分割区域的颜色表示
            vector<int> regionColor(3);
            for (int comm2DIndex: community.second) {
                realCommunity.insert(comms2D[comm2DIndex - 1].begin(), comms2D[comm2DIndex - 1].end());
            }

            //分割区域边界
            vector<int> boundaries;
            if (showBDY) boundaries = findBoundary(realCommunity, img.cols);

            if (randomColor) {
                regionColor = randomRegionColor();
            } else if (is2S) {
                regionColor = compute2SCommAve(community.second);
            } else {
                regionColor = computeCommAve(img, realCommunity);
            }

            //给每个分割区域中的像素点重新赋值
            for (int p: realCommunity) {
                int loc = p - 1;
                int y = loc / img.cols;
                int x = loc % img.cols;
                Vec3b pix;
                if (showBDY && CONTAINS(boundaries, p)) {
                    for (int c = 0; c < N_CHANNELS; ++c)
                        pix[c] = COLOR_RED[c];
                } else {
                    for (int c = 0; c < N_CHANNELS; ++c)
                        pix[c] = regionColor[c];
                }
                segMap[y][x] = ++label;
                outputImage.at<Vec3b>(y, x) = pix;
            }
        }

        if (segFile != "") {
            char sp = ' ';
            write2file(segMap, segFile, sp);
        }

        return outputImage;
    }

    vector<int> compute2SCommAve(const set<int> &community) {
        vector<int> color(3);
        for (int c: community) {
            for (int i = 0; i < N_CHANNELS; i++) {
                color[i] += commAverage[c - 1][i];
            }
        }

        for (int i = 0; i < N_CHANNELS; i++) {
            color[i] /= community.size();
        }

        return color;
    }

    /**
     * 根据原图像计算某个划分区域的像素RGB均值
     *
     * @param img
     * @param community
     * @return
     */
    vector<int> computeCommAve(const Mat &img, const set<int> &community) {
        vector<int> color(3);
        for (int pix: community) {
            int loc = pix - 1;
            int y = loc / img.cols;
            int x = loc % img.cols;
            Vec3b pc = img.at<Vec3b>(y, x);
            for (int i = 0; i < N_CHANNELS; i++) {
                color[i] += (int) pc[i];
            }
        }

        for (int &c: color) {
            c /= community.size();
        }

        return color;
    }

    /**
     * 随机生成一个RGB值
     *
     * @return
     */
    vector<int> randomRegionColor() {
        vector<int> rColor(3);
        for (int &c: rColor) {
            c = RANDOM(255);
        }
        return rColor;
    }


    vector<int> findBoundary(const set<int> &community, int width) {
        vector<int> boundaries;
        for (int p: community) {
            if (!(COMMUNITY_CONTAINS(p + 1) && COMMUNITY_CONTAINS(p - 1)
                  && COMMUNITY_CONTAINS(p + width) && COMMUNITY_CONTAINS(p - width))) {
                boundaries.push_back(p);
            }
        }

        return boundaries;
    }

    void write2file(vector<vector<int>> segMap, const string& file, char spl) {
        ofstream fs(file, ios::trunc);
        int m = segMap.size();
        int n = segMap[0].size();

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n - 1; ++j) {
                fs << segMap[i][j] << spl;
            }
            fs << segMap[i][n - 1] << "\n";
        }

        fs.close();
    }

    void clear() {
        commAverage.clear();
    }


};


#endif //DEIIMG2_IMGSEG_H
