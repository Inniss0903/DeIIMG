#include <cstdio>
#include <chrono>
#include <string>
#include <io.h>
#include <se.h>
#include <ctime>
#include <set>
#include <opencv2/opencv.hpp>
#include <imgseg.h>


using namespace std;
using namespace cv;


void timespec_subtract(struct timespec *result, struct timespec *x, struct timespec *y) {
/* Perform the carry for the later subtraction by updating y. */
    if (x->tv_nsec < y->tv_nsec) {
        long nsec = (y->tv_nsec - x->tv_nsec) / 1000000000 + 1;
        y->tv_nsec -= 1000000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_nsec - y->tv_nsec > 1000000000) {
        long nsec = (x->tv_nsec - y->tv_nsec) / 1000000000;
        y->tv_nsec += 1000000000 * nsec;
        y->tv_sec -= nsec;
    }

/* Compute the time remaining to wait.
   tv_nsec is certainly positive. */
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_nsec = x->tv_nsec - y->tv_nsec;
}

void printTimeLog(struct timespec *diff) {
    long second = diff->tv_sec;
    long ms = diff->tv_nsec / 1000000 + 1;
    printf("time cost is %lds and %ldms\n", second, ms);
}
void splitString(const std::string &s, std::vector<std::string> &v, const std::string &c) {
    std::string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while (std::string::npos != pos2) {
        v.push_back(s.substr(pos1, pos2 - pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1));
}

void segDir(const char *dirPath) {
    imgseg imgseg;
    vector<double> tSize = {1.0};
    DIR *dir = opendir(dirPath);
    struct dirent *ent;
    if (dir == nullptr) {
        cout << "error: can not open dir!!!" << endl;
    }
    /* print all the files and directories within directory */
    Mat img;
    while ((ent = readdir(dir)) != nullptr) {
        string name = ent->d_name;
        if (name != "." && name != ".." && name != "Thumbs.db") {
            vector<string> splites;
            splitString(name, splites, ".");
            string imgName = splites[0];
            cout << imgName << endl;
            img = imread(format(dirPath) + "/" + name);

            for (double t: tSize) {
                t1 = t;
                Graph g = imgseg.constructGraph(img);
                TwoDimSE se(g);
                se.min2dSE(true);
                imgseg.deIIMG_2D(se.communities, img, false, false, "/home/gemwang/segRes/test/" + imgName);
                g.clear();
                se.clear();
                imgseg.clear();
            }
        }
    }

    closedir(dir);
}



void imgsegm(const string &path) {
    imgseg imgseg;
    segDir(path.c_str());
//    cv::Mat img = cv::imread(path+"376086.jpg", IMREAD_COLOR);
//    Graph g = imgseg.constructGraph(img);
//    TwoDimSE se(g);
//    se.min2dSE("t2d", true);
//    cv::Mat res2D = imgseg.deIIMG_2D(se.communities, img, false, false,"/home/gemwang/segRes/test/33");

//    Graph g2 = imgseg.constructGraphBy2D(se.communities, 5);
//    TwoDimSE se2(g2);
//    se2.min2dSE("3d", true);
//    cv::Mat res3D = imgseg.deIIMG_3D(se.communities, se2.communities, img);
//    cv::imshow("2d", res2D);
//    imshow("3d", res3D);
//    cv::waitKeyEx();
}

int main(int argc, char **argv) {
    std::string path = "/Users/gem/PyProject/SE_image_seg/image_graph";
    //std::string path = "/Users/gem/PyProject/SE_image_seg/DataTxt/Lymph6Graph";
    struct timespec end, start, diff;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

//    Graph g = io::getUndirGraphFromFile(path);
//    TwoDimSE twoDimSe(g);
//    twoDimSe.min2dSE("2D", true);
    path = "/home/gemwang/superpixel-benchmark/data/BSDS500/images/test/";
//    imgsegm(path);
    segDir(path.c_str());


    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    timespec_subtract(&diff, &end, &start);
    printTimeLog(&diff);

    return 0;
}

