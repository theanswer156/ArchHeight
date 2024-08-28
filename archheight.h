#ifndef ARCHHEIGHT_H
#define ARCHHEIGHT_H
#include <math.h>
#include <vector>
#include <QPointF>
using namespace std;
struct Point {
    double x;
    double y;
    //      ????????????????????????????????????????????????????????????????
    Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
    operator QPointF()const{
        return QPointF(x,y);
    }
    Point operator-()const {
        return Point(-x, -y);
    }
    double dist(const Point& p1, const Point& p2) {
        return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
    }
    friend Point operator*(const Point& p, double scalar) {
        return { p.x * scalar,p.y * scalar };
    }
    friend Point operator*(double scalar, const Point& p) {
        return p * scalar;
    }
    friend Point operator+(double scalar, const Point& p) {
        return { p.x + scalar,p.y + scalar };
    }
    friend Point operator+(const Point& p, double scalar) {
        return scalar + p;
    }
    friend Point operator+(const Point& p1, const Point& p2) {
        return { p1.x + p2.x,p1.y + p2.y };
    }
    friend Point operator-(const Point& p1, const Point& p2) {
        return { p1.x - p2.x,p1.y - p2.y };
    }
    bool operator==(const Point& other) const {
        return (x == other.x) && (y == other.y);
    }
    Point& operator+=(const Point& p1) {
        x += p1.x;
        y += p1.y;
        return *this;
    }
    Point& operator=(const Point& other) {
        if (this != &other) {
            x = other.x;
            y = other.y;
        }
        return *this;
    }
};
class ArchHeight
{
public:
    ArchHeight();
    ArchHeight(const double height,const size_t index);
    ~ArchHeight();

    vector<vector<Point>> outPiecePoint() const;
    vector<vector<Point>> outDesPoint() const;
    vector<vector<Point>> outSrcPoint() const;
    vector<vector<Point>> outAdaptPoint() const;
    size_t outIndex() const;

private:
    void doInit();
    void setSrcPoint();
    void getBezier();
    void getBezierPrime();
    void getDesPoint();
    void doIsoHeightSeg();

    double dotProduct(const Point& p1, const Point& p2);
    double velocityAngleChange(const Point& point_1,const Point& point_2,const Point& point_3);
    //  速度的角度变化  cos(α)  用余弦能更好的表示变化情况  取反后+1  这样得出来的值在(0,2)之间  值越大  表示变化越大
    double circumRadius(const Point& point_1,const Point& point_2,const Point& point_3);          //  计算三个点的额外接圆半径
    double doComputeWeightedFunction(const size_t& index,const double& begintime,const double& endtime);
    void doAdaptiveSampling();      //  NOTE:曲线的自适应采样(非均匀采样) 
    void doComputeChangeTime();
    void doComputePiecePoint();     //  NOTE:根据拱高分段结果，计算拱高分段点
    void doComputeAdaptPoint();
    double doComputePointLineDist(const size_t& index, const double& begintime, const double& endtime, const double& heighttime);//  NOTE:计算点到直线的距离

    double doComputeArchHeight(const size_t& index, const double& begintime, const double& endtime);//  NOTE:计算两点的拱高
    double doComputeArchHeightTime(const size_t& index, const double& begintime, const double& endtime);
    double Distance(const Point&, const Point&);
    Point doComputePoint(const size_t& index, const double& t);    //  NOTE:根据时间计算点

public:
    vector<int> m_ivecPascalTri = { 0,1,3,3,1 };

private:
    size_t m_iIndex = 2;
    double m_dHeight = 5;
    double m_dPrecis = 1e-2;
    int m_dThresholdRatio = 60;
    vector<vector<Point>> m_vecBezier;
    vector<vector<Point>> m_vecPrime;
    vector<vector<Point>> m_vecSrcPoint;
    vector<vector<Point>> m_vecDesPoint;
    vector<vector<Point>> m_vecPiecePoint;
    vector<vector<Point>> m_vecAdaptPoint;
    vector<vector<double>> m_vecPieceTime;
    vector<vector<double>> m_vecAdaptTime;
    vector<vector<double>> m_vecChangeTime;
};



#endif // ARCHHEIGHT_H

