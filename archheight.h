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
    friend Point operator*(const Point& p1, const Point& p2) {
        return { p1.x * p2.x,p1.y * p2.y };
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
    void doAdaptiveSampling();      //  NOTE:���ߵ�����Ӧ����(�Ǿ��Ȳ���) 
    void doComputeChangeTime();
    void doComputePiecePoint();     //  NOTE:���ݹ��߷ֶν�������㹰�߷ֶε�
    void doComputeAdaptPoint();
    double doComputePointLineDist(const size_t& index, const double& begintime, const double& endtime, const double& heighttime);//  NOTE:����㵽ֱ�ߵľ���

    double doComputeArchHeight(const size_t& index, const double& begintime, const double& endtime);//  NOTE:��������Ĺ���
    double Distance(const Point&, const Point&);
    Point doComputePoint(const size_t& index, const double& t);    //  NOTE:����ʱ������

public:
    vector<int> m_ivecPascalTri = { 0,1,3,3,1 };

private:
    size_t m_iIndex = 1;
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

