#include "archheight.h"
#include <iostream>
#include <random>
#include <vector>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <cmath>

ArchHeight::ArchHeight(double height, size_t index):m_iIndex(index),m_dHeight(height)
{
    doInit();
    setSrcPoint();
    getDesPoint();
    getBezier();
    getBezierPrime();
    doIsoHeightSeg();
}

ArchHeight::~ArchHeight()
{
}

void ArchHeight::doInit()
{
    m_vecBezier.resize(m_iIndex);
    m_vecPrime.resize(m_iIndex);
    m_vecSrcPoint.resize(m_iIndex);
    m_vecDesPoint.resize(m_iIndex);
    m_vecPiecePoint.resize(m_iIndex);
    m_vecPieceTime.resize(m_iIndex);
    m_vecChangeTime.resize(m_iIndex);
}

void ArchHeight::setSrcPoint()
{
    srand(time(0));
    for (size_t i = 0; i < m_iIndex; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            m_vecSrcPoint[i].emplace_back(rand() / 4e1, rand() / 4e1);
        }
    }
}

void ArchHeight::getBezier()
{
    for (size_t i = 0; i < m_iIndex; ++i)
    {
        m_vecBezier[i].emplace_back(-m_vecSrcPoint[i][0] + 3 * m_vecSrcPoint[i][1] - 3 * m_vecSrcPoint[i][2] + m_vecSrcPoint[i][3]);
        m_vecBezier[i].emplace_back(3 * m_vecSrcPoint[i][0] - 6 * m_vecSrcPoint[i][1] + 3 * m_vecSrcPoint[i][2]);
        m_vecBezier[i].emplace_back(-3 * m_vecSrcPoint[i][0] + 3 * m_vecSrcPoint[i][1]);
        m_vecBezier[i].emplace_back(1 * m_vecSrcPoint[i][0]);
    }
}

void ArchHeight::getBezierPrime()
{
    for (size_t i = 0;i<m_iIndex;++i)
    {
        m_vecPrime[i].emplace_back(3 * m_vecBezier[i][0]);
        m_vecPrime[i].emplace_back(2 * m_vecBezier[i][1]);
        m_vecPrime[i].emplace_back(1 * m_vecBezier[i][2]);
    }
}

void ArchHeight::getDesPoint()
{
    size_t size = m_vecSrcPoint[0].size();
    for (double t = 0; t < 1.0000; t += m_dPrecis)
    {
        vector<Point> resultPoint(size, Point{ 0,0 });
        vector<double> bernstein1(size, 1), bernstein2(size, 1);
        for (size_t i = 1; i <= size - 1; ++i)
        {
            bernstein1[i] *= bernstein1[i - 1] * t;
            bernstein2[size - i - 1] *= bernstein2[size - i] * (1 - t);
        }
        for (size_t i = 0; i < m_iIndex; ++i)
        {
            for (size_t j = 0; j < size; ++j)
            {
                resultPoint[i] += m_vecSrcPoint[i][j] * m_ivecPascalTri[j + 1] * bernstein1[j] * bernstein2[j];
            }
        }
        for (size_t i = 0; i < m_iIndex; ++i)
        {
            m_vecDesPoint[i].emplace_back(resultPoint[i]);
        }
    }
}

double ArchHeight::doComputePointLineDist(const size_t& index, const double& begintime, const double& endtime, const double& heighttime)
{
    if (begintime >= heighttime || heighttime >= endtime)
    {
        printf("Got an error data ,we can't compute the ArchHeight \n");
        return 0;
    }
    double dDist = 0;
    Point beginPoint = doComputePoint(index, begintime);
    Point endPoint = doComputePoint(index, endtime);
    Point heightPoint = doComputePoint(index, heighttime);
    if (abs(beginPoint.x - endPoint.x) < 1e-6) {
        dDist = fabs(heightPoint.x - beginPoint.x);
    }
    else
    {
        double dA = endPoint.y - beginPoint.y;
        double dB = endPoint.x - beginPoint.x;
        double dC = endPoint.x * beginPoint.y - beginPoint.x * endPoint.y;
        dDist = fabs(dA * heightPoint.x + dB * heightPoint.y + dC) / sqrt(dA * dA + dB * dB);
    }
    return dDist;
}

void ArchHeight::doComputeChangeTime()
{
    for (size_t i = 0; i < m_iIndex; ++i) {
        m_vecChangeTime[i].emplace_back(0);
        if (true) {
            if (m_vecPrime[0].at(0).x != 0) {
                double delta = m_vecPrime[i].at(1).x * m_vecPrime[i].at(1).x - 4 * m_vecPrime[i].at(0).x * m_vecPrime[i].at(2).x;
                if (delta >= 0) {
                    double time1 = (-m_vecPrime[i].at(1).x + sqrt(delta)) / (2 * m_vecPrime[i].at(0).x);
                    double time2 = (-m_vecPrime[i].at(1).x - sqrt(delta)) / (2 * m_vecPrime[i].at(0).x);
                    if (time1 > 0 && time1 < 1) {
                        m_vecChangeTime[i].emplace_back(time1);
                    }
                    if (time2 > 0 && time2 < 1) {
                        m_vecChangeTime[i].emplace_back(time2);
                    }
                }
            }
            else if (m_vecPrime[i].at(1).x != 0) {
                double time1 = -m_vecPrime[i].at(2).x / m_vecPrime[i].at(1).x;
                if (time1 > 0 && time1 < 1) {
                    m_vecChangeTime[i].emplace_back(time1);
                }
            }
        }
        //  计算 y 参数方程单调性变化的两个时间点
        if (true) {
            if (m_vecPrime[i].at(0).y != 0) {
                double delta = m_vecPrime[i].at(1).y * m_vecPrime[i].at(1).y - 4 * m_vecPrime[i].at(0).y * m_vecPrime[i].at(2).y;
                if (delta >= 0) {
                    double time1 = (-m_vecPrime[i].at(1).y + sqrt(delta)) / (2 * m_vecPrime[i].at(0).y);
                    double time2 = (-m_vecPrime[i].at(1).y - sqrt(delta)) / (2 * m_vecPrime[i].at(0).y);
                    if (time1 > 0 && time1 < 1) {
                        m_vecChangeTime[i].emplace_back(time1);
                    }
                    if (time2 > 0 && time2 < 1) {
                        m_vecChangeTime[i].emplace_back(time2);
                    }
                }
            }
            else if (m_vecPrime[i].at(1).y != 0) {
                double time1 = -m_vecPrime[i].at(2).y / m_vecPrime[i].at(1).y;
                if (time1 > 0 && time1 < 1) {
                    m_vecChangeTime[i].emplace_back(time1);
                }
            }
        }
        //      对单调性变化的时间进行逆排序   从0到1的顺序排序
        m_vecChangeTime[i].emplace_back(1);
        std::sort(m_vecChangeTime[i].begin(), m_vecChangeTime[i].end());
    }
}

void ArchHeight::doComputePiecePoint()
{
    for (size_t i = 0; i < m_vecPieceTime.size(); ++i)
    {
        if (m_vecPieceTime[i].empty())
        {
            continue;
        }
        for(size_t j = 0; j < m_vecPieceTime[i].size(); ++j)
        {
            m_vecPiecePoint[i].emplace_back(doComputePoint(i, m_vecPieceTime[i][j]));
        }
    }
}

//  endTime要给的合适 使之不能有两个拱高点
double ArchHeight::doComputeArchHeight(const size_t& index, const double& begintime, const double& endtime)
{

    double dArchHeight = 0;
    double dHeightTime = 0;
    Point beginPoint = doComputePoint(index, begintime);
    Point endPoint = doComputePoint(index, endtime);

    double dKx = endPoint.x -beginPoint.x ;
    double dKy = endPoint.y - beginPoint.y;
    double dA = m_vecPrime[index][0].x;
    double dB = m_vecPrime[index][1].x;
    double dC = m_vecPrime[index][2].x;
    double dD = m_vecPrime[index][0].y;
    double dE = m_vecPrime[index][1].y;
    double dF = m_vecPrime[index][2].y;

    double dSecDegree = dD * dKx - dA * dKy;
    double dFirstDegree = dE*dKx-
    if (fabs(beginPoint.x - endPoint.x) < 1e-8)         //  NOTE:两个点形成的直线斜率不存在时
    {
        for (const auto& time : m_vecChangeTime[index])
        {   if (time >= begintime && time <= endtime)
            {
                dHeightTime = time;
                break;
            }
        }
    }
    else                                                //  FIXME:如果有变化率与斜率成倍数的点呢  NOTE:这样的点会不会存在
    {
        double dResDiscri = (dB * dD - dA * dE) * (dB * dD - dA * dE);
        dResDiscri -= (dA * dF - dC * dD) * (dC * dE - dB * dF);
        double dFirstDegree = dB * dD - dA * dE;
        double dConstTerm = (dA * dF - dC * dD);
        double dDel = dA * dE - dB * dD;
        //  NOTE:在确保存在公共解的情况下  这样的判断是否不需要了
        //  NOTE:结式判别在存在公共解的情况下确实是不需要了
        //  NOTE:如果要继续判断是否有两个拱高点  利用韦达定理计算就好了
        //  NOTE:最好避免使用除法  最好使用加法  X1+X2 = -B/A
        //  NOTE:已知解一定是存在的  是不是不用判断直接得解就行了
        //  NOTE:系数矩阵非奇异才能用高斯消元
        if(fabs(dDel)>1e-8)
        {
            dHeightTime = dConstTerm / dFirstDegree;
        }
        if (fabs(dResDiscri) < 1e-8 && fabs(dFirstDegree)>1e-8)
        {
            dHeightTime = dConstTerm / dFirstDegree;
        }
        if (fabs(dDel) < 1e-8)
        {
            //  NOTE:根据两个方程的韦达定理  如果dDel==0  那么方程有两个根
            //  NOTE:但是我们不知道上面的到的解是哪一个而我们要选择一个接近begin time的
            dHeightTime = min(dHeightTime, -(dB / dA) - dHeightTime);
        }
    }
    dArchHeight = doComputePointLineDist(index, begintime, endtime, dHeightTime);
    return dArchHeight;
}
//  NOTE:曲线和导矢都是连续变化的  但是拱高不清楚
void ArchHeight::doIsoHeightSeg()
{
    for (size_t i = 0; i < m_vecSrcPoint.size(); ++i)
    {
        m_vecPieceTime[i].emplace_back(0);
        double dBegin = 0.0;                    //  NOTE:每条曲线的拱高分段都从 0 开始
        double dEnd = 0.0;
        double dMaxRenge = 1.0;                 //  NOTE:每条曲线的拱高分段都以 1 结束
        double dSteps = 1e-2;                   //  NOTE:步长设置  步长也与拱高有关
        double dTol = 1e-2;                     //  NOTE:容忍误差与要求的拱高是相关的
        while (dBegin < dMaxRenge && dEnd < 1.0)
        {
            dEnd += dSteps;
            double dArchHeight = doComputeArchHeight(i, dBegin, dEnd);
            if (fabs(dArchHeight - m_dHeight) < 1e1)
            {
                m_vecPieceTime[i].emplace_back(dEnd);
                dBegin = dEnd;
            }
        }
    }

}







Point ArchHeight::doComputePoint(const size_t& index, const double& t)
{
    Point resultPoint{ 0,0 };
    size_t size = m_vecSrcPoint[index].size();
    vector<double> bernstein1(size, 1), bernstein2(size, 1);
    for (size_t i = 1; i <= size - 1; ++i) {
        bernstein1[i] *= bernstein1[i - 1] * t;
        bernstein2[size - i - 1] *= bernstein2[size - i] * (1 - t);
    }
    for (size_t i = 0; i < size; ++i) {
        resultPoint += m_vecSrcPoint[index][i] * m_ivecPascalTri[i + 1] * bernstein1[i] * bernstein2[i];
    }
    return resultPoint;
}

vector<vector<Point>> ArchHeight::outPiecePoint() const
{
    return m_vecPiecePoint;
}

vector<vector<Point> > ArchHeight::outDesPoint() const
{
    return m_vecDesPoint;
}

vector<vector<Point>> ArchHeight::outSrcPoint() const
{
    return m_vecSrcPoint;
}

size_t ArchHeight::outIndex() const
{
    return m_iIndex;
}
