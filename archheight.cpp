#include "archheight.h"
#include <iostream>
#include <random>
#include <vector>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <cmath>

ArchHeight::ArchHeight(const double height, const size_t index):m_iIndex(index),m_dHeight(height)
{
    doInit();
    setSrcPoint();
    getDesPoint();
    getBezier();
    getBezierPrime();
    doIsoHeightSeg();
    doAdaptiveSampling();
    doComputePiecePoint();
    doComputeAdaptPoint();
}

ArchHeight::ArchHeight()
{
    doInit();
    setSrcPoint();
    getDesPoint();
    getBezier();
    getBezierPrime();
    doIsoHeightSeg();
    doAdaptiveSampling();
    doComputePiecePoint();
    doComputeAdaptPoint();
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
    m_vecAdaptTime.resize(m_iIndex);
    m_vecAdaptPoint.resize(m_iIndex);
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
//    m_vecSrcPoint[0].emplace_back(Point(199.55,560.725));
//    m_vecSrcPoint[0].emplace_back(Point(104.85,374.45));
//    m_vecSrcPoint[0].emplace_back(Point(768.125,597.525));
//    m_vecSrcPoint[0].emplace_back(Point(323.15,172.375));

//    m_vecSrcPoint[1].emplace_back(Point(580.85,205.575));
//    m_vecSrcPoint[1].emplace_back(Point(371.525,812.425));
//    m_vecSrcPoint[1].emplace_back(Point(89.725,652.45));
//    m_vecSrcPoint[1].emplace_back(Point(213.25,52.55));
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
/***
* @brief        : 计算点到直线的距离，即拱高
* @in_param : index:第几条曲线   begintime:起始时间点     endtime:终止时间    heighttime:拱高点所在的时间
* @out_param : dDist:点到直线的距离
* @return : double类型的dDist
*/
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
        dDist = fabs(dA * heightPoint.x - dB * heightPoint.y + dC) / sqrt(dA * dA + dB * dB);
    }
    return dDist;
}
/***
* @brief        : 计算曲线变化率发生改变的时间点
* @in_param : none
* @out_param : none
* @return : none
*/
void ArchHeight::doComputeChangeTime()
{
    for (size_t i = 0; i < m_iIndex; ++i) {
        m_vecChangeTime[i].emplace_back(0);
        //  计算 X 参数方程单调性变化的两个时间点    
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
        //  计算 Y 参数方程单调性变化的两个时间点
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
/***
* @brief        : 根据拱高分段时间计算相应的拱高分段点
* @in_param : none
* @out_param : none
* @return : none
*/
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
/***
* @brief        : 根据自适应采样的时间计算数据点
* @in_param : none
* @out_param : none
* @return : none
*/
void ArchHeight::doComputeAdaptPoint()
{
    for (size_t i = 0; i < m_vecAdaptTime.size(); ++i)
    {
        if (m_vecAdaptTime[i].empty())
        {
            continue;
        }
        std::cout << m_vecAdaptTime[i].size() << std::endl;
        for(size_t j = 0; j < m_vecAdaptTime[i].size(); ++j)
        {
            m_vecAdaptPoint[i].emplace_back(doComputePoint(i, m_vecAdaptTime[i][j]));
        }
    }
}
/***
* @brief        : 计算两点之间的拱高
* @in_param : index:曲线指标    begintime:起始时间      endtime:终止时间
* @out_param : dArchHeight:拱高
* @return : double类型的拱高dArchHeight
*/
double ArchHeight::doComputeArchHeight(const size_t& index, const double& begintime, const double& endtime)
{

    double dArchHeight = 0;
    double dHeightTime = 0;

    double dMid = (begintime+endtime)/2;

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
    double dFirstDegree = dE * dKx - dB * dKy;
    double dConstTrem = dF * dKx - dC * dKy;
    double delta = dFirstDegree * dFirstDegree - 4 * dSecDegree * dConstTrem;

    if (delta >= 0)
    {
        double dX1 = (-dFirstDegree + sqrt(delta)) / (2 * dSecDegree);
        double dDistX1 = 0;
        double dX2 = (-dFirstDegree - sqrt(delta)) / (2 * dSecDegree);
        double dDistX2 = 0;

        //  NOTE:如果两个都有效则取较大的拱高
        if(dX1 > begintime && dX1 < endtime)
        {
            dDistX1 = doComputePointLineDist(index, begintime, endtime, dX1);
        }
        if(dX2 > begintime && dX2 < endtime)
        {
            dDistX2 = doComputePointLineDist(index, begintime, endtime, dX2);
        }
        dArchHeight = std::max(dDistX1,dDistX2);
    }
    return dArchHeight;
}

double ArchHeight::doComputeArchHeightTime(const size_t &index, const double &begintime, const double &endtime)
{
    double dArchHeightTime = -1;

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
    double dFirstDegree = dE * dKx - dB * dKy;
    double dConstTrem = dF * dKx - dC * dKy;
    double delta = dFirstDegree * dFirstDegree - 4 * dSecDegree * dConstTrem;

    if (delta >= 0)
    {
        double dX1 = (-dFirstDegree + sqrt(delta)) / (2 * dSecDegree);
        double dDistX1 = 0;
        bool bIsValidX1 = false;
        double dX2 = (-dFirstDegree - sqrt(delta)) / (2 * dSecDegree);
        double dDistX2 = 0;
        bool bIsValidX2 = false;
        //  NOTE:如果两个都有效则取较大的拱高
        if(dX1 > begintime && dX1 < endtime)
        {
            dDistX1 = doComputePointLineDist(index, begintime, endtime, dX1);
            bIsValidX1 = true;
        }
        if(dX2 > begintime && dX2 < endtime)
        {
            dDistX2 = doComputePointLineDist(index, begintime, endtime, dX2);
            bIsValidX2 = true;
        }
        if(bIsValidX1 && bIsValidX2)
        {
            dArchHeightTime =  std::fabs(std::min(dDistX1,dDistX2) - dDistX1)<1e-4?dX1:dX2;
        }
        else if(bIsValidX1)
        {
            dArchHeightTime = dX1;
        }
        else
        {
            dArchHeightTime = dX2;
        }
    }
    return dArchHeightTime;
}

double ArchHeight::Distance(const Point& point1, const Point& point2)
{
    return sqrt(pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2));
}
/***
* @brief        : 从零开始,逐步计算拱高,当拱高合适时,计算下一个点.
* @in_param : none
* @out_param : none
* @return : none
*/
void ArchHeight::doIsoHeightSeg()
{
    for (size_t i = 0; i < m_vecSrcPoint.size(); ++i)
    {
        m_vecPieceTime[i].emplace_back(0);
        double dBegin = 0.0;                    //  NOTE:每条曲线的拱高分段都从 0 开始
        double dEnd = 0.0;
        double dMaxRenge = 1.0;                 //  NOTE:每条曲线的拱高分段都以 1 结束
        double dSteps = 1e-2;                   //  NOTE:步长设置  步长也与拱高有关
        double dTol = 1e-1;                     //  NOTE:容忍误差与要求的拱高是相关的
        while (dBegin < dMaxRenge && dEnd < 1.0)
        {
            dEnd += dSteps;
            double dArchHeight = doComputeArchHeight(i, dBegin, dEnd);
            if ((m_dHeight - dArchHeight < dTol) || dArchHeight >= m_dHeight)   //  NOTE:当计算而来的拱高大于或接近设定拱高时 结束
            {
                m_vecPieceTime[i].emplace_back(dEnd);
                dBegin = dEnd;
            }
        }
        m_vecPieceTime[i].emplace_back(1.0);
    }

}

double ArchHeight::dotProduct(const Point &p1, const Point &p2)
{
    return p1.x * p2.x+p1.y * p2.y;
}

double ArchHeight::velocityAngleChange(const Point &point_1, const Point &point_2, const Point &point_3)
{
       Point p1p2 = point_2-point_1;
       Point p2p3 = point_3-point_2;
       double dDotProduct  = p1p2.x*p2p3.x+p1p2.y*p2p3.y;
       return 1-(dDotProduct/(Distance(point_1,point_2)+Distance(point_2,point_3)));

}



double ArchHeight::circumRadius(const Point &point_1, const Point &point_2, const Point &point_3)
{
    double r = 0;

    //1、利用向量叉乘计算三角形面积
    Point P21 = point_1-point_2;
    Point P23 = point_3-point_2;
    double dArea = std::fabs(P21.x*P23.y-P21.y*P23.x);

    //2、利用余弦定理计算中间点对应的叫的正弦值
    double dA = Distance(point_1,point_2);
    double dB = Distance(point_2,point_3);
    double dC = Distance(point_1,point_3);
    double dCosAlpha = (dA*dA+dB*dB-dC*dC)/(2*dA*dB);
    double dSinAlpha = sqrt(1-dCosAlpha*dCosAlpha);

    //3、利用三角形面积与两倍正弦值之比计算三角形外接圆半径
    r = dArea/(2*dSinAlpha);

    return r;
}

double ArchHeight::doComputeWeightedFunction(const size_t &index, const double &begintime, const double &endtime)
{
    //  NOTE:设置算法的一些常数
    double RADIUS = 10;
    double MAXLENGTHRATIO = 10;
    double dOmega1 = 0.2;
    double dOmega2 = 0.3;
    double dOmega3 = 1-dOmega1-dOmega2;

    double result = 0;

    //  NOTE:计算各个权重
    //  1、计算弧长的权重
    double dArchHeightTime = doComputeArchHeightTime(index,begintime,endtime);
    double dArchHeight = doComputePointLineDist(index,begintime,endtime,dArchHeightTime);
    std::cout << "BeginTime : " << begintime << "\t" << "EndTime : " << endtime << "\t" << "HeightTime : " << dArchHeightTime << "\t" << "ArchHeight : " << dArchHeight <<std::endl;
    double dDistance = Distance(doComputePoint(index,begintime),doComputePoint(index,endtime));
    double dLengthRatio = dDistance/(dArchHeight*MAXLENGTHRATIO);

    //  2、计算速度变化
    Point point1 = doComputePoint(index,begintime);
    Point point2 = doComputePoint(index,dArchHeightTime);
    Point point3 = doComputePoint(index,endtime);
    double dVelocityAngleChange = velocityAngleChange(point1,point2,point3);
    std::cout << "VelocityAngleChange : " << dVelocityAngleChange << std::endl;
    //  3、计算曲率半径
    double dCircumRadius = circumRadius(point1,point2,point3)-RADIUS;
    std::cout << "CircumRadius : " << dCircumRadius << std::endl;


    //  4、计算加权函数和
    result = dOmega1 * std::max(1.0,dLengthRatio)+dOmega2*dVelocityAngleChange+dOmega3*dCircumRadius;
    return result;
}
//  我们可以根据这个初步自适应采样的结果 对其再进行操作
//  使得在拐点和变化比较大的地方能够更好的描述曲线
//  自适应随机采样方法？？？
void ArchHeight::doAdaptiveSampling()
{
    double THREDHOLDS = 5e4;

    double dBegin = 0.0;
    double dSteps = 1e-2;
    double dEnd = dBegin+dSteps;
    for(size_t index = 0;index<m_iIndex;++index)
    {
        while(dEnd<1.0 )
        {
            double dWeightSum = doComputeWeightedFunction(index,dBegin,dEnd);
            std::cout << "WeightSum : " << dWeightSum << std::endl;
            if(dWeightSum>THREDHOLDS)
            {
                m_vecAdaptTime[index].emplace_back(dEnd);
                dBegin = dEnd;
            }
            dEnd+=dSteps;
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

vector<vector<Point>> ArchHeight::outAdaptPoint() const
{
    return m_vecAdaptPoint;
}

size_t ArchHeight::outIndex() const
{
    return m_iIndex;
}
