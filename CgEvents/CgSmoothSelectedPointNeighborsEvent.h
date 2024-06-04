#ifndef CGSMOOTHSELECTEDPOINTNEIGHBORSEVENT_H
#define CGSMOOTHSELECTEDPOINTNEIGHBORSEVENT_H

#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"

/**
 * @brief The CgSmoothSelectedPointNeighborsEvent class
 */
class CgSmoothSelectedPointNeighborsEvent : public CgBaseEvent {
public:
    CgSmoothSelectedPointNeighborsEvent() { }
    CgSmoothSelectedPointNeighborsEvent(size_t neighborCount, size_t bivariateFunctionDegree)
        : neighborCount(neighborCount)
        , bivariateFunctionDegree(bivariateFunctionDegree)
    {
    }

    Cg::EventType getType()
    {
        return Cg::CgSmoothSelectedPointNeighborsEvent;
    }

    CgBaseEvent* clone()
    {
        return new CgSmoothSelectedPointNeighborsEvent(neighborCount, bivariateFunctionDegree);
    }

    size_t getNeighborCount()
    {
        return neighborCount;
    }

    size_t getBivariateFunctionDegree()
    {
        return bivariateFunctionDegree;
    }

private:
    size_t neighborCount;
    size_t bivariateFunctionDegree;
};

#endif // CGSMOOTHSELECTEDPOINTNEIGHBORSEVENT_H
