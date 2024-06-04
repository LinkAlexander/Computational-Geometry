#ifndef CGSMOOTHSURFACEEVENT_H
#define CGSMOOTHSURFACEEVENT_H

#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"

/**
 * @brief The CgSmoothSurfaceEvent class
 */
class CgSmoothSurfaceEvent : public CgBaseEvent {

public:
    CgSmoothSurfaceEvent()
        : neighborCount(10)
        , bivariateFunctionDegree(2)
    {
    }
    CgSmoothSurfaceEvent(size_t neighborCount, size_t bivariateFunctionDegree)
        : neighborCount(neighborCount)
        , bivariateFunctionDegree(bivariateFunctionDegree)
    {
    }

    Cg::EventType getType()
    {
        return Cg::CgSmoothSurfaceEvent;
    }

    CgBaseEvent* clone()
    {
        return new CgSmoothSurfaceEvent(neighborCount, bivariateFunctionDegree);
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

#endif // CGSMOOTHSURFACEEVENT_H
