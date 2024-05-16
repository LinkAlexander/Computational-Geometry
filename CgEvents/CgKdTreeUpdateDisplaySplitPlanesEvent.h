#ifndef CGKDTREEUPDATEDISPLAYSPLITPLANESEVENT_H
#define CGKDTREEUPDATEDISPLAYSPLITPLANESEVENT_H

#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"

/**
 * @brief The CgKdTreeUpdateDisplaySplitPlanesEvent class
 */
class CgKdTreeUpdateDisplaySplitPlanesEvent : public CgBaseEvent {
public:
    CgKdTreeUpdateDisplaySplitPlanesEvent() { }
    CgKdTreeUpdateDisplaySplitPlanesEvent(bool showSplitPlanes, size_t displayDepth)
        : showSplitPlanes(showSplitPlanes)
        , displayDepth(displayDepth)
    {
    }

    Cg::EventType getType()
    {
        return Cg::CgKdTreeUpdateDisplaySplitPlanesEvent;
    }

    CgBaseEvent* clone()
    {
        return new CgKdTreeUpdateDisplaySplitPlanesEvent(showSplitPlanes, displayDepth);
    }

    bool getShowSplitPlanes() const
    {
        return showSplitPlanes;
    }

    size_t getDisplayDepth() const
    {
        return displayDepth;
    }

private:
    bool showSplitPlanes = false;
    size_t displayDepth = 1;
};
#endif // CGKDTREEUPDATEDISPLAYSPLITPLANESEVENT_H
