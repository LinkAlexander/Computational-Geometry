#ifndef CG_PICKRAY_EVENT
#define CG_PICKRAY_EVENT



#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"


#include <iostream>
#include <string>

class CgPickRayEvent : public CgBaseEvent
{
public:
    CgPickRayEvent();
    CgPickRayEvent(Cg::EventType type, bool use_spatting);
    ~CgPickRayEvent(){};

    //inherited
    Cg::EventType  getType();
    CgBaseEvent* clone();

    bool showPickRay(){return m_show_pickray;}
    friend std::ostream& operator <<(std::ostream& os, const CgPickRayEvent& e);

private:
     Cg::EventType m_type;
     bool m_show_pickray;



};








#endif // CGLOADOBJFILEEVENT

