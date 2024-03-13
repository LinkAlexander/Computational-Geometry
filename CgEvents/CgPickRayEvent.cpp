#include "CgPickRayEvent.h"





CgBaseEvent* CgPickRayEvent::clone()
{
    return new CgPickRayEvent(m_type,m_show_pickray);
}

Cg::EventType  CgPickRayEvent::getType()
{
    return m_type;
}


std::ostream& operator<<(std::ostream& os,const CgPickRayEvent& e)
{
    os << "CgSplatEvent: use splatting "<< e.m_show_pickray<< std::endl;
    return os;
}




CgPickRayEvent::CgPickRayEvent(Cg::EventType type, bool show_pickray)
{
    m_type=type;
    m_show_pickray=show_pickray;
}


CgPickRayEvent::CgPickRayEvent()
{
    m_type=Cg::CgPickRayEvent;
}


