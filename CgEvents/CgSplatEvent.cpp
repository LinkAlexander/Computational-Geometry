#include "CgSplatEvent.h"





CgBaseEvent* CgSplatEvent::clone()
{
    return new CgSplatEvent(m_type,m_use_spatting);
}

Cg::EventType  CgSplatEvent::getType()
{
    return m_type;
}


std::ostream& operator<<(std::ostream& os,const CgSplatEvent& e)
{
    os << "CgSplatEvent: use splatting "<< e.m_use_spatting << std::endl;
    return os;
}




CgSplatEvent::CgSplatEvent(Cg::EventType type, bool use_splatting)
{
    m_type=type;
    m_use_spatting=use_splatting;
}


CgSplatEvent::CgSplatEvent()
{
    m_type=Cg::CgSplatEvent;
}


