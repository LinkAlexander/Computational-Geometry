#include "CgLoadPointCloudEvent.h"





CgBaseEvent* CgLoadPointCloudEvent::clone()
{
    return new CgLoadPointCloudEvent(m_type,m_filename);
}

Cg::EventType  CgLoadPointCloudEvent::getType()
{
    return m_type;
}


std::ostream& operator<<(std::ostream& os,const CgLoadPointCloudEvent& e)
{
    os << "LoadPointCloudEvent: "<< e.m_filename<< std::endl;
    return os;
}




CgLoadPointCloudEvent::CgLoadPointCloudEvent(Cg::EventType type, std::string filename)
{
    m_type=type;
    m_filename=filename;
}


CgLoadPointCloudEvent::CgLoadPointCloudEvent()
{
    m_type=Cg::CgLoadPointCloudEvent;
}


