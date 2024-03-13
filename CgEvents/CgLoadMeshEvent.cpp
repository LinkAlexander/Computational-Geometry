#include "CgLoadMeshEvent.h"





CgBaseEvent* CgLoadMeshEvent::clone()
{
    return new CgLoadMeshEvent(m_type,m_filename);
}

Cg::EventType  CgLoadMeshEvent::getType()
{
    return m_type;
}


std::ostream& operator<<(std::ostream& os,const CgLoadMeshEvent& e)
{
    os << "LoadObjFileEvent: "<< e.m_filename<< std::endl;
    return os;
}




CgLoadMeshEvent::CgLoadMeshEvent(Cg::EventType type, std::string filename)
{
    m_type=type;
    m_filename=filename;
}


CgLoadMeshEvent::CgLoadMeshEvent()
{
    m_type=Cg::CgLoadMeshEvent;
}


