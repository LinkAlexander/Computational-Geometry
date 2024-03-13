#include "CgLoadHalfEdgeMeshEvent.h"





CgBaseEvent* CgLoadHalfEdgeMeshEvent::clone()
{
    return new CgLoadHalfEdgeMeshEvent(m_type,m_filename);
}

Cg::EventType  CgLoadHalfEdgeMeshEvent::getType()
{
    return m_type;
}


std::ostream& operator<<(std::ostream& os,const CgLoadHalfEdgeMeshEvent& e)
{
    os << "LoadObjFileEvent: "<< e.m_filename<< std::endl;
    return os;
}




CgLoadHalfEdgeMeshEvent::CgLoadHalfEdgeMeshEvent(Cg::EventType type, std::string filename)
{
    m_type=type;
    m_filename=filename;
}


CgLoadHalfEdgeMeshEvent::CgLoadHalfEdgeMeshEvent()
{
    m_type=Cg::CgLoadHalfEdgeMeshEvent;
}


