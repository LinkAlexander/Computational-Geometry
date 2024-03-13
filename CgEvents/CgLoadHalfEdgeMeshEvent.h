#ifndef CGLOAD_HE_OBJFILEEVENT
#define CGLOAD_HE_OBJFILEEVENT



#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"


#include <iostream>
#include <string>

class CgLoadHalfEdgeMeshEvent : public CgBaseEvent
{
public:
    CgLoadHalfEdgeMeshEvent();
    CgLoadHalfEdgeMeshEvent(Cg::EventType type, std::string filename);
    ~CgLoadHalfEdgeMeshEvent(){};

    //inherited
    Cg::EventType  getType();
    CgBaseEvent* clone();

    std::string FileName(){return m_filename;}


    friend std::ostream& operator <<(std::ostream& os, const CgLoadHalfEdgeMeshEvent& e);

private:
     Cg::EventType m_type;
     std::string m_filename;



};








#endif // CGLOADOBJFILEEVENT

