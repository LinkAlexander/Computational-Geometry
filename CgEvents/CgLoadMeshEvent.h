#ifndef CGLOADOBJFILEEVENT
#define CGLOADOBJFILEEVENT



#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"


#include <iostream>
#include <string>

class CgLoadMeshEvent : public CgBaseEvent
{
public:
    CgLoadMeshEvent();
    CgLoadMeshEvent(Cg::EventType type, std::string filename);
    ~CgLoadMeshEvent(){};

    //inherited
    Cg::EventType  getType();
    CgBaseEvent* clone();

    std::string FileName(){return m_filename;}


    friend std::ostream& operator <<(std::ostream& os, const CgLoadMeshEvent& e);

private:
     Cg::EventType m_type;
     std::string m_filename;



};








#endif // CGLOADOBJFILEEVENT

