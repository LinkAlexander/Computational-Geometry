#ifndef CGLOADPOINTCLOUDEVENT
#define CGLOADPOINTCLOUDEVENT



#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"


#include <iostream>
#include <string>

class CgLoadPointCloudEvent : public CgBaseEvent
{
public:
    CgLoadPointCloudEvent();
    CgLoadPointCloudEvent(Cg::EventType type, std::string filename);
    ~CgLoadPointCloudEvent(){};

    //inherited
    Cg::EventType  getType();
    CgBaseEvent* clone();

    std::string FileName(){return m_filename;}


    friend std::ostream& operator <<(std::ostream& os, const CgLoadPointCloudEvent& e);

private:
     Cg::EventType m_type;
     std::string m_filename;



};








#endif // CGLOADPOINTCLOUDEVENT

