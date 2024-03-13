#ifndef CG_SPLAT_EVENT
#define CG_SPLAT_EVENT



#include "../CgBase/CgBaseEvent.h"
#include "../CgBase/CgEnums.h"


#include <iostream>
#include <string>

class CgSplatEvent : public CgBaseEvent
{
public:
    CgSplatEvent();
    CgSplatEvent(Cg::EventType type, bool use_spatting);
    ~CgSplatEvent(){};

    //inherited
    Cg::EventType  getType();
    CgBaseEvent* clone();

    bool useSplatting(){return m_use_spatting;}


    friend std::ostream& operator <<(std::ostream& os, const CgSplatEvent& e);

private:
     Cg::EventType m_type;
     bool m_use_spatting;



};








#endif // CGLOADOBJFILEEVENT

