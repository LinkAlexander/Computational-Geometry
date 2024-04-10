#include "CgButtonPressedEvent.h"

CgButtonPressedEvent::CgButtonPressedEvent()
{

}
CgButtonPressedEvent::CgButtonPressedEvent(Cg::EventType type, int eventcode)
{
    m_type = type;
    this->m_eventcode = eventcode;
}

Cg::EventType CgButtonPressedEvent::getType()
{
    return m_type;
}

CgBaseEvent* CgButtonPressedEvent::clone()
{
    return new CgButtonPressedEvent(m_type, m_eventcode);
}

std::ostream& operator<<(std::ostream& os,const CgButtonPressedEvent& e)
{
    os << "Value Changed Event of Type: "<< e.m_type;
    return os;
}
