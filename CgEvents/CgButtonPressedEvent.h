#ifndef CGBUTTONPRESSEDEVENT_H
#define CGBUTTONPRESSEDEVENT_H

#include "../CgBase/CgBaseEvent.h"
#include "../CgUtils/CgEventEnums.h"
#include <iostream>
#include <string>

/**
 * @brief The CgButtonPressedEvent class
 */
class CgButtonPressedEvent : public CgBaseEvent {
public:
    CgButtonPressedEvent();
    CgButtonPressedEvent(Cg::EventType type, int eventcode);

    //inherited
    Cg::EventType getType();
    CgBaseEvent* clone();

    int getCode() const { return m_eventcode; }

    friend std::ostream& operator<<(std::ostream& os, const CgButtonPressedEvent& event);

private:
    Cg::EventType m_type;
    int m_eventcode;
};

#endif // CGBUTTONPRESSEDEVENT_H
