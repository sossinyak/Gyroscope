#pragma once
#ifndef OSTREAM_H
#define OSTREAM_H

#include <iostream>
#include <vector>
#include "euler_krylov.h"

template <class T>
inline std::ostream& operator << (std::ostream& o, const std::vector<T>& mvect) {
    o << '(';
    if (mvect.size() != 0) {
        int i = 0;
        for (; i < mvect.size() - 1; ++i)
            o << mvect[i] << ", ";
        o << mvect[i];
    }
    o << ')';
    return o;
}

inline std::ostream& operator << (std::ostream& o, const elem element) {
    o << "t = " << element.t << " lat = " << element.lat << " lon = " << element.lon << " h = " << element.height
        << " vxg = " << element.vxg << " vyg = " << element.vyg << " vzg = " << element.vzg
        << " pitch = " << element.pitch << " roll = " << element.roll << " thdg = " << element.thdg;
    return o;
}

#endif