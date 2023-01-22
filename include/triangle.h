#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "point.h"

namespace Liver
{
    Centroid triangle_centroid(const PointPointers&);
    Area triangle_area(const PointPointers&);
}

#endif
