#include "triangle.h"

namespace Liver
{
    Centroid triangle_centroid(const PointPointers& points)
    {
        return vertex_centroid(points);
    }

    Area triangle_area(const PointPointers& points)
    {
        return 0.5 * len(cross((points[1]->coordinate() - points[0]->coordinate()), (points[2]->coordinate() - points[0]->coordinate())));
    }
}
