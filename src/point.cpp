#include "point.h"

namespace Liver
{
    Point::Point(PointTag tag, Coordinate coordinate): tag_(tag), coordinate_(coordinate)
    {
    }

    const ParentCells& Point::parent_cells() const
    {
        return parent_cells_;
    }

    const ParentFaces& Point::parent_faces() const
    {
        return parent_faces_;
    }

    PointTag Point::tag() const
    {
        return tag_;
    }

    void Point::add_parent_cell(ParentCell parent_cell)
    {
        parent_cells_.push_back(parent_cell);
    }

    void Point::add_parent_face(ParentFace parent_face)
    {
        parent_faces_.push_back(parent_face);
    }

    const Coordinate& Point::coordinate() const
    {
        return coordinate_;
    }

    Centroid vertex_centroid(const PointPointers& points)
    {
        Centroid vc(0., 0., 0.);

        for (PointPointer p: points)
        {
            vc += p->coordinate();
        }

        vc /= points.size();

        return vc;
    }
}
