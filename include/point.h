#ifndef POINT_H
#define	POINT_H

#include "definitions.h"

namespace Liver
{
    class Point
    {
        public:

            Point(PointTag, Coordinate);

            void add_parent_face(ParentFace);
            void add_parent_cell(ParentCell);
            const ParentCells& parent_cells() const;
            const ParentFaces& parent_faces() const;
            PointTag tag() const;
            const Coordinate& coordinate() const;

        private:

            PointTag tag_;
            Coordinate coordinate_;
            ParentFaces parent_faces_;
            ParentCells parent_cells_;
    };

    Centroid vertex_centroid(const PointPointers&);
}

#endif
