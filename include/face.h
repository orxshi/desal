#ifndef FACE_H
#define	FACE_H

#include "definitions.h"
#include "point.h"
#include "triangle.h"

namespace Liver
{
    class Face
    {    
        public:

            Face(FaceTag, const PointPointers&);

            void centroid();
            void normal();
            void area();
            void set_tag(FaceTag);
            void add_parent_cell(ParentCell);
            FaceTag tag() const;
            const PointPointers& points() const;
            const ParentCells& parent_cells() const;


            FaceTag tag_;
            PointPointers points_;
            ParentCells parent_cells_;
            Centroid centroid_;
            Normal normal_;
            Area area_;
            Coefs aFl_;
            Coefs aFr_;
            Left l;
            Right r;
            double g;
            Scalar diff_flux_;
            Vector3 CF_;
            Scalar lenCF_;
    };
}

#endif
