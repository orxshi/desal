#ifndef CELL_H
#define CELL_H

#include "face.h"

namespace Liver
{
    enum Shape
    {
        tet,
        hex,
        pri,
        quad,
        tri,
        undef,
    };

    enum PropertyIndex
    {
        u = 0,
        v = 1,
        w = 2,
        p = 3,
        pp = 4,
        phi = 5,
    };

    enum PatchName
    {
        interior = 1,
        wall = 2,
        movingwall = 3,
        symmetry = 4,
        inlet = 5,
        outlet = 6,
    };

    class Cell
    {
        public:

            Cell(CellTag, const PointPointers&, Shape, PatchName);

            Velocity velocity() const;
            void volume();
            const Centroid& centroid();
            void add_neighbor(Neighbor);
            void add_face(FacePointer);
            void set_parent_cell_of_vertices();
            CellTag tag() const;
            Shape shape() const;
            PointPointers& points();
            FacePointers& faces();
            Neighbors& neighbors();


            CellTag tag_;
            Shape shape_;
            Neighbors neighbors_;
            PointPointers points_;
            FacePointers faces_;
            Volume volume_;
            Centroid centroid_;
            Properties var_;
            Properties oldvar_;
            PatchName patch_;
            SquareMatrixD<N_DIM> LU_;
            Coefs rhs_;
            Coefs aC_;
            //Gradient gradp_;
            Gradients grad_;
            Scalar summ_;

            void LS_coef();
    };

    CellPointer opposing_neighbor(FacePointer, CellPointer);
}


#endif
