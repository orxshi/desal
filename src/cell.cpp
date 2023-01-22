#include <algorithm>
#include "cell.h"

namespace Liver
{
    Cell::Cell(CellTag celltag, const PointPointers& points, Shape shape, PatchName patch): tag_(celltag), points_(points), shape_(shape), volume_(0.), centroid_(origin), patch_(patch)
    {
        var_[p] = 0.;
        var_[pp] = 0.;
        summ_ = 0.;
    }

    Shape Cell::shape() const
    {
        return shape_;
    }

    PointPointers& Cell::points()
    {
        return points_;
    }

    FacePointers& Cell::faces()
    {
        return faces_;
    }

    Neighbors& Cell::neighbors()
    {
        return neighbors_;
    }

    CellTag Cell::tag() const
    {
        return tag_;
    }

    void Cell::add_neighbor(Neighbor neighbor)
    {
        neighbors_.push_back(neighbor);
    }

    void Cell::add_face(FacePointer face)
    {
        faces_.push_back(face);
    }

    void Cell::set_parent_cell_of_vertices()
    {
        for (PointPointer p: points_)
        {
            auto it = std::find_if(p->parent_cells().begin(), p->parent_cells().end(), [&](ParentCell pc){return pc->tag() == tag_;});

            if (it == p->parent_cells().end())
            {
                p->add_parent_cell(this);
            }
        }
    }

    void Cell::volume()
    {
        if (shape_ == Shape::tri || shape_ == Shape::quad)
        {
            return;
        }

        volume_ = 0.;

        Centroid vc_cell = vertex_centroid(points_);

        for (FacePointer face: faces_)
        {
            const PointPointers& pts = face->points();

            Centroid vc_pyramid = vertex_centroid(pts);

            Centroid ac_pyramid = 0.75 * face->centroid_ + 0.25 * vc_pyramid;

            Vector3 dGF = face->centroid_ - vc_cell;

            Volume vol_pyramid = (face->area_ * dot(dGF, face->normal_)) / 3.;

            volume_ += std::abs(vol_pyramid);

        }

        assert(!std::isnan(volume_));
    }

    const Centroid& Cell::centroid()
    {
        if (shape_ == Shape::tri || shape_ == Shape::quad)
        {
            centroid_ = faces_[0]->centroid_;
            return centroid_;
        }

        centroid_ = origin;

        Centroid vc_cell = vertex_centroid(points_);

        for (FacePointer face: faces_)
        {
            const PointPointers& pts = face->points();

            Centroid vc_pyramid = vertex_centroid(pts);

            Centroid ac_pyramid = 0.75 * face->centroid_ + 0.25 * vc_pyramid;

            Vector3 dGF = face->centroid_ - vc_cell;

            Volume vol_pyramid = (face->area_ * dot(dGF, face->normal_)) / 3.;
            vol_pyramid = std::abs(vol_pyramid);

            centroid_ += vol_pyramid * ac_pyramid; 
        }

        centroid_ /= volume_;
        assert(!std::isnan(centroid_(0)));

        return centroid_;
    }

    CellPointer opposing_neighbor(FacePointer face, CellPointer cell)
    {
        assert(face->parent_cells().size() == 2);

        for (ParentCell parent: face->parent_cells())
        {
            if (parent->tag() != cell->tag())
            {
                return parent;
            }
        }
    }

    Velocity Cell::velocity() const
    {
        return Vector3(var_[u], var_[v], var_[w]);
    }
}
