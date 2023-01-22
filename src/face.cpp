#include "face.h"

namespace Liver
{
    Face::Face(FaceTag tag, const PointPointers& points): tag_(tag), points_(points), area_(0.), centroid_(origin)
    {
    }

    void Face::set_tag(FaceTag tag)
    {
        tag_ = tag;
    }

    const PointPointers& Face::points() const
    {
        return points_;
    }

    const ParentCells& Face::parent_cells() const
    {
        return parent_cells_;
    }

    FaceTag Face::tag() const
    {
        return tag_;
    }

    void Face::add_parent_cell(ParentCell parent_cell)
    {
        parent_cells_.push_back(parent_cell);
    }

    void Face::area()
    {
        if (points_.size() == 3)
        {
            area_ = triangle_area(points_);
            assert(area_ != 0.);
            return;
        }

        area_ = 0.;
        Centroid vc = vertex_centroid(points_);
        PointPointer pvc = new Point(-1, vc);

        for (int i=0; i<points_.size(); ++i)
        {
            PointPointers pts;
            pts.push_back(points_[i]);

            if (i == points_.size() - 1)
            {
                pts.push_back(points_[0]);
            }
            else
            {
                pts.push_back(points_[i+1]);
            }

            pts.push_back(pvc);

            Area tri_area = triangle_area(pts);

            area_ += tri_area;
        }

        assert(area_ != 0.);
        assert(!std::isnan(area_));
    }

    void Face::centroid()
    {
        Centroid vc = vertex_centroid(points_);

        if (points_.size() == 3)
        {
            centroid_ = vc;
            return;
        }

        centroid_ = origin;
        PointPointer pvc = new Point(-1, vc);

        for (int i=0; i<points_.size(); ++i)
        {
            PointPointers pts;
            pts.push_back(points_[i]);

            if (i == points_.size() - 1)
            {
                pts.push_back(points_[0]);
            }
            else
            {
                pts.push_back(points_[i+1]);
            }

            pts.push_back(pvc);

            Area tri_area = triangle_area(pts);
            Centroid tri_centroid = triangle_centroid(pts);

            centroid_ += tri_area * tri_centroid;
        }

        centroid_ /= area_;
        assert(!std::isnan(centroid_(0)));
    }

    void Face::normal()
    {
        normal_ = cross((points_.back()->coordinate() - points_[0]->coordinate()), (points_[1]->coordinate() - points_[0]->coordinate()));

        normal_ = normalize(normal_);
    }
}
