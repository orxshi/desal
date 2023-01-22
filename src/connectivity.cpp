#include "mesh.h"
#include "meshpoint.h"
#include "meshcell.h"

namespace Tailor
{
    void Mesh::add_interior_face(MeshFace& mf, MeshFace& cf, MeshCell& mc, MeshCell& cc)
    {
        assert(mf.parent_cells().size() == 1);
        assert(cf.parent_cells().size() == 1);

        Tag mctag = mc.tag();
        Tag cctag = cc.tag();

        mf.set_btype(BouType::interior);
        cf.set_btype(BouType::interior);

        auto ft = gen_face_tag(mf);

        mf.set_tag(ft);
        cf.set_tag(ft);

        auto ittt = std::find_if(mf.parent_cell().begin(), mf.parent_cell().end(), [&](const Tag& pcc){return pcc == cctag;});
        assert(ittt == mf.parent_cell().end());
        if (ittt == mf.parent_cell().end())
        {
            mf.add_parent_cell(cctag);
        }

        ittt = std::find_if(cf.parent_cell().begin(), cf.parent_cell().end(), [&](const Tag& pcc){return pcc == mctag;});
        if (ittt == cf.parent_cell().end())
        {
            cf.add_parent_cell(mctag);
            if (cf.parent_cell().size() == 2)
            {
                assert(cf.parent_cell(0) != cf.parent_cell(1));
            }
        }

        if (std::count(mc.pnei().begin(), mc.pnei().end(), cctag) == 0)
        {
            mc.add_pnei(cctag);
        }
        if (std::count(cc.pnei().begin(), cc.pnei().end(), mctag) == 0)
        {
            cc.add_pnei(mctag);
        }

        // for connect_cells this step is redundant.
        // for connect_after_exchange this step is necessary.
        for (const Tag& p: mf.mesh_point())
        {
            auto pp = point_p(p);
            assert(pp != nullptr);
            pp->add_parent_cell(cctag);
        }

        for (const Tag& p: cf.mesh_point())
        {
            auto pp = point_p(p);
            assert(pp != nullptr);
            pp->add_parent_cell(mctag);
        }

        assert(std::find_if(mc.face().begin(), mc.face().end(), [&](const MeshFace& ff){return ff.tag() == mf.tag();}) != mc.face().end());
        assert(std::find_if(cc.face().begin(), cc.face().end(), [&](const MeshFace& ff){return ff.tag() == mf.tag();}) != cc.face().end());

        assert(mf.parent_cell().size() <= 2);
        assert(cf.parent_cell().size() <= 2);
    }

    void Mesh::addface_bou(MeshFace& mf)
    {
        assert(mf.btype() != BouType::undefined);
        assert(mf.is_boundary());

        auto ft = gen_face_tag(mf);

        mf.set_tag(ft);
    }

    void Mesh::connect_cells(std::function<bool(const Vector3&)> is_resi, int rank, Profiler* profiler, std::string name)
    {
        for (MeshCell& mc: cell_)
        {
            for (MeshFace& mf: mc.face_p())
            {
                if (mf.is_boundary())
                {
                    continue;
                }

                if (!mf.parent_cell().empty())
                {
                    std::cout << "mf.btype(): " << static_cast<int>(mf.btype()) << std::endl;
                }
                assert(mf.parent_cell().empty());
                mf.add_parent_cell(mc.tag());
                assert(mf.parent_cell().size() == 1);
            }
        }

        for (MeshCell& mc: cell_)
        {
            assert(mc.tag().isvalid());

            for (MeshFace& mf: mc.face_p())
            {
                if (mf.is_boundary())
                {
                    assert(mf.parent_cell().size() == 2);
                    addface_bou(mf);
                    continue;
                }

                if (mf.parent_cell().size() == 2)
                {
                    assert(mf.btype() != BouType::undefined);
                    continue;
                }

                assert(!mf.mesh_point().empty());
                //MeshPoint& p0 = point_p(mf.mesh_point(0));
                auto pp = point_p(mf.mesh_point(0));
                assert(pp != nullptr);
                MeshPoint& p0 = *pp;
                assert(!p0.parent_cell().empty());

                for (auto& nei: p0.parent_cell())
                {
                    //assert(query(nei) != nullptr);
                    MeshCell& nc = cell_p(nei);
                    assert(nei.isvalid());
                    if (nc.tag() == mc.tag())
                    {
                        continue;
                    }
                    
                    if (is_potential_parent_cell(mf, nc.tag()))
                    {
                        MeshFace* cf = common_face(nc, mc.tag());
                        assert(cf != nullptr);

                        add_interior_face(mf, *cf, mc, nc);

                        break;
                    }
                }
            }
        }

        for (MeshCell& mc: cell_)
        {
            int nnei = 0;
            for (const MeshFace& mf: mc.face())
            {
                if (!mf.is_boundary())
                {
                    ++nnei;
                }
            }
            if (mc.pnei().size() != nnei)
            {
                mc.set_btype(BouType::partition);
            }
            for (MeshFace& mf: mc.face_p())
            {
                if (mf.parent_cell().size() == 1)
                {
                    assert(!mf.is_boundary());
                    add_partition_face(mf);
                    assert(mf.parent_cell().size() <= 1);
                }
            }
        }

        for (const MeshCell& mc: cell_)
        {
            for (const MeshFace& mf: mc.face())
            {
                assert(!mf.parent_cell().empty());

                if (mf.is_boundary())
                {
                    if (mf.parent_cell().size() != 2)
                    {
                        std::cout << mf.parent_cell().size() << std::endl;
                    }
                    assert(mf.parent_cell().size() == 2);
                }
                else
                {
                    if (mf.parent_cell().size() == 2)
                    {
                        assert(mf.btype() == BouType::interior);
                    }
                }
            }
        }

        for (MeshCell& mc: cell_)
        {
            if (!is_resi(mc.poly().centroid()))
            {
                mc.set_oga_cell_type(OGA_cell_type_t::non_resident);
            }
        }
        //if (profiler != nullptr) {profiler->stop(name + "-connect-cells-pre");}

        for (const auto& mc: cell_)
        {
            for (const auto& mf: mc.face_)
            {
                if (mf.is_boundary()) {
                    continue;
                }
                if (mf.parent_cell(0) != mc.tag())
                {
                    std::cout << "mc: " << mc.tag()() << std::endl;
                    std::cout << "pc 0: " << mf.parent_cell(0)() << std::endl;
                    std::cout << "pc 1: " << mf.parent_cell(1)() << std::endl;
                    std::cout << "mf btype: " << static_cast<int>(mf.btype()) << std::endl;
                }
                assert(mf.parent_cell(0) == mc.tag());
                assert(mf.right_cell() == mc.tag());
            }
        }
    }
}
