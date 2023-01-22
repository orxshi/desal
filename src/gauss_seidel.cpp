#include "gauss_seidel.h"
#include "blood_properties.h"
#include <fstream>

namespace Liver
{
    void solve_steady(Mesh& mesh, Residual& maxresu, Residual& maxresv, Residual& maxresw)
    {
        int max_nsweep = 1;
        //double lin_tol = 1e-7;

        for (int i=0; i<max_nsweep; ++i)
        {
            for (Cell& cell: mesh.cells())
            {
                cell.aC_[u] = 0.;
                cell.aC_[v] = 0.;
                cell.aC_[w] = 0.;

                cell.rhs_[u] = 0.;
                cell.rhs_[v] = 0.;
                cell.rhs_[w] = 0.;
            }

            spatial_discre(mesh);
            temporal_discre(mesh, 0., steady);
            jacobi_mom(mesh, maxresu, maxresv, maxresw);
            
            //if (maxresu < lin_tol && maxresv < lin_tol && maxresw < lin_tol)
            //{
            //    break;
            //}
        }
    }

    void lefts_rights(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            int i = 0;

            for (ParentCell parent: face.parent_cells())
            {
                if (parent->patch_ == interior)
                {
                    face.l = parent;
                    face.r = face.parent_cells()[1-i];
                    break;
                }

                ++i;
            }

            Vector3 CF = face.r->centroid() - face.l->centroid();
            Vector3 e = normalize(CF);
            Vector3 S = face.area_ * face.normal_;

            if (dot(e, S) < 0)
            {
                CellPointer a = face.l;
                face.l = face.r;
                face.r = a;
            }

            assert(face.l != face.r);
        }
    }

    void interfac(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            Centroid cntf = face.centroid_;

            double fl = len(cntf - face.l->centroid_);
            double fr = len(cntf - face.r->centroid_);

            face.g = fl / (fl + fr);
        }
    }

    void diff_flux(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            double gDiff = face.area_ / face.lenCF_; 
            face.diff_flux_ = MU * gDiff;
        }
    }

    void spatial_discre(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            Left l = face.l;
            Right r = face.r;

            Scalar diff_flux = face.diff_flux_;

            Vector3 n = face.normal_;

            Vector3& CF = face.CF_;
            Scalar area = face.area_;
            Vector3 e = normalize(CF);
            Vector3 S = area * n;

            double eS = dot(e, S);
            Vector3 E = eS * e;
            Vector3 T = S - E;

            double gl = face.g;
            double gr = 1. - gl;

            double ul = l->var_[u];
            double vl = l->var_[v];
            double wl = l->var_[w];
            double pl = l->var_[p];

            double ur = r->var_[u];
            double vr = r->var_[v];
            double wr = r->var_[w];
            double pr = r->var_[p];

            if (r->patch_ == symmetry)
            {
                double a = 2. * MU * area / dot(CF, n);
                double pb = pr * area;

                face.aFl_[u] = 0.;
                face.aFl_[v] = 0.;
                face.aFl_[w] = 0.;

                double anx = a * n(0);
                double any = a * n(1);
                double anz = a * n(2);

                l->aC_[u] += anx * n(0);
                l->aC_[v] += any * n(1);
                l->aC_[w] += anz * n(2);

                double udn = ul * n(0);
                double vdn = vl * n(1);
                double wdn = wl * n(2);

                l->rhs_[u] -= anx * (vdn + wdn) + pb * n(0);
                l->rhs_[v] -= any * (udn + wdn) + pb * n(1);
                l->rhs_[w] -= anz * (udn + vdn) + pb * n(2);

                continue;
            }
            else if (r->patch_ == wall || r->patch_ == movingwall)
            {
                //double a = MU * face.area() / dot(CF, n);
                double pb = pr * area;
                //double bx = 1. - n(0) * n(0);
                //double by = 1. - n(1) * n(1);
                //double bz = 1. - n(2) * n(2);

                face.aFl_[u] = 0.;
                face.aFl_[v] = 0.;
                face.aFl_[w] = 0.;

                //l->aC_[u] += a * bx;
                //l->aC_[v] += a * by;
                //l->aC_[w] += a * bz;

                l->aC_[u] += diff_flux;
                l->aC_[v] += diff_flux;
                l->aC_[w] += diff_flux;

                //l->aC_[u] += a * bx + 2. * a * n(0) * n(0);
                //l->aC_[v] += a * by + 2. * a * n(1) * n(1);
                //l->aC_[w] += a * bz + 2. * a * n(2) * n(2);

                //double udn = (ul - ur) * n(0);
                //double vdn = (vl - vr) * n(1);
                //double wdn = (wl - wr) * n(2);

                //l->rhs_[u] += a * (ur * bx + (vdn + wdn) * n(0)) - pb * n(0);
                //l->rhs_[v] += a * (vr * by + (udn + wdn) * n(1)) - pb * n(1);
                //l->rhs_[w] += a * (wr * bz + (udn + vdn) * n(2)) - pb * n(2);

                //l->rhs_[u] += a * (ur * bx + (vdn + wdn) * n(0)) - pb * n(0) - 2. * a * (vl * n(1) + wl * n(2)) * n(0);
                //l->rhs_[v] += a * (vr * by + (udn + wdn) * n(1)) - pb * n(1) - 2. * a * (ul * n(0) + wl * n(2)) * n(1);
                //l->rhs_[w] += a * (wr * bz + (udn + vdn) * n(2)) - pb * n(2) - 2. * a * (ul * n(0) + vl * n(1)) * n(2);

                l->rhs_[u] += diff_flux * ur - pb * n(0);
                l->rhs_[v] += diff_flux * vr - pb * n(1);
                l->rhs_[w] += diff_flux * wr - pb * n(2);

                continue;
            }
            else if (r->patch_ == inlet || r->patch_ == outlet)
            {
                double uf = ur;
                double vf = vr;
                double wf = wr;

                Vector3 Vf(uf, vf, wf);

                double conv_flux_l = std::min(dot(Vf, S), 0.);
                double pb = (pl + dot(l->grad_[p], CF)) * area;

                face.aFl_[u] = -diff_flux + conv_flux_l;
                face.aFl_[v] = -diff_flux + conv_flux_l;
                face.aFl_[w] = -diff_flux + conv_flux_l;

                double conv_flux_r = std::min(dot(-Vf, S), 0.);

                l->aC_[u] += diff_flux - conv_flux_r;
                l->aC_[v] += diff_flux - conv_flux_r;
                l->aC_[w] += diff_flux - conv_flux_r;

                l->rhs_[u] -= pb * n(0);
                l->rhs_[v] -= pb * n(1);
                l->rhs_[w] -= pb * n(2);

                continue;
            }

            double uf = gl * ur + gr * ul;
            double vf = gl * vr + gr * vl;
            double wf = gl * wr + gr * wl;

            Vector3 Vf(uf, vf, wf); // can be used in ppcoef.

            double Vfs = dot(Vf, S); // similarly.

            //double conv_flux = 0.5 * RHO * dot(Vf, S);
            double conv_flux_l = -RHO * std::min(-Vfs, 0.); // no need to remin.
            double conv_flux_r =  RHO * std::min( Vfs, 0.);

            double fluxl = -diff_flux + conv_flux_r;
            double fluxr = -diff_flux - conv_flux_l;

            double fluxlc = diff_flux + conv_flux_l;
            double fluxrc = diff_flux - conv_flux_r;

            face.aFl_[u] = fluxl;
            face.aFl_[v] = fluxl;
            face.aFl_[w] = fluxl;

            face.aFr_[u] = fluxr;
            face.aFr_[v] = fluxr;
            face.aFr_[w] = fluxr;

            l->aC_[u] += fluxlc;
            l->aC_[v] += fluxlc;
            l->aC_[w] += fluxlc;

            r->aC_[u] += fluxrc;
            r->aC_[v] += fluxrc;
            r->aC_[w] += fluxrc;
            
            double pb = (gl * pr + gr * pl) * area;
            double pbx = pb * n(0); // may use nx.
            double pby = pb * n(1);
            double pbz = pb * n(2);

            l->rhs_[u] -= pbx;
            l->rhs_[v] -= pby;
            l->rhs_[w] -= pbz;

            r->rhs_[u] += pbx;
            r->rhs_[v] += pby;
            r->rhs_[w] += pbz;

            continue;

            //if (r->patch_ == "interior")
            //{
            //    Gradient agu = gl * r->grad_[u] + gr * l->grad_[u];
            //    Gradient agv = gl * r->grad_[v] + gr * l->grad_[v];
            //    Gradient agw = gl * r->grad_[w] + gr * l->grad_[w];

            //    Scalar cu = MU * dot(agu, T);
            //    Scalar cv = MU * dot(agv, T);
            //    Scalar cw = MU * dot(agw, T);

            //    l->rhs_[u] -= cu;
            //    l->rhs_[v] -= cv;
            //    l->rhs_[w] -= cw;

            //    r->rhs_[u] += cu;
            //    r->rhs_[v] += cv;
            //    r->rhs_[w] += cw;
            //}

            double cu, cv, cw;
            Gradient grcu, grcv, grcw;
            double upu, upv, upw;
            double dwu, dwv, dww;
            double baru, barv, barw;
            double HRu, HRv, HRw;
            double dmu, dmv, dmw;
            Vector3 dcd;

            if (Vfs > 0.)
            {
                cu = ul;
                cv = vl;
                cw = wl;
                dwu = ur;
                dwv = vr;
                dww = wr;
                grcu = l->grad_[u];
                grcv = l->grad_[v];
                grcw = l->grad_[w];
                dcd = CF;
            }
            else
            {
                cu = ur;
                cv = vr;
                cw = wr;
                dwu = ul;
                dwv = vl;
                dww = wl;
                grcu = r->grad_[u];
                grcv = r->grad_[v];
                grcw = r->grad_[w];
                dcd = -CF;
            }

            upu = dwu - 2. * dot(grcu, dcd);
            upv = dwv - 2. * dot(grcv, dcd);
            upw = dww - 2. * dot(grcw, dcd);

            dmu = dwu - upu;
            dmv = dwv - upv;
            dmw = dww - upw;

            baru = (cu - upu) / dmu;
            barv = (cv - upv) / dmv;
            barw = (cw - upw) / dmw;

            if (baru >= 0. && baru < 0.5)
            {
                HRu = 1.5 * baru;
            }
            else if (baru >= 0.5 && baru < 1.0)
            {
                HRu = 0.5 * baru + 0.5;
            }
            else
            {
                HRu = baru;
            }

            if (barv >= 0. && barv < 0.5)
            {
                HRv = 1.5 * barv;
            }
            else if (barv >= 0.5 && barv < 1.0)
            {
                HRv = 0.5 * barv + 0.5;
            }
            else
            {
                HRv = barv;
            }

            if (barw >= 0. && barw < 0.5)
            {
                HRw = 1.5 * barw;
            }
            else if (barw >= 0.5 && barw < 1.0)
            {
                HRw = 0.5 * barw + 0.5;
            }
            else
            {
                HRw = barw;
            }

            HRu = Vfs * ((HRu - baru) * dmu);
            HRv = Vfs * ((HRv - barv) * dmv);
            HRw = Vfs * ((HRw - barw) * dmw);

            if (std::isnan(HRu))
            {
                HRu = 0.;
            }
            if (std::isnan(HRv))
            {
                HRv = 0.;
            }
            if (std::isnan(HRw))
            {
                HRw = 0.;
            }

            l->rhs_[u] -= HRu;
            l->rhs_[v] -= HRv;
            l->rhs_[w] -= HRw;

            r->rhs_[u] += HRu;
            r->rhs_[v] += HRv;
            r->rhs_[w] += HRw;
        }
    }

    void temporal_discre(Mesh& mesh, TimeStep deltat, TemporalScheme tscheme)
    {
        if (tscheme == steady)
        {
            for (Face& face: mesh.faces())
            {
                Left l = face.l;
                Right r = face.r;

                if (r->patch_ == symmetry || r->patch_ == wall || r->patch_ == movingwall)
                {
                    continue;
                }

                l->rhs_[u] -= face.aFl_[u] * r->var_[u];
                l->rhs_[v] -= face.aFl_[v] * r->var_[v];
                l->rhs_[w] -= face.aFl_[w] * r->var_[w];

                r->rhs_[u] -= face.aFr_[u] * l->var_[u];
                r->rhs_[v] -= face.aFr_[v] * l->var_[v];
                r->rhs_[w] -= face.aFr_[w] * l->var_[w];
            }
        }
        else if (tscheme == forward_euler)
        {
            for (Cell& cell: mesh.cells())
            {
                double aC_fc = RHO * cell.volume_ / deltat;
                double aC_ec = -aC_fc;

                cell.rhs_[u] -= (cell.aC_[u] + aC_ec) * cell.oldvar_[u];
                cell.rhs_[v] -= (cell.aC_[v] + aC_ec) * cell.oldvar_[v];
                cell.rhs_[w] -= (cell.aC_[w] + aC_ec) * cell.oldvar_[w];

                cell.aC_[u] = aC_fc;
                cell.aC_[v] = aC_fc;
                cell.aC_[w] = aC_fc;
            }

            for (Face& face: mesh.faces())
            {
                Left l = face.l;
                Right r = face.r;

                if (r->patch_ == symmetry || r->patch_ == wall || r->patch_ == movingwall)
                {
                    continue;
                }

                l->rhs_[u] -= face.aFl_[u] * r->oldvar_[u];
                l->rhs_[v] -= face.aFl_[v] * r->oldvar_[v];
                l->rhs_[w] -= face.aFl_[w] * r->oldvar_[w];

                r->rhs_[u] -= face.aFr_[u] * l->oldvar_[u];
                r->rhs_[v] -= face.aFr_[v] * l->oldvar_[v];
                r->rhs_[w] -= face.aFr_[w] * l->oldvar_[w];
            }
        }
        else if (tscheme == backward_euler)
        {
            for (Cell& cell: mesh.cells())
            {
                double aC_fc = RHO * cell.volume_ / deltat;
                double aC_ec = -aC_fc;

                cell.rhs_[u] -= aC_ec * cell.oldvar_[u];
                cell.rhs_[v] -= aC_ec * cell.oldvar_[v];
                cell.rhs_[w] -= aC_ec * cell.oldvar_[w];

                cell.aC_[u] += aC_fc;
                cell.aC_[v] += aC_fc;
                cell.aC_[w] += aC_fc;
            }

            for (Face& face: mesh.faces())
            {
                Left l = face.l;
                Right r = face.r;

                if (r->patch_ == symmetry || r->patch_ == wall || r->patch_ == movingwall)
                {
                    continue;
                }

                l->rhs_[u] -= face.aFl_[u] * r->var_[u];
                l->rhs_[v] -= face.aFl_[v] * r->var_[v];
                l->rhs_[w] -= face.aFl_[w] * r->var_[w];

                r->rhs_[u] -= face.aFr_[u] * l->var_[u];
                r->rhs_[v] -= face.aFr_[v] * l->var_[v];
                r->rhs_[w] -= face.aFr_[w] * l->var_[w];
            }
        }
        else
        {
            assert(false);
        }
    }

    //Residual jacobi(Mesh& mesh, PropertyIndex pindex)
    //{
    //    Scalar ref_res = BIG_NEG_NUM;
    //    Scalar max_res = BIG_NEG_NUM;

    //    for (Cell& cell: mesh.cells())
    //    {
    //        //double aCTu = cell.aC_[u] * cell.var_[u];
    //        //double aCTv = cell.aC_[v] * cell.var_[v];
    //        //double aCTw = cell.aC_[w] * cell.var_[w];

    //        //cell.u_ = cell.rhs_[u] / cell.aC_[u];
    //        //cell.v_ = cell.rhs_[v] / cell.aC_[v];
    //        //cell.w_ = cell.rhs_[w] / cell.aC_[w];

    //        //ref_res = std::max(ref_res, std::abs(aCT));
    //        //Scalar res = std::abs(aCT - cell.rhs_[pindex]) / ref_res;
    //        //max_res = std::max(max_res, res);
    //    }

    //    return max_res;
    //}

    Residual residual(Mesh& mesh, PropertyIndex pindex)
    {
        Scalar max_res = BIG_NEG_NUM;

        for (Cell& cell: mesh.cells())
        {
            max_res = std::max(max_res, std::abs(cell.rhs_[pp]));
        }

        return max_res;
    }

    void jacobi_mom(Mesh& mesh, Residual& maxresu, Residual& maxresv, Residual& maxresw)
    {
        maxresu = BIG_NEG_NUM;
        maxresv = BIG_NEG_NUM;
        maxresw = BIG_NEG_NUM;

        Scalar resu, resv, resw;
        Scalar oldu, oldv, oldw;
        Scalar acu, acv, acw;
        Scalar rhsu, rhsv, rhsw;

        int max_nsweep = 1;
        double lin_tol = 1e-5;

        for (int i=0; i<max_nsweep; ++i)
        {
            Scalar refresu = BIG_NEG_NUM;
            Scalar refresv = BIG_NEG_NUM;
            Scalar refresw = BIG_NEG_NUM;

            for (Cell& cell: mesh.cells())
            {
                oldu = cell.var_[u];
                oldv = cell.var_[v];
                oldw = cell.var_[w];

                acu = cell.aC_[u];
                acv = cell.aC_[v];
                acw = cell.aC_[w];

                rhsu = cell.rhs_[u];
                rhsv = cell.rhs_[v];
                rhsw = cell.rhs_[w];

                resu = std::abs(acu * oldu - rhsu);
                resv = std::abs(acv * oldv - rhsv);
                resw = std::abs(acw * oldw - rhsw);

                cell.var_[u] = rhsu / acu;
                cell.var_[v] = rhsv / acv;
                cell.var_[w] = rhsw / acw;

                assert(!std::isnan(cell.var_[u]));

                refresu = std::max(refresu, std::abs(acu * cell.var_[u]));
                refresv = std::max(refresv, std::abs(acv * cell.var_[v]));
                refresw = std::max(refresw, std::abs(acw * cell.var_[w]));

                if (std::abs(refresu) > ZERO)
                {
                    resu /= refresu;
                }

                if (std::abs(refresv) > ZERO)
                {
                    resv /= refresv;
                }

                if (std::abs(refresw) > ZERO)
                {
                    resw /= refresw;
                }
                
                maxresu = std::max(maxresu, resu);
                maxresv = std::max(maxresv, resv);
                maxresw = std::max(maxresw, resw);
            }

            if (maxresu < lin_tol && maxresv < lin_tol && maxresw < lin_tol)
            {
                break;
            }
        }
    }

    void gauss_seidel_mom(Mesh& mesh, Residual& maxresu, Residual& maxresv, Residual& maxresw)
    {
        maxresu = BIG_NEG_NUM;
        maxresv = BIG_NEG_NUM;
        maxresw = BIG_NEG_NUM;

        Scalar resu, resv, resw;
        Scalar oldu, oldv, oldw;
        Scalar acu, acv, acw;
        Scalar rhsu, rhsv, rhsw;

        int max_nsweep = 1;
        double lin_tol = 1e-5;

        for (int i=0; i<max_nsweep; ++i)
        {
            Scalar refresu = BIG_NEG_NUM;
            Scalar refresv = BIG_NEG_NUM;
            Scalar refresw = BIG_NEG_NUM;

            for (Cell& cell: mesh.cells())
            {
                oldu = cell.var_[u];
                oldv = cell.var_[v];
                oldw = cell.var_[w];

                acu = cell.aC_[u];
                acv = cell.aC_[v];
                acw = cell.aC_[w];

                rhsu = cell.rhs_[u];
                rhsv = cell.rhs_[v];
                rhsw = cell.rhs_[w];

                resu = std::abs(acu * oldu - rhsu);
                resv = std::abs(acv * oldv - rhsv);
                resw = std::abs(acw * oldw - rhsw);

                cell.var_[u] = rhsu / acu;
                cell.var_[v] = rhsv / acv;
                cell.var_[w] = rhsw / acw;

                assert(!std::isnan(cell.var_[u]));

                for (FacePointer face: cell.faces())
                {
                    Neighbor neighbor = opposing_neighbor(face, &cell);

                    Scalar afu, afv, afw;

                    if (neighbor == face->l)
                    {
                        afu = face->aFl_[u];
                        afv = face->aFl_[v];
                        afw = face->aFl_[w];
                    }
                    else
                    {
                        afu = face->aFr_[u];
                        afv = face->aFr_[v];
                        afw = face->aFr_[w];
                    }

                    neighbor->rhs_[u] -= afu * (-oldu + cell.var_[u]);
                    neighbor->rhs_[v] -= afv * (-oldv + cell.var_[v]);
                    neighbor->rhs_[w] -= afw * (-oldw + cell.var_[w]);
                }

                refresu = std::max(refresu, std::abs(acu * cell.var_[u]));
                refresv = std::max(refresv, std::abs(acv * cell.var_[v]));
                refresw = std::max(refresw, std::abs(acw * cell.var_[w]));

                if (std::abs(refresu) > ZERO)
                {
                    resu /= refresu;
                }

                if (std::abs(refresv) > ZERO)
                {
                    resv /= refresv;
                }

                if (std::abs(refresw) > ZERO)
                {
                    resw /= refresw;
                }
                
                maxresu = std::max(maxresu, resu);
                maxresv = std::max(maxresv, resv);
                maxresw = std::max(maxresw, resw);
            }

            if (maxresu < lin_tol && maxresv < lin_tol && maxresw < lin_tol)
            {
                break;
            }
        }
    }

    void jacobi(Mesh& mesh, PropertyIndex pindex)
    {
        int max_nsweep = 1;

        for (int i=0; i<max_nsweep; ++i)
        {
            for (Cell& cell: mesh.cells())
            {
                cell.var_[pindex] = cell.rhs_[pindex] / cell.aC_[pindex];
            }
        }
    }

    Residual gauss_seidel(Mesh& mesh, PropertyIndex pindex)
    {
        int max_nsweep = 1;
        double lin_tol = 1e-5;

        Scalar ref_res = BIG_NEG_NUM;
        Scalar max_res = BIG_NEG_NUM;

        for (int i=0; i<max_nsweep; ++i)
        {
            ref_res = BIG_NEG_NUM;
            max_res = BIG_NEG_NUM;

            for (Cell& cell: mesh.cells())
            {
                //if (&cell == &mesh.cells()[0])
                //{
                    //continue;
                //}

                double old = cell.var_[pindex];
                double aCTold = cell.aC_[pindex] * old;
                cell.var_[pindex] = cell.rhs_[pindex] / cell.aC_[pindex];

                for (FacePointer face: cell.faces())
                {
                    Neighbor neighbor = opposing_neighbor(face, &cell);

                    if (neighbor == face->l)
                    {
                        neighbor->rhs_[pindex] -= face->aFl_[pindex] * (-old + cell.var_[pindex]);
                    }
                    else
                    {
                        neighbor->rhs_[pindex] -= face->aFr_[pindex] * (-old + cell.var_[pindex]);
                    }
                }

                double aCT = cell.aC_[pindex] * cell.var_[pindex];
                ref_res = std::max(ref_res, std::abs(aCT));
                Scalar res = std::abs(aCTold - cell.rhs_[pindex]) / ref_res;
                //Scalar res = std::abs(aCTold - cell.rhs_[pindex]);
                max_res = std::max(max_res, res);
                
                //std::cout << "uu: " << cell.var_[u] << std::endl;
            }

            //std::cout << "max_reS: " << max_res << std::endl;
            //std::cout << "pp: " << mesh.cells()[0].var_[pp] << std::endl;

            if (max_res < lin_tol)
            {
                break;
            }
        }

        return max_res;
    }

    //void set_boundary_conditions(Mesh& mesh)
    //{
    //    {
    //        Patch* patch = mesh.patches().patch("inlet");

    //        if (patch != nullptr)
    //        {
    //            for (Cell& cell: patch->cells())
    //            {
    //                cell.properties_[0] = 0.0;
    //                cell.properties_[1] = 0.0;
    //                cell.properties_[2] = 0.0;
    //                cell.properties_[3] = 100.;
    //            }
    //        }
    //    }

    //    {
    //        Patch* patch = mesh.patches().patch("outlet");

    //        if (patch != nullptr)
    //        {
    //            for (Cell& cell: patch->cells())
    //            {
    //                cell.properties_[0] = 0.0;
    //                cell.properties_[1] = 0.0;
    //                cell.properties_[2] = 0.0;
    //                cell.properties_[3] = 500.;
    //            }
    //        }
    //    }

    //    {
    //        Patch* patch = mesh.patches().patch("wall");

    //        if (patch != nullptr)
    //        {
    //            for (Cell& cell: patch->cells())
    //            {
    //                cell.properties_[0] = 0.0;
    //                cell.properties_[1] = 0.0;
    //                cell.properties_[2] = 0.0;
    //            }
    //        }
    //    }
    //}

    //void init(Mesh& mesh)
    //{
    //    for (Cell& cell: mesh.cells())
    //    {
    //        cell.properties_[0] = 0.0;
    //        cell.properties_[1] = 0.0;
    //        cell.properties_[2] = 0.0;
    //        cell.properties_[3] = 0.0;
    //    }
    //}

    void ppcoef(Mesh& mesh)
    {
        for (Face& face: mesh.faces_)
        {
            Left l = face.l;
            Right r = face.r;

            Vector3& CF = face.CF_;
            Vector3 e = normalize(CF);
            Vector3 S = face.area_ * face.normal_;
            Vector3 E = dot(e, S) * e;
            //Vector3 T = S - E;

            double gl = face.g;
            double gr = 1. - gl;

            double& ul = l->var_[u];
            double& vl = l->var_[v];
            double& wl = l->var_[w];
            double& pl = l->var_[p];

            double& ur = r->var_[u];
            double& vr = r->var_[v];
            double& wr = r->var_[w];
            double& pr = r->var_[p];

            Coefs& lac = l->aC_;
            Coefs& rac = r->aC_;

            Coef& laf = face.aFl_[pp];
            Coef& raf = face.aFl_[pp];

            Vector3 D;
            if (r->patch_ == interior)
            {
                double gll = gl * r->volume_;
                double grr = gr * l->volume_;

                D(0) = gll / rac[u] + grr / lac[u]; 
                D(1) = gll / rac[v] + grr / lac[v]; 
                D(2) = gll / rac[w] + grr / lac[w]; 
            }
            else
            {
                Scalar lvol = l->volume_;

                D(0) = lvol / lac[u];
                D(1) = lvol / lac[v];
                D(2) = lvol / lac[w];
            }

            E(0) *= D(0);
            E(1) *= D(1);
            E(2) *= D(2);

            double Df = len(E) / face.lenCF_; 

            if (r->patch_ == interior)
            {
                laf = -Df;
                raf = -Df;
            }
            else
            {
                laf = 0.;
                raf = 0.;
            }
            
            lac[pp] -= laf; 
            rac[pp] -= raf; 

            double pb;
            Vector3 pdifb;
            double aveuf;
            double avevf;
            double avewf;

            if (r->patch_ == symmetry)
            {
                //pb = pr;
                //pdifb = l->grad_[p];

                //aveuf = ul;
                //avevf = vl;
                //avewf = wl;

                //aveuf = 0.;
                //avevf = 0.;
                //avewf = 0.;

                continue;
            }
            else if (r->patch_ == wall || r->patch_ == movingwall)
            {
                //pb = pr;
                //pdifb = l->grad_[p];

                //aveuf = ul;
                //avevf = vl;
                //avewf = wl;

                continue;
            }
            else if (r->patch_ == inlet)
            {
                pdifb = l->grad_[p];

                aveuf = ul;
                avevf = vl;
                avewf = wl;
            }
            else if (r->patch_ == outlet)
            {
                pdifb = l->grad_[p];

                aveuf = ul;
                avevf = vl;
                avewf = wl;
            }
            else
            {
                pb = pr;
                pdifb = gl * r->grad_[p] + gr * l->grad_[p];

                aveuf = gl * ur + gr * ul;
                avevf = gl * vr + gr * vl;
                avewf = gl * wr + gr * wl;
            }

            double pdifa = Df * (pb - pl);

            pdifb(0) *= D(0);
            pdifb(1) *= D(1);
            pdifb(2) *= D(2);

            Vector3 aveVf(aveuf, avevf, avewf);

            double a = dot(pdifb, S) - pdifa;
            double d = dot(aveVf, S);

            double b = a + d;
            double c = a - d;

            l->rhs_[pp] -= b;
            r->rhs_[pp] -= c;

            //l->summ_ += b;
            //r->summ_ += c;
        }
    }

    Scalar maxpp(Mesh& mesh)
    {
        double maxmf = BIG_NEG_NUM;

        for (Cell& cell: mesh.cells())
        {
            //maxmf = std::max(maxmf, std::abs(cell.summ_)); 
            maxmf = std::max(maxmf, std::abs(cell.rhs_[pp])); 
        }

        return maxmf;
    }

    void correct(Mesh& mesh)
    {
        mesh.LS_grad(pp);

        for (Cell& cell: mesh.cells())
        {
            Gradient& grad_pp = cell.grad_[pp];

            Scalar vol = cell.volume_;

            double up = -vol * grad_pp(0) / cell.aC_[u];
            double vp = -vol * grad_pp(1) / cell.aC_[v];
            double wp = -vol * grad_pp(2) / cell.aC_[w];

            cell.var_[u] += up; 
            cell.var_[v] += vp; 
            cell.var_[w] += wp; 
            cell.var_[p] += 0.01 * cell.var_[pp]; 
        }
    }

    void updatepp(Mesh& mesh)
    {
        {
            Patch* patch = mesh.patches().patch(symmetry);

            if (patch != nullptr)
            {
                for (Cell& cell: patch->cells())
                {
                    CellPointer interior = cell.neighbors()[0];

                    cell.var_[pp] = interior->var_[pp];
                }
            }
        }
        {
            Patch* patch = mesh.patches().patch(wall);

            if (patch != nullptr)
            {
                for (Cell& cell: patch->cells())
                {
                    CellPointer interior = cell.neighbors()[0];

                    cell.var_[pp] = interior->var_[pp];
                }
            }
        }
        {
            Patch* patch = mesh.patches().patch(movingwall);

            if (patch != nullptr)
            {
                for (Cell& cell: patch->cells())
                {
                    CellPointer interior = cell.neighbors()[0];

                    cell.var_[pp] = interior->var_[pp];
                }
            }
        }
    }

    void updatebc(Mesh& mesh)
    {
        {
            Patch* patch = mesh.patches().patch(symmetry);

            if (patch != nullptr)
            {
                for (Cell& cell: patch->cells())
                {
                    CellPointer interior = cell.neighbors()[0];

                    cell.var_[p] = interior->var_[p];
                }
            }
        }
        {
            Patch* patch = mesh.patches().patch(wall);

            if (patch != nullptr)
            {
                for (Cell& cell: patch->cells())
                {
                    CellPointer interior = cell.neighbors()[0];

                    cell.var_[p] = interior->var_[p];
                }
            }
        }
        {
            Patch* patch = mesh.patches().patch(movingwall);

            if (patch != nullptr)
            {
                for (Cell& cell: patch->cells())
                {
                    CellPointer interior = cell.neighbors()[0];

                    //cell.var_[p] = interior->var_[p] + dot(interior->grad_[p], CF);
                    cell.var_[p] = interior->var_[p];
                }
            }
        }
    }

    Scalar checkcon(Mesh& mesh)
    {
        for (Face& face: mesh.faces_)
        {
            Left l = face.l;
            Right r = face.r;

            Vector3& CF = face.CF_;
            Vector3 e = normalize(CF);
            Vector3 S = face.area_ * face.normal_;
            Vector3 E = dot(e, S) * e;

            double gl = face.g;
            double gr = 1. - gl;

            Vector3 D;
            if (r->patch_ == interior)
            {
                double gll = gl * r->volume_;
                double grr = gr * l->volume_;

                D(0) = gll / r->aC_[u] + grr / l->aC_[u]; 
                D(1) = gll / r->aC_[v] + grr / l->aC_[v]; 
                D(2) = gll / r->aC_[w] + grr / l->aC_[w]; 
            }
            else
            {
                Scalar lvol = l->volume_;

                D(0) = lvol / l->aC_[u];
                D(1) = lvol / l->aC_[v];
                D(2) = lvol / l->aC_[w];
            }

            E(0) *= D(0);
            E(1) *= D(1);
            E(2) *= D(2);

            double Df = len(E) / face.lenCF_; 

            //double pdifa = Df * (pb - pl);
            Scalar mfp = Df * (r->var_[pp] - l->var_[pp]);

            //if (l->tag_ == 1219)
            //{
            //    std::cout << "r: " << r->tag_ << std::endl;
            //}
            //if (r->tag_ == 1219)
            //{
            //    std::cout << "l: " << l->tag_ << std::endl;
            //}

            l->summ_ -= mfp;
            r->summ_ += mfp;

            //if (l->tag_ == 1219)
            //{
            //    std::cout << "mfpl: " << mfp << std::endl;
            //    //std::cout << "updated: " << l->summ_ << std::endl;
            //}
            //if (r->tag_ == 1219)
            //{
            //    std::cout << "mfpr: " << mfp << std::endl;
            //    //std::cout << "updated: " << r->summ_ << std::endl;
            //}
        }

        double maxmf = BIG_NEG_NUM;

        int jj;
        for (int ii = 0; ii<mesh.cells().size(); ++ii)
        //for (Cell& cell: mesh.cells())
        {
            Cell& cell = mesh.cells()[ii];
            if (maxmf < std::abs(cell.summ_))
            {
                maxmf = std::max(maxmf, std::abs(cell.summ_)); 
                jj = ii;
            }
            //assert(std::abs(cell.summ_) < 1e-6);
            maxmf = std::max(maxmf, std::abs(cell.summ_)); 
        }

        //std::cout << "jj: " << mesh.cells()[jj].tag_ << std::endl;

        return maxmf;
    }

    void simple(Mesh& mesh, TemporalScheme tscheme)
    {
        std::ofstream out;
        out.open("res.dat", std::ios_base::app);

        int max_nsweep = 10000;
        double lin_tol = 1e-10;

        diff_flux(mesh);
        mesh.LS_grad_uvwp();

        Residual maxresu, maxresv, maxresw, maxrespp, continuity;

        for (int i=0; i<max_nsweep; ++i)
        {
            solve_steady(mesh, maxresu, maxresv, maxresw);

            for (Cell& cell: mesh.cells())
            {
                cell.aC_[pp] = 0.;
                cell.rhs_[pp] = 0.;
                cell.var_[pp] = 0.;
                cell.summ_ = 0.;
            }

            ppcoef(mesh);
            continuity = maxpp(mesh);

            std::cout << i << " " << continuity << std::endl;

            if (continuity < lin_tol)
            {
                break;
            }

            out << i;
            out << " ";
            out << continuity;
            out << std::endl;

            //maxrespp = gauss_seidel(mesh, pp);
            jacobi(mesh, pp);

            updatepp(mesh);
            correct(mesh);

            mesh.LS_grad_uvwp();

            updatebc(mesh);
        }

        out.close();
    }
}
