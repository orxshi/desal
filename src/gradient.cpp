#include "cell.h"
#include "mesh.h"

namespace Liver
{
    void Mesh::GG_grad_uvwp()
    {
        Vector3 S, a;
        Centroid rf, rl, rr;
        Cells& cells = this->cells();
        Face& first = faces_[0];
        std::vector<Properties> phifs(faces_.size());

        for (int k=0; k<2; ++k)
        {
            for (Face& face: faces_)
            {
                Left l = face.l;
                Right r = face.r;

                int j = &face - &first;

                PatchName rpatch = r->patch_;
                Gradients& lgrad = l->grad_;
                Properties& lvar = l->var_;
                Properties& rvar = r->var_;

                rf = face.centroid_;
                rl = l->centroid_;
                rr = r->centroid_;

                if (k == 0)
                {
                    phifs[j][u] = 0.5 * (lvar[u] + rvar[u]);
                    phifs[j][v] = 0.5 * (lvar[v] + rvar[v]);
                    phifs[j][w] = 0.5 * (lvar[w] + rvar[w]);
                    phifs[j][p] = 0.5 * (lvar[p] + rvar[p]);
                }
                else
                {
                    if (rpatch == interior)
                    {
                        Gradients& rgrad = r->grad_;
                        a = rf - 0.5 * (rl + rr);

                        phifs[j][u] += 0.5 * dot((lgrad[u] + rgrad[u]), a);
                        phifs[j][v] += 0.5 * dot((lgrad[v] + rgrad[v]), a);
                        phifs[j][w] += 0.5 * dot((lgrad[w] + rgrad[w]), a);
                        phifs[j][p] += 0.5 * dot((lgrad[p] + rgrad[p]), a);
                    }
                    else
                    {
                        phifs[j][u] += dot((lgrad[u]), a);
                        phifs[j][v] += dot((lgrad[v]), a);
                        phifs[j][w] += dot((lgrad[w]), a);
                        phifs[j][p] += dot((lgrad[p]), a);
                    }
                }
            }

            for (Cell& cell: cells)
            {
                Gradients& grad  = cell.grad_;

                grad[u] = 0.;
                grad[v] = 0.;
                grad[w] = 0.;
                grad[p] = 0.;
            }

            for (Cell& cell: cells)
            {
                Gradients& grad = cell.grad_;

                for (FacePointer face: cell.faces_)
                {
                    int j = face - &first;

                    Properties& phif = phifs[j];

                    S = face->area_ * face->normal_;

                    if (cell.tag_ == face->l->tag_)
                    {
                        grad[u] += phif[u] * S;
                        grad[v] += phif[v] * S;
                        grad[w] += phif[w] * S;
                        grad[p] += phif[p] * S;
                    }
                    else
                    {
                        grad[u] -= phif[u] * S;
                        grad[v] -= phif[v] * S;
                        grad[w] -= phif[w] * S;
                        grad[p] -= phif[p] * S;
                    }
                }

                Scalar vol = cell.volume_;

                grad[u] /= vol;
                grad[v] /= vol;
                grad[w] /= vol;
                grad[p] /= vol;
            }
        }
    }

    void Mesh::GG_grad(PropertyIndex i)
    {
        Vector3 S, a;
        Centroid rf, rl, rr;
        Cells& cells = this->cells();
        Face& first = faces_[0];
        std::vector<Property> phif(faces_.size());

        for (int k=0; k<2; ++k)
        {
            for (Face& face: faces_)
            {
                Left l = face.l;
                Right r = face.r;

                int j = &face - &first;

                PatchName rpatch = r->patch_;
                Gradients& lgrad = l->grad_;
                Properties& lvar = l->var_;
                Properties& rvar = r->var_;

                rf = face.centroid_;
                rl = l->centroid_;
                rr = r->centroid_;

                if (k == 0)
                {
                    phif[j] = 0.5 * (lvar[u] + rvar[u]);
                }
                else
                {
                    if (rpatch == interior)
                    {
                        Gradients& rgrad = r->grad_;
                        a = rf - 0.5 * (rl + rr);

                        phif[j] += 0.5 * dot((lgrad[i] + rgrad[i]), a);
                    }
                    else
                    {
                        phif[j] += dot((lgrad[i]), a);
                    }
                }
            }

            for (Cell& cell: cells)
            {
                Gradients& grad  = cell.grad_;

                grad[i] = 0.;
            }

            for (Cell& cell: cells)
            {
                Gradients& grad = cell.grad_;

                for (FacePointer face: cell.faces_)
                {
                    int j = face - &first;

                    Property& phf = phif[j];

                    S = face->area_ * face->normal_;

                    if (cell.tag_ == face->l->tag_)
                    {
                        grad[i] += phf * S;
                    }
                    else
                    {
                        grad[i] -= phf * S;
                    }
                }

                Scalar vol = cell.volume_;

                grad[i] /= vol;
            }
        }
    }
    void Mesh::LS_grad_uvwp()
    {    
        double dq, a, b, c;
        Cells& cells = this->cells();

        for (Cell& cell: cells)
        {
            Gradients& grad = cell.grad_;

            grad[u] = 0.;
            grad[v] = 0.;
            grad[w] = 0.;
            grad[p] = 0.;
        }

        for (Face& face: faces_)
        {
            Left l = face.l;
            Right r = face.r;

            Vector3& CF = face.CF_;
            double CFx = CF(0);
            double CFy = CF(1);
            double CFz = CF(2);
            double lenCF = face.lenCF_;

            Properties& lvar = l->var_;
            Properties& rvar = r->var_;
            PatchName rpatch = r->patch_;

            auto aaa = [&](PropertyIndex i)
            {
                Gradient& lgrad = l->grad_[i];

                dq = (rvar[i] - lvar[i]) / lenCF;

                a = CFx * dq;
                b = CFy * dq;
                c = CFz * dq;

                lgrad(0) += a;
                lgrad(1) += b;
                lgrad(2) += c;

                if (rpatch == interior)
                {
                    Gradient& rgrad = r->grad_[i];

                    rgrad(0) += a;
                    rgrad(1) += b;
                    rgrad(2) += c;
                }
            };

            aaa(u);
            aaa(v);
            aaa(w);
            aaa(p);
        }

        for (Cell& cell: cells)
        {
            Gradients& grad = cell.grad_;
            auto& LU = cell.LU_;

            grad[u] = LU_solve(LU, grad[u]);
            grad[v] = LU_solve(LU, grad[v]);
            grad[w] = LU_solve(LU, grad[w]);
            grad[p] = LU_solve(LU, grad[p]);
        }
    }
    void Mesh::LS_grad(PropertyIndex i)
    {    
        double dphi;
        Cells& cells = this->cells();

        for (Cell& cell: cells)
        {
            cell.grad_[i] = 0.;
        }

        for (Face& face: faces_)
        {
            Left l = face.l;
            Right r = face.r;

            Vector3& CF = face.CF_;

            dphi = (r->var_[i] - l->var_[i]) / face.lenCF_;

            double a = CF(0) * dphi;
            double b = CF(1) * dphi;
            double c = CF(2) * dphi;

            Gradient& lgrad = l->grad_[i];

            lgrad(0) += a;
            lgrad(1) += b;
            lgrad(2) += c;

            if (r->patch_ == interior)
            {
                Gradient& rgrad = r->grad_[i];

                rgrad(0) += a;
                rgrad(1) += b;
                rgrad(2) += c;
            }
        }

        for (Cell& cell: cells)
        {
            cell.grad_[i] = LU_solve(cell.LU_, cell.grad_[i]);
        }
    }
    //Gradient Cell::LS_grad(PropertyIndex i)
    //{    
    //    double dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz, w, dphi;
    //    Vector3 d, dif(0., 0., 0.);

    //    for (FacePointer face: faces_)
    //    {
    //        Neighbor neighbor = opposing_neighbor(face, this);

    //        auto d = neighbor->centroid_ - this->centroid_;

    //        dx = d(0);
    //        dy = d(1);
    //        dz = d(2);

    //        dxx = dx * dx;
    //        dxy = dx * dy;
    //        dxz = dx * dz;

    //        dyy = dy * dy;
    //        dyz = dy * dz;

    //        dzz = dz * dz;

    //        dphi = (neighbor->var_[i] - var_[i]) / std::sqrt(dx*dx + dy*dy + dz*dz);

    //        dif(0) += dx * dphi;
    //        dif(1) += dy * dphi;
    //        dif(2) += dz * dphi;
    //    }

    //    return LU_solve(LU_, dif);
    //}

    void Cell::LS_coef()
    {
        double dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz, w;
        Vector3 d;
        SquareMatrixD<N_DIM> A;

        for (FacePointer face: faces_)
        {
            Neighbor neighbor = opposing_neighbor(face, this);

            //if (neighbor->patch_ == "symmetry")
            //{
            //    continue;
            //}

            auto d = neighbor->centroid() - this->centroid();

            dx = d(0);
            dy = d(1);
            dz = d(2);

            dxx = dx * dx;
            dxy = dx * dy;
            dxz = dx * dz;

            dyy = dy * dy;
            dyz = dy * dz;

            dzz = dz * dz;

            w = 1. / std::sqrt(dx*dx + dy*dy + dz*dz);

            A(0,0) += w * dxx;
            A(0,1) += w * dxy;
            A(0,2) += w * dxz;
            A(1,1) += w * dyy;
            A(1,2) += w * dyz;
            A(2,2) += w * dzz;
        }

        A(1,0) = A(0,1);
        A(2,0) = A(0,2);
        A(2,1) = A(1,2);

        auto [L, U] = LU(A);
        LU_ = LU_combine(L, U);
    }
}
