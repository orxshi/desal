#ifndef GAUSS_SEIDEL_H
#define	GAUSS_SEIDEL_H

#include "mesh.h"

namespace Liver
{
    enum TemporalScheme
    {
        steady,
        forward_euler,
        backward_euler,
    };

    void ppcoef(Mesh&);
    void spatial_discre(Mesh&);
    void temporal_discre(Mesh&, TimeStep, TemporalScheme, PropertyIndex);
    void temporal_discre(Mesh&, TimeStep, TemporalScheme);
    void solve_unsteady(Mesh&, TimeStep, Time, TemporalScheme, PropertyIndex);
    void solve_steady(Mesh&, Residual&, Residual&, Residual&);
    void gauss_seidel_mom(Mesh& mesh, Residual&, Residual&, Residual&);
    void jacobi_mom(Mesh&, Residual&, Residual&, Residual&);
    Residual gauss_seidel(Mesh&, PropertyIndex);
    //Residual jacobi(Mesh&, PropertyIndex);
    //void set_boundary_conditions(Mesh&);
    //void init(Mesh&);
    void lefts_rights(Mesh&);
    void interfac(Mesh&);
    void simple(Mesh&, TemporalScheme);
    void diff_flux(Mesh&);
    Scalar maxpp(Mesh&);
}

#endif
