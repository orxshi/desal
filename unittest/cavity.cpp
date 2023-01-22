#include <cassert>
#include <fstream>
#include <iostream>
#include "gauss_seidel.h"
#include "blood_properties.h"

using namespace Liver;

int main()
{
    RHO = 10.;
    MU = 0.1;

    double L = 0.1;
    double Re = 10.0;
    double U = Re * MU / (RHO * L);

    std::cout << "Re: " << RHO * U * L / MU << std::endl;
    std::cout << "U: " << U << std::endl;

    Mesh mesh;
    mesh.read_mesh_gmsh("/home/orhan/liver/unittest/cavity.msh");
    
    mesh.set_parent_cells_of_vertices();
    mesh.create_faces();

    set_face_areas(mesh);
    set_face_centroids(mesh);
    set_face_normals(mesh);
    set_volumes(mesh);
    set_centroids(mesh);

    lefts_rights(mesh);
    set_CF(mesh);
    interfac(mesh);

    for (Cell& cell: mesh.cells())
    {
        cell.LS_coef();
    }

    for (Cell& cell: mesh.cells())
    {
        cell.var_[u] = 0.0;
        cell.var_[v] = 0.0;
        cell.var_[w] = 0.0;
        cell.var_[p] = 0.0;
    }
    {
        Patch* patch = mesh.patches().patch(symmetry);

        if (patch != nullptr)
        {
            for (Cell& cell: patch->cells())
            {
                cell.var_[u] = 0.0;
                cell.var_[v] = 0.0;
                cell.var_[w] = 0.0;
                cell.var_[p] = 0.0;
            }
        }
    }
    {
        Patch* patch = mesh.patches().patch(wall);

        if (patch != nullptr)
        {
            for (Cell& cell: patch->cells())
            {
                cell.var_[u] = 0.0;
                cell.var_[v] = 0.0;
                cell.var_[w] = 0.0;
                cell.var_[p] = 0.0;
            }
        }
    }
    {
        Patch* patch = mesh.patches().patch(movingwall);

        if (patch != nullptr)
        {
            for (Cell& cell: patch->cells())
            {
                cell.var_[u] = U;
                cell.var_[v] = 0.0;
                cell.var_[w] = 0.0;
                cell.var_[p] = 0.0;
            }
        }
    }

    simple(mesh, forward_euler);

    mesh.export_vtk("cavity.vtk");

    return 0;
}
