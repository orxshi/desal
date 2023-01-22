#ifndef MESH_H
#define	MESH_H

#include "point.h"
#include "face.h"
#include "patch.h"

namespace Liver
{
    class Mesh
    {    
        public:

            void LS_grad(PropertyIndex);
            void LS_grad_uvwp();
            void GG_grad(PropertyIndex);
            void GG_grad_uvwp();
            void export_vtk(FileName file_name);
            void create_faces();
            void read_mesh_gmsh(FileName);
            void add_point(Point, NumberOfPoints=0);
            void add_cell(Cell, NumberOfCells=0);
            void set_parent_cells_of_vertices();
            Points& points();
            Faces& faces();
            Patches& patches();
            Cells& cells();


            Points points_;
            Faces faces_;
            Patches patches_;
    };

    void set_centroids(Mesh&);
    void set_CF(Mesh&);
    void set_face_areas(Mesh&);
    void set_face_centroids(Mesh&);
    void set_face_normals(Mesh&);
    void set_volumes(Mesh&);
}

#endif
