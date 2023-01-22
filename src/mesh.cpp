#include <iostream>
#include <algorithm>
#include <fstream>
#include "mesh.h"

namespace Liver
{
    void set_CF(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            Left l = face.l;
            Right r = face.r;

            face.CF_ = r->centroid_ - l->centroid_;
            face.lenCF_ = len(face.CF_);
        }
    }

    void set_volumes(Mesh& mesh)
    {
        for (Patch& patch: mesh.patches_)
        {
            for (Cell& cell: patch.cells())
            {
                cell.volume();
            }
        }
    }

    void set_centroids(Mesh& mesh)
    {
        for (Patch& patch: mesh.patches_)
        {
            for (Cell& cell: patch.cells())
            {
                cell.centroid();
            }
        }
    }

    void set_face_normals(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            face.normal();
        }
    }

    void set_face_centroids(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            face.centroid();
        }
    }

    void set_face_areas(Mesh& mesh)
    {
        for (Face& face: mesh.faces())
        {
            face.area();
        }
    }

    void Mesh::set_parent_cells_of_vertices()
    {
        for (Patch& patch: patches_)
        {
            for (Cell& cell: patch.cells())
            {
                cell.set_parent_cell_of_vertices();
            }
        }
    }

    void Mesh::add_point(Point point, NumberOfPoints npoints)
    {
        if (points_.capacity() < npoints)
        {
            points_.reserve(npoints);
        }

        points_.push_back(point);
    }

    void Mesh::add_cell(Cell cell, NumberOfCells ncells)
    {
        patches_.add_cell(cell, ncells);
    }

    Points& Mesh::points()
    {
        return points_;
    }

    Faces& Mesh::faces()
    {
        return faces_;
    }

    Patches& Mesh::patches()
    {
        return patches_;
    }

    Faces make_faces(Cell& cell)
    {
        Shape shape = cell.shape();
        PointPointers points = cell.points();
        Faces faces;

        if (shape == Shape::quad)
        {
            auto pts0 = {points[0], points[1], points[2], points[3]};
            faces.push_back(Face("dummy", pts0));
        }
        else if (shape == Shape::tri)
        {
            auto pts0 = {points[0], points[1], points[2]};
            faces.push_back(Face("dummy", pts0));
        }
        else if (shape == Shape::tet)
        {
            auto pts0 = {points[1], points[2], points[0]};
            auto pts1 = {points[0], points[3], points[1]};
            auto pts2 = {points[0], points[2], points[3]};
            auto pts3 = {points[2], points[1], points[3]};

            faces.push_back(Face("dummy", pts0));
            faces.push_back(Face("dummy", pts1));
            faces.push_back(Face("dummy", pts2));
            faces.push_back(Face("dummy", pts3));
        }
        else if (shape == Shape::pri)
        {
            auto pts0 = {points[0], points[1], points[2]};
            auto pts1 = {points[5], points[4], points[3]};
            auto pts2 = {points[1], points[0], points[3], points[4]};
            auto pts3 = {points[3], points[0], points[2], points[5]};
            auto pts4 = {points[2], points[1], points[4], points[5]};

            faces.push_back(Face("dummy", pts0));
            faces.push_back(Face("dummy", pts1));
            faces.push_back(Face("dummy", pts2));
            faces.push_back(Face("dummy", pts3));
            faces.push_back(Face("dummy", pts4));
        }
        else if (shape == Shape::hex)
        {
            auto pts0 = {points[0], points[1], points[2], points[3]};
            auto pts1 = {points[6], points[5], points[4], points[7]};
            auto pts2 = {points[2], points[1], points[5], points[6]};
            auto pts3 = {points[4], points[0], points[3], points[7]};
            auto pts4 = {points[1], points[0], points[4], points[5]};
            auto pts5 = {points[3], points[2], points[6], points[7]};

            faces.push_back(Face("dummy", pts0));
            faces.push_back(Face("dummy", pts1));
            faces.push_back(Face("dummy", pts2));
            faces.push_back(Face("dummy", pts3));
            faces.push_back(Face("dummy", pts4));
            faces.push_back(Face("dummy", pts5));
        }
        else
        {
            assert(false);
        }

        return faces;
    }

    //void make_edges(const Face& face)
    //{
    //    PointPointers pts = face.points();

    //    EdgePointers edges;

    //    if (pts.size() == 3)
    //    {
    //        edges.push_back(Edge(PointPointers{points[0], points[1])});
    //        edges.push_back(Edge(PointPointers{points[1], points[2])});
    //        edges.push_back(Edge(PointPointers{points[2], points[0])});
    //    }
    //    else if (pts.size() == 4)
    //    {
    //        edges.push_back(Edge(PointPointers{points[0], points[1])});
    //        edges.push_back(Edge(PointPointers{points[0], points[1])});
    //        edges.push_back(Edge(PointPointers{points[0], points[1])});
    //    }

    //    return edges;
    //}

    //void Mesh::create_edges()
    //{
    //    for (Face& face: faces_)
    //    {
    //        const PointPointers& points = face.points();

    //        for (int i=0; i<points.size(); ++i)
    //        {
    //            if (i == points.size() - 1)
    //            {
    //                break;
    //            }

    //            edges_.push_back(Edge(PointPointers{points[i], points[i+1]}));
    //            face.add_edge(&edges_.back());
    //        }
    //    }
    //}

    void Mesh::create_faces()
    {
        int n_faces = 0;

        for (Cell& cell: cells())
        {
            if (cell.shape() == Shape::tet)
            {
                n_faces += 4;
            }
            else if (cell.shape() == Shape::hex)
            {
                n_faces += 6;
            }
            else if (cell.shape() == Shape::pri)
            {
                n_faces += 5;
            }
            else
            {
                assert(false);
            }
        }

        for (Patch& patch: patches_)
        {
            if (patch.name() == interior)
            {
                continue;
            }

            n_faces += patch.cells().size();
        }

        n_faces /= 2; 

        faces_.reserve(n_faces);

        for (Cell& cell: cells())
        {
            Faces faces = make_faces(cell);

            for (Face& face: faces)
            {
                bool this_face_exists = false;

                const PointPointers& fpoints = face.points();
                PointPointer p = fpoints.front();

                FaceTag ftag;

                std::vector<int> tags;
                for (PointPointer q: fpoints)
                {
                    tags.push_back(q->tag());
                }

                std::sort(tags.begin(), tags.end(), std::less<int>());

                for (int t: tags)
                {
                    ftag.append(",");
                    ftag.append(std::to_string(t));
                }

                face.set_tag(ftag);

                for (ParentFace pf: p->parent_faces())
                {
                    if (pf->tag() == ftag)
                    {
                        this_face_exists = true;
                        break;
                    }
                }

                if (this_face_exists)
                {
                    continue;
                }

                for (ParentCell pc: p->parent_cells())
                {
                    if (pc->tag() == cell.tag())
                    {
                        continue;
                    }

                    int count = 0;

                    for (PointPointer q: fpoints)
                    {
                        for (ParentCell pc2: q->parent_cells())
                        {
                            if (pc2->tag() == pc->tag())
                            {
                                ++count;
                            }
                        }
                    }

                    if (count == face.points().size())
                    {
                        face.add_parent_cell(&cell);
                        face.add_parent_cell(pc);

                        faces_.push_back(face);
                        assert(faces_.size() <= n_faces);

                        cell.add_face(&faces_.back());
                        pc->add_face(&faces_.back());

                        for (PointPointer r: fpoints)
                        {
                            r->add_parent_face(&faces_.back());
                        }

                        cell.add_neighbor(pc);
                        pc->add_neighbor(&cell);

                        break;
                    }
                }
            }
        }
    }

    Cells& Mesh::cells()
    {
        Patch* patch = patches_.patch(interior);

        if (patch == nullptr)
        {
            assert(false);
        }

        return patch->cells();
    }

    void Mesh::export_vtk(FileName file_name)
    {
        int cell_list_size = 0;
        std::ofstream out;

        out.open (file_name);

        out << "# vtk DataFile Version 3.0" << std::endl;
        out << "All in VTK format" << std::endl;
        out << "ASCII" << std::endl;
        out << "DATASET UNSTRUCTURED_GRID" << std::endl;
        out << "POINTS " << points_.size() << " float" << std::endl;

        for (const Point& p: points_)
        {
            const Coordinate& c = p.coordinate();

            out << c(0);
            out << " ";
            out << c(1);
            out << " ";
            out << c(2);
            out << std::endl;
        }

        Cells& cells = this->cells();

        // get cell list size.
        for (Cell& cell: cells)
        {
            cell_list_size += (cell.points().size() + 1);
        }

        out << std::endl;    
        out << "CELLS " << cells.size() << " " << cell_list_size << std::endl;

        for (Cell& cell: cells)
        {        
            out << cell.points().size();
            out << " ";

            for (PointPointer p: cell.points())
            {
                out << p - &points_[0];
                out << " ";
            }

            out << std::endl;
        }    

        out << "CELL_TYPES " << cells.size() << std::endl;
        for (Cell& cell: cells)
        {
            //https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

            if (cell.shape() == Shape::tet)
            {
                out << 10;
            }
            else if (cell.shape() == Shape::hex)
            {
                out << 12;
            }        
            else if (cell.shape() == Shape::pri)
            {
                out << 13;
            }        
            else if (cell.shape() == Shape::tri)
            {
                out << 5;
            }        
            else if (cell.shape() == Shape::quad)
            {
                out << 9;
            }        
            else
            {
                assert(false);
            }

            out << std::endl;
        }

        out << "CELL_DATA " << cells.size() << std::endl;

        out << "SCALARS " << "u " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.var_[u] << std::endl;
            assert(!std::isnan(cell.var_[u]));
        } 
        out << "SCALARS " << "v " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.var_[v] << std::endl;
        } 
        out << "SCALARS " << "w " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.var_[w] << std::endl;
        } 
        out << "SCALARS " << "phi " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.var_[phi] << std::endl;
        } 
        out << "SCALARS " << "p " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.var_[p] << std::endl;
        } 
        out << "SCALARS " << "pp " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.var_[pp] << std::endl;
        } 
        out << "VECTORS " << "V " << "float" << std::endl;
        for (Cell& cell: cells)
        {
            out << cell.var_[u];
            out << " ";
            out << cell.var_[v];
            out << " ";
            out << cell.var_[w] << std::endl;
        } 

        out << "SCALARS " << "x " << "float " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.centroid()(0) << std::endl;
        } 

        out << "SCALARS " << "tag " << "int " << "1" << std::endl;
        out << "LOOKUP_TABLE default" << std::endl;    
        for (Cell& cell: cells)
        {
            out << cell.tag() << std::endl;
        } 

        out.close();
    }
}
