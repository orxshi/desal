#include <fstream>
#include <iostream>
#include <cassert>
#include <map>
#include "mesh.h"

namespace Liver
{
    std::ifstream& go_to_beg_of_line(std::ifstream& file, int num)
    {
        file.seekg(std::ios::beg);
        for(int i=0; i<num-1; ++i)
        {
            file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        }

        return file;
    }

    void Mesh::read_mesh_gmsh(FileName file_name)
    {
        enum class gmsh
        {
            tri = 2,
            quad = 3,
            tet = 4,
            hex = 5,
            pri = 6,
        };

        std::string temps;
        std::string phys_name;
        int tempi;
        int line_number = 1;
        int tag, igs, n_tags, geo, n_part, n_total, n_point, part_tag;
        int phys_tag, n_phys;
        gmsh gs;
        //std::map<int, std::string> phys_tag_name;

        std::ifstream in;
        in.open(file_name);
        if (!in.is_open())
        {
            std::cout << "file " << file_name << " could not be opened." << std::endl;
            return;
        }
        assert(in.is_open());

        in >> temps; // mesh format.
        in >> temps; // version.
        assert(temps == "2.2"); // msh2 format.
        in >> temps; // version.
        in >> temps; // version.
        in >> temps; // end mesh format.
        //in >> temps; // $PhysicalNames
        //in >> n_phys; // number of physical names
        //for (int i=0; i<n_phys; ++i)
        //{
        //    in >> temps; // dimension
        //    in >> phys_tag; // physical tag
        //    //in >> phys_name; // physical name
        //    //phys_name.erase(0, 1);
        //    //phys_name.erase(phys_name.size() - 1);
        //    //phys_tag_name.insert(std::make_pair(phys_tag, phys_name));
        //}
        //in >> temps; // $EndPhysicalNames
        in >> temps; // nodes.
        in >> n_point; // number under "$Nodes" is the number of points.
        line_number += 5;    
        //line_number += n_phys + 3;

        for (int i=0; i<n_point; ++i)
        {
            int ptag;
            double x, y, z;
            in >> ptag;
            in >> x;
            in >> y;
            in >> z;
            ++line_number;

            add_point(Point(ptag, Coordinate(x, y, z)), n_point);
        }

        in >> temps; // end of nodes.
        in >> temps; // elements.
        in >> n_total; // the number under "$Elements" is total number of elements which includes boundary faces and cells.
        line_number += 3;

        // read elements.
        int n_bface = 0;
        for (int e=0; e<n_total; ++e)
        {
            in >> tag; // tag of boundary face or cell.
            in >> igs; // geometric shape of boundary face or cell.        
            gs = static_cast<gmsh>(igs);
            in >> n_tags; // number of GMSH tags.

            int tag_count = 0;

            // read GMSH tags.
            while (n_tags > 0)
            {
                // read physical number.
                in >> tempi;
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> geo; // geometrical number.
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> n_part; // number of partitions to which element belongs.
                for (int i=0; i<n_part; ++i)
                {
                    in >> temps;
                }            

                break;
            }

            if (gs == gmsh::tri)
            {
                in >> tempi;
                in >> tempi;
                in >> tempi;
            }
            else if (gs == gmsh::quad)
            {
                in >> tempi;
                in >> tempi;
                in >> tempi;
                in >> tempi;
            }
            else
            {
                break;
            }        

            ++n_bface;
        }

        go_to_beg_of_line(in, line_number);

        for (int e=0; e<n_total; ++e)
        {
            in >> tag; // read tag of boundary face or cell.
            in >> igs; // read geometric shape of boundary face or cell.        
            gs = static_cast<gmsh>(igs);
            in >> n_tags; // read number of GMSH tags.

            int tag_count = 0;

            while (n_tags > 0)
            {
                // read physical number.
                in >> phys_tag;
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> geo; // read geometrical number.
                ++ tag_count;
                if (tag_count == n_tags) break;
                in >> n_part; // read number of partitions to which element belongs.
                //assert(n_part == 1);
                for (int i=0; i<n_part; ++i)
                {
                    in >> part_tag;
                }            

                break;
            }

            if (e < n_bface)
            {
                Shape shape = Shape::undef;
                int n_vertex;
                switch (gs)
                {
                    case gmsh::tri:					    
                        n_vertex = 3;
                        shape = Shape::tri;
                        break;
                    case gmsh::quad:			
                        n_vertex = 4;
                        shape = Shape::quad;
                        break;
                    case gmsh::tet:
                        n_vertex = 4;						
                        shape = Shape::tet;
                        break;
                    case gmsh::hex:
                        n_vertex = 8;						
                        shape = Shape::hex;
                        break;
                    case gmsh::pri:
                        n_vertex = 6;						
                        shape = Shape::pri;
                        break;
                    default:
                        std::cout << "file: " << file_name << std::endl;
                        std::cout << "igs: " << igs << std::endl;
                        std::cout << "invalid n_vertex" << std::endl;
                        std::cout << "igs = " << igs << std::endl;
                        std::cout << "tag = " << tag << std::endl;
                        std::cout << "n_tags = " << n_tags << std::endl;
                        assert(false);
                }

                PointPointers vtx;
                vtx.reserve(n_vertex);

                for (int i=0; i<n_vertex; ++i)
                {
                    in >> tempi;
                    vtx.push_back(&points_[tempi-1]);
                }

                //add_cell(Cell(tag-1, vtx, shape, phys_tag_name[phys_tag]));
                add_cell(Cell(tag-1, vtx, shape, PatchName(phys_tag)));
            }
            else
            {
                Shape shape;
                int n_vertex;
                switch (gs)
                {
                    case gmsh::tri:					    
                        n_vertex = 3;
                        shape = Shape::tri;
                        break;
                    case gmsh::quad:			
                        n_vertex = 4;
                        shape = Shape::quad;
                        break;
                    case gmsh::tet:
                        n_vertex = 4;						
                        shape = Shape::tet;
                        break;
                    case gmsh::hex:
                        n_vertex = 8;						
                        shape = Shape::hex;
                        break;
                    case gmsh::pri:
                        n_vertex = 6;						
                        shape = Shape::pri;
                        break;
                    default:
                        std::cout << "file: " << file_name << std::endl;
                        std::cout << "igs: " << igs << std::endl;
                        std::cout << "invalid n_vertex" << std::endl;
                        std::cout << "igs = " << igs << std::endl;
                        std::cout << "tag = " << tag << std::endl;
                        std::cout << "n_tags = " << n_tags << std::endl;
                        assert(false);
                }

                PointPointers vtx;
                vtx.reserve(n_vertex);

                for (int i=0; i<n_vertex; ++i)
                {
                    in >> tempi;
                    vtx.push_back(&points_[tempi-1]);
                }

                //add_cell(Cell(tag-1, vtx, shape, phys_tag_name[phys_tag]), n_total - n_bface);
                add_cell(Cell(tag-1, vtx, shape, PatchName(phys_tag)), n_total - n_bface);
            }        
        }

        in.close();
    }
}

