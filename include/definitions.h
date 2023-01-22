#ifndef DEFINITIONS_H
#define	DEFINITIONS_H

#include <vector>
#include "matrix.h"

namespace Liver
{
    using Scalar = double;

    using PointTag = int;
    using FaceTag = std::string;
    using CellTag = int;

    using Coordinate = Vector3;
    using Coordinates = std::vector<Coordinate>;
    using Centroid = Coordinate;
    using Normal = Vector3;

    class Point;
    class Edge;
    class Face;
    class Cell;

    using Points = std::vector<Point>;
    using PointPointer = Point*;
    using PointPointers = std::vector<PointPointer>;

    using Edges = std::vector<Edge>;
    using EdgePointer = Edge*;
    using EdgePointers = std::vector<EdgePointer>;

    using Faces = std::vector<Face>;
    using FacePointer = Face*;
    using FacePointers = std::vector<FacePointer>;

    using Cells = std::vector<Cell>;
    using CellPointer = Cell*;
    using CellPointers = std::vector<CellPointer>;

    using ParentFace = FacePointer;
    using ParentFaces = FacePointers;

    using ParentCell = CellPointer;
    using ParentCells = CellPointers;

    using Neighbor = CellPointer;
    using Neighbors = CellPointers;

    using NumberOfPoints = int;
    using NumberOfCells = int;

    //using PatchName = std::string;
    using FileName = std::string;

    using Area = Scalar;
    using Volume = Scalar;

    const Coordinate origin = Coordinate(0., 0., 0.);
    constexpr int NVAR = 6;
    
    using LSCoef = std::vector<double>;
    using Gradient = Vector3;
    using Gradients = std::array<Gradient, NVAR>;

    using Property = double;
    using Properties = std::array<double, NVAR>;
    //using PropertyIndex = int;

    using Velocity = Vector3;
    
    using Residual = Scalar;

    using TimeStep = double;

    using IsSteady = bool;

    using Time = double;

    using Coef= double;
    using Coefs = std::array<Coef, NVAR>;

    using Left = CellPointer;
    using Right = CellPointer;
}

#endif
