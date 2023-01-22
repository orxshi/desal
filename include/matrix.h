#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <array>
#include <cassert>
#include <iostream>
#include <tuple>

namespace Liver
{
    template<int nrow, int ncol, class T>
        class Matrix
        {
            public:

                Matrix();
                Matrix(T a, T b, T c);
                Matrix(T s);

                Matrix transpose();
                T& operator()(int row, int col);
                T operator()(int row, int col) const;
                T& operator()(int row);
                T operator()(int row) const;
                bool isnan() const;
                template<int nrowo, int ncolo>
                    Matrix<nrow, ncolo, T> operator*(const Matrix<nrowo, ncolo, T>& other) const
                    {
                        // matrix-matrix product

                        assert(nrowo == ncol);

                        Matrix<nrow, ncolo, T> mat;

                        for (int r=0; r<nrow; ++r)
                        {
                            for (int c=0; c<ncolo; ++c)
                            {
                                mat(r, c) = 0.; // take a row

                                for (int j=0; j<ncol; ++j)
                                {
                                    mat(r,c) += (*this)(r,j) * other(j,c);
                                }
                            }
                        }

                        return mat;
                    }
                Matrix operator-();
                Matrix operator*(T s) const;
                Matrix operator/(T s) const;
                Matrix& operator+=(const Matrix<nrow, ncol, T>& M);
                Matrix& operator-=(const Matrix<nrow, ncol, T>& M);
                Matrix& add_diag(T d);
                Matrix& operator+=(T d);
                Matrix& operator*=(T d);
                Matrix& operator/=(T s);
                Matrix operator+(const Matrix<nrow, ncol, T>& M) const;
                Matrix operator-(const Matrix<nrow, ncol, T>& M) const;
                Matrix& operator=(T s);
                bool operator==(const Matrix<nrow, ncol, T>& M) const;
                bool operator!=(const Matrix<nrow, ncol, T>& M) const;
                double len();
                int nelm() const;

            private:

                std::array<T, nrow * ncol> data_;
        };

    using Matrix5 = Matrix<5, 5, double>;
    using Matrix3 = Matrix<3, 3, double>;

    template<int nrow, class T>
        using SquareMatrix = Matrix<nrow, nrow, T>;

    template<int nrow>
        using SquareMatrixD = SquareMatrix<nrow, double>;

    template<int nrow, class T>
        class LowerTriMat: public SquareMatrix<nrow, T>
        {
            public:

            LowerTriMat();
        };

    template<int nrow, class T>
        class UnitLowerTriMat: public LowerTriMat<nrow, T>
        {
            public:

            UnitLowerTriMat();
        };

    template<int nrow, class T>
        class UpperTriMat: public SquareMatrix<nrow, T>
        {
            public:

            UpperTriMat();
        };

    template<int nrow>
        using LowerTriMatD = LowerTriMat<nrow, double>;

    template<int nrow>
        using UnitLowerTriMatD = UnitLowerTriMat<nrow, double>;

    template<int nrow>
        using UpperTriMatD = UpperTriMat<nrow, double>;

    template<int nrow, class T>
        using Vector = Matrix<nrow, 1, T>;

    template<int nrow>
        using VectorD = Vector<nrow, double>;

    template<int nrow>
        using SolVectorD = VectorD<nrow>;

    template<int nrow>
        using RhsVectorD = VectorD<nrow>;

    using Vector5 = Vector<5, double>;
    using Vector3 = Vector<3, double>;
    using Vector3Int = Vector<3, int>;

    template<int nrow, int ncol, class T>
        std::ostream& operator<<(std::ostream& os, const Matrix<nrow, ncol, T>& mat);

    template<int nrow, class T>
        double len(const Vector<nrow, T> a);

    template<int nrow, class T>
        Vector<nrow, T> cross(const Vector<nrow, T>& a, const Vector<nrow, T>& b);

    template<int nrow, class T>
        double dot(const Vector<nrow, T>& a, const Vector<nrow, T>& b);

    template<int nrow, class T>
        Vector<nrow, T> normals(const Vector<nrow, T>& a, const Vector<nrow, T>& b, const Vector<nrow, T>& c);

    template<int nrow, class T>
        Vector<nrow, T> normalize(const Vector<nrow, T>& a);

    Vector5 operator/(const Vector5& f, const Matrix5& A);

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> operator*(T s, const Matrix<nrow, ncol, T>& m);

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> operator/(T s, const Matrix<nrow, ncol, T>& m);

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> abs(const Matrix<nrow, ncol, T>& m);

    template<int nrow, int ncol, class T>
        T max(const Matrix<nrow, ncol, T>& m);
    
    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> unit_matrix();

    template<int nrow>
        VectorD<nrow> Lc(const UnitLowerTriMatD<nrow>&, const VectorD<nrow>&);

    template<int nrow>
        VectorD<nrow> Lc(const SquareMatrixD<nrow>&, const VectorD<nrow>&);

    template<int nrow>
        VectorD<nrow> Uc(const UpperTriMatD<nrow>&, const VectorD<nrow>&);

    template<int nrow>
        VectorD<nrow> Uc(const SquareMatrixD<nrow>&, const VectorD<nrow>&);

    template<int nrow>
        std::tuple<UnitLowerTriMatD<nrow>, UpperTriMatD<nrow>> LU(SquareMatrixD<nrow>&);

    template<int nrow>
        SolVectorD<nrow> LU_solve(const UnitLowerTriMatD<nrow>&, const UpperTriMatD<nrow>&, const RhsVectorD<nrow>&);

    template<int nrow>
        SolVectorD<nrow> LU_solve(const SquareMatrixD<nrow>&, const RhsVectorD<nrow>&);

    template<int nrow>
        SquareMatrixD<nrow> LU_combine(const UnitLowerTriMatD<nrow>&, const UpperTriMatD<nrow>&);
}

#include "matrix.hpp"

#endif
