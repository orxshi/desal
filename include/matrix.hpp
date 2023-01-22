#ifndef MATRIX_HPP
#define	MATRIX_HPP

#include "numeric_constants.h"

namespace Liver
{
    template<int nrow, class T>
        LowerTriMat<nrow, T>::LowerTriMat()
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=r; c<nrow; ++c)
                {
                    (*this)(r,c) = 0.;
                }
            }
        };

    template<int nrow, class T>
        UnitLowerTriMat<nrow, T>::UnitLowerTriMat()
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=r; c<nrow; ++c)
                {
                    (*this)(r,c) = 0.;
                }
            }

            for (int r=0; r<nrow; ++r)
            {
                (*this)(r,r) = 1.;
            }
        };

    template<int nrow, class T>
        UpperTriMat<nrow, T>::UpperTriMat()
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<r; ++c)
                {
                    (*this)(r,c) = 0.;
                }
            }
        };

    template<int nrow, int ncol, class T>
        int Matrix<nrow, ncol, T>::nelm() const
        {
            return nrow * ncol;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>::Matrix()
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    (*this)(r,c) = 0.;
                }
            }
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>::Matrix(T a, T b, T c)
        {
            assert(nrow == 3);
            assert(ncol == 1);

            (*this)(0) = a;
            (*this)(1) = b;
            (*this)(2) = c;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>::Matrix(T s)
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    (*this)(r,c) = s;
                }
            }
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> Matrix<nrow, ncol, T>::transpose()
        {
            auto m = *this;

            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    if (c == r) {
                        continue;
                    }

                    m(r,c) = (*this)(c,r);
                }
            }

            return m;
        }



    template<int nrow, int ncol, class T>
        T& Matrix<nrow, ncol, T>::operator()(int row, int col)
        {
            return data_[row * ncol + col]; 
        }

    template<int nrow, int ncol, class T>
        T Matrix<nrow, ncol, T>::operator()(int row, int col) const
        {
            return data_[row * ncol + col]; 
        }

    template<int nrow, int ncol, class T>
        T& Matrix<nrow, ncol, T>::operator()(int row)
        {
            assert(ncol == 1);
            int col = 0;
            return (*this)(row, col);
        }

    template<int nrow, int ncol, class T>
        T Matrix<nrow, ncol, T>::operator()(int row) const
        {
            assert(ncol == 1);
            int col = 0;
            return (*this)(row, col);
        }

    template<int nrow, int ncol, class T>
        bool Matrix<nrow, ncol, T>::isnan() const
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    if (std::isnan((*this)(r,c)))
                    {
                        return true;
                    }
                }
            }

            return false;
        }


    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> Matrix<nrow, ncol, T>::operator-()
        {
            Matrix<nrow, ncol, T> a;
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //data_[r * ncol + c] *= -1.;
                    //(*this)(r,c) *= -1.;
                    a(r,c) = (*this)(r,c) * -1.;
                }
            }

            return a;
        }

    /*Matrix operator*(const Vector<ncol>& v) const
      {
    // matrix-diag matrix product
    // v is ncol x ncol diagonal square matrix

    Matrix<nrow, ncol> mat;

    for (int r=0; r<nrow; ++r)
    {
    for (int c=0; c<ncol; ++c)
    {
    //mat(r, c) = data_(r * ncol + c) * v(c);
    mat(r,c) = (*this)(r,c) * v(c);
    }
    }

    return mat;
    }*/

    //Vector<nrow> operator*(const Vector<ncol>& v) const
    //{
    //    // matrix-vector product

    //    Vector<nrow> vec;
    //    vec = 0.;

    //    for (int r=0; r<nrow; ++r)
    //    {
    //        for (int c=0; c<ncol; ++c)
    //        {
    //            //vec[r] += data_(r * ncol + c) * v(c);
    //            vec[r] += (*this)(r,c) * v[c];
    //        }
    //    }

    //    return vec;
    //}

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> Matrix<nrow, ncol, T>::operator*(T s) const
        {
            Matrix<nrow, ncol, T> mat;

            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //mat(r, c) = s * data_(r * ncol + c);
                    mat(r,c) = s * (*this)(r,c);
                }
            }

            return mat;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> Matrix<nrow, ncol, T>::operator/(T s) const
        {
            Matrix<nrow, ncol, T> mat;

            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //mat(r, c) = s * data_(r * ncol + c);
                    mat(r,c) = (*this)(r,c) / s;
                }
            }

            return mat;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>& Matrix<nrow, ncol, T>::operator+=(const Matrix<nrow, ncol, T>& M)
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //data_[r * ncol + c] += M(r, c);
                    (*this)(r,c) += M(r, c);
                }
            }

            return *this;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>& Matrix<nrow, ncol, T>::operator/=(T s)
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    (*this)(r,c) /= s;
                }
            }

            return *this;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>& Matrix<nrow, ncol, T>::operator-=(const Matrix<nrow, ncol, T>& M)
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    (*this)(r,c) -= M(r, c);
                }
            }

            return *this;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>& Matrix<nrow, ncol, T>::add_diag(T d)
        {
            for (int r=0; r<nrow; ++r)
            {
                (*this)(r,r) += d;
            }

            return *this;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>& Matrix<nrow, ncol, T>::operator+=(T d)
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //data_[r * ncol + c] += d;
                    (*this)(r,c) += d;
                }
            }

            return *this;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>& Matrix<nrow, ncol, T>::operator*=(T d)
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    (*this)(r,c) *= d;
                }
            }

            return *this;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> Matrix<nrow, ncol, T>::operator+(const Matrix<nrow, ncol, T>& M) const
        {
            Matrix<nrow, ncol, T> mat;

            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //mat(r, c) = data_(r * ncol + c) + M(r, c);
                    mat(r,c) = (*this)(r,c) + M(r,c);
                }
            }

            return mat;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> Matrix<nrow, ncol, T>::operator-(const Matrix<nrow, ncol, T>& M) const
        {
            Matrix<nrow, ncol, T> mat;

            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //mat(r, c) = data_(r * ncol + c) - M(r, c);
                    mat(r,c) = (*this)(r,c) - M(r,c);
                }
            }

            return mat;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T>& Matrix<nrow, ncol, T>::operator=(T s)
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    //data_[r * ncol + c] = s;
                    (*this)(r,c) = s;
                }
            }

            return *this;
        }

    template<int nrow, int ncol, class T>
        bool Matrix<nrow, ncol, T>::operator==(const Matrix<nrow, ncol, T>& M) const
        {
            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    if ((*this)(r,c) != M(r,c))
                    {
                        return false;
                    }
                }
            }

            return true;
        }

    template<int nrow, int ncol, class T>
        bool Matrix<nrow, ncol, T>::operator!=(const Matrix<nrow, ncol, T>& M) const
        {
            return !(M == *this);
        }

    template<int nrow, int ncol, class T>
        double Matrix<nrow, ncol, T>::len()
        {
            return std::sqrt(std::pow(data_[0], 2.) + std::pow(data_[1], 2.) + std::pow(data_[2], 2.));
        }

    template<int nrow, int ncol, class T>
        std::ostream& operator<<(std::ostream& os, const Matrix<nrow, ncol, T>& mat)
        {
            for (int i=0; i<nrow; ++i)
            {
                os << "MAT(" << i << ",:): ";
                for (int j=0; j<ncol; ++j)
                {
                    os << mat(i,j) << ", "; 
                }
                os << "\n"; 
            }

            return os;
        }

    template<int nrow, class T>
        double len(const Vector<nrow, T> a)
        {
            return std::sqrt(std::pow(a(0), 2.) + std::pow(a(1), 2.) + std::pow(a(2), 2.));
        }

    template<int nrow, class T>
        Vector<nrow, T> cross(const Vector<nrow, T>& a, const Vector<nrow, T>& b)
        {
            Vector<nrow, T> c;

            c(0) = a(1) * b(2) - a(2) * b(1);
            c(1) = a(2) * b(0) - a(0) * b(2);
            c(2) = a(0) * b(1) - a(1) * b(0);

            return c;
        }

    template<int nrow, class T>
        double dot(const Vector<nrow, T>& a, const Vector<nrow, T>& b)
        {
            auto c = a(0) * b(0) + a(1) * b(1) + a(2) * b(2);

            return c;
        }

    template<int nrow, class T>
        double angle(const Vector<nrow, T>& a, const Vector<nrow, T>& b)
        {
            return std::acos(dot(a, b) / (len(a) * len(b))); // radian
        }

    template<int nrow, class T>
        Vector<nrow, T> normalize(const Vector<nrow, T>& a)
        {
            return a / len(a);
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> operator*(T s, const Matrix<nrow, ncol, T>& m)
        {
            return m * s;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> operator/(T s, const Matrix<nrow, ncol, T>& m)
        {
            return m / s;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> abs(const Matrix<nrow, ncol, T>& m)
        {
            Matrix<nrow, ncol, T> mm;

            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    mm(r,c) = std::abs(m(r,c));
                }
            }

            return mm;
        }

    template<int nrow, int ncol, class T>
        T max(const Matrix<nrow, ncol, T>& m)
        {
            T mx = BIG_NEG_NUM;

            for (int r=0; r<nrow; ++r)
            {
                for (int c=0; c<ncol; ++c)
                {
                    mx = std::max(mx, m(r,c));
                }
            }

            return mx;
        }

    template<int nrow, int ncol, class T>
        Matrix<nrow, ncol, T> unit_matrix()
        {
            assert(nrow == ncol);

            Matrix<nrow, ncol, T> m = 0.;

            for (int r=0; r<nrow; ++r)
            {
                m(r,r) = 1.;
            }

            return m;
        }

    template<int nrow>
        VectorD<nrow> Lc(const UnitLowerTriMatD<nrow>& L, const VectorD<nrow>& b)
        {
            // Solve Lc = b.

            VectorD<nrow> c;

            c(0) = b(0);

            for (int i=1; i<nrow; ++i)
            {
                c(i) = b(i);

                for (int j=0; j<i; ++j)
                {
                    c(i) -= L(i,j) * c(j);
                }
            }

            return c;
        }

    template<int nrow>
        VectorD<nrow> Lc(const SquareMatrixD<nrow>& L, const VectorD<nrow>& b)
        {
            // Solve Lc = b.

            VectorD<nrow> c;

            c(0) = b(0);

            for (int i=1; i<nrow; ++i)
            {
                c(i) = b(i);

                for (int j=0; j<i; ++j)
                {
                    c(i) -= L(i,j) * c(j);
                }
            }

            return c;
        }

    template<int nrow>
        VectorD<nrow> Uc(const UpperTriMatD<nrow>& U, const VectorD<nrow>& b)
        {
            // Solve Uc = b.

            VectorD<nrow> c;

            int last_row = nrow - 1;

            c(last_row) = b(last_row) / U(last_row,last_row);

            for (int i=last_row-1; i>=0; --i)
            {
                c(i) = b(i);

                for (int j=nrow; j>i; --j)
                {
                    c(i) -= U(i,j) * c(j);
                }

                c(i) /= U(i,i);
            }

            return c;
        }

    template<int nrow>
        VectorD<nrow> Uc(const SquareMatrixD<nrow>& U, const VectorD<nrow>& b)
        {
            // Solve Uc = b.

            VectorD<nrow> c;

            int last_row = nrow - 1;

            c(last_row) = b(last_row) / U(last_row,last_row);

            for (int i=last_row-1; i>=0; --i)
            {
                c(i) = b(i);

                for (int j=last_row; j>i; --j)
                {
                    c(i) -= U(i,j) * c(j);
                }

                c(i) /= U(i,i);
            }

            return c;
        }

    template<int nrow>
        std::tuple<UnitLowerTriMatD<nrow>, UpperTriMatD<nrow>> LU(SquareMatrixD<nrow>& A)
        {
            // Decompose A into unit lower triangular matrix and upper triangular matrix.

            UnitLowerTriMatD<nrow> L;
            UpperTriMatD<nrow> U;

            for (int j=0; j<nrow; ++j)
            {
                U(0,j) = A(0,j);
            }

            for (int i=1; i<nrow; ++i)
            {
                L(i,0) = A(i,0) / U(0,0);
            }

            for (int i=1; i<nrow-1; ++i)
            {
                for (int j=i; j<nrow; ++j)
                {
                    U(i,j) = A(i,j);
                    //assert(U(j,j) != 0.);

                    for (int k=0; k<i; ++k)
                    {
                        U(i,j) -= L(i,k) * U(k,j);

                        //if (U(i,j) == 0.)
                        //{
                            //std::cout << "L(i,k): " << L(i,k) << std::endl;
                            //std::cout << "U(k,j): " << U(k,j) << std::endl;
                            //std::cout << "U(i,j): " << U(i,j) << std::endl;
                            //std::cout << "A(i,j): " << A(i,j) << std::endl;
                            //std::cout << "pro: " << L(i,k) * U(k,j) << std::endl;
                        //}
                    }
                    if (j == i)
                    {
                        if (U(j,j) == 0.)
                        {
                            std::cout << "j: " << i << std::endl;
                            for (int k=0; k<i; ++k)
                            {
                                //U(i,j) -= L(i,k) * U(k,j);

                                //if (U(i,j) == 0.)
                                //{
                                std::cout << "L(i,k): " << L(i,k) << std::endl;
                                std::cout << "U(k,j): " << U(k,j) << std::endl;
                                std::cout << "U(i,j): " << U(i,j) << std::endl;
                                std::cout << "A(i,j): " << A(i,j) << std::endl;
                                std::cout << "pro: " << L(i,k) * U(k,j) << std::endl;
                                //}
                            }

                        }
                        assert(U(j,j) != 0.);
                    }
                }

                for (int k=i+1; k<nrow; ++k)
                {
                    L(k,i) = A(k,i);

                    for (int j=0; j<i; ++j)
                    {
                        L(k,i) -= L(k,j) * U(j,i);
                    }

                    if (U(i,i) == 0.)
                    {
                        std::cout << "i: " << i << std::endl;

                    }
                    assert(U(i,i) != 0.);
                    L(k,i) /= U(i,i);
                }
            }

            U(nrow-1,nrow-1) = A(nrow-1,nrow-1);

            for (int i=0; i<nrow-1; ++i)
            {
                U(nrow-1,nrow-1) -= L(nrow-1,i) * U(i,nrow-1);
            }

            return std::make_tuple(L, U);
        }

    template<int nrow>
        SolVectorD<nrow> LU_solve(const UnitLowerTriMatD<nrow>& L, const UpperTriMatD<nrow>& U, const RhsVectorD<nrow>& b)
        {
            VectorD<nrow> c = Lc(L, b);
            VectorD<nrow> x = Uc(U, c);

            return x;
        }

    template<int nrow>
        SolVectorD<nrow> LU_solve(const SquareMatrixD<nrow>& A, const RhsVectorD<nrow>& b)
        {
            VectorD<nrow> c = Lc(A, b);
            VectorD<nrow> x = Uc(A, c);

            return x;
        }

    template<int nrow>
        SquareMatrixD<nrow> LU_combine(const UnitLowerTriMatD<nrow>& L, const UpperTriMatD<nrow>& U)
        {
            SquareMatrixD<nrow> A;

            for (int i=0; i<nrow; ++i)
            {
                for (int j=0; j<i; ++j)
                {
                    A(i,j) = L(i,j);
                }

                for (int j=i; j<nrow; ++j)
                {
                    A(i,j) = U(i,j);
                }
            }

            return A;
        }
}

#endif
