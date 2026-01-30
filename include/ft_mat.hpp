#ifndef MATRIX_H
#define MATRIX_H
#include "ft_vec.hpp"
#include "utils.hpp"
#include <iostream>
#include <stdexcept>
#include "iterators.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>

namespace ft
{
   //--------------------------------------Iterator--------------------------------------------//
    template <typename T>
        class matrix_iterator : public iterator<bidirectional_iterator_tag, T>
        {
        public:
          typedef typename iterator<bidirectional_iterator_tag, T>::iterator_category             iterator_category;
          typedef typename iterator<bidirectional_iterator_tag, T>::value_type                    value_type;
          typedef typename iterator<bidirectional_iterator_tag, T>::difference_type               difference_type;
          typedef  T*                                            pointer;
          typedef  T&                                            reference;
          typedef  size_t                                          size_type;

        private :
          vector< vector<T> >          *mat;
          unsigned int                  m;
          unsigned int                  n;

        public :
          explicit matrix_iterator(vector< vector<T> >  &mat): mat(&mat.mat),m(0),n(0) {;}
          matrix_iterator(size_type _n, size_type m, vector< vector<T> >  &mat): mat(&mat),m(m),n(_n) {}
          matrix_iterator(const matrix_iterator &cpy):mat(0) {mat = cpy.mat; m = cpy.m; n = cpy.n;}
          matrix_iterator &operator=(const matrix_iterator &cpy) {mat = cpy.mat; m = cpy.m; n = cpy.n; return *this;}

          T &operator  *() {return  (*mat)[n][m];}
          vector<T> *operator ->()  {return &(*mat)[n];}


          matrix_iterator  operator ++() {
            if (m >= (*mat)[n].size() - 1 && n == (*mat).size() - 1) {
                m = (*mat)[n].size();
                n = (*mat).size();
               return (*this);
            }
            if (m < (*mat)[n].size() - 1)
                ++m;
            else if (n < (*mat).size() - 1) {
                m = 0;
                    ++n;
            }
            return (*this);
          }

          matrix_iterator operator --() {
            if (m == (*mat)[0].size() && n == (*mat).size()) {
                m = (*mat)[0].size() - 1; n = (*mat).size() - 1;
                return (*this);
            }
            if (m > 0)
                --m;
            else if (n > 0) {
                --n;
                m = (*mat)[n].size() - 1;
            }

            return (*this);
          }

          matrix_iterator operator ++(int) {
            matrix_iterator out(*this);
            ++(*this);
            return out;
          }

          matrix_iterator operator --(int) {
            matrix_iterator out(*this);
            --(*this);
            return out;
          }

          matrix_iterator operator +=(const int& d) {
            unsigned int i,j;
            i = d / (*mat)[0].size();
            j = d % (*mat)[0].size();
            n += i;
            m += j;
            if (n + i >= (*mat)[0].size()) {
                n = (*mat).size();
                m = (*mat)[0].size();
                return *this;
            }
            return *this;
          }

          matrix_iterator operator -=(const int& d) {
            unsigned int i,j;
            i = d / (*mat)[0].size();
            j = d % (*mat)[0].size();
            n += i;
            m += j;
            if (n - i >= (*mat)[0].size()) {
                n = (*mat).size();
                m = (*mat)[0].size();
                return *this;
            }
            return *this;
          }

          matrix_iterator operator +(const int &d) {
              matrix_iterator out(*this);
              out += d;
              return out;
          }

          difference_type operator+(const matrix_iterator r) {
              return *this + r.n + (*mat)[0].size()  + r.m;
          }

          friend matrix_iterator operator +(const int &d, matrix_iterator &r) {
              matrix_iterator out(r);
              out += d;
              return out;
          }
          matrix_iterator operator -(const int &d) {
              matrix_iterator out(*this);
              out -= d;
              return out;
          }

          friend matrix_iterator operator -(const int &d, matrix_iterator &r) {
              matrix_iterator out(r);
              out -= d;
              return out;
          }

          difference_type operator-(const matrix_iterator &r) {
              return n * (*mat)[0].size() + m - r.n * (*r.mat)[0].size() + r.m;
          }

        friend difference_type operator-(const matrix_iterator &a, const matrix_iterator &b) {
                return a.n * (*a.mat)[0].size() + a.m - b.n * (*b.mat)[0].size() + b.m;
          }

        friend bool operator==(const matrix_iterator &a, const matrix_iterator &b) {return *(a.mat) == *(b.mat) &&  a.n * (*a.mat)[0].size() + a.m == b.n * (*b.mat)[0].size() + b.m;}
        friend bool operator!=(const matrix_iterator &a, const matrix_iterator &b) {return !(a == b);}
        };


   //--------------------------------------matrix--------------------------------------------//
	template <typename T>
    class matrix {

	public:
        typedef vector<vector<T> > Container;
		typedef T value_type;
		typedef typename Container::size_type size_type;
		typedef typename Container::reference reference;
		typedef typename Container::const_reference const_reference;
        typedef typename Container::pointer    pointer;
        typedef matrix_iterator<T>        iterator;
        typedef matrix_iterator<const T>  const_iterator;
        typedef ft::reverse_iterator<iterator>              reverse_iterator;
        typedef ft::reverse_iterator<const_iterator>              const_reverse_iterator;

    protected:
		Container mat;

        size_type m;
        size_type n;
    public:



    matrix(size_type n, size_type m, const value_type &val = value_type()) : mat(n, vector<T>(m, val)), m(m), n(n) {if (n == 0 && m > 0)this->n = 1;}
    matrix(size_type n) : mat(n, vector<T>(n, value_type())) , m(n), n(n) {
        for (size_type i = 0; i < n; i++)
            mat[i][i] = 1;
    }
    matrix(std::initializer_list<std::initializer_list<T> > l) : mat(), m(l.begin()->size()), n(l.size()) {
        for (auto it = l.begin(); it != l.end(); it++)
            if (it->size() != m)
                throw std::invalid_argument("matrix variable row size");
        for (auto it = l.begin(); it != l.end(); it++)
            mat.push_back(vector<T>(it->begin(), it->end()));
    }
    ~matrix() {}
    matrix(const matrix &other) : mat(other.mat), m(other.m), n(other.n) {}
    matrix(const ft::vector<T> &v){
        auto tmp = matrix(1, v.size());
        for (size_type i = 0; i < v.size(); i++)
            tmp[0][i] = v[i];
        *this = tmp;
    }

    matrix &operator%=(const value_type &val) {
        if (val == 0)
            throw std::invalid_argument("matrix modulus by zero");
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                mat[i][j] %= val;
        return *this;
    }
    matrix &operator/=(const value_type &val) {
        if (val == (T)0)
            throw std::invalid_argument("matrix division by zero");
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                mat[i][j] /= val;
        return *this;
    }

    matrix &operator*=(const value_type &val) {
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                mat[i][j] *= val;
        return *this;
    }


    matrix& operator*=(const matrix& other) {
        if (other.n != m) {
            throw std::invalid_argument("Matrix multiplication of " + std::to_string(n) + "x" + std::to_string(m) + " and " + std::to_string(other.n) + "x" + std::to_string(other.m) + " matrices");
        }

        matrix result(n, other.m);  // Result matrix dimensions: n x other.m

        for (size_type i = 0; i < result.n; i++) {
            for (size_type j = 0; j < result.m; j++) {
                result.mat[i][j] = 0;
                for (size_type k = 0; k < m; k++) {
                    result.mat[i][j] += mat[i][k] * other.mat[k][j];
                }
            }
        }

        *this = result;
        return *this;
    }

    matrix &operator-=(const matrix &other) {
        if ((n != other.n) || (m != other.m))
            throw std::invalid_argument("matrix subtraction of " + std::to_string(n) + "x" + std::to_string(m) + " and " + std::to_string(other.n) + "x" + std::to_string(other.m) + " matrix");
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                mat[i][j] -= other.mat[i][j];
        return *this;
    }

    matrix &operator+=(const matrix &other) {
        if ((n != other.n) || (m != other.m))
            throw std::invalid_argument("matrix addition of " + std::to_string(n) + "x" + std::to_string(m) + " and " + std::to_string(other.n) + "x" + std::to_string(other.m) + " matrix");
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                mat[i][j] += other.mat[i][j];
        return *this;
    }

    ft::vector<T> operator*(const ft::vector<T> &v) {
        if (m != v.size() && (m != 4 && v.size() != 3))
            throw std::invalid_argument("matrix multiplication of " + std::to_string(n) + "x" + std::to_string(m) + " and " + std::to_string(v.size()) + "x1 vector");
        ft::vector<T> res(n);
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                res[i] += mat[i][j] * v[j];
        return res;
    }

    matrix &operator=(const matrix &other) {
        mat = other.mat;
        m = other.m;
        n = other.n;
        return *this;
    }

    matrix &operator=(const vector<T> &other) {
        *this = matrix(other);
        return *this;
    }

    reference at (size_type  i, size_type j) {
        if ((n <= i) || (m <= j))
            throw std::out_of_range("matrix");
        return mat[i][j];
    }

    reference at (size_type  i, size_type j) const {
        if ((n <= i) || (m <= j))
            throw std::out_of_range("matrix");
        return mat[i][j];
    }

    reference operator[] (size_type  i) {
        if (i >= n)
            throw std::out_of_range("matrix");
        return mat[i];
    }

    const_reference operator[] (size_type i) const  {
        if (i >= n)
            throw std::out_of_range("matrix");
        return mat[i];
    }

    pointer data() {return mat.data();}
    const pointer data() const {return mat.data();}

    size_type rows() const {return n;}
    size_type cols() const {return m;}
    void set(size_type row, size_type col, const matrix &block) {
        if (row + block.n > n || col + block.m > m)
            throw std::out_of_range("matrix set out of range");
        for (size_type i = 0; i < block.n; ++i)
            for (size_type j = 0; j < block.m; ++j)
                mat[row + i][col + j] = block.mat[i][j];
    }
    void set_block(size_type row, size_type col, const matrix &block) {
        set(row, col, block);
    }
    void set_diag(const ft::vector<T> &d, size_type offset = 0) {
        if (offset + d.size() > n || offset + d.size() > m)
            throw std::out_of_range("matrix set_diag out of range");
        for (size_type i = 0; i < d.size(); ++i)
            mat[offset + i][offset + i] = d[i];
    }
    iterator begin() {return iterator(0, 0, mat);}
    iterator end()   {return iterator(n, m, mat);}
    reverse_iterator rbegin() {return   end();}
    reverse_iterator rend() {return begin();}
    const_iterator begin() const {return const_iterator(0, 0, mat);}
    const_iterator end() const   {return const_iterator(n, m, mat);}
    const_reverse_iterator rbegin() const {return   end();}
    const_reverse_iterator rend() const {return begin();}

    T &operator()(size_type i, size_type j) {
        if ((n <= i) || (m <= j))
            throw std::out_of_range("matrix");
        return mat[i][j];
    }

    T max() const {
        T max = mat[0][0];
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                if (mat[i][j] > max)
                    max = mat[i][j];
        return max;
    }

    T trace() {
        if (n != m)
            throw std::invalid_argument("matrix must be square");
        T sum = 0;
        for (size_type i = 0; i < n; i++)
            sum += mat[i][i];
        return sum;
    }

    matrix operator+(const matrix &r) {
        matrix tmp(*this);
        tmp += r;
        return tmp;
    }

    matrix transposed() const {
        matrix tmp(m, n);
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                tmp.mat[j][i] = mat[i][j];
        return tmp;
    }

    matrix c_transposed() {
        matrix tmp(m, n);
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                tmp.mat[j][i] = std::conj(mat[i][j]);
        return tmp;
    }

    size_type i_pivot(const vector<T> &v) {
        size_type i = 0;
        for (i = 0; i < v.size(); i++)
            if (v[i] != (T)0)
                break;
        return i;
    }
    matrix row_echelon() {
        matrix tmp = *this;
        for (size_type i = 0; i < tmp.n; i++) {
            size_type pivot = i_pivot(tmp.mat[i]);
            if (pivot == tmp.m)
                continue;
            if (pivot != i)
                std::swap(tmp.mat[i], tmp.mat[pivot]);
            for (size_type j = i + 1; j < tmp.n; j++) {
                T coef = (std::fabs(tmp.mat[i][i]) != static_cast<T>(0)) ? tmp.mat[j][i] / tmp.mat[i][i] : static_cast<T>(1);
                for (size_type k = i; k < tmp.m; k++)
                    tmp.mat[j][k] -= coef * tmp.mat[i][k];
            }
        }
        return tmp;
    }

    matrix r_row_echelon() {
        matrix tmp = *this;
        int lead = 0;
        while (lead < n) {
            for (int r = 0; r < n; r++) {
                T div = (tmp[lead][lead] != static_cast<T>(0)) ? tmp[lead][lead] : static_cast<T>(1);
                T mult = (tmp[lead][lead] != static_cast<T>(0)) ? tmp[r][lead] / tmp[lead][lead] : static_cast<T>(0);
                for (int c = 0; c < m; c++) {
                    if (r == lead)
                        tmp[r][c] /= div;
                    else
                        tmp[r][c] -= tmp[lead][c] * mult;
                }
            }
            lead++;
        }
        for (int i = 0; i < n; i++)
            ft::reverse(tmp[i].begin(), tmp[i].end());
        ft::reverse(tmp.mat.begin(), tmp.mat.end());
        return tmp;
    }

    T det() {
        if (n != m)
            throw std::invalid_argument("matrix must be square");
        matrix tmp = row_echelon();
        T d = 1;
        for (size_type i = 0; i < tmp.n; i++)
            d *= tmp.mat[i][i];
        return d;
    }

    matrix adj() {
        if (n != m)
            throw std::invalid_argument("matrix must be square");
        matrix tmp(n, n);
        for (size_type i = 0; i < n; i++) {
            for (size_type j = 0; j < n; j++) {
                matrix minor(n - 1, n - 1);
                for (size_type k = 0; k < n; k++) {
                    if (k != i) {
                        for (size_type l = 0; l < n; l++) {
                            if (l != j) {
                                minor.mat[k > i ? k - 1 : k][l > j ? l - 1 : l] = mat[k][l];
                            }
                        }
                    }
                }
                tmp.mat[i][j] = minor.det() * static_cast<double>((i + j) % 2 ? -1 : 1);
            }
        }
        return tmp;
    }


    matrix<T> inv() {
        if (n != m)
            throw std::invalid_argument("Matrix must be square");

        matrix<T> augmented(n, 2 * n);

        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < n; j++)
                augmented[i][j] = mat[i][j];

        for (size_type i = 0; i < n; i++)
            augmented[i][i + n] = 1;


        for (size_type i = 0; i < n; i++) {

            size_type pivot = i;
            for (size_type j = i + 1; j < n; j++)
                if (std::fabs(augmented[j][i]) > std::fabs(augmented[pivot][i]))
                    pivot = j;

            if (pivot != i)
                std::swap(augmented[i], augmented[pivot]);

            T pivotValue = augmented[i][i];

            if (pivotValue == 0)
                throw std::invalid_argument("Matrix is not invertible");

            for (size_type j = 0; j < 2 * n; j++)
                augmented[i][j] /= pivotValue;

            for (size_type j = 0; j < n; j++) {
                if (j != i) {
                    T factor = augmented[j][i];
                    for (size_type k = 0; k < 2 * n; k++)
                        augmented[j][k] -= factor * augmented[i][k];
                }
            }
        }

        matrix<T> inverse(n, n);
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < n; j++)
                inverse[i][j] = augmented[i][j + n];

        return inverse;
    }

    T sum() {
        T sum = 0;
        for (size_type i = 0; i < n; i++)
            for (size_type j = 0; j < m; j++)
                sum += mat[i][j];
        return sum;
    }

    friend T frob_norm(const matrix &m) {
        T sum = 0;
        for (size_type i = 0; i < m.n; i++)
            for (size_type j = 0; j < m.m; j++)
                sum += m.mat[i][j] * m.mat[i][j];
        return sqrt(sum);
    }

    friend T norm(const matrix &m) {
        T sum = 0;
        for (size_type i = 0; i < m.n; i++)
            for (size_type j = 0; j < m.m; j++)
                sum += m.mat[i][j] * m.mat[i][j];
        return sqrt(sum);
    }

    vector<T> to_vec() const {
        if (n != 1 && m != 1)
            throw std::invalid_argument("matrix must have either 1 row or 1 column but is " + std::to_string(n) + "x" + std::to_string(m));
        matrix cp = *this;
        if (m == 1)
            cp = transposed();
        vector<T> tmp;
        for (size_type i = 0; i < n; i++)
            tmp.push_back(mat[i][0]);
        return tmp;
    }


    size_type rank() {
        matrix tmp = r_row_echelon();
        size_type r = 0;
        for (size_type i = 0; i < tmp.n; i++)
            if (tmp.mat[i] != vector<T>(tmp.m, 0))
                r++;
        return r;
    }

    friend matrix operator +(const matrix &lhs, const matrix &rhs) {
        matrix tmp(rhs);
        tmp += lhs;
        return tmp;
    }


    matrix operator-(const matrix &r) {
        matrix tmp(*this);
        tmp -= r;
        return tmp;
    }

    friend matrix operator -(const matrix &lhs, const matrix &rhs) {
        matrix tmp = rhs;
        tmp -= lhs;
        return tmp;
    }


    matrix operator*(const matrix &r) {
        matrix tmp(*this);
        tmp *= r;
        return tmp;
    }

    friend matrix operator *(const matrix &lhs, const matrix &rhs) {
        matrix tmp = lhs;
        tmp *= rhs;
        return tmp;
    }

    matrix operator*(const value_type &val) {
        matrix tmp(*this);
        tmp *= val;
        return tmp;
    }

    friend matrix operator *(const matrix &lhs, const value_type &val) {
        matrix tmp = lhs * val;
        return tmp;
    }

    friend matrix operator *(const value_type &val, const matrix &rhs) {
        matrix tmp = rhs;
        tmp *= val;
        return tmp;
    }

    matrix &operator/(const value_type &val) {
        *this /= val;
        return *this;
    }

    friend matrix operator /(const matrix &lhs, const value_type &val) {
        matrix tmp = lhs / val;
        return tmp;
    }

    friend matrix operator /(const value_type &val, const matrix &rhs) {
        matrix tmp = rhs / val;
        return tmp;
    }

    matrix &operator%(const value_type &val) {
        *this %= val;
        return *this;
    }

    friend matrix operator %(const matrix &lhs, const value_type &val) {
        matrix tmp = lhs % val;
        return tmp;
    }

    Container to_vecs() {
        return mat;
    }

    vector<T> to_buffer() {
        vector<T> tmp(mat.begin(), mat.end());
        return tmp;
    }


    matrix row(size_type i) {
        if (i >= n)
            throw std::out_of_range("matrix");
        return this->mat[i];;
    }

    matrix col(size_type i) {
        if (i >= m)
            throw std::out_of_range("matrix");
        matrix tmp(n, 1);
        for (size_type j = 0; j < n; j++)
            tmp.mat[j][0] = mat[j][i];
        return tmp;
    }

    friend std::ostream &operator<<(std::ostream &os, matrix mat) {
        size_type width = 0;
        for (auto it = mat.begin(); it != mat.end(); it++)
            width = std::max(sstr(*it).length(), width);
        for (size_type i = 0; i < mat.n; i++) {
            os << "[ " << std::setw(width);
            for (size_type j = 0; j < mat.m - 1; j++) {
                os  << std::setw(width) << mat[i][j]  << " ,";
            }
            if (mat.m > 0) os << std::setw(width) << mat[i][mat.m - 1] ;
             os << "]";
            if (i < mat.n - 1) os << std::endl;
        }
        return os;
    }

        friend matrix linear_interpolation(const matrix &a, const matrix &b, T t) {
            return a + (b - a) * t;
        }

        friend bool operator==(const matrix& lhs, const matrix& rhs) {
                return (lhs.mat == rhs.mat);
        }

        friend bool operator!=( const matrix& lhs,
                        const matrix& rhs ) {
                            return !(lhs.mat == rhs.mat);
        }

        matrix<double> &rot_z(double theta) {
            *this *=  matrix<double> ({
                {cos(theta), sin(theta), 0, 0},
                {-sin(theta), cos(theta), 0, 0},
                {        0,          0, 1, 0},
                {        0,          0, 0, 1}
                 });
            return *this;
        }

        matrix<double> &rot_y(double theta) {
            *this *=  (matrix<double>) {
                {cos(theta), 0, -sin(theta), 0},
                {        0, 1,           0, 0},
                {sin(theta), 0,  cos(theta), 0},
                {        0, 0,           0, 1}
            };
            return *this;
        }

        matrix<double> &rot_x(double theta) {
            *this *=  (matrix<double>) {
                {1, 0, 0, 0},
                {0, cos(theta), sin(theta), 0},
                {0, -sin(theta), cos(theta), 0},
                {0, 0, 0, 1}
            };
            return *this;
        }


        matrix<double> friend sqrtm(const matrix<double> &m) {
            if (m.n != m.m)
                throw std::invalid_argument("matrix must be square");
            matrix<double> tmp(m.n, m.m);
            for (size_type i = 0; i < m.n; i++)
                for (size_type j = 0; j < m.m; j++)
                    tmp[i][j] = sqrt(m[i][j]);
            return tmp;
        }


    };

    matrix<double> eulerToRotation(double roll, double pitch, double yaw) {
        return matrix<double>({
                {cos(yaw) * cos(pitch), cos(yaw) * sin(pitch) * sin(roll) - sin(yaw) * cos(roll), cos(yaw) * sin(pitch) * cos(roll) + sin(yaw) * sin(roll)},
                {sin(yaw) * cos(pitch), sin(yaw) * sin(pitch) * sin(roll) + cos(yaw) * cos(roll), sin(yaw) * sin(pitch) * cos(roll) - cos(yaw) * sin(roll)},
                {-sin(pitch)          , cos(pitch) * sin(roll)                                  , cos(pitch) * cos(roll)},
            });
    }

    matrix<double> eulerToRotation(const vector<double>& angles) {
        if (angles.size() != 3)
            throw std::invalid_argument("Euler angles must be a vector of size 3");
        return eulerToRotation(angles[0], angles[1], angles[2]);
    }

    matrix<double> rotate(matrix<double> m, double theta, vector<double> axis) {
        theta = -theta;
        axis = axis / norm(axis);
        return matrix<double>({
                {cos(theta) + axis[0] * axis[0] * (1 - cos(theta)), axis[0] * axis[1] * (1 - cos(theta)) + axis[2] * sin(theta), axis[0] * axis[2] * (1 - cos(theta)) - axis[1] * sin(theta), 0},
                {axis[1] * axis[0] * (1 - cos(theta)) - axis[2] * sin(theta), cos(theta) + axis[1] * axis[1] * (1 - cos(theta)), axis[1] * axis[2] * (1 - cos(theta)) + axis[0] * sin(theta), 0},
                {axis[2] * axis[0] * (1 - cos(theta)) + axis[1] * sin(theta), axis[2] * axis[1] * (1 - cos(theta)) - axis[0] * sin(theta), cos(theta) + axis[2] * axis[2] * (1 - cos(theta)), 0},
                {0, 0, 0, 1}
            }) * m;
    }

    matrix<double> rotate(double theta, vector<double> axis) {
        theta = -theta;
        axis = axis / norm(axis);
        return matrix<double>({
                {cos(theta) + axis[0] * axis[0] * (1 - cos(theta)), axis[0] * axis[1] * (1 - cos(theta)) + axis[2] * sin(theta), axis[0] * axis[2] * (1 - cos(theta)) - axis[1] * sin(theta), 0},
                {axis[1] * axis[0] * (1 - cos(theta)) - axis[2] * sin(theta), cos(theta) + axis[1] * axis[1] * (1 - cos(theta)), axis[1] * axis[2] * (1 - cos(theta)) + axis[0] * sin(theta), 0},
                {axis[2] * axis[0] * (1 - cos(theta)) + axis[1] * sin(theta), axis[2] * axis[1] * (1 - cos(theta)) - axis[0] * sin(theta), cos(theta) + axis[2] * axis[2] * (1 - cos(theta)), 0},
                {0, 0, 0, 1}
            });
    }

    matrix<double> perspective(double fov, double ratio, double near, double far) {
        double f = 1 / tanf(fov / 2);
        return matrix<double>({
            {f / ratio, 0, 0, 0},
            {0, f, 0, 0},
            {0, 0, (far + near) / (near - far), -1},
            {0, 0, 2 * far * near / (near - far), 0}
        });
    }


    matrix<double> translation(vector<double> v) {
        return matrix<double>({
        {1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{v[0], v[1], v[2], 1}
        });
    }


    matrix<double> scale(vector<double> v) {
        return matrix<double>({
        {v[0], 0, 0, 0},
		{0, v[1], 0, 0},
		{0, 0, v[2], 0},
		{0, 0, 0, 1}
        });
    }

    matrix<double> scale(double s) {
        return matrix<double>({
        {s, 0, 0, 0},
		{0, s, 0, 0},
		{0, 0, s, 0},
		{0, 0, 0, 1}
        });
    }

    matrix<double> camera(vector<double> const &eye, vector<double> const &center, vector<double> const &up) {
        vector<double> f = normalize(center - eye);
        vector<double> u = normalize(up);
        vector<double> s = normalize(cross(f, u));
        u = cross(s, f);

       return matrix<double> ({
            {s[0], u[0], -f[0], 0},
            {s[1], u[1], -f[1], 0},
            {s[2], u[2], -f[2], 0},
            {-dot(s, eye), -dot(u, eye), dot(f, eye), 1.0f}
        });
    }

    



}
#endif
