#ifndef VECTOR_H
#define VECTOR_H
#include "iterators.hpp"
#include "utils.hpp"
#include <memory>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <bits/c++allocator.h>
#include <iostream>

namespace ft
{
	template < class T, typename Alloc = std::allocator<T> > 
	class vector
	{ 
    public:
        typedef T                                           value_type;
        typedef Alloc                                       allocator_type;
        typedef std::size_t                                 size_type;
        typedef random_access_iterator<T>                   iterator;
        typedef random_access_iterator<const T>             const_iterator;
        typedef ft::reverse_iterator<iterator>              reverse_iterator;
        typedef ft::reverse_iterator<const_iterator>        const_reverse_iterator;
        typedef typename allocator_type::reference          reference;
        typedef typename allocator_type::const_reference    const_reference;
        typedef	typename allocator_type::pointer            pointer;
        typedef typename allocator_type::const_pointer      const_pointer;

    private:
        pointer arr;
        size_type _size;
        Alloc alloc;
       	size_type _capacity;

    public:
        //member functions
        // vector():arr(NULL),_size(0), alloc(Alloc()), _capacity(0) {;}
        explicit vector(const Alloc& alloc = Alloc()):arr(NULL),_size(0), alloc(alloc), _capacity(0) {;}
		
		vector(std::initializer_list<T> init, const allocator_type& alloc = Alloc()):arr(NULL),_size(init.size()), alloc(alloc), _capacity(init.size()) {
			arr = this->alloc.allocate(init.size());
			int i = 0;
			for (auto it = init.begin(); it != init.end(); it++)
				this->alloc.construct(arr + i++, *it);
		}


        template<class InputIt>
        vector(InputIt first, InputIt last, const allocator_type& alloc = Alloc(),
        typename enable_if<!is_integral<InputIt>::value, bool>::type is = 0)
        :arr(NULL),_size(is), alloc(alloc), _capacity(0)
        {
		    _size = ft::distance(first, last);
			_capacity = ft::distance(first, last);
			arr = this->alloc.allocate(_capacity);
			for (int i = 0; first != last; i++)
				this->alloc.construct(arr + i, *first++);
		}

        vector(size_type count, const T& value = T(), const allocator_type& alloc = Alloc())
        :arr(NULL),_size(count), alloc(alloc), _capacity(count) {
			arr = this->alloc.allocate(count);
            for (size_type i = 0; i < _size; i++) {
				this->alloc.construct(arr + i, value);
			}
        }

        vector( const vector& other )
        : _size(other._size), alloc(other.alloc), _capacity(other._capacity){ 
        	arr = this->alloc.allocate(_capacity);
            for (size_type i = 0; i < _size; i++)
				this->alloc.construct(arr + i, *(other.arr + i)); 
        }
        ~vector() {
            clear();
            if (_capacity)
                alloc.deallocate(arr, _capacity);
        }


		value_type& x() {
			if (_size < 1)
				throw std::length_error("vector x: vector is not 1D");
			return (arr[0]);
		}

		value_type& y() {
			if (_size < 2)
				throw std::length_error("vector y: vector is not 2D");
			return (arr[1]);
		}

		value_type& z() {
			if (_size < 3)
				throw std::length_error("vector z: vector is not 3D");
			return (arr[2]);
		}

		vector& operator+=(const vector& other) {
			if (other._size > _capacity)
				reserve(other._size);
			for (size_type i = 0; i < other._size; i++)
				arr[i] += other.arr[i];
			return (*this);
		}
		
		vector& operator-=(const vector& other) {
			if (other._size > _capacity)
				reserve(other._size);
			for (size_type i = 0; i < other._size; i++)
				arr[i] -= other.arr[i];
			return (*this);
		}

		vector& operator*=(const T& val) {
			for (size_type i = 0; i < _size; i++)
				arr[i] *= val;
			return (*this);
		}
		
		vector& operator/=(const T& val) {
			for (size_type i = 0; i < _size; i++)
				arr[i] /= val;
			return (*this);
		}

		friend T dot(const vector& a, const vector& b) {
			if (a._size != b._size)
				throw std::length_error("vector dot: vectors are not the same size");
			T res = 0;
			for (size_type i = 0; i < a._size; i++)
				res += a.arr[i] * b.arr[i];
			return (res);
		}

		vector cross(const vector &other) {
			if ((*this)._size != other._size)
				throw std::length_error("vector cross: vectors (*this)re not the s(*this)me size");
			if ((*this)._size != 3 && (*this)._size != 2)
				throw std::length_error("vector cross: only vectors 2D and 3D vectors supported");
			else if ((*this)._size == 2) {
				return vector({(*this).arr[0] * other.arr[1] - (*this).arr[1] * other.arr[0]});
			}
			else if ((*this)._size == 3) {
				return (vector({(*this).arr[1] * other.arr[2] - (*this).arr[2] * other.arr[1],
							 	(*this).arr[2] * other.arr[0] - (*this).arr[0] * other.arr[2],
								(*this).arr[0] * other.arr[1] - (*this).arr[1] * other.arr[0]}));
			} 	
			return 0;
		}

		friend vector cross(const vector& a, const vector& b) {
			if (a._size != b._size)
				throw std::length_error("vector cross: vectors are not the same size");
			if (a._size != 3 && a._size != 2)
				throw std::length_error("vector cross: only vectors 2D and 3D vectors supported");
			else if (a._size == 2) {
				return vector({a.arr[0] * b.arr[1] - a.arr[1] * b.arr[0]});
			}
			else if (a._size == 3) {
				return (vector({a.arr[1] * b.arr[2] - a.arr[2] * b.arr[1],
							 	a.arr[2] * b.arr[0] - a.arr[0] * b.arr[2],
								a.arr[0] * b.arr[1] - a.arr[1] * b.arr[0]}));
			} 	
			return 0;
		}

		vector& operator=( const vector& other ) {
            clear();
            alloc.deallocate(arr, _capacity);
            alloc = other.alloc; _capacity = other._capacity;  _size = other._size;
            arr = alloc.allocate(_capacity);
            for (size_type i = 0; i < _size; i++)
				alloc.construct(arr + i, *(other.arr + i));
            return (*this);
        }
        
		template <class InputIterator>
  		void assign (InputIterator first, InputIterator last,
		typename ft::enable_if<!ft::is_integral<InputIterator>::value, bool>::type tmp = 0) { // if is_integer<InputIt>::val == 0 , if_is::exist will not exist :)
			tmp += 1; // for the flag
			clear();
			int sz = ft::distance(first, last);
			if (sz > (int)_capacity)
				reserve(sz);
			while (first != last)
				alloc.construct(arr + _size++, *first++);
		}

		void assign (size_type n, const value_type &val) {clear(), resize(n, val);}
        
        allocator_type getallocator() const    {
            return alloc;
        }

        //element access
        reference operator[] (size_type n) {return arr[n];}
        const_reference operator[] (size_type n) const {return arr[n];}
		reference at (size_type n) {
			if (n >= _size) throw std::out_of_range("vector");
			return arr[n];
		}
		const_reference at (size_type n) const {
			if (n >= _size) throw std::out_of_range("vector");
			return arr[n];
		}
		reference front() {return arr[0];}
		const_reference front() const {return arr[0];}
		reference back() {return arr[_size-1];}
		const_reference back() const {return arr[_size-1];}
        pointer data() {return arr;}
        const_pointer data() const {return arr;}
    
        //iterators
		iterator            	begin() 		{return iterator(arr);}
		iterator            	end() 			{return iterator(arr + _size);}
		const_iterator      	begin() const 	{return const_iterator(arr);}
		const_iterator      	end() const 	{return const_iterator(arr + _size);}
		reverse_iterator        rbegin() 		{return reverse_iterator((arr + _size));}
		reverse_iterator        rend() 			{return reverse_iterator(arr);}
		const_reverse_iterator  rbegin() const 	{return const_reverse_iterator((arr + _size));}
		const_reverse_iterator  rend() const 	{return const_reverse_iterator(arr);}

        //_capacity
        bool empty() 				{return !_size;}
        size_type size() const 		{return _size;}
        size_type max_size() 		{return alloc.max_size();}
		size_t capacity() const 	{return _capacity;}

		void reserve (size_type n)
		{
			// if (n > alloc.max_size())
			// 	throw std::length_error("vector::reserve"); //profiling reasons
			if (n > _capacity)
			{
				T * tmp = alloc.allocate(n + 1);

				for (size_type i = 0 ; i < _size && i < n ; i++)
				{
					alloc.construct(tmp + i, arr[i]);
					alloc.destroy(arr + i);
				}
				alloc.deallocate(arr, _capacity);

				_capacity = n;
				arr = tmp;
			}
		}

        //modifiers
        void clear() {
            for (size_type i = 0; i < _size; i++) 
				alloc.destroy(arr + i);
            _size = 0;
        }
	
		void resize (size_type n, value_type val = value_type()) {
			while (_size > n)
				alloc.destroy(arr + --_size);
			if (n > _capacity) {
				int sz = ft::max(_capacity * 2, n);
				pointer tmp = alloc.allocate(sz);
				for (int i = 0; i < (int)_size; i++) {
					alloc.construct(tmp + i, *(arr + i));
					alloc.destroy(arr + i);
				}
				alloc.deallocate(arr,_capacity);
				_capacity = sz;
				arr = tmp;
			}
			while (_size < n)
				alloc.construct(arr + _size++, val);
		}

		iterator insert (iterator position, const value_type &val) {
			size_t gap = ft::distance(begin(), position);
			if (_size + 1 >= _capacity)
				this->reserve(ft::max(_capacity * 2, _size + 1));
			for (size_t i = _size ; i > gap; i--) {
				alloc.construct(arr + i, *(arr + i - 1));
				alloc.destroy(arr + i - 1);
			}
			alloc.destroy(arr + gap);
			alloc.construct(arr + gap, val);
			_size++;
			return begin() + gap;
		}

   	    void insert (iterator pos, size_type n, const value_type& val) {
			if (!n) return;
			size_t dist = pos - begin();
			if (_size + n >= _capacity)
				this->reserve(ft::max(_capacity * 2, _size + n));
			for (size_t i = _size - 1 + n ; i > dist + n - 1; i--) {
				alloc.construct(arr + i, *(arr + i - n));
				alloc.destroy(arr + i - n);
			}
			for (size_t i =  dist ; i < dist + n; i++) {
				alloc.destroy(arr + i);
				alloc.construct(arr + i, val);
			}
			_size += n;
		}

		template <class InputIterator>
  		void insert (iterator position, InputIterator first, InputIterator last,
			typename ft::enable_if<!ft::is_integral<InputIterator>::value, bool>::type tmp = 0) {
			tmp += 1; 
			int gap = distance(begin(), position);
			while (first != last)
				insert(iterator(begin()+gap++), *first), ++first;
  		}
		
		iterator erase(iterator position)
		{
			return (this->erase(position, position + 1));
		}

		iterator erase (iterator first, iterator last) {
			int gap = last - first;
			int beg = first - iterator(arr);

			for (int i = beg; first < last || i < (int)_size; i++) {
				alloc.destroy(&(*first++));
				if (i + gap < (int)_size)
					alloc.construct(arr + i, *(arr + i + gap));
				first++;
			}
			_size -= gap;
			return iterator(arr + beg);
		}			
		
		void push_back (const value_type& val) {
				if (_size == _capacity) {
					pointer _tmp = arr;
					arr = alloc.allocate(_capacity ? _capacity *= 2 : _capacity = 1);
					for (int i = 0 ; i < (int)_size; i++) {
						alloc.construct(arr + i, _tmp[i]);
						alloc.destroy(_tmp + i);
					}
					if (_size)
						alloc.deallocate(_tmp, _size);
				}
				alloc.construct(arr + _size++, val);
		}

        void pop_back() {
   			if (_size) {
				alloc.destroy(arr + _size - 1);
				_size--;
			}
        }

		void swap (vector& x) {
			ft::swap(arr, x.arr);
			ft::swap(alloc, x.alloc);
			ft::swap(_capacity, x._capacity);
			ft::swap(_size, x._size);
		}
            
        friend bool operator==( const vector<T,Alloc>& lhs,
                    const vector<T,Alloc>& rhs ) {
                if (lhs.size() != rhs.size())  return 0;    
                for (size_t i = 0; i < lhs.size(); i++) 
                    if (lhs[i] != rhs[i]) return 0;
                return 1;
        }
        friend bool operator!=( const vector<T,Alloc>& lhs,
                    const vector<T,Alloc>& rhs ) {
                        return !(lhs == rhs);
        }

        friend bool operator<( const vector<T,Alloc>& lhs,
                    const vector<T,Alloc>& rhs ) {
            return lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
        }
        friend bool operator<=( const vector<T,Alloc>& lhs,
                    const vector<T,Alloc>& rhs ) {
            return (lhs == rhs || lhs < rhs);
        }
        friend bool operator>( const vector<T,Alloc>& lhs,
                    const vector<T,Alloc>& rhs ) {
            return lexicographical_compare(rhs.begin(), rhs.end(), lhs.begin(), lhs.end());
        }
        friend bool operator>=( const vector<T,Alloc>& lhs,
                    const vector<T,Alloc>& rhs ) {
            return (lhs == rhs || lhs > rhs);
        }

		friend vector operator *(const T &v, const vector &rhs) {
			vector res(rhs);
			res *= v;
			return res;
		}

		vector operator+(const vector &v) const {
			vector res;
			if (v.size() != _size)
				throw std::invalid_argument("operator + called on vectors of different sizes");
			for (size_t i = 0; i < _size; i++)
				res.push_back(arr[i] + v[i]);
			return res;
		}

		// friend vector operator +(const vector &rhs, const vector &lhs) {
		// 	vector res;
		// 	if (rhs.size() != lhs.size())
		// 		throw std::exception();
		// 	for (size_t i = 0; i < rhs.size(); i++)
		// 		res.push_back(rhs[i] + lhs[i]);
		// 	return res;
		// }

		vector operator-(const vector &v) const {
			vector res(*this);
			res -= v;
			return res;
		}

		vector operator-()  {
			for (size_t i = 0; i < _size; i++)
				arr[i] = -arr[i];
			return *this;
		}
		// friend vector operator -(const vector &rhs, const vector &lhs) {
		// 	vector res;
		// 	if (rhs.size() != lhs.size())
		// 		throw std::exception();
		// 	for (size_t i = 0; i < rhs.size(); i++)
		// 		res.push_back(rhs[i] - lhs[i]);
		// 	return res;
		// }

		friend vector direction(T roll, T pitch, T yaw) {
			vector res(3);
			res[0] = cos(pitch) * cos(yaw);
			res[1] = cos(roll) * sin(yaw) + sin(roll) * sin(pitch) * cos(yaw);
			res[2] = sin(roll) * sin(yaw) - cos(roll) * sin(pitch) * cos(yaw);
			return res;
		}

		friend vector direction(vector eurler_angles) {
			vector res(3);
			res[0] = cos(eurler_angles[1]) * cos(eurler_angles[2]);
			res[1] = cos(eurler_angles[0]) * sin(eurler_angles[2]) + sin(eurler_angles[0]) * sin(eurler_angles[1]) * cos(eurler_angles[2]);
			res[2] = sin(eurler_angles[0]) * sin(eurler_angles[2]) - cos(eurler_angles[0]) * sin(eurler_angles[1]) * cos(eurler_angles[2]);
			return res;
		}

		vector operator*(const T &v) const {
			vector res(*this);
			res *= v;
			return res;
		}

		vector operator*(const T &v)  {
		*this *= v;
		return *this;
		}

		vector operator/(const T &v)  {
		*this /= v;
		return *this;
		}

		friend vector normalize(const vector v)  {
			vector res(v);
			res /= norm(v);
			return res;
		}

		vector &normalize() {
			*this /= norm(*this);
			return *this;
		}

		friend double norm(const vector& a) {
			double res = 0;
			for (size_type i = 0; i < a._size; i++)
				res += fabs(a.arr[i] * a.arr[i]);
			return (sqrt(res));
		}

		friend double norm_1(const vector& a) {
			double res = 0;
			for (size_type i = 0; i < a._size; i++)
				res += fabs(a.arr[i]);
			return (res);
		}

		friend double norm_inf(const vector& a) {
			double res = 0;
			for (size_type i = 0; i < a._size; i++)
				res = max(res, fabs(a.arr[i]));
			return (res);
		}

		friend double angle_cos(const vector& a, const vector& b) {
			return fabs(dot(a, b)) / (norm(a) * norm(b));
		}
	
		friend vector linear_interpolation(const vector &a, const vector &b, T t) {
			return a + (b - a) * t;
		}

	};
    
    template< class T, class Alloc >
    void swap(vector<T,Alloc>& lhs, vector<T,Alloc>& rhs) {
        lhs.swap(rhs);
    }

    template< class T, class Alloc >
	std::ostream &operator<<(std::ostream &out, const vector<T,Alloc> &v) {
		std::cout << "[";
		for (size_t i = 0; i < v.size() - 1; i++)
			out << v[i] << ", ";
		if (v.size())
			out << v[v.size() - 1];
		std::cout << "]";
		return out;
	}


	template<class T>
	ft::vector<T> linear_combination(ft::vector<ft::vector<T>> v, ft::vector<T> a)  {
		ft::vector<T> res(v[0].size(), 0);
		if (v.size() != a.size())
			throw std::invalid_argument("vectors must be of same size");
		for	(int i = 0; i < (int)v.size(); i++) {
			res += v[i] * a[i];
		}
		return res;
	}


}


#endif 