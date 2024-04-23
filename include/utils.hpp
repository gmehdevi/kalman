#ifndef UTILS_H
#define UTILS_H
#include <complex>
#include <sstream>

template < typename... Args >
std::string sstr( Args &&... args )
{
    std::ostringstream sstr;
    ( sstr << std::dec << ... << args );
    return sstr.str();
}

namespace ft
{
    //misc
    template <class T>
    void swap(T &a, T &b) {
        T tmp = a;
        a = b;
        b = tmp;
    }
    
    template <class T>
    T max (T a, T b) {
        return (a > b ? a : b);
    }

    //enable if
	template <class T, T v>
	struct integral_constant {
  		typedef T value_type;
  		typedef integral_constant<T,v> type;
  		static const value_type value = v;  
  		operator value_type() {return v;}
	};  

	typedef integral_constant<bool,false> false_type;
	typedef integral_constant<bool,true> true_type;

	template <typename T> struct is_integral : public false_type {};
	template <typename T> struct is_integral<const T> : public false_type {};
	template <> struct is_integral<bool> : public true_type {};
	template <> struct is_integral<char> : public true_type {};
	template <> struct is_integral<wchar_t> : public true_type {};
	template <> struct is_integral<signed char> : public true_type {};
	template <> struct is_integral<short> : public true_type {};
	template <> struct is_integral<int> : public true_type {};
	template <> struct is_integral<long> : public true_type {};
	template <> struct is_integral<long long> : public true_type {};
	template <> struct is_integral<unsigned char> : public true_type {};
	template <> struct is_integral<unsigned short> : public true_type {};
	template <> struct is_integral<unsigned int> : public true_type {};
	template <> struct is_integral<unsigned long> : public true_type {};
	template <> struct is_integral<unsigned long long> : public true_type {};

    template<bool B, class T = void>
    struct enable_if {};
    
    template<class T>
    struct enable_if<true, T> { typedef T type; };

    //compare

    template< class InputIt1, class InputIt2, class Compare >
    bool lexicographical_compare( InputIt1 first1, InputIt1 last1,
                                    InputIt2 first2, InputIt2 last2,
                                    Compare comp )

    {
        while (first1 != last1 && first2 != last2)
        {
            if      (comp(*first1, *first2))  return true;
            else if (comp(*first2, *first1)) return false;
            ++first1; ++first2;
        }
        return (first1 == last1 && first2 != last2);
    }

    template< class InputIt1, class InputIt2 >
    bool lexicographical_compare( InputIt1 first1, InputIt1 last1,
                                    InputIt2 first2, InputIt2 last2)
    {
        while (first1 != last1 && first2 != last2)
        {
            if      (*first1 < *first2) return true;
            else if (*first2 < *first1) return false;
            ++first1; ++first2;
        }
        return (first1 == last1 && first2 != last2);
    }

    template< class InputIt1, class InputIt2 >
    bool equal( InputIt1 first1, InputIt1 last1,
                InputIt2 first2 );

    template< class InputIt1,
            class InputIt2,
            class BinaryPredicate >
    bool equal( InputIt1 first1,
                InputIt1 last1,
                InputIt2 first2,
                BinaryPredicate p);


    template <bool flag, class IsTrue, class IsFalse>
    struct choose;
    
    template <class IsTrue, class IsFalse>
    struct choose<true, IsTrue, IsFalse> {
    typedef IsTrue type;
    };
    
    template <class IsTrue, class IsFalse>
    struct choose<false, IsTrue, IsFalse> {
    typedef IsFalse type;
    };

    template<typename T>
	int distance(T first, T end) {
		int cnt=0;
		while (first != end) {
			first++;
			cnt++;
		}
		return cnt;
	}

    //this is pretty cool actually even if it's ugly
    template<bool b, typename F, typename T>
    struct isConst {};

    template<typename T, typename F>
    struct isConst<0, F, T> {
        typedef F type;
    };

    template<typename T, typename F>
    struct isConst<1, F, T> {
        typedef T type;
    };
    
	template<typename T = void>
	struct less {
		typedef bool result_type;
		typedef T first_argument_type;
		typedef T second_argument_type;

		bool operator()(const T &lhs, const T &rhs) const {
			return (lhs < rhs);
		}
	};

    template <typename T>
    bool operator> (const std::complex<T> & lhs, const std::complex<T> &rhs) {
        return (lhs.real() > rhs.real()) || (lhs.real() == rhs.real() && lhs.imag() > rhs.imag());
    }

    template <class It>
    void reverse(It first, It last) {
        while ((first != last) && (first != --last)) {
            ft::swap(*first, *last);
            ++first;
        }
    }
}


#endif 