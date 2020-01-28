#ifndef QMC_DIGITALSEQ_B2G_HPP
#define QMC_DIGITALSEQ_B2G_HPP

////
// (C) Dirk Nuyens, KU Leuven, 2014,2015,2016,2017,2018,...

#include <cstdint>
#include <limits>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <istream>
#include <algorithm>
#include <iterator>
#include <functional>
#include <cassert>

namespace qmc {

    /** Count the number of trailing zero bits.
      *
      * @see Count the consecutive zero bits (trailing) on the right in parallel
      * http://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel
      */
    inline unsigned count_trailing_zero_bits(std::uint32_t v)
    {
        unsigned c = 32;
        v &= -signed(v);
        if (v) c--;
        if (v & 0x0000FFFF) c -= 16;
        if (v & 0x00FF00FF) c -= 8;
        if (v & 0x0F0F0F0F) c -= 4;
        if (v & 0x33333333) c -= 2;
        if (v & 0x55555555) c -= 1;
        return c;
    }

#ifndef QMC_BITREVERSE // want to share this with latticeseq_b2.hpp without having a common include
#define QMC_BITREVERSE
    /** Reverse the bits in a std::uint32_t.
      *
      * @see Bit reverse code from Stanford Bit hacks page:
      * https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
      */
    inline std::uint32_t bitreverse(std::uint32_t k)
    {
        std::uint32_t v = k;
        v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);  // swap odd and even bits
        v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);  // swap consecutive pairs
        v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);  // swap nibbles ...
        v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);  // swap bytes
        v = ( v >> 16             ) | ( v               << 16); // swap 2-byte long pairs
        return v;
    }
#endif

    /** Reverse the bits in a std::uint64_t.
      *
      * Calls twice the version for std::uint32_t.
      * \see bitreverse(std::uint32_t)
      */
    inline std::uint64_t bitreverse(std::uint64_t k)
    {
        return (std::uint64_t(bitreverse(std::uint32_t(k))) << 32)
              | std::uint64_t(bitreverse(std::uint32_t(k >> 32)));
    }

    struct generating_matrices_range_error : std::range_error
    {
        unsigned s, m, s_, m_; // current values of parsing
        generating_matrices_range_error(const char* what,
                unsigned s, unsigned m, unsigned s_, unsigned m_)
        : std::range_error(what), s(s), m(m), s_(s_), m_(m_)
        {
        }
    };

    /** Load generating matrices form an input stream (e.g., a file).
      *
      * @param is       the input stream to read from, see the description of
      *                 digitalseq_b2g for a description of how the generating
      *                 matrices should be represented
      * @param s        number of dimensions to read, if set to 0, set from input
      * @param m        number of columns to read, is set to 0, set from input
      *
      * @throws generating_matrices_range_error when s or m are non-zero and the
      * input stream does not have enough dimensions or columns to fulfill this
      * request.
      */
    template <typename UINT_T=std::uint64_t>
    std::vector<UINT_T>
    load_generating_matrices(std::istream& is, unsigned& s, unsigned& m)
    // throw(generating_matrices_range_error)
    {
        std::vector<UINT_T> Cs;
        // parse generating matrices from input stream `is`
        unsigned m_ = 0, s_ = 0; // unsigned should be big enough (modified because clang warns)
        std::string line;
        while(std::getline(is, line) && (s == 0 || s_ < s)) {
            s_++;
            std::istringstream str(line);
            if(s_ == 1) {
                // read first line completely
                std::copy(std::istream_iterator<UINT_T>(str),
                          std::istream_iterator<UINT_T>(), std::back_inserter(Cs));
                m_ = static_cast<unsigned>(Cs.size());
                if(m == 0) m = m_;
                if(m > m_) throw generating_matrices_range_error("load_generating_matrices::"
                       "m > m_::You want more columns than are present", s, m, s_, m_);
                // now throw away those columns we don't need
                Cs.resize(m);
            } else {
                unsigned n = 0;
                std::copy_if(std::istream_iterator<UINT_T>(str),
                             std::istream_iterator<UINT_T>(),
                             std::back_inserter(Cs),
                             [&](const UINT_T) { return n++ < m; });
                if(n == 0) { s_--; break; } // empty line terminates input (e.g. from stdin)
                if(n != m_) throw generating_matrices_range_error("load_generating_matrices::"
                        "incorrect number of integers on line", s, m, s_, m_);
            }
        }
        if(s == 0) s = s_;
        if(s > s_) throw generating_matrices_range_error("load_generating_matrices::"
                "s > s_::You want more dimensions than are present", s, m, s_, m_);
        return Cs;
    }

    /// Determine maximum bit depth needed to represent stream of integers.
    template <typename ForwardIterator,
              typename UINT_T=typename std::iterator_traits<ForwardIterator>::value_type>
    unsigned bit_depth(ForwardIterator Cs_begin, ForwardIterator Cs_end)
    {
        UINT_T mx = *std::max_element(Cs_begin, Cs_end);
        unsigned t = std::numeric_limits<UINT_T>::digits - 1;
        while(t && !((UINT_T(1) << (t-1)) & mx)) { t--; }
        return t;
    }

    /// Return 1 if bit `i` in `a` is set and 0 otherwise.
    /// If i is negative we index from the MSB.
    template <typename T>
    T bitget(T a, int i)
    {
        unsigned t = std::numeric_limits<T>::digits - 1;
        if(i < 0) return (a >> (t+i+1)) & T(1);
        return (a >> i) & T(1);
    }

    /** Binary matrix-matrix multiplication, C = A * B over Z_2.
      *
      * A is t2-by-t; represented as t integers
      * B is t-by-m; represented as m integers
      * C = A * B is t2-by-m; represented as m integers
      */
    template <typename RAIterator, typename RAIterator2, typename RAIterator3>
    void matrix_matrix_multiplication_over_Z2(
            unsigned t /* nb of cols of A */,
            unsigned m /* nb of cols of B and C */,
            RAIterator A, RAIterator2 B, RAIterator3 C)
    {
        // C = A * B where A is n-times-p, B is p-times-m and C is n-times-m
        for(unsigned k = 0; k < m; ++k) {    // we calculate the result column by column
            C[k] = 0;                        // column k of the result
            for(unsigned i = 0; i < t; ++i)  // iterate over the columns of A
                C[k] ^= bitget(B[k], i) * A[i];
        }
    }

    /** Produce list of scrambled generating matrices, Cout_j = M_j * C_j over Z_2.
      *
      * Scramble a list of s binary matrix of t-by-m by multiplying it from the
      * left with a binary matrix of t2-by-t.
      * The matrices in M are t2-by-t represented by t integers (up to t2 bits
      * precision).
      * The matrices in C are t-by-m represented by m integers.
      */
    template <typename FWIterator, typename FWIterator2, typename FWIterator3>
    void scramble(unsigned s, unsigned t, unsigned m, FWIterator M, FWIterator2 C, FWIterator3 Cout)
    {
        // M is t2-by-t; t integers
        // C is t-by-m; m integers
        // Cout is t2-by-m; m integers
        for(unsigned j = 0; j < s; ++j, M += t, C += m, Cout += m)
            matrix_matrix_multiplication_over_Z2(t, m, M, C, Cout);
    }

    /** Generate a random linear scramble.
      *
      * Generate s random linear binary scramble matrices of t2-by-t. Each
      * matrix is represented as a list of t integers with bits set up to
      * precision t2. These matrices are generated to be lower triangular with
      * all below diagonal entries uniformly random chosen from 0 and 1, and
      * guarantee the diagonal entries to be all 1's.
      *
      * @param s        number of matrices to generate (number of dimensions)
      * @param t2       bit depth of scramble, this is the number of rows of
      *                 the scramble matrices
      * @param t        bit depth of digital net, this is the number of columns
      *                 of the scramble matrices
      * @param[out] M   iterator to store the s matrices as lists of t
      *                 unsigned integers, needs have size s*t, the integers
      *                 should have at least t2 bits
      * @param gen      random number generator to use, should produce integers
      *                 with at least t2 bits, it should generate uniform
      *                 random integers on a call like `gen()`, e.g., the
      *                 generators defined in the C++11 <random> header
      */
    template <typename FWIterator, typename RNG>
    void generate_random_linear_scramble(unsigned s, unsigned t2, unsigned t, FWIterator M, RNG gen)
    {
        typedef typename std::iterator_traits<FWIterator>::value_type UINT_T;
        assert(t2 <= std::numeric_limits<UINT_T>::digits); // assuming this is base 2 digits
        assert(t <= t2);
        for(unsigned j = 0; j < s; ++j)
            for(unsigned k = 0; k < t; ++k, ++M) // iterate over columns of M_j
            {
                auto r = gen(); // generate a random integer
                //r = bitreverse(r); // bit reverse to not have trouble with LS bits? or not needed?
                r &= ((UINT_T(1) << (t2-k-1)) - 1) | (UINT_T(1) << (t2-k-1)); // truncate to t2-k LS bits
                r |= 1;  // set the diagonal bit (we will shift this in place next)
                // we now have the t2-k bits of the column but completely
                // shifted to the right, so do left shift:
                r <<= k; // we now have r = (x ... x 1 0 ... 0)
                         //                  ^       ^       ^
                         //           pos = t2-1     k       0
                         // which represents a binary column vector (of length t2):
                         //      0
                         //     ...
                         //      0
                         //      1
                         //      x
                         //     ...
                         //      x
                *M = r;
            }
    }

    /// Serve as an iterator or generator returning a constant stream of
    /// `zero`'s.
    template <typename T=int>
    struct cst_zero : std::iterator<std::input_iterator_tag, T>
    {
        const T zero;
        int n;

        explicit cst_zero(const T zero=0, int n=0) : zero(zero), n(n)
        {
        }

        cst_zero& operator++()
        {
            ++n;
            return *this;
        }

        cst_zero operator++(int)
        {
            cst_zero w(*this);
            ++(*this);
            return w;
        }

        cst_zero operator+(int m) const
        {
            cst_zero w(*this);
            w.n += m;
            return w;
        }

        bool operator==(const cst_zero& other) const
        {
            return (zero == other.zero) && (n == other.n);
        }

        bool operator!=(const cst_zero& other) const
        {
            return !(*this == other);
        }

        const T operator*() const
        {
            return zero;
        }

        const T operator()() const
        {
            return zero;
        }

    };


    /// Base 2 digital sequence point generator (in gray coded radical inverse ordering).
    template <typename FLOAT_T=long double, typename UINT_T=std::uint64_t>
    struct digitalseq_b2g : std::iterator<std::input_iterator_tag, std::vector<FLOAT_T>>
    {
        typedef FLOAT_T float_t;
        typedef UINT_T  uint_t;

        // TODO: add static_assert for base 2 types?
        static const int uint_t_digits2
            = std::numeric_limits<uint_t>::digits;     // assuming this is base 2 digits
        static const int float_t_digits2
            = std::numeric_limits<float_t>::digits;    // assuming this is base 2 digits
        static constexpr const float_t recipd
            = 1 / float_t(uint_t(1) << (uint_t_digits2-1)) / 2; // assuming base 2 types

        std::uint32_t k;         ///< index of the current point
        unsigned m;              ///< this base-2-digital-sequence has `n=2^m` points, nb of columns
        unsigned s;              ///< number of generating matrices, i.e., number of dimensions
        unsigned t;              ///< bit depth of generating matrices (nb of rows)
        uint_t n;                ///< maximum number of points `n=2^m`
        std::vector<uint_t> Csr; ///< radical inverse of `s` generating matrices consisting of
                                 /// `m` integers encoding the columns with lsb in the top row
        std::vector<uint_t> shfr;///< radical inverse of `s` dimensional integer digital shift vector
        std::vector<uint_t> cur; ///< `s` dimensional vector of current point as binary number
                                 /// completely shifted to the left, `cur_j = C_j k over Z_2`
        std::vector<float_t> x;  ///< the current point converted to float, `x = (cur ^ shfr) * recipd`

        std::size_t Csr_hash;
        std::size_t shfr_hash;

        /** Constructor from iterator.
          *
          * @param s            number of dimensions
          * @param m            number of points is `2^m`, asserted to fit the
          *                     `uint_t` type, so `m` is at least one less than
          *                     the number of bits available in the `uint_t`
          *                     type
          * @param Cs_begin     iterator to the `s` generating matrices,
          *                     they are copied into this `digitalseq_b2g` object
          * @param shift_begin  iterator to the `s` digital shifts of type uint_t
          *
          * The generating matrices are stored as a flat list of $s \times m$
          * integers which represent the columns of the generating matrices
          * stored as integers with the LSB in the top row. The constructor
          * will reverse the bits with the old LSB shifted completely to the
          * left of the `uint_t` type, this is stored in the member variable
          * `Csr` (`r` for radical inverse). Similarly for the digital shift
          * vector which is stored in the variable `shfr`.
          *
          * The first output point of this sequence will have been calculated
          * and is immediately available as a floating point vector of type
          * `std::vector<float_t>` as the member variable `x`, or by using the
          * dereference operator on this object. The un-shifted vector can be
          * obtained as an integer aligned as far to the left as possible as
          * the member variable `cur`.
          *
          * The sequence can be brought back into the initial state by calling
          * `reset()`, which will calculate the first point again and set the
          * appropriate state.
          *
          * Skipping to an arbitrary point in the sequence is possible by
          * calling the `set_state(std::uint32_t)` method.
          *
          * Advancing to the next point in the sequence is done by the `next()`
          * method or `operator++`.
          *
          * @see x, operator*(), cur, reset(), set_state(std::uint32_t), next()
          */
        template <typename InputIterator, typename InputIterator2 = cst_zero<uint_t>>
        digitalseq_b2g(const unsigned& s, const unsigned& m,
                       InputIterator Cs_begin, InputIterator2 shift_begin = cst_zero<uint_t>{})
        : k(0), m(m), s(s), t(bit_depth(Cs_begin, Cs_begin+m*s)), n(uint_t(1) << m),
          Csr(Cs_begin, Cs_begin+m*s), shfr(shift_begin, shift_begin+s), cur(s), x(s)
        {
            assert(this->m < uint_t_digits2); // assuming the type is base 2
            bitreverse_and_hash();
            reset();
        }

        digitalseq_b2g() {};

        /** Constructor from input stream (e.g., file).
          *
          * @param s            the number of dimensions to generate
          *                     if set to 0 then the number of available
          *                     dimensions from the input stream will be used
          * @param m            number of points 2^m to generate
          *                     if set to 0 the number of columns will be used for m
          * @param is           an input stream, e.g., std::ifstream("filename")
          *
          * Note that `s` and `m` are passed by reference such that if
          * they were passed in as 0 they could be read out by the caller.
          * There is another version of this constructor who takes those
          * parameters as rvalue references (but then they both need to be).
          *
          * @see The other constructor for more explanation.
          */
        digitalseq_b2g(unsigned& s, unsigned& m,
                       std::istream& is)
        : digitalseq_b2g(s, m, std::begin(load_generating_matrices(is, s, m)))
        {
        }

        /** Constructor from input stream (e.g., file) alternative form.
          *
          * This version allows to call the constructor with two literals for
          * `s` and `m` in case of parsing from input stream.
          *
          * @see The other constructor for more explanation.
          */
        digitalseq_b2g(unsigned&& s, unsigned&& m,
                       std::istream& is)
        : digitalseq_b2g(s, m, std::begin(load_generating_matrices(is, s, m)))
        {
        }

        /** Advance to the next point in the sequence.
          *
          * No checks are made if the sequence is already past the end.
          * @see past_end to check if all points have been generated
          */
        std::vector<float_t> next()
        {
            // figure out which bit changed in gray code ordering
            k += 1;
            unsigned ctz = count_trailing_zero_bits(k) % m;
            // then update the point by adding in that column of the generating matrix
            for(unsigned j = 0; j < s; ++j) {
                cur[j] ^= Csr[j*m+ctz];
                x[j] = recipd * (cur[j] ^ shfr[j]);
            }
            return x;
        }

        /** Advance to the next point in the sequence.
          *
          * This one is here to be able to use this class as a forward iterator.
          */
        digitalseq_b2g& operator++()
        {
            next();
            return *this;
        }

        digitalseq_b2g operator++(int)
        {
            digitalseq_b2g w(*this);
            ++(*this);
            return w;
        }

        /** Return the current point as a floating point vector.
          *
          * This one is here to be able to use this class as a forward iterator.
          */
        const std::vector<float_t>& operator*() const
        {
            return x;
        }

        /// Skip to the k-th point in the sequence.
        void set_state(std::uint32_t new_k)
        {
            if(new_k == 0) {
                reset();
                return;
            }
            std::fill(cur.begin(), cur.end(), 0);
            k = new_k;
            std::uint32_t gk = (k >> 1) ^ k; // we use Gray code ordering
            for(unsigned j = 0; j < s; ++j) {
                for(unsigned i = 0; i < m; ++i) {
                    if((uint_t(1) << i) & gk)
                        cur[j] ^= Csr[j*m+i];
                }
                x[j] = recipd * (cur[j] ^ shfr[j]);
            }
        }

        /// Reset the sequence to the initial state of the object `k = 0`.
        void reset()
        {
            k = 0;
            std::fill(cur.begin(), cur.end(), 0);
            for(unsigned j = 0; j < s; ++j)
                x[j] = recipd * (cur[j] ^ shfr[j]); // cur[j] is zero, but we don't care
        }

        /// Bitreverse the generating matrices and shift vector. This function
        /// is its own inverse. This function also calculates a hash value over
        /// the generating matrices and shift vector such that `operator==` can
        /// use this to compare to an other class.
        void bitreverse_and_hash()
        {
            const std::size_t seed = s;
            const std::size_t m1 = std::size_t(0x9e3779b9);
            const std::size_t m2 = std::size_t(bitreverse(uint64_t(m1)) + bitreverse(uint64_t(m1))/2);
            const std::size_t m3 = (seed << 6) + (seed >> 2);
            const std::size_t m4 = std::size_t(bitreverse(uint64_t(m3)) + bitreverse(uint64_t(m3))/2);
            const std::size_t v1 = m1 + m2 + m3 + m4;
            const std::size_t v2 = std::size_t(bitreverse(uint64_t(v1)));
            std::size_t w1;
            std::size_t w2;
            w1 = v1;
            w2 = v2;
            Csr_hash  = 0;
            for(auto& a : Csr) {
                a = bitreverse(a);
                Csr_hash ^= a + w1 + w2;
                w1 <<= 1; w2 >>= 1;
            }
            w1 = v1;
            w2 = v2;
            shfr_hash = 0;
            for(auto& a : shfr) {
                a = bitreverse(a);
                shfr_hash ^= a & (w1 + w2); // modified such that zero shift gives zero hash
                w1 <<= 1; w2 >>= 1;
            }
        }

        /// Check if we have passed the end of the sequence.
        bool past_end() const
        {
            return k >= n;
        }

        struct digitalseq_b2g_dummy_end
        {
        };

        digitalseq_b2g_dummy_end end() const
        {
            return digitalseq_b2g_dummy_end();
        }

        /// Test if two digital sequences are the same and at the same position.
        /// Note: this is using `Csr_hash` and `shfr_hash` to check for
        /// equality of the generating matrices and the shift vector.
        bool operator==(const digitalseq_b2g& other) const
        {
            return (Csr_hash == other.Csr_hash) && (shfr_hash == other.shfr_hash) && (k == other.k);
        }

        // This is not the opposite of the above `operator==`, this
        // `operator!=` compares to the special type `digitalseq_b2g_dummy_end`
        // to use this class directly as an iterator. Might be removed in
        // favour of real iterators.
        bool operator!=(const digitalseq_b2g_dummy_end&) const
        {
            return !past_end();
        }

    };

    /// Helper to generate a vector of scrambled and shifted digital sequences.
    /// (Such that we have RAII in the constructor of `multiplexed_scrambled_digitalseq_b2g`.)
    template <typename GEN, typename InputIterator, typename InputIterator2, typename InputIterator3>
    std::vector<GEN> make_scrambled_generators(unsigned M, unsigned s, unsigned m, unsigned t,
                                               InputIterator scramble_matrices_begin,
                                               InputIterator2 Cs_begin, InputIterator3 shifts_begin)
    {
        typedef typename GEN::uint_t uint_t;
        std::vector<GEN> seq;
        seq.reserve(M);
        for(unsigned v = 0; v < M; ++v)
        {
            std::vector<uint_t> scrambled_Cs(s*m);
            scramble(s, t, m, scramble_matrices_begin+v*s*t, Cs_begin, std::begin(scrambled_Cs));
            seq.emplace_back(GEN(s, m, std::begin(scrambled_Cs), shifts_begin+v*s));
        }
        return seq;
    }

    /// Base 2 digital sequence point generator (in gray coded radical inverse ordering).
    template <typename FLOAT_T=long double, typename UINT_T=std::uint64_t>
    struct multiplexed_scrambled_digitalseq_b2g
    {
        typedef digitalseq_b2g<FLOAT_T, UINT_T> seq_t;
        std::vector<digitalseq_b2g<FLOAT_T, UINT_T>> seq;
        std::vector<UINT_T> scramble_matrices;
        std::vector<UINT_T> shifts;

        template <typename InputIterator, typename InputIterator2, typename InputIterator3>
        multiplexed_scrambled_digitalseq_b2g(unsigned M, unsigned s, unsigned m, unsigned t,
                                             InputIterator scramble_matrices_begin,
                                             InputIterator2 Cs_begin, InputIterator3 shifts_begin)
        : seq(make_scrambled_generators<seq_t>(M, s, m, t, scramble_matrices_begin, Cs_begin, shifts_begin)),
          scramble_matrices(scramble_matrices_begin, scramble_matrices_begin + M*s*t),
          shifts(shifts_begin, shifts_begin + M*s)
        {
        }

        multiplexed_scrambled_digitalseq_b2g operator++()
        {
            for(auto& g : seq) ++g;
            return *this;
        }

    };

    /** Some default Sobol sequence with a maximum of 250 dimensions and 2^32
      * points.
      */
    template <typename ForwardIterator>
    void JK2008_sobol_matrices(ForwardIterator Cs, unsigned s=250, unsigned m=32)
    {
        const int mmax = 32;
        typedef uint32_t UINT_T;
        UINT_T Cs_all[][mmax] = { // we know these are 32 bit numbers
            // and we know this has 32 columns (we cannot deduce that here)
{ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648 },
{ 1, 3, 5, 15, 17, 51, 85, 255, 257, 771, 1285, 3855, 4369, 13107, 21845, 65535, 65537, 196611, 327685, 983055, 1114129, 3342387, 5570645, 16711935, 16843009, 50529027, 84215045, 252645135, 286331153, 858993459, 1431655765, 4294967295 },
{ 1, 3, 6, 9, 23, 58, 113, 163, 278, 825, 1655, 2474, 5633, 14595, 30470, 43529, 65815, 197434, 394865, 592291, 1512982, 3815737, 7436151, 10726058, 18284545, 54132739, 108068870, 161677321, 370540567, 960036922, 2004287601, 2863268003 },
{ 1, 3, 4, 10, 31, 46, 69, 201, 283, 676, 1946, 2919, 4126, 12333, 16449, 41155, 127236, 189066, 284639, 826286, 1155333, 2781833, 7989211, 11983780, 16777498, 50332327, 67110814, 167775085, 520097793, 771764227, 1157644292, 3372261386 },
{ 1, 2, 4, 13, 31, 59, 94, 185, 346, 1012, 1669, 3343, 4443, 9206, 18049, 56578, 123204, 234445, 370399, 740795, 1327134, 3985465, 7204954, 14397620, 16900421, 33788879, 67479259, 218844598, 521420801, 993841154, 1584263172, 3118182413 },
{ 1, 2, 6, 12, 19, 36, 106, 223, 263, 526, 1557, 3112, 4985, 9467, 27501, 56785, 65554, 131110, 393324, 786643, 1245460, 2359850, 6948479, 14617847, 17240702, 34481909, 102067576, 204001785, 326636395, 620307677, 1802201857, 3722304770 },
{ 1, 3, 5, 11, 26, 41, 124, 199, 381, 964, 1144, 2255, 7778, 8678, 25118, 58913, 73315, 205285, 352795, 779818, 1777273, 2761164, 7955047, 13428461, 23469828, 65937928, 67200543, 134475298, 520356198, 570950638, 1713308673, 3995207683 },
{ 1, 2, 5, 10, 17, 36, 72, 180, 366, 633, 1040, 2086, 5197, 10430, 17791, 37469, 74840, 186514, 372003, 642759, 1048943, 2097787, 5243925, 10487852, 17831004, 37759130, 75515191, 188781289, 383853878, 663935723, 1090890035, 2187970273 },
{ 1, 2, 5, 10, 20, 43, 86, 142, 284, 538, 1111, 2188, 5401, 10768, 21571, 42151, 85327, 135838, 271711, 577213, 1053976, 2107922, 5264454, 10527917, 21056859, 45224629, 90449161, 149474867, 298849284, 562046984, 1159745553, 2288036897 },
{ 1, 2, 7, 13, 25, 41, 81, 218, 460, 923, 1102, 2300, 7555, 14181, 25034, 44948, 89168, 218328, 477643, 982934, 1055831, 2110677, 7366098, 13674431, 26296326, 43212815, 85404702, 229529636, 483482696, 965956851, 1150765469, 2401173313 },
{ 1, 2, 5, 8, 16, 54, 121, 196, 490, 949, 1534, 2953, 4498, 11123, 20497, 61492, 110716, 256204, 396794, 860035, 1159559, 2353997, 5638264, 9245894, 17932783, 54772669, 130440686, 214822847, 531628523, 940573623, 1478493691, 3025144705 },
{ 1, 2, 4, 12, 26, 53, 105, 212, 299, 656, 1351, 2634, 5234, 14563, 31046, 58952, 116854, 232687, 334172, 598653, 1170463, 2339899, 4546679, 13223149, 28317016, 57686641, 114326533, 226552846, 319937566, 710128697, 1390764147, 2850669793 },
{ 1, 3, 5, 10, 31, 49, 71, 204, 496, 644, 1449, 3706, 4123, 9272, 26717, 60663, 90536, 259705, 448542, 601138, 1142850, 3395782, 5682671, 11128501, 31457774, 52429494, 69207531, 207621820, 519049716, 726673037, 1591766451, 3948995137 },
{ 1, 3, 6, 9, 28, 35, 66, 197, 399, 597, 1855, 2209, 4103, 12298, 24602, 36906, 114782, 143590, 270797, 807568, 1636016, 2448116, 7595832, 9054379, 16777245, 50331680, 100663364, 150995148, 469762451, 587203190, 1107298173, 3305113700 },
{ 1, 2, 4, 15, 21, 42, 89, 185, 376, 826, 1214, 2225, 4388, 9102, 17879, 64430, 83414, 166828, 378322, 740259, 1476039, 3267465, 5207454, 8583984, 16864486, 33728522, 67469400, 252424379, 353878396, 708006709, 1498068139, 3112810651 },
{ 1, 3, 4, 11, 27, 35, 98, 161, 421, 878, 1461, 3414, 5556, 15701, 21936, 36190, 124331, 179581, 509385, 568796, 1938540, 4083378, 4637145, 15410148, 18715757, 54415025, 71746013, 199959535, 438146166, 540954258, 1715913151, 2867536718 },
{ 1, 2, 4, 15, 28, 40, 111, 230, 497, 522, 1361, 4080, 5448, 12245, 21823, 57108, 103613, 155088, 483442, 798924, 1597850, 2929379, 5059772, 13848018, 18374774, 36483267, 72167814, 265503435, 488142035, 707566900, 1800233351, 3922957001 },
{ 1, 3, 4, 15, 22, 38, 108, 182, 386, 519, 1867, 4057, 5936, 16202, 22234, 52724, 77829, 172044, 516114, 651305, 2064506, 2715792, 6967790, 14058161, 18847433, 53038558, 74055803, 265699475, 387887594, 623628990, 1752020703, 3118431736 },
{ 1, 2, 5, 10, 25, 47, 117, 221, 487, 937, 1080, 3083, 5595, 12202, 16895, 33668, 115784, 154844, 447525, 986156, 1664434, 3657554, 6037620, 16273631, 18391522, 37065635, 89846817, 183983140, 402769326, 755130231, 1879494680, 3608084525 },
{ 1, 3, 7, 13, 29, 60, 115, 130, 388, 906, 1680, 3745, 7759, 14705, 16390, 49166, 114714, 213041, 475246, 983230, 1884663, 2130696, 6358804, 14847275, 27531487, 61372368, 127114825, 240957823, 268435484, 805306431, 1879048308, 3489661071 },
{ 1, 3, 7, 11, 22, 60, 81, 138, 405, 955, 1498, 2972, 8105, 11114, 16464, 49289, 115090, 181168, 361932, 986016, 1335288, 2272224, 6619589, 15598386, 24593480, 48613420, 133159525, 181772490, 268443560, 805317481, 1879064663, 2952839298 },
{ 1, 2, 6, 11, 28, 49, 126, 142, 275, 800, 1496, 3805, 6482, 15811, 16485, 32950, 98679, 181140, 459945, 806210, 2072039, 2338992, 4496892, 13152392, 24632728, 62206012, 106061033, 259978147, 268441692, 536886992, 1610630597, 2952826606 },
{ 1, 3, 5, 9, 16, 38, 86, 145, 421, 723, 1048, 2485, 4597, 12110, 16804, 49872, 82973, 149948, 266725, 634728, 1425906, 2392641, 6817208, 11963247, 16913917, 40167133, 73961479, 200748303, 268518428, 805456319, 1342444000, 2416553825 },
{ 1, 3, 4, 11, 18, 49, 107, 156, 427, 588, 1357, 2126, 7114, 13121, 16851, 49890, 66697, 182933, 302910, 817522, 1744703, 2619761, 7053115, 9517434, 22454057, 35518795, 116236098, 217659863, 268508905, 805501851, 1074060196, 2953622485 },
{ 1, 3, 4, 10, 27, 47, 124, 169, 497, 700, 1231, 3956, 5045, 12368, 17617, 53074, 70614, 176348, 427335, 724328, 1971860, 2686397, 8049612, 10760432, 19760447, 66391499, 76994430, 213024558, 270407551, 807993133, 1081792379, 2695118631 },
{ 1, 2, 5, 13, 25, 37, 94, 171, 320, 518, 1930, 3601, 5169, 9186, 18384, 36533, 87405, 221644, 393245, 639018, 1491010, 2719875, 5652743, 8962696, 32835220, 61770940, 89936635, 142159861, 269931627, 539599686, 1347815052, 3498593435 },
{ 1, 3, 5, 12, 24, 44, 81, 185, 486, 577, 1794, 3846, 4745, 8852, 18100, 52733, 87144, 208863, 414503, 778307, 1392644, 2981903, 8306717, 9977888, 30589001, 65323157, 80708023, 135987960, 269747940, 808357191, 1350096267, 3230887314 },
{ 1, 2, 7, 11, 16, 50, 64, 173, 364, 835, 1064, 2784, 8152, 11274, 17682, 35125, 121547, 189757, 286686, 797699, 1148165, 2689292, 6216347, 14591394, 18767602, 44003309, 129794113, 201015471, 269500779, 539738952, 1884996664, 2966555346 },
{ 1, 3, 7, 10, 22, 50, 110, 239, 364, 747, 2017, 3319, 5573, 8875, 21956, 58024, 103875, 156322, 267733, 975504, 2069947, 3482239, 4494551, 9429140, 29152054, 62008419, 80578291, 186175176, 297587511, 867314784, 1959626484, 2870529730 },
{ 1, 2, 6, 9, 19, 46, 74, 192, 471, 758, 1574, 3931, 6632, 16020, 22949, 48735, 120935, 195726, 335384, 619323, 1541741, 3755929, 6428729, 8682603, 32137104, 55603498, 127205829, 217875162, 300566122, 592465810, 1737840684, 2633776204 },
{ 1, 3, 5, 11, 29, 32, 118, 233, 335, 951, 1736, 2362, 6875, 10124, 23440, 58419, 73024, 142626, 411632, 658023, 1809753, 3190026, 4802462, 14633637, 33365117, 43773556, 122683503, 189300929, 301728033, 848937845, 1465272300, 3141700164 },
{ 1, 3, 7, 12, 22, 55, 68, 212, 364, 697, 1110, 3688, 6066, 12492, 22217, 62017, 103132, 261361, 263188, 921267, 1398479, 3928654, 5231309, 9862346, 22938694, 50499152, 78010343, 252293155, 291477495, 855732251, 1956796322, 3473387764 },
{ 1, 3, 4, 12, 20, 43, 81, 255, 289, 974, 1224, 3395, 7391, 14444, 22530, 62727, 83080, 249240, 269759, 654202, 1121710, 3537758, 6118895, 13613958, 25042059, 61503644, 102095155, 226309870, 292487173, 866091023, 1179156496, 3451027495 },
{ 1, 2, 5, 10, 29, 33, 88, 245, 445, 1019, 1832, 2754, 7259, 14450, 23346, 45676, 67220, 133435, 400246, 717498, 1116020, 3508799, 7871998, 12663714, 27417695, 38189178, 114023210, 210989639, 296961745, 577964527, 1454611091, 2891116980 },
{ 1, 2, 7, 14, 16, 47, 111, 209, 403, 682, 1254, 4047, 7724, 8938, 23384, 44941, 98709, 197286, 378097, 544750, 2088531, 4170260, 8346276, 10382838, 22856928, 56250307, 112680507, 188818123, 289577767, 593948531, 1993195625, 3956736221 },
{ 1, 2, 7, 9, 22, 47, 70, 251, 414, 699, 1902, 3892, 7667, 14218, 23315, 47649, 98777, 180802, 280311, 626054, 1317515, 3487889, 7751334, 9473360, 29219668, 56736472, 111312960, 237289712, 296239503, 592855709, 1993377598, 2663239648 },
{ 1, 3, 6, 10, 24, 59, 66, 236, 260, 655, 1556, 2345, 4577, 15381, 22186, 63335, 114719, 147506, 491612, 803037, 1474910, 3523160, 5343058, 10013514, 29628145, 48871347, 100467039, 223330907, 299059028, 853412672, 1705727721, 2916857736 },
{ 1, 3, 4, 15, 31, 44, 70, 175, 277, 831, 1137, 4042, 7907, 11945, 16413, 42795, 65869, 197567, 263462, 986197, 2038680, 2892400, 4610761, 11503335, 18112166, 54298882, 74727431, 264677131, 519635216, 780731187, 1077937258, 2812743657 },
{ 1, 3, 5, 15, 31, 55, 126, 134, 314, 871, 1467, 3668, 7369, 13256, 29643, 40654, 67009, 200414, 335337, 995479, 2061841, 3641387, 8329292, 8732407, 20933283, 57189226, 94726818, 241214313, 481520295, 860005222, 1927965368, 2643632977 },
{ 1, 3, 4, 13, 26, 52, 89, 159, 316, 845, 1205, 3164, 6545, 12578, 22388, 33270, 66807, 199924, 268528, 865021, 1725159, 3444435, 5773706, 10236181, 20469545, 54559076, 78513617, 208342925, 432806940, 833633854, 1448233546, 2190458044 },
{ 1, 3, 4, 13, 27, 53, 113, 144, 377, 901, 1363, 3788, 7849, 14864, 27940, 44118, 73410, 211638, 289832, 896846, 1710835, 3553226, 7698339, 10328834, 23291911, 62104073, 86442262, 240255790, 533660740, 967642593, 1746737897, 2722052092 },
{ 1, 2, 7, 15, 21, 52, 69, 180, 347, 671, 1591, 3392, 4796, 14913, 21182, 34374, 70577, 145491, 480645, 1018134, 1312305, 3542349, 4330158, 12270202, 21910254, 44795590, 108204719, 234240632, 333468393, 953255625, 1419580090, 2335193676 },
{ 1, 3, 7, 12, 19, 62, 65, 242, 381, 951, 1668, 3797, 5662, 13870, 22904, 54207, 69788, 211188, 479090, 845219, 1198518, 4057991, 4677330, 16637458, 23593277, 59114310, 114362366, 236260718, 389091465, 894777541, 1589661479, 3727680867 },
{ 1, 3, 4, 8, 25, 52, 96, 179, 337, 984, 1535, 2741, 8030, 15309, 29914, 38136, 72377, 209231, 289760, 569230, 1601323, 3647067, 6649030, 12096711, 21830340, 65483200, 92433096, 171617489, 511142885, 949071749, 1880190774, 2637014631 },
{ 1, 2, 5, 9, 25, 42, 92, 185, 382, 766, 1256, 2773, 8118, 9577, 19393, 34702, 72451, 143111, 349452, 630800, 1594419, 2663030, 5788133, 11709127, 23284608, 48203542, 78191165, 172830051, 512250079, 655206568, 1309888079, 2402586509 },
{ 1, 3, 7, 13, 16, 33, 77, 157, 381, 1010, 1706, 3871, 5173, 11110, 20685, 44497, 70627, 206728, 475733, 887205, 1136395, 2269470, 4915510, 9634657, 23726016, 64883905, 116855490, 262350021, 352667080, 671712728, 1461404921, 2705295028 },
{ 1, 3, 6, 12, 30, 36, 121, 226, 486, 751, 1275, 2253, 5774, 12849, 19020, 43943, 84565, 240524, 477748, 1026886, 1657781, 2861423, 6311633, 13214639, 25582914, 36084668, 92079483, 183175911, 334631916, 953759997, 1505253825, 2473759632 },
{ 1, 3, 7, 13, 30, 57, 119, 216, 467, 707, 1507, 2208, 5649, 11811, 23108, 34236, 88620, 230750, 416655, 952136, 1597345, 3640594, 7293988, 15706953, 29098402, 36438549, 73270313, 177146455, 306713499, 614540898, 1209159153, 2703534468 },
{ 1, 2, 6, 8, 26, 62, 67, 135, 432, 1003, 1364, 3247, 5593, 12035, 25348, 62478, 91666, 194596, 488061, 709828, 1905719, 3419227, 6252479, 11664379, 32012150, 55630042, 77594631, 264830986, 287375388, 544276534, 1943011417, 3608871097 },
{ 1, 2, 4, 12, 29, 53, 78, 141, 448, 843, 1667, 2521, 4722, 11984, 29281, 53500, 90134, 196130, 332398, 929511, 1667627, 3788925, 5398987, 12254556, 25892003, 63449011, 123363737, 179043814, 346232105, 666051705, 1637443271, 4259772225 },
{ 1, 3, 7, 14, 17, 34, 82, 226, 443, 578, 1219, 2286, 7598, 13929, 31374, 43327, 92005, 237467, 400404, 865320, 1474122, 2925023, 4356826, 12822992, 32466632, 44932597, 91166868, 176372742, 422142477, 1045009686, 1725498668, 2648115779 },
{ 1, 3, 4, 10, 27, 63, 111, 171, 509, 564, 1911, 3728, 4248, 9860, 22197, 51403, 82244, 254185, 329497, 593464, 2036065, 3482810, 7623129, 9750116, 26431923, 43828934, 110055505, 212714695, 401692498, 673368771, 1180452440, 4006801112 },
{ 1, 2, 6, 10, 26, 53, 86, 161, 479, 789, 1057, 3189, 4306, 11783, 24840, 61724, 95295, 184396, 459156, 525193, 1836468, 4132862, 4986208, 9708019, 25785458, 62181594, 100597787, 258080823, 351600720, 578289835, 1908408773, 3744334624 },
{ 1, 3, 5, 10, 23, 34, 122, 181, 503, 534, 1569, 2687, 5311, 12256, 17972, 62043, 90314, 260936, 262646, 786965, 1312292, 2624117, 6034600, 8925122, 31999566, 47510254, 131948861, 140245342, 411043799, 703596650, 1390940827, 3215467925 },
{ 1, 3, 6, 13, 24, 32, 91, 144, 324, 535, 1597, 3184, 6878, 12780, 17192, 45793, 76697, 171517, 262424, 787081, 1574759, 3411530, 6298784, 8404268, 23877865, 37781898, 85028808, 140386683, 418461282, 835625704, 1801740169, 3347615694 },
{ 1, 2, 4, 10, 17, 57, 116, 160, 469, 568, 1142, 2212, 5599, 8233, 30287, 57552, 87423, 232444, 263287, 526502, 1054171, 2629667, 4486750, 14999785, 30496011, 42175324, 122684834, 149424798, 298324397, 582494343, 1463313281, 2151727296 },
{ 1, 3, 4, 10, 19, 60, 124, 230, 415, 563, 1638, 2253, 5589, 9402, 32547, 62237, 122214, 202126, 263688, 788500, 1053746, 2631269, 5011657, 15792607, 32624297, 60504351, 108555617, 146857344, 428453905, 589285435, 1461565042, 2451030271 },
{ 1, 2, 4, 13, 26, 34, 126, 150, 474, 571, 1114, 2273, 7003, 14297, 16957, 61523, 81654, 236387, 263557, 527061, 1056543, 3424023, 6838530, 8963902, 33127791, 39545245, 124600051, 149396332, 293340571, 599263994, 1830035323, 3755491235 },
{ 1, 2, 5, 13, 18, 46, 67, 231, 429, 605, 1241, 3013, 6388, 8577, 21531, 38963, 123002, 203911, 268571, 533305, 1331562, 3445158, 4848709, 12269293, 17323447, 60070508, 111273126, 155391823, 321389053, 799277718, 1657275698, 2242644338 },
{ 1, 2, 5, 15, 25, 45, 74, 224, 509, 602, 1244, 2953, 7371, 14271, 21161, 36634, 129965, 243856, 269665, 537893, 1333214, 3966475, 6677012, 12027441, 19682920, 59278515, 132245296, 156733421, 319835750, 781760680, 1913760024, 3712169896 },
{ 1, 3, 7, 14, 31, 50, 101, 145, 374, 600, 1768, 4050, 8005, 14910, 27774, 55466, 65794, 177380, 268745, 800126, 1864266, 3727554, 8219531, 13263790, 26281963, 38452021, 96470667, 156503886, 465575467, 1069037658, 2096917228, 3933525979 },
{ 1, 2, 4, 12, 29, 60, 123, 251, 404, 588, 1175, 2374, 7054, 15482, 29945, 59280, 120384, 224906, 270202, 538101, 1076865, 3201389, 7710684, 15958743, 32102348, 65498356, 105126799, 151534712, 305689853, 634972060, 1872877149, 4076302006 },
{ 1, 2, 5, 11, 31, 60, 118, 152, 266, 595, 1246, 3061, 5609, 14802, 30122, 63270, 70671, 150038, 266278, 536641, 1336561, 2934190, 8221999, 15894037, 31290401, 40535119, 70873317, 157978509, 320169317, 788707067, 1443959219, 3858645782 },
{ 1, 2, 6, 11, 19, 61, 114, 234, 257, 592, 1215, 3505, 5631, 8549, 28808, 64466, 125742, 151069, 266273, 536145, 1597629, 2938295, 5082100, 16161654, 30251701, 60995488, 68812740, 157707548, 313287793, 907335918, 1449238284, 2263267400 },
{ 1, 2, 4, 13, 21, 43, 95, 159, 329, 664, 1346, 2692, 8049, 8421, 18842, 40812, 95444, 134143, 280970, 565064, 1143962, 3541830, 5261449, 10788196, 23949518, 38405061, 91509747, 180663709, 342172519, 677033672, 2024894924, 2303853035 },
{ 1, 2, 7, 13, 26, 56, 78, 242, 389, 750, 1458, 3261, 8053, 14700, 27985, 41237, 100336, 241751, 289999, 566209, 1933835, 3652629, 7107109, 15243387, 19262118, 67047759, 108653349, 189993902, 401206922, 822288166, 2075259307, 4007077504 },
{ 1, 2, 5, 10, 17, 44, 69, 192, 386, 751, 1472, 2080, 4185, 10482, 19952, 43607, 122091, 230856, 280628, 565887, 1428646, 2862364, 4195728, 12060868, 16781706, 52964091, 105662950, 185376884, 370526389, 595830581, 1198002655, 2826709525 },
{ 1, 2, 7, 11, 29, 56, 65, 223, 439, 737, 1535, 3108, 4731, 13465, 26467, 47947, 117030, 255386, 287360, 571221, 1956635, 3133911, 7341129, 15208135, 18356099, 55589046, 111703301, 191139818, 385461764, 871325198, 1328917009, 3656077358 },
{ 1, 3, 5, 8, 18, 53, 96, 169, 416, 659, 2012, 2064, 5683, 11885, 31411, 58759, 80582, 231189, 294169, 847104, 1387810, 2337663, 4981156, 14156440, 26478539, 42207277, 113776193, 177745649, 519076723, 585426877, 1582120122, 3019607446 },
{ 1, 2, 6, 11, 31, 44, 88, 255, 463, 713, 1460, 3708, 4796, 13073, 19786, 35249, 109683, 209070, 280873, 558369, 1686330, 3083547, 7870295, 12074395, 24665632, 64010318, 121758922, 194492336, 397801073, 973571240, 1305590562, 3343322430 },
{ 1, 3, 6, 8, 20, 63, 77, 218, 508, 766, 1937, 3652, 5837, 9669, 27835, 46431, 116103, 226942, 289422, 833795, 1693036, 2318291, 5510785, 15738140, 19688769, 55097257, 128041474, 185038853, 527727630, 980720668, 1558034475, 2630871154 },
{ 1, 2, 6, 12, 16, 50, 111, 254, 366, 664, 1416, 4062, 8017, 10981, 30018, 59603, 107815, 149553, 289387, 574196, 1692018, 3318458, 4574165, 12748623, 27369665, 63063827, 98858066, 177269381, 398827944, 1014646149, 2030334412, 2731337061 },
{ 1, 2, 5, 14, 29, 62, 82, 207, 265, 641, 1530, 2527, 7098, 14129, 27866, 36131, 118011, 190794, 292456, 571537, 1418718, 3829176, 7446836, 15886036, 22658878, 51271877, 71420184, 181330087, 381971864, 612929375, 1810221122, 3577668331 },
{ 1, 2, 5, 13, 17, 43, 87, 237, 431, 677, 1407, 2430, 7548, 12153, 16756, 39781, 112974, 234265, 286196, 570971, 1441022, 3605377, 4324607, 10945411, 21888250, 59442062, 115997931, 168231845, 351665340, 652149576, 1917189395, 3025601005 },
{ 1, 2, 4, 10, 16, 44, 88, 139, 309, 746, 1474, 2840, 4784, 11597, 17959, 32842, 91811, 169831, 282737, 566491, 1128362, 2765305, 4547424, 12229757, 24396485, 34376092, 75854237, 187871647, 371983771, 787745681, 1331958657, 3164090285 },
{ 1, 2, 6, 15, 16, 39, 103, 224, 264, 702, 1466, 3856, 6276, 10732, 24511, 63259, 107677, 165844, 279535, 578491, 1697554, 4084866, 4564451, 9864623, 25184060, 62185210, 75103028, 173699303, 389146373, 1062622376, 1739695506, 2935915879 },
{ 1, 3, 4, 13, 28, 53, 69, 243, 490, 1012, 1478, 3466, 5992, 11947, 18223, 51295, 74444, 189370, 336675, 974912, 1386237, 3862002, 6631884, 15647647, 24213324, 52562119, 114839980, 217113859, 298147873, 987113058, 1282414235, 2151689572 },
{ 1, 3, 4, 8, 30, 54, 126, 255, 280, 771, 1342, 3418, 7579, 8419, 20265, 46961, 83927, 224874, 345808, 1011036, 1393556, 3105009, 7167753, 16204583, 26547552, 51116531, 99896884, 251966585, 285819635, 979591950, 1735482155, 3116584310 },
{ 1, 3, 5, 12, 27, 47, 75, 213, 281, 777, 1335, 3905, 6119, 8276, 31989, 60748, 86015, 221822, 347826, 1008002, 1133257, 4156732, 8214864, 8602074, 22099495, 64454751, 87299812, 267395953, 300710284, 855171796, 1315264794, 3013296396 },
{ 1, 2, 5, 12, 16, 48, 78, 159, 338, 865, 1832, 3518, 5360, 14257, 20711, 49544, 125621, 228151, 388505, 751239, 1177980, 3883274, 5351877, 16482349, 23609964, 42268900, 69577103, 252412092, 412369707, 969472953, 1258663673, 3891760557 },
{ 1, 2, 6, 10, 20, 39, 120, 176, 498, 960, 1949, 2899, 7707, 13373, 32327, 33521, 128777, 145629, 389479, 667724, 1960167, 2232104, 5081263, 12499395, 28049817, 41681759, 106956293, 225708046, 403445784, 537934905, 1613784651, 4034445039 },
{ 1, 2, 4, 13, 28, 52, 75, 177, 270, 808, 1869, 3978, 5231, 11984, 17348, 57041, 106950, 232149, 370635, 755401, 1422847, 4170882, 6677838, 15686540, 21500518, 36184769, 91509228, 243334318, 408023356, 815756650, 1234719213, 2332359852 },
{ 1, 3, 6, 10, 27, 55, 79, 166, 481, 928, 1298, 2259, 7962, 11469, 19233, 51347, 110507, 141578, 372962, 918879, 1964656, 2543351, 7360892, 16313375, 21625400, 40816216, 114315392, 214495106, 317325182, 796837402, 1875660340, 2653470793 },
{ 1, 3, 5, 9, 26, 35, 102, 220, 412, 1003, 1288, 3805, 8095, 12270, 25857, 45255, 77244, 156552, 340957, 937819, 1135191, 2814592, 7890176, 10793156, 31796665, 64643969, 126170055, 207507320, 333009457, 931328604, 1817470108, 2344923951 },
{ 1, 2, 7, 8, 29, 45, 95, 226, 395, 1016, 1863, 2089, 6738, 9456, 29107, 56727, 75735, 156959, 342723, 674244, 1659173, 2797739, 7079772, 10225678, 18094615, 52175927, 119828602, 215801504, 441262404, 751726124, 1920804445, 3163449061 },
{ 1, 2, 6, 10, 29, 40, 75, 182, 420, 974, 1834, 2588, 7210, 8781, 32444, 64953, 109030, 177505, 373418, 695694, 1954179, 2245014, 6816165, 8389580, 20973356, 38799894, 130554935, 222306917, 454328055, 630258959, 1787930690, 2927539887 },
{ 1, 2, 7, 15, 28, 52, 97, 158, 327, 795, 2037, 2051, 5125, 10248, 19987, 37416, 104021, 225535, 361945, 741468, 1695470, 3664374, 6815750, 14942221, 31981595, 43253819, 82313341, 247726250, 457965862, 773325701, 1255933618, 2353007384 },
{ 1, 3, 7, 12, 18, 44, 124, 216, 263, 839, 1483, 2659, 4325, 13136, 26607, 32774, 71695, 217109, 330784, 1010286, 1670900, 4182907, 5247391, 9187532, 27550500, 65046830, 94702899, 248473866, 314932054, 671078880, 1354498067, 3924557871 },
{ 1, 3, 4, 15, 25, 57, 114, 154, 480, 955, 1403, 3160, 4292, 9516, 23685, 41937, 109019, 134094, 376298, 926126, 1439071, 3319836, 7345205, 13904495, 29116588, 45123979, 103670040, 207782633, 329276793, 1038655583, 1366142159, 2869036858 },
{ 1, 2, 4, 8, 16, 33, 73, 137, 445, 817, 1670, 3491, 6924, 14070, 28498, 51792, 111802, 132049, 383993, 706978, 1339150, 2673394, 5354842, 10879040, 23190683, 44746648, 130315120, 248050207, 482199103, 988738164, 2001214201, 4026283340 },
{ 1, 3, 4, 15, 30, 53, 92, 176, 399, 822, 1199, 4025, 5486, 9744, 22568, 50285, 129763, 141601, 369292, 978373, 1333670, 3196760, 6946376, 15345303, 17676284, 35310560, 123208146, 237503365, 350749988, 852757632, 1175196127, 2763008412 },
{ 1, 2, 7, 12, 25, 54, 85, 193, 491, 771, 1669, 2828, 5273, 12081, 21209, 56287, 71505, 140872, 344830, 718252, 1591265, 3733268, 7946925, 15398243, 20203028, 64287789, 115688548, 255524504, 500716339, 552392926, 1228476883, 2478976328 },
{ 1, 3, 6, 12, 20, 49, 116, 195, 373, 926, 1209, 2983, 5840, 14155, 17912, 38471, 85170, 203708, 330995, 1017090, 1856347, 3973568, 4786221, 16223832, 25445528, 63215083, 69027961, 222761172, 392502082, 539996646, 1214780014, 4034936035 },
{ 1, 3, 6, 14, 18, 54, 114, 220, 297, 899, 1149, 2765, 4377, 15359, 19635, 33746, 85736, 251230, 332631, 1011016, 1839466, 3423008, 5262748, 14712927, 25250439, 61595557, 90767932, 225662571, 357326064, 621591409, 1438713609, 3291176908 },
{ 1, 2, 6, 9, 26, 34, 88, 143, 370, 943, 1612, 2236, 7956, 9555, 30707, 65228, 112093, 250499, 384357, 723352, 1842725, 2893395, 7368851, 11335001, 17411565, 45612257, 78492566, 210491956, 484819565, 744697056, 1695836564, 3071858738 },
{ 1, 3, 7, 15, 19, 39, 66, 189, 342, 928, 1108, 2193, 4367, 15142, 28020, 49641, 125682, 255996, 366297, 990634, 1577036, 3160237, 5269878, 10537453, 21363450, 40624096, 72744685, 209442767, 336719539, 754446150, 1429833600, 3463652889 },
{ 1, 2, 7, 8, 18, 59, 73, 202, 310, 929, 1644, 2727, 8083, 13877, 22097, 60132, 71500, 253744, 352171, 721017, 1579156, 2637768, 5525188, 13423918, 23943567, 63375894, 100140598, 228854356, 494149355, 607663446, 1739870489, 4070892505 },
{ 1, 3, 6, 11, 25, 54, 116, 142, 318, 983, 1099, 2790, 7051, 9871, 22845, 35793, 124480, 219391, 370109, 1046267, 1843123, 2108655, 8143767, 16035508, 28687195, 48486768, 86792508, 224955346, 365815366, 619403508, 2137105828, 2977167053 },
{ 1, 3, 5, 12, 29, 53, 96, 191, 503, 799, 1073, 3183, 4775, 9166, 21346, 44731, 111608, 164103, 390664, 917522, 1049645, 3935321, 6296258, 14951293, 32789128, 36220305, 106805161, 243437015, 387835224, 1036918979, 1514689406, 3199774349 },
{ 1, 2, 5, 10, 26, 40, 90, 174, 310, 993, 1595, 3695, 7891, 8658, 26187, 57502, 123713, 206623, 368057, 780432, 1050451, 2100520, 8133057, 8924286, 16809699, 35180453, 95525557, 222291228, 472915902, 909124255, 1672765763, 3202149658 },
{ 1, 3, 6, 14, 23, 41, 66, 222, 402, 828, 1033, 2586, 4152, 12901, 31883, 55141, 80620, 172011, 349057, 946974, 1835092, 3408116, 4456918, 10224620, 22283660, 64751887, 121640051, 260585633, 376204577, 614000188, 1355104878, 3791536788 },
{ 1, 3, 4, 15, 20, 40, 82, 199, 502, 846, 1038, 3095, 4140, 15453, 20691, 41438, 84764, 202953, 513505, 863074, 1048659, 3145924, 4194802, 15729473, 20972570, 41946175, 85987454, 208682138, 526405925, 887136912, 1088507666, 3245544670 },
{ 1, 2, 7, 10, 20, 57, 126, 255, 450, 926, 1053, 2088, 7242, 10399, 20752, 58907, 130818, 258312, 466494, 936790, 1049018, 2098025, 7341516, 10488709, 20977725, 59777137, 132140262, 267439589, 471968717, 971185317, 1104607595, 2190318313 },
{ 1, 3, 4, 14, 18, 56, 117, 159, 310, 937, 1086, 3196, 4231, 14610, 19446, 58604, 121248, 159367, 321459, 975901, 1049641, 3148873, 4198628, 14694843, 18893479, 58779618, 122805469, 166879693, 325382716, 983487194, 1139825510, 3354451403 },
{ 1, 2, 6, 11, 18, 61, 112, 135, 429, 956, 1073, 2155, 6318, 11755, 19289, 61734, 117258, 132004, 431132, 960544, 1049682, 2099411, 6297973, 11546331, 18894551, 64027340, 117561061, 141698723, 450283078, 1003386705, 1124172080, 2257777210 },
{ 1, 3, 7, 11, 25, 44, 72, 194, 288, 541, 1122, 3202, 7676, 12058, 25134, 47118, 77846, 211001, 313462, 565439, 1055110, 3154871, 7372701, 11573213, 26301185, 46367238, 75889738, 204182726, 303353132, 568897039, 1178599511, 3365929190 },
{ 1, 3, 5, 10, 25, 55, 112, 148, 505, 661, 1098, 3313, 5452, 11064, 25109, 54673, 121397, 158150, 501498, 720102, 1053035, 3155811, 5273292, 10550421, 26328570, 57848464, 117861440, 155758824, 528991611, 694823752, 1155649153, 3481419880 },
{ 1, 2, 5, 12, 31, 63, 120, 212, 486, 675, 1237, 2532, 5798, 13529, 30203, 60057, 119969, 206127, 488319, 743938, 1078778, 2157211, 5362852, 12789027, 32994144, 66804285, 126907778, 224455247, 506582338, 720578432, 1286043573, 2646299609 },
{ 1, 2, 7, 9, 22, 57, 96, 244, 287, 596, 1278, 2316, 7779, 8321, 20951, 64293, 107073, 230594, 298338, 584328, 1069489, 2161626, 7447375, 9668262, 23366122, 60351305, 101732013, 258005499, 293730169, 634615515, 1317351734, 2472234502 },
{ 1, 3, 4, 14, 21, 34, 67, 200, 499, 583, 1166, 3382, 5029, 16101, 24046, 38521, 73968, 213380, 459498, 566776, 1068639, 3188925, 4291929, 14931771, 22422524, 36431385, 69683320, 213385403, 526459216, 614679332, 1240926149, 3517491821 },
{ 1, 3, 7, 8, 24, 62, 113, 246, 330, 1016, 1178, 3470, 7892, 10170, 26633, 59419, 129081, 235641, 372974, 974196, 1079177, 3198060, 7443652, 8617260, 25450272, 65888647, 117811919, 257881987, 348127344, 1060110581, 1213202765, 3682599920 },
{ 1, 2, 5, 13, 16, 47, 101, 238, 422, 971, 1186, 2343, 5787, 13053, 19988, 42938, 112748, 250097, 411036, 947084, 1072185, 2137152, 5365907, 13845882, 17281591, 50137037, 105358504, 248329535, 438823593, 1011716775, 1264997040, 2424193687 },
{ 1, 3, 6, 11, 18, 47, 100, 134, 488, 784, 1245, 3360, 6679, 10938, 17208, 44209, 108978, 157659, 490972, 863052, 1072157, 3180595, 6419540, 11726036, 19303735, 50317877, 106351303, 140471171, 519534866, 831333952, 1283475048, 3497042436 },
{ 1, 2, 7, 14, 30, 46, 126, 250, 267, 686, 1214, 2484, 7960, 16289, 30235, 43293, 120466, 241911, 318736, 761481, 1076432, 2137470, 7453251, 14888331, 31811412, 48804620, 133499802, 260583004, 282374567, 707461933, 1247926264, 2573447823 },
{ 1, 3, 4, 9, 19, 36, 113, 156, 342, 518, 1370, 3615, 5481, 10819, 22935, 48045, 106166, 183405, 314533, 707869, 1154773, 3329238, 4509062, 10145678, 18979530, 41074923, 123002308, 156953422, 339859298, 571626243, 1391796159, 3903361681 },
{ 1, 3, 5, 11, 31, 61, 83, 156, 456, 682, 1501, 3723, 4534, 8783, 27901, 54544, 75737, 175834, 399660, 538505, 1123907, 3321071, 5643580, 12073891, 31535671, 65175574, 90595364, 158911584, 492861659, 692124022, 1488018188, 3787611064 },
{ 1, 2, 6, 12, 31, 39, 95, 181, 365, 913, 1353, 3018, 7670, 15028, 25183, 41861, 71039, 158635, 301367, 805676, 1118257, 2252911, 6587608, 13398451, 31552046, 38978424, 93687957, 194910524, 389410621, 995000335, 1360304155, 3072084013 },
{ 1, 3, 7, 8, 22, 39, 97, 185, 391, 781, 1463, 3912, 6504, 11902, 17396, 46307, 118122, 149114, 518139, 924925, 1165659, 3295804, 7863075, 9321923, 24220625, 37933190, 107378140, 186437603, 419936455, 849250572, 1560305355, 4288714346 },
{ 1, 2, 4, 14, 21, 44, 75, 154, 363, 655, 1369, 2801, 5557, 13161, 17045, 34146, 90770, 193898, 321165, 566621, 1140479, 2293152, 4519749, 15245022, 21050872, 48382969, 83216357, 155899860, 401926060, 706258762, 1366385353, 2792551888 },
{ 1, 2, 5, 9, 17, 51, 90, 238, 370, 708, 1323, 2605, 4181, 11505, 20828, 60085, 76221, 229336, 308055, 557698, 1128940, 2315042, 5562911, 9988109, 16781338, 51391527, 99635296, 242280613, 372321711, 796096492, 1460971265, 2889393780 },
{ 1, 3, 7, 14, 22, 45, 109, 144, 258, 669, 1404, 3626, 6156, 13330, 17444, 35957, 126123, 193858, 376416, 536814, 1177013, 3326957, 7733060, 15249922, 24117355, 48234653, 111149331, 165675710, 289408263, 721424023, 1360009571, 3953145375 },
{ 1, 3, 7, 8, 29, 39, 97, 177, 466, 881, 1438, 4024, 6357, 11631, 27313, 47837, 129661, 151450, 456888, 1000908, 1173329, 3304951, 7768852, 9425184, 29389439, 37787550, 107057335, 194180569, 485950315, 896338353, 1608981444, 4037878851 },
{ 1, 2, 6, 13, 17, 61, 93, 214, 477, 740, 1292, 2870, 8191, 14858, 20613, 53534, 73487, 252841, 442048, 613748, 1134490, 2296456, 6791572, 14321158, 17297558, 66737445, 94825311, 228165486, 481316640, 761326541, 1415663168, 3171092579 },
{ 1, 3, 6, 10, 26, 42, 106, 204, 298, 807, 1512, 3608, 7923, 9120, 30845, 38133, 127350, 252833, 377982, 958707, 1149308, 3362747, 6566996, 11368601, 28592560, 43901585, 114665331, 207493489, 326115240, 809508962, 1507884243, 4013988108 },
{ 1, 2, 5, 12, 30, 59, 83, 194, 455, 569, 1309, 3022, 4960, 15928, 24863, 50123, 71532, 238118, 449828, 680856, 1144750, 2319329, 5754653, 13224581, 33011808, 60374145, 91675941, 214590362, 500266923, 549676013, 1431817987, 3015297726 },
{ 1, 3, 5, 9, 16, 58, 92, 244, 309, 809, 1530, 3651, 4642, 11909, 24408, 53607, 69570, 216297, 337173, 965466, 1135969, 3380174, 5519600, 10344767, 18100028, 60684745, 92561935, 267192044, 303379436, 832477317, 1528896936, 3899876008 },
{ 1, 3, 7, 14, 20, 35, 92, 217, 389, 713, 1285, 3999, 6740, 12394, 16549, 43271, 72603, 245341, 516208, 598162, 1130872, 3352350, 7781121, 15349564, 22526814, 34192354, 94379733, 222313791, 425746419, 775984892, 1422999934, 4132644627 },
{ 1, 3, 6, 14, 20, 37, 127, 188, 321, 717, 1440, 3906, 7810, 12613, 18118, 45500, 105341, 142028, 286115, 556868, 1164940, 3319121, 6659811, 15390147, 22300609, 40387469, 125934446, 185740004, 355753422, 778600403, 1554106278, 4253198085 },
{ 1, 3, 7, 14, 17, 58, 87, 194, 258, 575, 1465, 3992, 6665, 12795, 22377, 52035, 76589, 242623, 365169, 704861, 1139425, 3359895, 7664068, 15304500, 17143705, 60476938, 85031420, 217282407, 289710930, 549012247, 1608820712, 4130950835 },
{ 1, 2, 6, 11, 26, 57, 84, 243, 444, 562, 1328, 2830, 7725, 8454, 31565, 55030, 71885, 217569, 488157, 660670, 1145202, 2338775, 6687712, 12172193, 28799868, 62527164, 87129098, 262382616, 436606015, 537499743, 1475808489, 3161115013 },
{ 1, 2, 6, 9, 23, 48, 125, 244, 462, 584, 1628, 3688, 5658, 13040, 28441, 44316, 90057, 163995, 371975, 834531, 1420494, 2931104, 7711363, 12367857, 17153267, 62686661, 114676313, 217852517, 439650818, 679388894, 2146514755, 3471974821 },
{ 1, 3, 4, 9, 16, 46, 87, 168, 315, 770, 1815, 2867, 7024, 16383, 32387, 50439, 106364, 210917, 509625, 784759, 1558419, 3930712, 5750999, 11268554, 22534796, 38522137, 69706575, 147522437, 395716216, 954730653, 1726060885, 2345877470 },
{ 1, 3, 4, 11, 21, 57, 95, 203, 475, 813, 1809, 2886, 7068, 13941, 24850, 33521, 128051, 180297, 367862, 606607, 1416179, 3754739, 5603380, 9010246, 17162472, 52034979, 83104661, 262878823, 423728420, 1030898352, 1749313748, 2301329908 },
{ 1, 2, 7, 12, 17, 63, 103, 154, 280, 643, 1705, 3832, 4639, 9151, 25709, 40069, 67889, 229079, 469573, 879429, 1519075, 2975503, 6764836, 15551201, 24596024, 52295677, 117843122, 180230478, 377492072, 578822936, 2080400661, 3441467024 },
{ 1, 3, 7, 11, 26, 33, 83, 214, 504, 584, 1561, 2763, 5951, 14694, 22340, 55805, 98887, 188943, 349947, 830295, 1393090, 3969593, 6624923, 9194478, 29660309, 43952406, 82953157, 267518175, 459572710, 726571621, 2056356443, 2272453158 },
{ 1, 2, 5, 10, 22, 39, 126, 225, 347, 564, 1606, 3772, 7129, 12712, 28480, 64543, 72758, 175190, 500867, 1022347, 1556278, 3123437, 4706624, 9393674, 19319299, 50198034, 113340986, 216158808, 325500552, 797784970, 2000683297, 3300917981 },
{ 1, 2, 6, 14, 27, 63, 93, 244, 267, 798, 1877, 3969, 5856, 11586, 18380, 46650, 118840, 140369, 414953, 697658, 1459032, 2788254, 7733975, 13250826, 29868828, 54015827, 68274063, 208800507, 347477373, 1035607953, 1637319374, 3311963443 },
{ 1, 2, 6, 8, 29, 46, 64, 253, 280, 928, 2010, 3926, 5851, 13526, 17748, 63310, 125689, 144517, 408978, 536101, 1453514, 2627322, 7764097, 11065756, 27772464, 38412793, 94388884, 227604540, 338820217, 934412424, 1777760643, 3260602904 },
{ 1, 2, 6, 15, 19, 59, 72, 166, 438, 1001, 1923, 3918, 5687, 10401, 25019, 34812, 113591, 202517, 367274, 716127, 1411828, 2804169, 7731984, 14314145, 21246277, 57269971, 96493993, 151029701, 514964473, 934483900, 1726323471, 3282759309 },
{ 1, 3, 4, 14, 28, 50, 102, 223, 291, 818, 2018, 2782, 6618, 9808, 23629, 35974, 78210, 136946, 400777, 786147, 1458595, 3895993, 5730566, 14247799, 26307441, 66229112, 132592494, 238667598, 360718086, 1036005261, 1781547543, 2479923417 },
{ 1, 3, 6, 12, 19, 59, 68, 210, 277, 768, 2045, 2712, 4379, 10006, 25548, 39619, 106977, 168572, 420954, 753913, 1444204, 3871642, 7847465, 16362678, 21078434, 54694561, 85355866, 240878542, 339087046, 1051400683, 1790422627, 2488895602 },
{ 1, 2, 4, 11, 19, 34, 65, 180, 503, 742, 1554, 3923, 7238, 13497, 29163, 60124, 116322, 217988, 267588, 959260, 1341672, 3046779, 5623653, 10339368, 23075920, 45102226, 90206653, 164686401, 441567175, 654529012, 2051281632, 3491671581 },
{ 1, 2, 7, 14, 27, 35, 99, 133, 342, 941, 1875, 3750, 4151, 11353, 18626, 55742, 93797, 201099, 444987, 647495, 1479582, 2692865, 6795894, 13291967, 30410343, 49286540, 133174837, 176169308, 316694461, 816905058, 1717647091, 3638813929 },
{ 1, 3, 5, 8, 29, 40, 106, 148, 446, 619, 1671, 2966, 7185, 14398, 16474, 62686, 97629, 220009, 329016, 855022, 1398974, 4048338, 4481778, 12477228, 26214856, 36700884, 119539542, 188746081, 479206197, 762327531, 1966103206, 2526070754 },
{ 1, 3, 6, 14, 26, 47, 114, 234, 509, 779, 1967, 2802, 4271, 11599, 19046, 59815, 70558, 180900, 369701, 527459, 1439959, 3719604, 7649190, 14339813, 30432412, 35704092, 115429092, 212995231, 424996122, 1021971178, 1645282437, 2432886069 },
{ 1, 2, 6, 15, 22, 56, 92, 224, 319, 871, 1611, 3359, 4887, 8905, 26771, 43449, 112354, 168188, 511258, 584476, 1473238, 2783410, 7451083, 13803110, 17291647, 58247074, 79038223, 229310181, 373325041, 950043914, 1984023339, 3841070748 },
{ 1, 3, 7, 8, 20, 62, 122, 185, 379, 702, 2015, 2368, 4808, 15229, 30745, 47142, 76895, 228577, 500114, 919352, 1477849, 4028899, 6368146, 11761044, 20427571, 57556170, 116848085, 139314172, 275775763, 549454386, 1802503885, 3162508145 },
{ 1, 2, 6, 14, 16, 53, 72, 201, 412, 745, 1857, 3554, 4716, 9756, 30345, 33761, 121026, 159110, 385732, 810786, 1326406, 3030852, 7458280, 12732020, 23430711, 62716644, 93688661, 267178463, 490075706, 583825136, 1822711656, 4079916425 },
{ 1, 3, 6, 11, 20, 32, 119, 203, 455, 607, 1958, 2489, 5765, 12876, 30606, 49628, 91749, 168936, 312634, 1014730, 1418568, 4055820, 7434384, 8560999, 19182433, 46119010, 103073000, 244283830, 460694175, 701360757, 1730242502, 3168964948 },
{ 1, 2, 6, 10, 22, 50, 96, 207, 420, 843, 1674, 3328, 4660, 13393, 30891, 35175, 90867, 188901, 356312, 664475, 1318670, 2684437, 7407629, 8543262, 17069090, 59263064, 124785849, 267401561, 471886479, 954252558, 1937840686, 3885131903 },
{ 1, 3, 7, 14, 16, 50, 109, 159, 333, 713, 2017, 2474, 4899, 8209, 29745, 37994, 83089, 231773, 492283, 1024908, 1496373, 3934830, 6365912, 13844432, 24556992, 63748018, 130023756, 182452938, 309331942, 629148068, 1807749939, 3037732899 },
{ 1, 2, 4, 9, 16, 38, 91, 163, 465, 870, 1785, 3345, 6774, 12429, 28056, 52199, 120698, 241348, 347512, 757409, 1343924, 2626474, 5351410, 11754327, 21293707, 49862132, 77622075, 141608535, 417451225, 1019452710, 1938105865, 3843770481 },
{ 1, 3, 4, 15, 29, 39, 110, 225, 506, 840, 1743, 2399, 7719, 10416, 17758, 64036, 72884, 148817, 320057, 739475, 1433919, 3745496, 5325161, 12754551, 26584087, 42589238, 119561296, 206624935, 419512680, 979541620, 1921360915, 2735332409 },
{ 1, 2, 5, 15, 23, 42, 114, 191, 509, 746, 1868, 3436, 8040, 11531, 27640, 58507, 87482, 228880, 367182, 902795, 1365972, 3053809, 4274456, 13866975, 19189995, 39752992, 104933256, 158578737, 413449288, 632192237, 1785752872, 4209044373 },
{ 1, 3, 6, 10, 25, 45, 94, 233, 454, 687, 1806, 2347, 5990, 13748, 16923, 52775, 126536, 174814, 275388, 819208, 1428508, 4068385, 7449677, 9594077, 32890293, 41822744, 73530913, 208836162, 434398919, 617380753, 1829003350, 3131863285 },
{ 1, 3, 5, 11, 18, 47, 107, 177, 354, 512, 2012, 2212, 6464, 14967, 25435, 64912, 68533, 218131, 468012, 878702, 1461434, 4112752, 4219439, 8454071, 24185877, 37968930, 127351927, 166547591, 277227828, 724449013, 1805743716, 2864946039 },
{ 1, 3, 4, 13, 27, 38, 85, 150, 345, 815, 1785, 2222, 7460, 10227, 24349, 65174, 119874, 205992, 488750, 596961, 1507117, 3851006, 5265575, 14742834, 32624590, 45288302, 82261589, 185188749, 282555102, 1000935672, 2052515254, 3103452839 },
{ 1, 3, 6, 13, 25, 56, 103, 251, 373, 848, 1790, 2183, 5555, 11849, 21825, 35629, 90683, 150936, 412161, 548346, 1538801, 3838107, 7356800, 14722618, 32584091, 55758343, 132595191, 205143784, 288449699, 981616103, 2027314881, 3181932782 },
{ 1, 2, 7, 9, 21, 57, 111, 168, 401, 999, 1710, 3161, 5337, 14610, 29225, 46491, 87031, 198297, 314410, 628829, 1439957, 2756873, 6296076, 11548120, 18900745, 50368463, 128001780, 152288389, 535137733, 888843047, 1944146846, 4223862330 },
{ 1, 2, 7, 11, 30, 37, 122, 207, 489, 694, 2002, 3096, 5164, 12387, 21729, 60813, 97884, 184897, 300641, 880151, 1421005, 3108639, 6291958, 9437841, 26216367, 48237788, 104863195, 245379824, 422597449, 663806298, 1672571801, 3940737172 },
{ 1, 2, 4, 14, 28, 48, 120, 198, 307, 973, 1635, 3085, 6170, 10298, 18538, 45290, 74107, 195443, 512918, 784115, 1484212, 2726484, 5243005, 12583114, 25166123, 65012723, 104859143, 257953019, 347085137, 816851761, 1962954554, 4227907370 },
{ 1, 3, 6, 15, 27, 35, 71, 222, 284, 776, 1372, 2055, 6156, 12317, 30764, 55388, 71933, 145755, 455638, 582720, 1592079, 2817360, 4194330, 12582944, 25165889, 62914769, 113246471, 146801451, 297796891, 931137753, 1191188752, 3254792981 },
{ 1, 3, 7, 13, 24, 58, 117, 149, 476, 798, 1816, 2077, 6195, 14439, 26810, 49547, 119748, 241188, 303690, 982505, 1623164, 3719052, 4194754, 12583722, 29362034, 54528191, 100669826, 243284950, 490761739, 625003037, 1996606771, 3347291456 },
{ 1, 2, 5, 15, 28, 34, 91, 167, 410, 684, 1338, 2094, 4160, 10383, 31186, 57909, 70873, 184958, 339015, 848791, 1379844, 2693038, 4195685, 8390788, 20976067, 62925326, 117472431, 142666398, 381748568, 700643578, 1720013109, 2869697312 },
{ 1, 2, 5, 11, 26, 34, 76, 225, 403, 721, 1229, 2093, 4191, 10460, 23030, 53862, 71069, 158499, 466594, 834991, 1493845, 2524712, 4195482, 8390892, 20976005, 46148327, 109075634, 142661816, 318832963, 943868723, 1690785432, 3024897476 },
{ 1, 3, 4, 8, 16, 47, 97, 235, 378, 567, 1753, 2151, 6372, 8550, 16896, 34439, 94402, 200980, 490176, 759700, 1126926, 3649176, 4200689, 12591426, 16794229, 33588811, 67203529, 197333991, 407338242, 986422575, 1586572617, 2381809070 },
{ 1, 2, 7, 11, 17, 33, 72, 187, 363, 620, 1931, 2154, 4342, 14812, 23325, 36311, 67464, 153711, 389370, 735686, 1291053, 4054462, 4200315, 8401343, 29379581, 46183457, 71392460, 138597666, 302315089, 784950142, 1523991863, 2604068478 },
{ 1, 2, 6, 13, 27, 48, 120, 186, 372, 539, 1217, 2127, 4301, 12755, 27474, 57017, 101380, 252302, 390079, 786228, 1075135, 2576049, 4199446, 8403374, 25197537, 54587357, 113370681, 201524508, 503610114, 780676687, 1561691577, 2263646797 },
{ 1, 3, 7, 14, 19, 37, 79, 196, 349, 890, 1667, 2168, 6317, 14801, 29421, 40119, 75172, 159257, 395712, 703180, 1770737, 3528061, 4202354, 12595920, 29387972, 58764560, 79777783, 155372378, 331833030, 822626470, 1465657800, 3736928978 },
{ 1, 3, 6, 9, 26, 49, 90, 179, 286, 576, 1760, 2170, 6343, 12745, 19364, 54566, 99759, 181492, 376377, 563598, 1225516, 3548286, 4202340, 12598090, 25187457, 37809901, 109170796, 205662435, 377776561, 751483750, 1201063136, 2419893784 },
{ 1, 3, 4, 14, 24, 51, 112, 141, 344, 975, 1232, 2182, 6477, 9184, 29853, 51299, 98475, 237831, 269071, 656749, 2091695, 2432890, 4224244, 12634334, 16875911, 58958411, 100933279, 214568471, 471856962, 593835153, 1447041148, 4076871826 },
{ 1, 2, 6, 10, 27, 32, 71, 231, 425, 957, 1494, 2299, 4481, 13287, 21787, 53518, 70395, 135135, 485050, 888048, 1884762, 2948755, 4215976, 8442513, 25237356, 42076105, 113735305, 135117962, 299668176, 971777921, 1778472571, 4022529300 },
{ 1, 3, 6, 15, 19, 46, 80, 255, 418, 891, 1983, 2279, 6549, 13087, 32546, 37335, 93133, 177878, 494172, 890307, 1742180, 4133926, 4226756, 12620401, 25259413, 63094164, 80188823, 193834385, 337311134, 1073724813, 1757509027, 3724709363 },
{ 1, 3, 4, 8, 22, 36, 91, 131, 395, 1000, 1088, 2222, 6597, 9035, 17816, 48083, 80938, 192598, 278686, 838077, 1975199, 2124942, 4211074, 12628989, 16852064, 33739005, 92537168, 151812836, 383714859, 551730960, 1652835599, 4207723131 },
{ 1, 2, 7, 9, 24, 56, 92, 148, 450, 907, 1590, 2229, 4516, 15172, 20329, 51940, 119841, 185755, 325356, 979398, 1964935, 3401248, 4212892, 8438232, 29475764, 37938787, 100968505, 235798622, 387735699, 624032203, 1891986323, 3796690446 },
{ 1, 2, 7, 9, 19, 35, 111, 155, 284, 711, 1880, 2297, 4504, 15333, 19913, 38757, 65706, 217426, 305716, 538326, 1553201, 3716132, 4217899, 8432691, 29445124, 37995585, 80063654, 147549442, 466816694, 653308850, 1196930319, 2977480361 },
{ 1, 2, 7, 12, 29, 60, 106, 161, 310, 906, 1294, 2263, 4521, 15057, 26153, 58610, 129736, 228892, 351374, 598596, 1931105, 2791411, 4224695, 8445355, 29465064, 50565685, 121944050, 252339893, 446188972, 677597668, 1306327592, 3794352078 },
{ 1, 2, 6, 12, 22, 43, 126, 222, 375, 909, 1717, 2186, 4565, 12834, 25667, 48548, 83905, 252675, 440170, 742066, 1963435, 3641299, 4223778, 8423182, 25244241, 50549897, 92771582, 180892734, 530187646, 935265980, 1570451604, 3824966114 },
{ 1, 2, 6, 15, 28, 46, 105, 238, 346, 759, 1912, 2220, 4566, 13142, 31976, 60792, 93867, 222193, 470315, 685587, 1460760, 3751456, 4221533, 8443553, 25232358, 63159568, 117849671, 193591981, 441454547, 1001398647, 1456544410, 3178473375 },
{ 1, 2, 4, 13, 27, 38, 91, 233, 417, 792, 1389, 2202, 4366, 8831, 28569, 54843, 77093, 188481, 463053, 891390, 1655804, 2645207, 4225956, 8453176, 16863756, 54662966, 113661788, 160409553, 383741664, 979829990, 1746813888, 3332760744 },
{ 1, 3, 7, 15, 20, 40, 86, 137, 416, 685, 1501, 2277, 6491, 15297, 32382, 44663, 87663, 171604, 269861, 915182, 1368047, 3009716, 4219725, 12623829, 29439174, 63046529, 84225755, 168820617, 361884752, 576856597, 1750199315, 2879529655 },
{ 1, 3, 7, 8, 18, 53, 113, 167, 370, 1001, 1703, 2428, 7152, 16007, 18736, 35636, 104295, 246746, 374493, 657870, 1874623, 3190029, 4298501, 12829515, 29734873, 34212488, 77372792, 225486765, 469862947, 688124951, 1564861547, 4165647570 },
{ 1, 3, 6, 11, 19, 32, 80, 216, 471, 746, 2032, 2502, 6863, 14253, 20742, 33579, 79159, 185982, 412642, 1021125, 1422537, 3773650, 4273383, 12768412, 25578503, 47160616, 81112645, 137995143, 339829862, 918705541, 1950842647, 3092070522 },
{ 1, 2, 4, 8, 29, 41, 72, 187, 353, 519, 1355, 2381, 4677, 9705, 18445, 63511, 92208, 131177, 340206, 672243, 1211182, 3114999, 4287847, 8521999, 17113831, 34235561, 122827359, 175078862, 306184276, 775946384, 1497366829, 2210398900 },
{ 1, 2, 6, 12, 31, 46, 102, 200, 473, 593, 1135, 2554, 4650, 13455, 26697, 59564, 82199, 217988, 441889, 1029563, 1156022, 2432621, 4277535, 8608436, 25603540, 51348329, 131190671, 195428117, 432019012, 847265143, 2009102949, 2537610287 },
{ 1, 3, 7, 10, 20, 59, 68, 228, 439, 997, 1277, 2448, 7050, 15467, 22681, 47379, 123579, 161457, 495271, 816792, 1910481, 2143787, 4319155, 12745769, 29855860, 42746039, 85811567, 249575964, 289425384, 968920289, 1812046267, 4223857656 },
{ 1, 2, 7, 8, 16, 49, 122, 216, 347, 838, 1587, 2420, 4887, 16055, 18612, 37304, 111216, 234900, 413302, 739788, 1690290, 3613929, 4306291, 8625997, 29776500, 34307481, 68785769, 209158629, 515922665, 914402426, 1459747978, 3548581256 },
{ 1, 2, 7, 12, 22, 58, 106, 225, 453, 836, 1292, 2473, 5035, 15576, 27355, 42715, 124330, 207071, 496038, 1024952, 1680617, 3015339, 4322838, 8590399, 29883152, 51316110, 93864080, 246229258, 449209046, 952942287, 1889422743, 3525392569 },
{ 1, 3, 7, 12, 18, 38, 109, 152, 363, 876, 1026, 2315, 7141, 15688, 27142, 36717, 66830, 203490, 292543, 671592, 1633542, 2448123, 4267655, 12789553, 29629909, 51032933, 77196312, 161892660, 453342128, 650716553, 1528601534, 3626984863 },
{ 1, 2, 6, 10, 22, 46, 119, 147, 337, 740, 1791, 2353, 4699, 14303, 23372, 42467, 85028, 259999, 263070, 756834, 1328967, 3336713, 4284207, 8635046, 25452075, 42672314, 93657390, 196160103, 495257482, 608758708, 1439781913, 3149545416 },
{ 1, 3, 5, 12, 23, 39, 88, 136, 458, 797, 1793, 2443, 7075, 11965, 27229, 42448, 73118, 177458, 312485, 997180, 1755616, 4037457, 4261803, 12752917, 21276705, 51353681, 98289811, 167709178, 373666658, 583870417, 1943998665, 3295680884 },
{ 1, 3, 6, 10, 23, 42, 125, 216, 510, 917, 1789, 2492, 6943, 14250, 23314, 42932, 92969, 239593, 410523, 978661, 2081165, 3158869, 4290354, 12819043, 25560442, 42891111, 98495171, 179504208, 520488209, 919500694, 2115955526, 3891463959 },
{ 1, 3, 7, 15, 18, 40, 84, 226, 376, 527, 1129, 2320, 6813, 15794, 29346, 36293, 94734, 181354, 504087, 713362, 1209760, 2452106, 4287889, 12776172, 29875474, 63593240, 76771067, 170335408, 357029911, 944394275, 1607692366, 2150050007 },
{ 1, 3, 5, 11, 30, 39, 66, 186, 423, 795, 2008, 2512, 7127, 11952, 21414, 61042, 72411, 152637, 329499, 805038, 1706397, 3824507, 4261675, 12742867, 21317964, 46957105, 127622420, 167298514, 281367719, 776775411, 1754993439, 3292008400 },
{ 1, 2, 4, 15, 17, 53, 127, 253, 388, 564, 1121, 2426, 5046, 9822, 29957, 37608, 105890, 240601, 462484, 863482, 1282094, 2596950, 4319394, 8599868, 17208099, 63866649, 72654548, 225086564, 529693047, 1071946659, 1615298148, 2387457387 },
{ 1, 2, 4, 15, 29, 56, 84, 227, 420, 658, 1683, 2383, 4898, 9263, 30626, 62207, 128550, 192524, 497691, 796723, 1542214, 3262662, 4293064, 8520229, 17199060, 63887993, 122689315, 238939123, 357965263, 947119602, 1749519787, 2815165777 },
{ 1, 3, 7, 11, 30, 57, 84, 169, 260, 669, 2020, 2461, 7127, 15845, 22332, 57433, 125119, 180516, 385787, 597767, 1542153, 3723290, 4306997, 12724289, 29626510, 46915945, 127023712, 242296393, 358072836, 716703406, 1082137500, 2772450678 },
{ 1, 3, 7, 13, 29, 44, 83, 154, 337, 897, 1818, 2540, 6800, 15415, 26581, 63529, 82010, 182413, 282989, 721891, 2095103, 4100212, 4268250, 12728829, 29690530, 55178319, 123455249, 188047871, 346000038, 645813317, 1438648065, 3821013454 },
{ 1, 3, 7, 15, 25, 54, 85, 211, 315, 536, 1847, 2505, 7002, 15768, 30540, 55674, 103136, 192103, 402241, 572842, 1306385, 4055487, 4318196, 12729372, 29863152, 63595915, 106276712, 230011144, 353628742, 876468213, 1329694930, 2231554360 },
{ 1, 3, 4, 12, 23, 38, 84, 217, 416, 1012, 1304, 2347, 6790, 10067, 28666, 44752, 66645, 193282, 404724, 836147, 1898048, 2952470, 4288827, 12745384, 17254172, 51289874, 98295554, 161693484, 358330211, 916701067, 1767562841, 4276208928 },
{ 1, 2, 7, 10, 27, 52, 69, 166, 295, 789, 1867, 2516, 4793, 15631, 23907, 50620, 99428, 155342, 367593, 540401, 1809333, 3964516, 4311773, 8562124, 29639813, 42643272, 115315153, 221653684, 287278366, 691412300, 1214238157, 3346459783 },
{ 1, 3, 4, 10, 19, 40, 65, 131, 296, 945, 1414, 2503, 6738, 10123, 24314, 36336, 94515, 156557, 325089, 573724, 1983426, 2835835, 4272143, 12763162, 17154111, 42600547, 81500369, 169961859, 270201560, 539176628, 1264475497, 3995822116 },
{ 1, 3, 4, 10, 26, 60, 71, 178, 435, 949, 1873, 2406, 6850, 9654, 23727, 50834, 119359, 152746, 343707, 846369, 2006172, 3980998, 4297391, 12762472, 17083841, 42882125, 110770107, 254992466, 295418013, 745566678, 1839158134, 4013868557 },
{ 1, 2, 6, 9, 26, 43, 119, 213, 472, 594, 1532, 2643, 5630, 14933, 24055, 59983, 67036, 217656, 503049, 889824, 1392475, 2258460, 4697352, 9278434, 26558301, 40007189, 113749266, 172856265, 475348778, 925005504, 1892134090, 2664272283 },
{ 1, 2, 4, 8, 26, 57, 71, 209, 330, 874, 1618, 2925, 5724, 11131, 22135, 64287, 106200, 180931, 463913, 562865, 1854076, 3808722, 4658230, 8951426, 18631205, 37362994, 113710369, 248026368, 282873210, 913970660, 1413940325, 3561535100 },
{ 1, 3, 7, 13, 18, 34, 83, 176, 276, 650, 1307, 2714, 7485, 12995, 30098, 41959, 92009, 146975, 322647, 688316, 1083653, 2237103, 4517189, 13271608, 30444555, 56764442, 80015417, 139112553, 353433822, 794986951, 1103467359, 2865268326 },
{ 1, 3, 4, 8, 21, 53, 103, 144, 284, 839, 1290, 2903, 7476, 11041, 23961, 33359, 128635, 245310, 360187, 739100, 2057674, 2314939, 4554739, 13321325, 18833571, 35871092, 84239310, 227241022, 417282135, 639881468, 1107593629, 3695770214 },
{ 1, 2, 7, 12, 18, 38, 87, 163, 498, 942, 1756, 2949, 5787, 13063, 30514, 40984, 81978, 135268, 280798, 899348, 1661657, 3919063, 4476500, 9286192, 31018423, 54262898, 79970216, 151853451, 337124161, 737796018, 2017460662, 3800040233 },
{ 1, 3, 4, 9, 26, 42, 78, 243, 335, 665, 1523, 2783, 7449, 11198, 22471, 65188, 68733, 168298, 441027, 570685, 1217511, 2494221, 4630421, 13147763, 18007403, 40301369, 113749505, 172623340, 311714559, 1051618633, 1427467118, 2888627932 },
{ 1, 2, 7, 14, 27, 59, 70, 251, 308, 680, 1661, 2772, 5826, 13084, 25936, 61394, 111042, 181122, 412903, 734458, 1279181, 3769495, 4611161, 9112880, 30668632, 62478423, 109545680, 239673691, 273943099, 1030795220, 1267876302, 2759562142 },
{ 1, 2, 5, 12, 25, 37, 90, 204, 507, 850, 1685, 2853, 5652, 8747, 29902, 61283, 86124, 139425, 459096, 893455, 1905814, 3704746, 4657563, 9276394, 22902760, 54020464, 101069390, 164580474, 400242183, 808968732, 2027701409, 3715803077 },
{ 1, 3, 4, 11, 16, 33, 87, 241, 313, 884, 1121, 2840, 7413, 10914, 17941, 44143, 82688, 132288, 455391, 767651, 1936774, 2701594, 4647716, 13341855, 18731578, 48871344, 71825077, 143203755, 349845842, 1050706927, 1243450816, 3584637385 },
{ 1, 2, 7, 14, 20, 52, 65, 134, 317, 634, 1501, 2577, 5425, 13250, 26271, 38335, 117479, 165058, 380953, 729125, 1497212, 2414813, 4571547, 9116466, 30846816, 61109879, 88374375, 227094868, 294830737, 573406208, 1243007403, 2468668280 },
{ 1, 3, 4, 8, 21, 49, 125, 130, 396, 928, 1881, 2855, 7902, 10390, 22974, 35801, 118739, 215742, 361807, 1025970, 1655323, 3554443, 4578533, 13577305, 18386353, 37058239, 83908659, 209750138, 507627662, 579027345, 1715837828, 3830419221 },
{ 1, 2, 5, 9, 24, 41, 125, 229, 412, 598, 1590, 2730, 6017, 8324, 22875, 60330, 66725, 206181, 425905, 997586, 1436102, 4059866, 4597623, 9425196, 22342480, 41868572, 104877919, 180406611, 503405488, 990102665, 1627876625, 2680962892 },
{ 1, 3, 6, 11, 16, 41, 108, 189, 276, 917, 1286, 2868, 7203, 14578, 16854, 45587, 71648, 257290, 308012, 781340, 1646780, 2117911, 4518803, 13327629, 26751780, 48477194, 71325854, 176193899, 444680967, 763585653, 1090904076, 4010326040 },
{ 1, 2, 6, 13, 16, 47, 99, 151, 363, 591, 1120, 2802, 5473, 14385, 30808, 45264, 70018, 236428, 378632, 543785, 1497615, 2653318, 4569895, 8975915, 26730137, 57401230, 71304559, 188745765, 423628924, 645935278, 1589668138, 2550170283 },
{ 1, 2, 4, 14, 19, 63, 99, 165, 350, 611, 1905, 2804, 5721, 10247, 24584, 47129, 100386, 229455, 315641, 649624, 1455000, 3243084, 4503526, 9001948, 18328746, 62051926, 75497494, 255852595, 432013428, 717226132, 1396703534, 2533360377 },
{ 1, 2, 5, 13, 26, 53, 118, 130, 452, 866, 1092, 3020, 5572, 8274, 30929, 63767, 115312, 222777, 374255, 886814, 2046010, 2900073, 4573370, 9271720, 22997973, 57609718, 113250348, 230697028, 482371836, 599839057, 2004986505, 3586373477 },
{ 1, 2, 5, 10, 18, 48, 83, 223, 397, 890, 1683, 2947, 6005, 8290, 16526, 47447, 119549, 134651, 488736, 978085, 1908680, 3707233, 4666019, 9329977, 22862981, 45791150, 79696359, 209726435, 360732670, 893425983, 1740739169, 3531781276 },
{ 1, 2, 4, 15, 21, 63, 85, 198, 280, 733, 1893, 2654, 5769, 10719, 27591, 36284, 101684, 138434, 522089, 666189, 1064637, 3336592, 4688683, 9080014, 17948538, 66089593, 84010738, 256020860, 340163161, 843593344, 1129826753, 3099116514 },
{ 1, 2, 5, 11, 27, 36, 79, 212, 510, 845, 1060, 3008, 5437, 8461, 19195, 63411, 95754, 179599, 503935, 911533, 1859935, 3109383, 4714007, 9256359, 22894629, 49098829, 109131985, 159531509, 310825814, 936342528, 2036046735, 3663883753 },
{ 1, 3, 4, 10, 27, 50, 68, 251, 356, 881, 1805, 3063, 7886, 10543, 19352, 65092, 116986, 188775, 420725, 612103, 1868780, 3327740, 4602219, 13208419, 18685728, 45151115, 109172330, 205703298, 268837321, 1011460803, 1608338788, 3505401647 },
{ 1, 2, 5, 12, 22, 55, 95, 223, 390, 913, 1179, 2860, 5504, 8652, 29566, 40267, 118904, 153788, 485657, 961253, 1754866, 2710673, 4655921, 9313704, 22770110, 53197791, 96584826, 222460655, 377954031, 885931193, 1679491907, 3920307465 },
{ 1, 3, 6, 15, 22, 54, 70, 150, 323, 847, 1658, 3006, 8156, 14692, 25395, 34490, 115282, 171392, 343281, 542706, 1988865, 3681235, 4540198, 13168862, 27255232, 66716220, 96585084, 222467396, 268772631, 642282929, 1432233147, 3732445044 },
{ 1, 3, 5, 15, 29, 53, 98, 214, 507, 893, 1528, 3052, 7321, 8673, 25409, 50573, 121629, 255272, 473629, 831478, 1938481, 2218234, 4669856, 13386705, 22820079, 65196309, 117562111, 226748012, 432488351, 910992630, 2040367404, 3525442181 },
{ 1, 3, 5, 13, 21, 48, 68, 228, 287, 934, 1586, 2747, 7570, 9345, 32561, 37174, 121796, 183975, 420617, 782298, 1675322, 3938394, 4698313, 13181299, 23067427, 58193853, 85860514, 217397757, 270027441, 985585036, 1132584108, 3858941789 },
{ 1, 2, 7, 11, 30, 62, 71, 136, 464, 894, 1695, 2724, 5629, 13314, 18197, 53609, 99065, 175429, 333191, 894152, 1859137, 3963858, 4495160, 9373959, 31017489, 49405039, 131649528, 255031375, 274581651, 545186274, 1912718124, 3527570999 },
{ 1, 2, 7, 15, 24, 60, 124, 153, 459, 968, 1916, 2586, 5258, 13597, 26118, 57412, 110831, 221506, 365096, 826604, 1607103, 3313598, 4465421, 9436772, 31265004, 66984263, 106738208, 264062203, 513365403, 623970302, 1967133672, 4265623350 },
{ 1, 2, 7, 9, 25, 38, 116, 213, 319, 871, 1311, 2660, 5679, 14292, 21937, 60316, 87119, 201207, 461426, 761371, 2049944, 2415935, 4610587, 8947428, 30957252, 40861371, 102325872, 153913184, 480183402, 933495174, 1208441513, 3532357428 },
{ 1, 2, 5, 9, 21, 54, 79, 251, 322, 936, 1258, 2758, 5920, 9911, 21562, 39914, 115727, 179106, 433400, 611068, 1642355, 2649711, 4676865, 9167094, 23036252, 39916430, 85857434, 220412519, 316235439, 1009318930, 1422297987, 3846796420 },
{ 1, 2, 5, 10, 18, 50, 107, 198, 505, 765, 1389, 2926, 5972, 10239, 20173, 42452, 127598, 238630, 475267, 905596, 1289100, 2100988, 4629944, 9327271, 22324679, 44923486, 81067080, 204339279, 466960535, 867105091, 2077488122, 2984283664 },
{ 1, 2, 6, 12, 25, 38, 96, 166, 426, 621, 1768, 2984, 5182, 15558, 32142, 59017, 97036, 232850, 300711, 1037183, 1344775, 4018002, 4585970, 9200259, 26242971, 53639395, 101712232, 154315627, 433132858, 707917586, 1830064440, 2440917910 },
{ 1, 2, 6, 8, 22, 58, 85, 186, 393, 556, 1312, 3055, 6040, 16175, 24252, 40114, 113123, 162511, 312331, 1034350, 1534189, 2347317, 4549549, 9277210, 26474229, 36473890, 97757197, 237768730, 345116706, 741238897, 1712771267, 2267130204 },
{ 1, 2, 4, 15, 19, 38, 89, 191, 357, 823, 1121, 2576, 6101, 12147, 26158, 42353, 94687, 137823, 266049, 540248, 2047375, 3020863, 4548705, 9046231, 18567594, 65440423, 77280938, 153496698, 387596844, 748961714, 1558077359, 3341559704 },
{ 1, 2, 6, 9, 17, 46, 118, 244, 510, 562, 1407, 2963, 6084, 16369, 22460, 42783, 71291, 220540, 408471, 907211, 1314793, 2209667, 4669255, 9355001, 26375286, 40695899, 68317318, 187494661, 475161556, 1072482162, 2048009908, 2246857941 },
{ 1, 3, 6, 9, 23, 63, 72, 230, 303, 855, 1661, 2562, 7583, 15761, 23948, 40354, 103900, 169317, 400531, 661155, 2024303, 4113983, 4690677, 13153414, 26787481, 41186088, 93896391, 255095794, 332971942, 993292362, 1339603171, 3635706144 },
{ 1, 3, 7, 15, 28, 62, 77, 245, 285, 810, 1618, 2641, 7662, 13699, 27905, 52396, 108157, 186938, 451936, 777269, 1961953, 4096186, 4676076, 13187647, 31161707, 66432033, 123430860, 259363029, 316346730, 1047880602, 1104857001, 3325509700 },
{ 1, 3, 7, 11, 19, 40, 113, 243, 446, 986, 1226, 2706, 7841, 13725, 18633, 47588, 74591, 216460, 411877, 991641, 1768372, 2909193, 4662118, 13424081, 31297643, 48236748, 77304299, 183472960, 499142076, 1052815535, 1786853683, 4265833096 },
{ 1, 3, 6, 12, 21, 38, 105, 196, 496, 522, 1555, 2864, 7673, 16133, 28842, 33083, 91113, 230451, 505183, 912939, 1383269, 3765601, 4630130, 13560808, 26463280, 53943641, 85200423, 159021936, 461401415, 851492379, 2025919276, 2281899456 },
{ 1, 2, 4, 9, 26, 32, 65, 160, 416, 586, 1858, 2899, 5262, 11799, 20524, 57424, 92318, 168392, 316081, 921186, 1396984, 3277094, 4582310, 9207227, 17886603, 41850358, 114588986, 145775720, 285214488, 717229009, 1870665067, 2562731252 },
{ 1, 2, 5, 12, 25, 44, 68, 212, 445, 670, 2036, 3050, 5365, 9944, 30751, 61479, 75860, 166125, 461260, 961058, 1079001, 3465245, 4610082, 9228376, 22434036, 54053344, 101790310, 196525581, 268435872, 914358972, 1749026732, 2944404254 },
{ 1, 3, 5, 11, 22, 36, 85, 201, 303, 803, 1457, 2797, 8179, 9354, 18924, 33305, 87731, 151298, 511425, 707186, 1938960, 2258573, 4642664, 13110625, 22746104, 48944306, 98353544, 165690020, 339812107, 809684433, 1325832792, 3272115800 },
{ 1, 3, 4, 8, 20, 51, 82, 129, 390, 594, 1915, 2932, 7361, 11926, 20966, 39554, 130680, 161960, 338428, 989866, 1330689, 3360852, 4489358, 13496734, 17859197, 37621518, 89500564, 204543380, 365028803, 578968905, 1703240927, 2588913334 },
{ 1, 2, 6, 15, 25, 39, 110, 173, 455, 629, 1151, 2857, 6032, 15832, 25066, 59913, 95461, 248497, 298264, 804996, 1431986, 2642681, 4562393, 9396712, 26227215, 65529066, 102021800, 158049599, 440496362, 671334687, 1980013374, 2538370988 },
{ 1, 3, 7, 8, 27, 36, 65, 137, 325, 756, 1176, 2928, 8120, 13434, 21229, 64696, 86846, 167714, 361728, 610428, 1120488, 2881922, 4496171, 13309718, 30776697, 35707032, 110186878, 166410910, 298157123, 545870467, 1448156268, 3056333534 },
{ 1, 3, 4, 15, 31, 58, 100, 160, 335, 649, 1730, 2832, 7558, 11841, 27145, 51167, 98630, 244382, 265964, 588638, 1215847, 3708912, 4585739, 13347451, 18313042, 66263107, 127066361, 242768256, 407135094, 721958177, 1427314437, 2939881616 },
{ 1, 2, 7, 11, 26, 60, 65, 202, 314, 971, 2022, 2669, 5618, 14106, 17166, 65224, 100748, 168568, 482778, 683889, 1756072, 3243871, 4621258, 8927204, 31218282, 49903097, 114986752, 263293746, 298262153, 818581830, 1239106370, 4281593361 },
{ 1, 3, 5, 14, 28, 59, 73, 240, 489, 1022, 1143, 2706, 7889, 9237, 25098, 55217, 102057, 183485, 426801, 826842, 1739141, 2761594, 4683052, 13543527, 22935702, 61112700, 123376263, 237637374, 323423344, 1062097562, 2082213574, 4075205692 },
{ 1, 3, 5, 9, 28, 42, 90, 240, 355, 518, 1362, 3048, 7819, 9500, 23341, 55197, 81849, 159728, 397114, 585276, 1438886, 2607127, 4681776, 13269116, 22184127, 40745381, 123071253, 174583806, 394735401, 1049273862, 1586734275, 2305816715 }
        };
        assert(m <= mmax);
        assert(s <= (sizeof(Cs_all) / (sizeof(UINT_T) * mmax)));
        // we need to make an explicit copy, can't get around that here:
        for(unsigned j = 0; j < s; ++j, Cs += m)
            std::copy_n(std::begin(Cs_all[j]), m, Cs);
    }
    /** Some default Sobol sequence with a maximum of 250 dimensions and 2^32
      * points.
      */
    template <typename FLOAT_T=double, typename UINT_T=std::uint32_t>
    digitalseq_b2g<FLOAT_T, UINT_T> JK2008_sobolseq(unsigned s=250, unsigned m=32)
    {
        // we need to make an explicit copy, can't get around that here:
        std::vector<UINT_T> Cs(m*s);
        JK2008_sobol_matrices(Cs.begin(), s, m);
        auto g = digitalseq_b2g<FLOAT_T, UINT_T>(s, m, std::begin(Cs));
        return g;
    }

}

#endif // QMC_DIGITALSEQ_B2G_HPP
