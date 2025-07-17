#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

#include <cstdint>
#include <vector>
#include <fstream>

namespace sketching {

class HyperLogLog
{
    private:
        using register_t = uint8_t;
        using hash_t = __uint128_t;
        
    public:
        HyperLogLog();
        HyperLogLog(const uint8_t kmer_length, const uint8_t msb_length);
        HyperLogLog(const uint8_t kmer_length, const double error_rate);
        HyperLogLog(std::istream& istrm);
        void naive_add(char const * const seq, const std::size_t length) noexcept;
        void buffered_add(char const * const seq, const std::size_t length, std::vector<uint64_t>& buffer) noexcept;
        void clear() noexcept;
        std::size_t size() const noexcept;
        std::size_t count() const noexcept;
        double standard_error() const noexcept;
        HyperLogLog operator+(const HyperLogLog& other) const;
        HyperLogLog& operator+=(const HyperLogLog& other);
        void store(std::ostream& ostrm) const;
        void store(std::string const& sketch_file) const;
        static HyperLogLog load(std::string const& sketch_filename);

    private:
        friend HyperLogLog load_hll(std::istream& istrm);
        void init();
        void sanitize_endianness() const;
        void sanitize_kmer_length(const std::size_t kmer_length) const;
        void sanitize_b(const std::size_t bval) const;
        bool compatible(const HyperLogLog& other) const noexcept;
        double harmonic_mean() const noexcept;
        double bias_correction(const double raw_estimate) const noexcept;
        int clz(const uint32_t x) const noexcept;
        int clz(const uint64_t x) const noexcept;
        int clz(const __uint128_t x) const noexcept;
        uint8_t k;
        uint8_t b;
        std::vector<register_t> registers;
        std::size_t shift; // optimization
        hash_t mask;
        std::size_t total_seen_kmers; // with repetitions = L1 norm
        double alpha_m;
};

} // namespace sketching

#endif // HYPERLOGLOG_HPP