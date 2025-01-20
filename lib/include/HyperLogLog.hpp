#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

#include <cstdint>
#include <vector>
#include <fstream>

namespace sketching {

class HyperLogLog
{
    public:
        HyperLogLog();
        HyperLogLog(uint8_t kmer_length, uint8_t msb_length);
        HyperLogLog(uint8_t kmer_length, double error_rate);
        HyperLogLog(std::istream& istrm);
        void add(char const * const seq, std::size_t length) noexcept;
        void clear() noexcept;
        std::size_t count() const noexcept;
        double standard_error() const noexcept;
        HyperLogLog operator+(const HyperLogLog& other) const;
        HyperLogLog& operator+=(const HyperLogLog& other);
        void store(std::ostream& ostrm) const;
        void store(std::string const& sketch_file) const;
        static HyperLogLog load(std::string const& sketch_filename);

    private:
        using register_t = uint8_t;
        using hash_t = uint64_t;
        friend HyperLogLog load_hll(std::istream& istrm);
        void init();
        void sanitize_kmer_length(std::size_t kmer_length) const;
        void sanitize_b(std::size_t bval) const;
        bool compatible(const HyperLogLog& other) const noexcept;
        double harmonic_mean() const noexcept;
        double bias_correction(double raw_estimate) const noexcept;
        uint8_t k;
        uint8_t b;
        std::vector<register_t> registers;
        std::size_t shift; // optimization
        std::size_t mask;
        double alpha_m;
};

} // namespace sketching

#endif // HYPERLOGLOG_HPP