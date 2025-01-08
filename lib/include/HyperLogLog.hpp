#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

#include <cstdint>
#include <vector>

class HyperLogLog
{
    public:
        HyperLogLog(uint8_t kmer_length, uint8_t msb_length);
        void add(char const * const seq, std::size_t length) noexcept;
        std::size_t count() const noexcept;
        double standard_error() const noexcept;
        HyperLogLog operator+(const HyperLogLog& other) const;

    private:
        using register_t = uint8_t;
        using hash_t = uint64_t;
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

#endif