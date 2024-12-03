#ifndef HYPERLOGLOG_HPP
#define HYPERLOGLOG_HPP

#include <vector>

class HyperLogLog
{
    public:
        HyperLogLog(uint8_t kmer_length, uint8_t msb_length);
        void add(char const * const seq, std::size_t length) noexcept;
        std::size_t count() const noexcept;
        HyperLogLog operator+(const HyperLogLog& other) const;

    private:
        bool compatible(const HyperLogLog& other) const noexcept;
        double harmonic_mean() const noexcept;
        double bias_correction(double raw_estimate) const noexcept;
        uint8_t b;
        uint8_t k;
        std::vector<uint8_t> table;
};

#endif