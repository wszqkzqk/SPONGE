#ifndef SPONGE_BIAS_META_IO_TRANSACTION_H
#define SPONGE_BIAS_META_IO_TRANSACTION_H

#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace sponge_meta_io
{

inline std::mutex& Process_Mutex()
{
    static std::mutex process_mutex;
    return process_mutex;
}

template <typename HillVector, typename SinkVector>
void Publish_History(HillVector* hills, SinkVector* sink,
                     HillVector* parsed_hills, const SinkVector& parsed_sink)
{
    using Hill = typename HillVector::value_type;
    using SinkValue = typename SinkVector::value_type;
    static_assert(std::is_nothrow_move_constructible<Hill>::value,
                  "Transactional hills publication requires noexcept move");
    static_assert(std::is_nothrow_copy_constructible<SinkValue>::value,
                  "Transactional sink publication requires noexcept copy");

    if (parsed_hills->size() > hills->max_size() - hills->size() ||
        parsed_sink.size() > sink->max_size() - sink->size())
    {
        throw std::length_error("metadynamics history capacity overflow");
    }

    // Capacity changes are not logical publication.  Reserve both vectors
    // before appending either one so failure of the second allocation leaves
    // both observable histories unchanged.
    hills->reserve(hills->size() + parsed_hills->size());
    sink->reserve(sink->size() + parsed_sink.size());
    for (Hill& parsed_hill : *parsed_hills)
    {
        hills->emplace_back(std::move(parsed_hill));
    }
    for (const SinkValue& parsed_value : parsed_sink)
    {
        sink->emplace_back(parsed_value);
    }
}

}  // namespace sponge_meta_io

#endif
