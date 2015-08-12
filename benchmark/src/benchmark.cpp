// copyright 2015 Peter Bright. See LICENSE.txt for licensing details.

#define NOMINMAX

#include <SDKDDKVer.h>

#include "encoder.hpp"

#include <chrono>
#include <iostream>
#include <vector>
#include <random>

struct high_priority_observer : tbb::task_scheduler_observer
{
	virtual void on_scheduler_entry(bool is_worker)
	{
		if(is_worker)
		{
			::SetThreadPriority(::GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
		}
	}
};

int main()
{
	high_priority_observer observer;
	observer.observe(true);

	static constexpr uint8_t DATA_COUNT = 16;
	static constexpr uint8_t PARITY_COUNT = 4;
	static constexpr uint8_t TOTAL_COUNT = DATA_COUNT + PARITY_COUNT;
	static constexpr size_t BUFFER_SIZE = 16 * 1024 * 1024;
	static constexpr size_t NUMBER_OF_BUFFER_SETS = 1;
	static constexpr std::chrono::seconds MEASUREMENT_DURATION{ 10 };

	struct buffer_set
	{
		buffer_set(size_t buffer_size, size_t shard_count) : all_shards{ new unsigned char[shard_count * buffer_size] }, shards{ new unsigned char*[shard_count] }
		{
			static std::independent_bits_engine<std::default_random_engine, 16, unsigned short> engine(0);
			static std::uniform_int_distribution<unsigned short> distribution(0, std::numeric_limits<unsigned char>::max());

			for(size_t i = 0; i < shard_count; ++i)
			{
				shards[i] = &all_shards[i * buffer_size];
			}
			for(size_t i = 0; i < shard_count * buffer_size; ++i)
			{
				all_shards[i] = static_cast<unsigned char>(distribution(engine));
			}
		}
		std::unique_ptr<unsigned char[]> all_shards;
		std::unique_ptr<unsigned char*[]> shards;
	};

	std::vector<buffer_set> buffers;
	for(size_t i = 0; i < NUMBER_OF_BUFFER_SETS; ++i)
	{
		buffers.push_back(buffer_set{ BUFFER_SIZE, TOTAL_COUNT });
	}

	reed_solomon rs{ DATA_COUNT, PARITY_COUNT };

	size_t passes_completed = 0;
	size_t bytes_encoded = 0;
	size_t current_buffer = 0;
	std::chrono::nanoseconds encoding_time{ 0 };
	std::cout << "starting..." << std::endl;
	while(encoding_time < MEASUREMENT_DURATION)
	{
		auto start = std::chrono::high_resolution_clock::now();
		rs.encode_parity(buffers[current_buffer % NUMBER_OF_BUFFER_SETS].shards.get(), 0, BUFFER_SIZE);
		auto end = std::chrono::high_resolution_clock::now();
		encoding_time += (end - start);
		bytes_encoded += BUFFER_SIZE * DATA_COUNT;
		++passes_completed;
	}
	std::cout << "done" << std::endl;
	size_t megabytes = bytes_encoded / (1024 * 1024);
	float seconds = std::chrono::duration_cast<std::chrono::duration<float>>(encoding_time).count();
	std::cout << megabytes << " MiB in " << seconds << " seconds = " << (megabytes / seconds) << " MiB/s in " << passes_completed << " iterations" << std::endl;
}

