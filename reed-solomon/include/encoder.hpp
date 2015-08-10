// simple utility wrapper. copyright 2015 Peter Bright, Backblaze. See LICENSE.txt for licensing details.

#pragma once

#include "reed-solomon.hpp"

struct encoder
{
	encoder(uint8_t data_shard_count_, uint8_t parity_shard_count_) : rs{data_shard_count_, parity_shard_count_}
	{
	}

	struct buffer
	{
		buffer(size_t buffer_size_, size_t shard_size_, size_t shard_count_) : buffer_size(buffer_size_),
		                                                                       shard_size(shard_size_),
		                                                                       shard_count(shard_count_),
		                                                                       data{new unsigned char[buffer_size] },
		                                                                       shards { new unsigned char*[shard_count] }
		{
			std::memset(data.get(), 0, buffer_size);
			for(size_t i = 0; i < shard_count; ++i)
			{
				shards[i] = &data[i * shard_size];
			}
		}
	
		size_t buffer_size;
		size_t shard_size;
		size_t shard_count;

		std::unique_ptr<unsigned char[]> data;
		std::unique_ptr<unsigned char*[]> shards;
	};

	buffer allocate_buffers_from_object_size(size_t object_size) const
	{
		const size_t shard_size = (((object_size + rs.get_data_shard_count() - 1) / rs.get_data_shard_count()) + rs.alignment) & ~(rs.alignment - 1);
		const size_t buffer_size = shard_size * get_shard_count();

		return buffer{ buffer_size, shard_size, get_shard_count() };
	}

	buffer allocate_buffers_from_shard_size(size_t shard_size) const
	{
		const size_t buffer_size = shard_size * get_shard_count();
		return buffer{ buffer_size, shard_size, get_shard_count() };
	}

	uint8_t get_shard_count() const
	{
		return rs.get_total_shard_count();
	}

	void encode(buffer& b) const
	{
		rs.encode_parity(b.shards.get(), 0, b.shard_size);
	}

	bool verify(buffer& b) const
	{
		return rs.is_parity_correct(const_cast<const uint8_t**>(b.shards.get()), 0, b.shard_size);
	}

	bool repair(buffer& b, bool* present)
	{
		return rs.decode_missing(b.shards.get(), present, 0, b.shard_size);
	}
private:
	reed_solomon rs;
};
