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
		buffer(size_t buffer_size_, size_t shard_size_, size_t shard_count_, size_t padding_size_) : buffer_size(buffer_size_),
		                                                                                             shard_size(shard_size_),
		                                                                                             shard_count(shard_count_),
		                                                                                             padding_size(padding_size_),
		                                                                                             data{new unsigned char[buffer_size] },
		                                                                                             shards { new unsigned char*[shard_count] }
		{
			std::memset(data.get(), 0, buffer_size);
			for(size_t i = 0; i < shard_count; ++i)
			{
				shards[i] = &data[i * shard_size];
			}
		}
	
		const size_t buffer_size;
		const size_t shard_size;
		const size_t shard_count;
		const size_t padding_size;

		std::unique_ptr<unsigned char[]> data;
		std::unique_ptr<unsigned char*[]> shards;
	};

	buffer allocate_buffers_from_object_size(const size_t object_size) const
	{
		const size_t shard_size = (((object_size + rs.get_data_shard_count() - 1) / rs.get_data_shard_count()) + reed_solomon::alignment) & ~(reed_solomon::alignment - 1);
		const size_t buffer_size = shard_size * get_shard_count();

		return buffer{ buffer_size, shard_size, get_shard_count(), 0 };
	}

	buffer allocate_buffers_from_object_size(const size_t object_size, const size_t minimum_padding) const
	{
		const size_t padding_size = (minimum_padding + reed_solomon::alignment) & ~(reed_solomon::alignment - 1);
		const size_t shard_size   = padding_size + ((((object_size + rs.get_data_shard_count() - 1) / rs.get_data_shard_count()) + reed_solomon::alignment) & ~(reed_solomon::alignment - 1));
		const size_t buffer_size  = shard_size * get_shard_count();

		return buffer{ buffer_size, shard_size, get_shard_count(), padding_size };
	}

	buffer allocate_buffers_from_shard_size(const size_t shard_size) const
	{
		const size_t buffer_size = shard_size * get_shard_count();
		return buffer{ buffer_size, shard_size, get_shard_count(), 0 };
	}

	buffer allocate_buffers_from_shard_size(const size_t shard_size, const size_t minimum_padding) const
	{
		const size_t padding_size = (minimum_padding + reed_solomon::alignment) & ~(reed_solomon::alignment - 1);
		const size_t buffer_size = shard_size * get_shard_count();
		return buffer{ buffer_size, shard_size, get_shard_count(), padding_size };
	}

	uint8_t get_data_shard_count() const
	{
		return rs.get_data_shard_count();
	}

	uint8_t get_shard_count() const
	{
		return rs.get_total_shard_count();
	}

	void encode(buffer& b) const
	{
		rs.encode_parity(b.shards.get(), b.padding_size, b.shard_size - b.padding_size);
	}

	bool verify(buffer& b) const
	{
		return rs.is_parity_correct(const_cast<const uint8_t**>(b.shards.get()), b.padding_size, b.shard_size - b.padding_size);
	}

	bool repair(buffer& b, bool* present)
	{
		return rs.decode_missing(b.shards.get(), present, b.padding_size, b.shard_size - b.padding_size);
	}
private:
	reed_solomon rs;
};
