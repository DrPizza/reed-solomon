// copyright 2015 Peter Bright. See LICENSE.txt for licensing details.

#include <SDKDDKVer.h>

#include "encoder.hpp"

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	const char* const filename = argc > 1 ? argv[1] : argv[0];
	encoder e{ 17, 3 };

	struct header_t
	{
		uint64_t file_size;
		uint8_t data_shard_count;
		uint8_t parity_shard_count;
	};

	{
		std::ifstream fin(filename, std::ifstream::ate | std::ifstream::binary);
		const size_t file_size = static_cast<size_t>(fin.tellg().seekpos());
		fin.seekg(0, std::ios::beg);

		const uint64_t stored_size = static_cast<uint64_t>(file_size);
		encoder::buffer buf = e.allocate_buffers_from_object_size(file_size, sizeof(uint64_t));
		for(size_t i = 0; i < buf.shard_count; ++i)
		{
			*reinterpret_cast<uint64_t*>(buf.shards[i]) = stored_size;
			fin.read(reinterpret_cast<char*>(buf.shards[i]) + buf.padding_size, buf.shard_size - buf.padding_size);
		}
		fin.close();
		e.encode(buf);
		for(size_t i = 0; i < e.get_shard_count(); ++i)
		{
			std::ofstream fout(std::string(filename) + "." + std::to_string(i), std::ofstream::binary | std::ofstream::trunc);
			fout.write(reinterpret_cast<const char*>(buf.shards[i]), buf.shard_size);
		}
	}

	{
		// obviously in real situations you can't assume that file 0 is always available
		std::ifstream first_file(std::string(filename) + "." + std::to_string(0), std::ifstream::ate | std::ifstream::binary);
		const uint64_t shard_file_size = static_cast<size_t>(first_file.tellg().seekpos());
		first_file.close();

		encoder::buffer buf = e.allocate_buffers_from_shard_size(static_cast<size_t>(shard_file_size), sizeof(uint64_t));
		for(size_t i = 0; i < e.get_shard_count(); ++i)
		{
			std::ifstream fin(std::string(filename) + "." + std::to_string(i), std::ifstream::binary);
			fin.read(reinterpret_cast<char*>(buf.shards[i]), buf.shard_size);
		}
		std::cout << std::boolalpha;
	
		std::cout << "Does reading all parts verify? " << e.verify(buf) << std::endl;
		std::unique_ptr<bool[]> present{ new bool[e.get_shard_count()] };
		for(size_t i = 0; i < e.get_shard_count(); ++i)
		{
			present[i] = true;
		}
		// clobber a data shard
		present[3] = false;
		std::memset(buf.shards[3], 0, buf.shard_size);
		std::cout << "Does clobbering a data shard verify? " << e.verify(buf) << std::endl;
		e.repair(buf, present.get());
		std::cout << "Does reading a repaired data shard verify? " << e.verify(buf) << std::endl;
		// clobber a parity shard
		present[3] = false;
		present[5] = false;
		std::memset(buf.shards[3], 0, buf.shard_size);
		std::memset(buf.shards[5], 0, buf.shard_size);
		std::cout << "Does clobbering a parity shard verify? " << e.verify(buf) << std::endl;
		e.repair(buf, present.get());
		std::cout << "Does reading a repaired parity shard verify? " << e.verify(buf) << std::endl;

		std::ofstream fout(std::string(filename) + ".recovered", std::ofstream::binary | std::ofstream::trunc);
		uint64_t bytes_remaining = *reinterpret_cast<const uint64_t*>(buf.shards[0]);
		for(size_t i = 0; i < e.get_data_shard_count(); ++i)
		{
			fout.write(reinterpret_cast<const char*>(buf.shards[i]) + buf.padding_size, std::min(buf.shard_size - buf.padding_size, bytes_remaining));
			bytes_remaining -= buf.shard_size - buf.padding_size;
		}
	}

	return 0;
}

