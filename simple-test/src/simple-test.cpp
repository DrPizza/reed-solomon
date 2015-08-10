// copyright 2015 Peter Bright. See LICENSE.txt for licensing details.

#include <SDKDDKVer.h>

#include "encoder.hpp"

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	const char* const filename = argc > 1 ? argv[1] : argv[0];
	encoder e{ 4, 2 };

	{
		std::ifstream fin(filename, std::ifstream::ate | std::ifstream::binary);
		const uint64_t file_size = static_cast<size_t>(fin.tellg().seekpos());
		fin.seekg(0, std::ios::beg);

		const size_t stored_size = static_cast<size_t>(file_size + sizeof(uint64_t));
		encoder::buffer buf = e.allocate_buffers_from_object_size(stored_size);
		*reinterpret_cast<uint64_t*>(buf.data.get()) = static_cast<uint64_t>(file_size);
		fin.read(reinterpret_cast<char*>(buf.data.get()) + sizeof(uint64_t), file_size);
		fin.close();
		e.encode(buf);
		for(size_t i = 0; i < e.get_shard_count(); ++i)
		{
			std::ofstream fout(std::string(filename) + "." + std::to_string(i), std::ofstream::binary | std::ofstream::trunc);
			fout.write(reinterpret_cast<char*>(buf.shards[i]), buf.shard_size);
		}
	}

	{
		// obviously in real situations you can't assume that file 0 is always available
		std::ifstream first_file(std::string(filename) + "." + std::to_string(0), std::ifstream::ate | std::ifstream::binary);
		const uint64_t file_size = static_cast<size_t>(first_file.tellg().seekpos());
		first_file.close();

		encoder::buffer buf = e.allocate_buffers_from_shard_size(static_cast<size_t>(file_size));
		for(size_t i = 0; i < e.get_shard_count(); ++i)
		{
			std::ifstream fin(std::string(filename) + "." + std::to_string(i), std::ifstream::binary);
			fin.read(reinterpret_cast<char*>(buf.shards[i]), buf.shard_size);
		}

		std::cout << e.verify(buf) << std::endl;
		std::unique_ptr<bool[]> present{ new bool[e.get_shard_count()] };
		for(size_t i = 0; i < e.get_shard_count(); ++i)
		{
			present[i] = true;
		}
		// clobber a data shard
		present[3] = false;
		std::memset(buf.shards[3], 0, buf.shard_size);
		std::cout << e.verify(buf) << std::endl;
		e.repair(buf, present.get());
		std::cout << e.verify(buf) << std::endl;
		// clobber a parity shard
		present[5] = false;
		std::memset(buf.shards[5], 0, buf.shard_size);
		std::cout << e.verify(buf) << std::endl;
		e.repair(buf, present.get());
		std::cout << e.verify(buf) << std::endl;

		std::ofstream fout(std::string(filename) + ".recovered", std::ofstream::binary | std::ofstream::trunc);
		const uint64_t original_size = *reinterpret_cast<const uint64_t*>(buf.data.get());
		fout.write(reinterpret_cast<char*>(buf.data.get()) + sizeof(uint64_t), original_size);
	}

	return 0;
}

