// reed solomon erasure correction. copyright 2015 Peter Bright, Backblaze. See LICENSE.txt for licensing details.

#pragma once

#include "galois.hpp"
#include "matrix.hpp"

#include <memory>

#include <ppl.h>

#include <xmmintrin.h>
#include <tmmintrin.h>

#define _mm_srli_epi8(_A, _Imm) (_mm_and_si128(_mm_set1_epi8(static_cast<int8_t>(                             0xFF >> _Imm  )), _mm_srli_epi32(_A, _Imm)))
#define _mm_slli_epi8(_A, _Imm) (_mm_and_si128(_mm_set1_epi8(static_cast<int8_t>(static_cast<uint8_t>(0xFF & (0xFF << _Imm)))), _mm_slli_epi32(_A, _Imm)))

struct reed_solomon
{
	static constexpr size_t alignment = 16;

	reed_solomon(uint8_t dsc, uint8_t psc) : data_shard_count(dsc), parity_shard_count(psc), total_shard_count(dsc + psc), m(build_matrix(dsc, dsc + psc)), parity_rows(new const uint8_t*[psc])
	{
		if(data_shard_count + parity_shard_count > 255)
		{
			throw - 9;
		}

		for(size_t i = 0; i < parity_shard_count; ++i)
		{
			parity_rows[i] = m.get_row(data_shard_count + i);
		}
	}

	~reed_solomon()
	{
		delete[] parity_rows;
	}

	uint8_t get_data_shard_count() const
	{
		return data_shard_count;
	}

	uint8_t get_parity_shard_count() const
	{
		return parity_shard_count;
	}

	uint8_t get_total_shard_count() const
	{
		return total_shard_count;
	}

	void encode_parity(uint8_t* __restrict* __restrict shards, size_t offset, size_t shard_size) const
	{
		// shards[0] through shards[data_shard_count - 1] contain the file data
		// shards[data_shard_count] through shards[data_shard_count + parity_shard_count - 1] are where parity data should be written to
		const uint8_t** inputs = const_cast<const uint8_t**>(&shards[0]);
		uint8_t** outputs = &shards[data_shard_count];
		code_some_shards(parity_rows, inputs, data_shard_count, outputs, parity_shard_count, offset, shard_size);
	}

	bool is_parity_correct(const uint8_t* __restrict* __restrict shards, size_t offset, size_t shard_size) const
	{
		const uint8_t** inputs = &shards[0];
		const uint8_t** parities = &shards[data_shard_count];
		return check_some_shards(parity_rows, inputs, data_shard_count, parities, parity_shard_count, offset, shard_size);
	}

	bool decode_missing(uint8_t* __restrict* __restrict shards, bool* shard_present, size_t offset, size_t shard_size) const
	{
		size_t number_present = 0;
		for(size_t i = 0; i < total_shard_count; ++i)
		{
			if(shard_present[i])
			{
				++number_present;
			}
		}
		if(number_present == total_shard_count)
		{
			return true;
		}
		if(number_present < data_shard_count)
		{
			return false;
		}

		matrix sub_matrix{ data_shard_count, data_shard_count };
		std::unique_ptr<const uint8_t*[]> sub_shards{ new const uint8_t*[data_shard_count] };
		{
			size_t sub_matrix_row = 0;
			for(size_t matrix_row = 0; matrix_row < total_shard_count && sub_matrix_row < data_shard_count; ++matrix_row)
			{
				if(shard_present[matrix_row])
				{
					for(size_t c = 0; c < data_shard_count; ++c)
					{
						sub_matrix.set(sub_matrix_row, c, m.get(matrix_row, c));
					}
					sub_shards[sub_matrix_row] = shards[matrix_row];
					++sub_matrix_row;
				}
			}
		}
		matrix data_decode_matrix = sub_matrix.invert();

		std::unique_ptr<uint8_t*[]> outputs{ new uint8_t*[parity_shard_count] };
		std::unique_ptr<const uint8_t*[]> matrix_rows{ new const uint8_t*[parity_shard_count] };
		uint8_t output_count = 0;
		for(int shard = 0; shard < data_shard_count; ++shard)
		{
			if(!shard_present[shard])
			{
				outputs[output_count] = shards[shard];
				matrix_rows[output_count] = data_decode_matrix.get_row(shard);
				++output_count;
			}
		}
		code_some_shards(matrix_rows.get(), sub_shards.get(), data_shard_count, outputs.get(), output_count, offset, shard_size);
		output_count = 0;
		for(int shard = data_shard_count; shard < total_shard_count; shard++)
		{
			if(!shard_present[shard])
			{
				outputs[output_count] = shards[shard];
				matrix_rows[output_count] = parity_rows[shard - data_shard_count];
				++output_count;
			}
		}
		code_some_shards(matrix_rows.get(), const_cast<const uint8_t**>(shards), data_shard_count, outputs.get(), output_count, offset, shard_size);
		return true;
	}

private:
	// http://www.snia.org/sites/default/files2/SDC2013/presentations/NewThinking/EthanMiller_Screaming_Fast_Galois_Field%20Arithmetic_SIMD%20Instructions.pdf
	void do_multiply(uint8_t matrix_value, const uint8_t* __restrict inputs, uint8_t* __restrict outputs, size_t offset, size_t byte_count) const
	{
		if((reinterpret_cast<size_t>(&inputs[offset]) & (alignment - 1)) != (reinterpret_cast<size_t>(&outputs[offset]) & (alignment - 1)))
		{
			// align on output, leave input unaligned. Rationale: make code more similar to xor case!
			size_t head = (alignment - (reinterpret_cast<size_t>(&outputs[offset]) & (alignment - 1))) % alignment;
			size_t body = (byte_count - head) & (~(alignment - 1));
			size_t tail = byte_count - body - head;
			for(size_t i = offset; i < offset + head; ++i)
			{
				outputs[i] = galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
			const __m128i low_table  = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_LOW [matrix_value].data()));
			const __m128i high_table = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_HIGH[matrix_value].data()));
			const __m128i mask       = _mm_set1_epi8(0x0f);
			const __m128i* __restrict input_ptr  = reinterpret_cast<const __m128i*>(&inputs [offset + head]);
			      __m128i* __restrict output_ptr = reinterpret_cast<      __m128i*>(&outputs[offset + head]);
			for(size_t i = offset + head; i < offset + head + body; i += 16, ++input_ptr, ++output_ptr)
			{
				__m128i input        = _mm_loadu_si128(input_ptr);
				__m128i low_indices  = _mm_and_si128(input, mask);
				__m128i high_indices = _mm_srli_epi8(input, 4);
				__m128i low_parts    = _mm_shuffle_epi8(low_table, low_indices);
				__m128i high_parts   = _mm_shuffle_epi8(high_table, high_indices);
				__m128i output       = _mm_xor_si128(low_parts, high_parts);
				                       _mm_store_si128(output_ptr, output);
			}
			for(size_t i = offset + head + body; i < offset + head + body + tail; ++i)
			{
				outputs[i] = galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
		}
		else
		{
			size_t head = (alignment - (reinterpret_cast<size_t>(&outputs[offset]) & (alignment - 1))) % alignment;
			size_t body = (byte_count - head) & (~(alignment - 1));
			size_t tail = byte_count - body - head;
			for(size_t i = offset; i < offset + head; ++i)
			{
				outputs[i] = galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
			const __m128i low_table  = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_LOW [matrix_value].data()));
			const __m128i high_table = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_HIGH[matrix_value].data()));
			const __m128i mask       = _mm_set1_epi8(0x0f);
			const __m128i* __restrict input_ptr  = reinterpret_cast<const __m128i*>(&inputs [offset + head]);
			      __m128i* __restrict output_ptr = reinterpret_cast<      __m128i*>(&outputs[offset + head]);
			for(size_t i = offset + head; i < offset + head + body; i += 16, ++input_ptr, ++output_ptr)
			{
				__m128i input        = _mm_load_si128(input_ptr);
				__m128i low_indices  = _mm_and_si128(input, mask);
				__m128i high_indices = _mm_srli_epi8(input, 4);
				__m128i low_parts    = _mm_shuffle_epi8(low_table, low_indices);
				__m128i high_parts   = _mm_shuffle_epi8(high_table, high_indices);
				__m128i output       = _mm_xor_si128(low_parts, high_parts);
				                       _mm_store_si128(output_ptr, output);
			}
			for(size_t i = offset + head + body; i < offset + head + body + tail; ++i)
			{
				outputs[i] = galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
		}
	}

	void do_multiply_xor(uint8_t matrix_value, const uint8_t* __restrict inputs, uint8_t* __restrict outputs, size_t offset, size_t byte_count) const
	{
		if((reinterpret_cast<size_t>(&inputs[offset]) & (alignment - 1)) != (reinterpret_cast<size_t>(&outputs[offset]) & (alignment - 1)))
		{
			// align on output, leave input unaligned. Rationale: input has one load; output has one load and one store.
			size_t head = (alignment - (reinterpret_cast<size_t>(&outputs[offset]) & (alignment - 1))) % alignment;
			size_t body = (byte_count - head) & (~(alignment - 1));
			size_t tail = byte_count - body - head;
			for(size_t i = offset; i < offset + head; ++i)
			{
				outputs[i] ^= galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
			const __m128i low_table  = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_LOW[matrix_value].data()));
			const __m128i high_table = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_HIGH[matrix_value].data()));
			const __m128i mask       = _mm_set1_epi8(0x0f);
			const __m128i* __restrict input_ptr  = reinterpret_cast<const __m128i*>(&inputs [offset + head]);
			      __m128i* __restrict output_ptr = reinterpret_cast<      __m128i*>(&outputs[offset + head]);
			for(size_t i = offset + head; i < offset + head + body; i += 16, ++input_ptr, ++output_ptr)
			{
				__m128i input          = _mm_loadu_si128(input_ptr);
				__m128i initial_output = _mm_load_si128(output_ptr);
				__m128i low_indices    = _mm_and_si128(input, mask);
				__m128i high_indices   = _mm_srli_epi8(input, 4);
				__m128i low_parts      = _mm_shuffle_epi8(low_table, low_indices);
				__m128i high_parts     = _mm_shuffle_epi8(high_table, high_indices);
				__m128i output         = _mm_xor_si128(low_parts, high_parts);
				        output         = _mm_xor_si128(initial_output, output);
				                         _mm_store_si128(output_ptr, output);
			}
			for(size_t i = offset + head + body; i < offset + head + body + tail; ++i)
			{
				outputs[i] ^= galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
		}
		else
		{
			size_t head = (alignment - (reinterpret_cast<size_t>(&outputs[offset]) & (alignment - 1))) % alignment;
			size_t body = (byte_count - head) & (~(alignment - 1));
			size_t tail = byte_count - body - head;
			for(size_t i = offset; i < offset + head; ++i)
			{
				outputs[i] ^= galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
			const __m128i low_table  = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_LOW [matrix_value].data()));
			const __m128i high_table = _mm_loadu_si128(reinterpret_cast<const __m128i*>(galois.MULTIPLICATION_TABLE_HIGH[matrix_value].data()));
			const __m128i mask       = _mm_set1_epi8(0x0f);
			const __m128i* __restrict input_ptr  = reinterpret_cast<const __m128i*>(&inputs [offset + head]);
			      __m128i* __restrict output_ptr = reinterpret_cast<      __m128i*>(&outputs[offset + head]);
			for(size_t i = offset + head; i < offset + head + body; i += 16, ++input_ptr, ++output_ptr)
			{
				__m128i input          = _mm_load_si128(input_ptr);
				__m128i initial_output = _mm_load_si128(output_ptr);
				__m128i low_indices    = _mm_and_si128(input, mask);
				__m128i high_indices   = _mm_srli_epi8(input, 4);
				__m128i low_parts      = _mm_shuffle_epi8(low_table, low_indices);
				__m128i high_parts     = _mm_shuffle_epi8(high_table, high_indices);
				__m128i output         = _mm_xor_si128(low_parts, high_parts);
				        output         = _mm_xor_si128(initial_output, output);
				                         _mm_store_si128(output_ptr, output);
			}
			for(size_t i = offset + head + body; i < offset + head + body + tail; ++i)
			{
				outputs[i] ^= galois.MULTIPLICATION_TABLE[matrix_value][inputs[i]];
			}
		}
	}

	void code_some_shards(const uint8_t* __restrict* __restrict matrix_rows, const uint8_t* __restrict* __restrict inputs, uint8_t input_count, uint8_t* __restrict* __restrict outputs, uint8_t output_count, size_t offset, size_t byte_count) const
	{
		static constexpr size_t chunk_size = 4096;
		static const size_t chunks = byte_count / chunk_size;
		concurrency::parallel_for(static_cast<size_t>(0), chunks, [&](size_t chunk)
		{
			for(int output_shard = 0; output_shard < output_count; ++output_shard)
			{
				{
					int input_shard = 0;
					do_multiply(matrix_rows[output_shard][input_shard], inputs[input_shard], outputs[output_shard], offset + (chunk * chunk_size), chunk_size);
				}
				for(int input_shard = 1; input_shard < input_count; ++input_shard)
				{
					do_multiply_xor(matrix_rows[output_shard][input_shard], inputs[input_shard], outputs[output_shard], offset + (chunk * chunk_size), chunk_size);
				}
			}
		});
		if(chunks * chunk_size < byte_count)
		{
			for(int output_shard = 0; output_shard < output_count; ++output_shard)
			{
				{
					int input_shard = 0;
					do_multiply(matrix_rows[output_shard][input_shard], inputs[input_shard], outputs[output_shard], offset + (chunks * chunk_size), byte_count - (chunks * chunk_size));
				}
				for(int input_shard = 1; input_shard < input_count; ++input_shard)
				{
					do_multiply_xor(matrix_rows[output_shard][input_shard], inputs[input_shard], outputs[output_shard], offset + (chunks * chunk_size), byte_count - (chunks * chunk_size));
				}
			}
		}
	}

	bool check_some_shards(const uint8_t* __restrict* __restrict matrix_rows, const uint8_t* __restrict* __restrict datas, uint8_t data_count, const uint8_t* __restrict* __restrict parities, uint8_t parity_count, size_t offset, size_t byte_count) const
	{
		static constexpr size_t chunk_size = 4096;
		static const size_t chunks = byte_count / chunk_size;
		std::unique_ptr<unsigned char[]> buffer(new unsigned char[(offset * parity_count) + (parity_count * byte_count)]);
	
		concurrency::combinable<bool> ok;
		concurrency::parallel_for(static_cast<size_t>(0), chunks, [&](size_t chunk)
		{
			ok.local() = true;
			for(int output_shard = 0; output_shard < parity_count; ++output_shard)
			{
				{
					int input_shard = 0;
					do_multiply(matrix_rows[output_shard][input_shard], datas[input_shard], buffer.get(), offset + (chunk * chunk_size), chunk_size);
				}
				for(int input_shard = 1; input_shard < data_count; ++input_shard)
				{
					do_multiply_xor(matrix_rows[output_shard][input_shard], datas[input_shard], buffer.get(), offset + (chunk * chunk_size), chunk_size);
				}
				if(0 != std::memcmp(buffer.get() + offset + (chunk * chunk_size), parities[output_shard] + offset + (chunk * chunk_size), chunk_size))
				{
					ok.local() = false;
					break;
				}
			}
		});
		bool all_ok = true;
		ok.combine_each([&](bool val)
		{
			all_ok = all_ok && val;
		});
		if(all_ok && (chunks * chunk_size) < byte_count)
		{
			for(int output_shard = 0; output_shard < parity_count; ++output_shard)
			{
				{
					int input_shard = 0;
					do_multiply(matrix_rows[output_shard][input_shard], datas[input_shard], buffer.get(), offset + (chunks * chunk_size), byte_count - (chunks * chunk_size));
				}
				for(int input_shard = 1; input_shard < data_count; ++input_shard)
				{
					do_multiply_xor(matrix_rows[output_shard][input_shard], datas[input_shard], buffer.get(), offset + (chunks * chunk_size), byte_count - (chunks * chunk_size));
				}
				if(0 != std::memcmp(buffer.get() + offset + (chunks * chunk_size), parities[output_shard] + offset + (chunks * chunk_size), byte_count - (chunks * chunk_size)))
				{
					return false;
				}
			}
		}
		return all_ok;
	}

	static matrix build_matrix(uint8_t data_shards, uint8_t total_shards)
	{
		matrix v = vandermonde(total_shards, data_shards);
		matrix top = v.submatrix(0, 0, data_shards, data_shards);
		return v.times(top.invert());
	}

	static matrix vandermonde(uint8_t rows, uint8_t cols)
	{
		matrix result{ rows, cols };
		for(uint8_t r = 0; r < rows; ++r)
		{
			for(uint8_t c = 0; c < cols; ++c)
			{
				result.set(r, c, galois.exp(r, c));
			}
		}
		return result;
	}

	uint8_t data_shard_count;
	uint8_t parity_shard_count;
	uint8_t total_shard_count;

	matrix m;

	const uint8_t* __restrict* __restrict parity_rows;
};
