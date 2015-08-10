// GF(2^8) galois field. copyright 2015 Peter Bright, Backblaze. See LICENSE.txt for licensing details.

#pragma once

#include <cstdint>
#include <array>

struct galois_t
{
	static constexpr size_t FIELD_SIZE = 256;
	static constexpr size_t GENERATING_POLYNOMIAL = 29;

	galois_t() : MULTIPLICATION_TABLE     (generate_multiplication_table()     ),
	             MULTIPLICATION_TABLE_HIGH(generate_multiplication_table_high()),
	             MULTIPLICATION_TABLE_LOW (generate_multiplication_table_low() )
	{
	}

	uint8_t add(uint8_t a, uint8_t b) { return a ^ b; }

	uint8_t subtract(uint8_t a, uint8_t b) { return a ^ b; }

	uint8_t multiply(uint8_t a, uint8_t b)
	{
		if(a == 0 || b == 0)
		{
			return 0;
		}
		return EXP_TABLE[LOG_TABLE[a] + LOG_TABLE[b]];
	}

	uint8_t divide(uint8_t a, uint8_t b)
	{
		if(a == 0)
		{
			return 0;
		}
		if(b == 0)
		{
			// generate a divide by zero while ensuring the compiler won't complain about a static divide by zero
			volatile int one = 1, zero = 0;
			volatile int error = one / zero;
			error = error;
		}
		int16_t logA = LOG_TABLE[a];
		int16_t logB = LOG_TABLE[b];
		int16_t logResult = logA - logB;
		if(logResult < 0)
		{
			logResult += 255;
		}
		return EXP_TABLE[logResult];
	}

	uint8_t exp(uint8_t a, size_t n)
	{
		if(n == 0)
		{
			return 1;
		}
		if(a == 0)
		{
			return 0;
		}
		size_t logResult = LOG_TABLE[a] * n;
		while(255 <= logResult)
		{
			logResult -= 255;
		}
		return EXP_TABLE[logResult];
	}

	std::array<uint8_t, FIELD_SIZE> generate_log_table(uint8_t polynomial)
	{
		std::array<uint8_t, FIELD_SIZE> result;
		result.assign(0);
		int16_t b = 2;
		for(uint16_t log = 0; log < FIELD_SIZE - 1; ++log)
		{
			result[b] = static_cast<uint8_t>(log);
			b <<= 1;
			if(FIELD_SIZE <= b)
			{
				b = (b - FIELD_SIZE) ^ polynomial;
			}
		}
		return result;
	}

	std::array<uint8_t, (2 * FIELD_SIZE) - 2> generate_exp_table(const std::array<uint8_t, FIELD_SIZE>& logTable)
	{
		std::array<uint8_t, (2 * FIELD_SIZE) - 2> result = { 0 };
		for(size_t i = 1; i < FIELD_SIZE; ++i)
		{
			uint8_t log = logTable[i];
			result[log] = static_cast<uint8_t>(i);
			result[log + FIELD_SIZE - 1] = static_cast<uint8_t>(i);
		}
		return result;
	}

	std::array<std::array<uint8_t, FIELD_SIZE>, FIELD_SIZE> generate_multiplication_table()
	{
		std::array<std::array<uint8_t, FIELD_SIZE>, FIELD_SIZE> result = { 0 };
		for(size_t a = 0; a < FIELD_SIZE; a++)
		{
			for(size_t b = 0; b < FIELD_SIZE; b++)
			{
				result[a][b] = multiply(static_cast<uint8_t>(a), static_cast<uint8_t>(b));
			}
		}
		return result;
	}

	// http://jerasure.org/jerasure/gf-complete/blob/master/src/gf_w8.c#L1194 scase 3
	std::array<std::array<uint8_t, 16>, FIELD_SIZE> generate_multiplication_table_high()
	{
		std::array<std::array<uint8_t, 16>, FIELD_SIZE> result = { 0 };
		for(size_t a = 0; a < FIELD_SIZE; a++)
		{
			for(size_t b = 0; b < FIELD_SIZE; b++)
			{
				if(b == (b & 0xf0))
				{
					result[a][b >> 4] = multiply(static_cast<uint8_t>(a), static_cast<uint8_t>(b));
				}
			}
		}
		return result;
	}

	// http://jerasure.org/jerasure/gf-complete/blob/master/src/gf_w8.c#L1194 scase 3
	std::array<std::array<uint8_t, 16>, FIELD_SIZE> generate_multiplication_table_low()
	{
		std::array<std::array<uint8_t, 16>, FIELD_SIZE> result = { 0 };
		for(size_t a = 0; a < FIELD_SIZE; a++)
		{
			for(size_t b = 0; b < FIELD_SIZE; b++)
			{
				if(b == (b & 0xf))
				{
					result[a][b] = multiply(static_cast<uint8_t>(a), static_cast<uint8_t>(b));
				}
			}
		}
		return result;
	}


	std::array<uint8_t, FIELD_SIZE> LOG_TABLE = {
		0,    0,    1,   25,    2,   50,   26,  198,
		3,  223,   51,  238,   27,  104,  199,   75,
		4,  100,  224,   14,   52,  141,  239,  129,
		28,  193,  105,  248,  200,    8,   76,  113,
		5,  138,  101,   47,  225,   36,   15,   33,
		53,  147,  142,  218,  240,   18,  130,   69,
		29,  181,  194,  125,  106,   39,  249,  185,
		201,  154,    9,  120,   77,  228,  114,  166,
		6,  191,  139,   98,  102,  221,   48,  253,
		226,  152,   37,  179,   16,  145,   34,  136,
		54,  208,  148,  206,  143,  150,  219,  189,
		241,  210,   19,   92,  131,   56,   70,   64,
		30,   66,  182,  163,  195,   72,  126,  110,
		107,   58,   40,   84,  250,  133,  186,   61,
		202,   94,  155,  159,   10,   21,  121,   43,
		78,  212,  229,  172,  115,  243,  167,   87,
		7,  112,  192,  247,  140,  128,   99,   13,
		103,   74,  222,  237,   49,  197,  254,   24,
		227,  165,  153,  119,   38,  184,  180,  124,
		17,   68,  146,  217,   35,   32,  137,   46,
		55,   63,  209,   91,  149,  188,  207,  205,
		144,  135,  151,  178,  220,  252,  190,   97,
		242,   86,  211,  171,   20,   42,   93,  158,
		132,   60,   57,   83,   71,  109,   65,  162,
		31,   45,   67,  216,  183,  123,  164,  118,
		196,   23,   73,  236,  127,   12,  111,  246,
		108,  161,   59,   82,   41,  157,   85,  170,
		251,   96,  134,  177,  187,  204,   62,   90,
		203,   89,   95,  176,  156,  169,  160,   81,
		11,  245,   22,  235,  122,  117,   44,  215,
		79,  174,  213,  233,  230,  231,  173,  232,
		116,  214,  244,  234,  168,   80,   88,  175
	};

	std::array<uint8_t, (2 * FIELD_SIZE) - 2> EXP_TABLE = {
		1, 2, 4, 8, 16, 32, 64, 128, 29,
		58, 116, 232, 205, 135, 19, 38, 76,
		152, 45, 90, 180, 117, 234, 201, 143,
		3, 6, 12, 24, 48, 96, 192, 157,
		39, 78, 156, 37, 74, 148, 53, 106,
		212, 181, 119, 238, 193, 159, 35, 70,
		140, 5, 10, 20, 40, 80, 160, 93,
		186, 105, 210, 185, 111, 222, 161, 95,
		190, 97, 194, 153, 47, 94, 188, 101,
		202, 137, 15, 30, 60, 120, 240, 253,
		231, 211, 187, 107, 214, 177, 127, 254,
		225, 223, 163, 91, 182, 113, 226, 217,
		175, 67, 134, 17, 34, 68, 136, 13,
		26, 52, 104, 208, 189, 103, 206, 129,
		31, 62, 124, 248, 237, 199, 147, 59,
		118, 236, 197, 151, 51, 102, 204, 133,
		23, 46, 92, 184, 109, 218, 169, 79,
		158, 33, 66, 132, 21, 42, 84, 168,
		77, 154, 41, 82, 164, 85, 170, 73,
		146, 57, 114, 228, 213, 183, 115, 230,
		209, 191, 99, 198, 145, 63, 126, 252,
		229, 215, 179, 123, 246, 241, 255, 227,
		219, 171, 75, 150, 49, 98, 196, 149,
		55, 110, 220, 165, 87, 174, 65, 130,
		25, 50, 100, 200, 141, 7, 14, 28,
		56, 112, 224, 221, 167, 83, 166, 81,
		162, 89, 178, 121, 242, 249, 239, 195,
		155, 43, 86, 172, 69, 138, 9, 18,
		36, 72, 144, 61, 122, 244, 245, 247,
		243, 251, 235, 203, 139, 11, 22, 44,
		88, 176, 125, 250, 233, 207, 131, 27,
		54, 108, 216, 173, 71, 142,
		// Repeat the table a second time, so multiply()
		// does not have to check bounds.
		1, 2, 4, 8, 16, 32, 64, 128, 29,
		58, 116, 232, 205, 135, 19, 38, 76,
		152, 45, 90, 180, 117, 234, 201, 143,
		3, 6, 12, 24, 48, 96, 192, 157,
		39, 78, 156, 37, 74, 148, 53, 106,
		212, 181, 119, 238, 193, 159, 35, 70,
		140, 5, 10, 20, 40, 80, 160, 93,
		186, 105, 210, 185, 111, 222, 161, 95,
		190, 97, 194, 153, 47, 94, 188, 101,
		202, 137, 15, 30, 60, 120, 240, 253,
		231, 211, 187, 107, 214, 177, 127, 254,
		225, 223, 163, 91, 182, 113, 226, 217,
		175, 67, 134, 17, 34, 68, 136, 13,
		26, 52, 104, 208, 189, 103, 206, 129,
		31, 62, 124, 248, 237, 199, 147, 59,
		118, 236, 197, 151, 51, 102, 204, 133,
		23, 46, 92, 184, 109, 218, 169, 79,
		158, 33, 66, 132, 21, 42, 84, 168,
		77, 154, 41, 82, 164, 85, 170, 73,
		146, 57, 114, 228, 213, 183, 115, 230,
		209, 191, 99, 198, 145, 63, 126, 252,
		229, 215, 179, 123, 246, 241, 255, 227,
		219, 171, 75, 150, 49, 98, 196, 149,
		55, 110, 220, 165, 87, 174, 65, 130,
		25, 50, 100, 200, 141, 7, 14, 28,
		56, 112, 224, 221, 167, 83, 166, 81,
		162, 89, 178, 121, 242, 249, 239, 195,
		155, 43, 86, 172, 69, 138, 9, 18,
		36, 72, 144, 61, 122, 244, 245, 247,
		243, 251, 235, 203, 139, 11, 22, 44,
		88, 176, 125, 250, 233, 207, 131, 27,
		54, 108, 216, 173, 71, 142
	};

	std::array<std::array<uint8_t, FIELD_SIZE>, FIELD_SIZE> MULTIPLICATION_TABLE;
	std::array<std::array<uint8_t, 16>, FIELD_SIZE> MULTIPLICATION_TABLE_HIGH;
	std::array<std::array<uint8_t, 16>, FIELD_SIZE> MULTIPLICATION_TABLE_LOW;
};

extern galois_t galois;
