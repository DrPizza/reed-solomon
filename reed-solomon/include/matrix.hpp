// simple galois matrix. copyright 2015 Peter Bright, Backblaze. See LICENSE.txt for licensing details.

#pragma once

#include "galois.hpp"

#include <stdexcept>
#include <iostream>

struct matrix;

std::ostream& operator<<(std::ostream& os, const matrix& rhs);

struct matrix
{
	matrix(size_t rows_, size_t columns_) : rows(rows_), columns(columns_), data(new uint8_t[rows_ * columns_]), row_pointers(new const uint8_t*[rows_])
	{
		std::memset(data, 0, rows * columns);
		build_row_pointers();
	}

	matrix(const matrix& rhs) : rows(rhs.rows), columns(rhs.columns), data(new uint8_t[rhs.rows * rhs.columns]), row_pointers(new const uint8_t*[rhs.rows])
	{
		std::memcpy(data, rhs.data, rows * columns);
		build_row_pointers();
	}

	matrix(matrix&& rhs) : rows(rhs.rows), columns(rhs.columns), data(rhs.data), row_pointers(rhs.row_pointers)
	{
		rhs.data = nullptr;
		rhs.row_pointers = nullptr;
	}

	~matrix()
	{
		delete[] row_pointers;
		delete[] data;
	}

	static matrix identity(size_t size)
	{
		matrix m{ size, size };
		for(size_t i = 0; i < size; ++i)
		{
			m.set(i, i, 1);
		}
		return m;
	}

	const uint8_t* get_row(size_t r) const
	{
		if(r >= rows)
		{
			throw std::out_of_range("no such row");
		}
		return row_pointers[r];
	}

	size_t get_rows() const
	{
		return rows;
	}

	size_t get_columns() const
	{
		return columns;
	}

	uint8_t get(size_t r, size_t c) const
	{
		if(r >= rows || c >= columns)
		{
			throw std::out_of_range("no such row or column");
		}
		return data[(r * columns) + c];
	}

	void set(size_t r, size_t c, uint8_t value)
	{
		if(r >= rows || c >= columns)
		{
			throw std::out_of_range("no such row or column");
		}
		data[(r * columns) + c] = value;
	}

	bool operator==(const matrix& rhs) const
	{
		return rows == rhs.rows && columns == rhs.columns && std::memcmp(data, rhs.data, rows * columns) == 0;
	}

	bool operator!=(const matrix& rhs) const
	{
		return !(*this == rhs);
	}

	matrix times(const matrix& rhs) const
	{
		if(columns != rhs.rows)
		{
			throw std::out_of_range("left.columns != right.rows");
		}
		matrix result{ rows, rhs.columns };
		for(size_t r = 0; r < rows; ++r)
		{
			for(size_t c = 0; c < rhs.columns; ++c)
			{
				uint8_t value = 0;
				for(size_t i = 0; i < columns; ++i)
				{
					value ^= galois.multiply(get(r, i), rhs.get(i, c));
				}
				result.set(r, c, value);
			}
		}
		return result;
	}

	matrix augment(const matrix& rhs) const
	{
		if(rows != rhs.rows)
		{
			throw std::out_of_range("left.rows != right.rows");
		}
		matrix result{ rows, columns + rhs.columns };
		for(size_t r = 0; r < rows; ++r)
		{
			for(size_t c = 0; c < columns; ++c)
			{
				result.set(r, c, get(r, c));
			}
			for(size_t c = 0; c < rhs.columns; ++c)
			{
				result.set(r, columns + c, rhs.get(r, c));
			}
		}
		return result;
	}

	matrix submatrix(size_t rmin, size_t cmin, size_t rmax, size_t cmax) const
	{
		matrix result{ rmax - rmin, cmax - cmin };
		for(size_t r = rmin; r < rmax; ++r)
		{
			for(size_t c = cmin; c < cmax; ++c)
			{
				result.set(r - rmin, c - cmin, get(r, c));
			}
		}
		return result;
	}

	void swap_rows(size_t r1, size_t r2)
	{
		if(r1 >= rows || r2 >= rows) {
			throw std::out_of_range("no such row");
		}
		for(size_t i = 0; i < columns; ++i) {
			uint8_t tmp = get(r1, i);
			set(r1, i, get(r2, i));
			set(r2, i, tmp);
		}
	}

	void multiply_row(size_t r, uint8_t m) {
		if(r >= rows) {
			throw std::out_of_range("no such row");
		}
		for(size_t c = 0; c < columns; ++c) {
			set(r, c, galois.multiply(get(r, c), m));
		}
	}

	// row dst <- dst + mul * scale
	void row_linear_combination(size_t dst, size_t src, uint8_t scale) {
		if(dst >= rows || src >= rows) {
			throw std::out_of_range("no such row");
		}
		for(size_t c = 0; c < columns; ++c) {
			uint8_t d = get(dst, c);
			uint8_t s = get(src, c);
			set(dst, c, d ^ galois.multiply(s, scale));
		}
	}

	matrix invert() const
	{
		if(rows != columns) {
			throw std::out_of_range("matrix not square");
		}
		matrix work = augment(identity(rows));
		// work = { M | I }
		work.gaussian_elimination();
		// work = { I | M^-1 }
		return work.submatrix(0, rows, columns, columns * 2);
	}

private:
	void gaussian_elimination()
	{
		for(size_t pivot = 0; pivot < rows; ++pivot) {
			if(get(pivot, pivot) == 0) {
				for(size_t row_below = pivot + 1; row_below < rows; ++row_below) {
					if(get(row_below, pivot) != 0) {
						swap_rows(pivot, row_below);
					}
				}
			}
			if(get(pivot, pivot) == 0) {
				throw std::runtime_error("matrix is singular");
			}
			// row pivot = row pivot * get(pivot, pivot)^-1
			if(get(pivot, pivot) != 1) {
				uint8_t scale = galois.divide(1, get(pivot, pivot));
				multiply_row(pivot, scale);
			}
			// if any other row has a non-zero element in this column, subtract this row from it until it's zero
			for(size_t d = 0; d < rows; ++d) {
				if(d == pivot) {
					continue;
				}
				if(get(d, pivot) != 0) {
					row_linear_combination(d, pivot, get(d, pivot));
				}
			}
		}
	}

	void build_row_pointers()
	{
		for(size_t r = 0; r < rows; ++r)
		{
			row_pointers[r] = data + (r * columns);
		}
	}

	const size_t rows;
	const size_t columns;

	uint8_t* data;
	const uint8_t** row_pointers;
};

std::ostream& operator<<(std::ostream& os, const matrix& rhs) {
	os << "{\n";
	for(size_t r = 0; r < rhs.get_rows(); ++r) {
		os << "\t";
		os << std::hex << static_cast<unsigned int>(rhs.get(r, 0));
		
		for(size_t c = 1; c < rhs.get_columns(); ++c) {
			os << ", " << std::hex << static_cast<unsigned int>(rhs.get(r, c));
		}
		os << "\n";
	}
	os << "}\n";
	return os << std::dec;
}
