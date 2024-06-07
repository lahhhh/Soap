#pragma once

#include "Identifier.h"

#include <vector>

#include "Custom.h"

template <typename T, typename Index = int>
requires std::is_arithmetic_v<T>
class SparseVector
{
public:

	G_CLASS_FUNCTION_DEFAULT(SparseVector);
	
	std::vector<T> data_;
	std::vector<Index> index_;

	void set_from_pairs(std::vector< std::pair<Index, T> > pairs) {

		if (pairs.empty()) {
			return;
		}
		else if (pairs.size() == 1) {

			this->index_.push_back(pairs.front().first);
			this->data_.push_back(static_cast<T>(pairs.front().second));
			return;
		}

		auto compare_func = [](const std::pair<Index, T>& left, const std::pair<Index, T>& right) {
			return left.first < right.first;
		};

		std::sort(pairs.begin(), pairs.end(), compare_func);

		Index now_loc = pairs.front().first;
		T now_value = pairs.front().second;

		for (auto iter = pairs.cbegin() + 1; iter != pairs.cend(); ++iter) {
			if (iter->first == now_loc) {
				now_value += iter->second;
			}
			else {

				this->index_.push_back(now_loc);
				this->data_.push_back(now_value);

				now_loc = iter->first;
				now_value = iter->second;
			}
		}

		this->index_.push_back(now_loc);
		this->data_.push_back(now_value);
	}

	// arithmetic, return value directly instead of const reference
	// no check
	T operator[](const Index loc) const{

		auto it = std::ranges::lower_bound(this->index_, static_cast<Index>(loc));

		if (it != this->index_.end() && *it == loc) {

			std::size_t pos = std::distance(this->index_.begin(), it);

			return this->data_[pos];
		}

		return 0;
	}

	[[nodiscard]] std::vector<T> segment(const Index start, const Index size) const {

		std::vector<T> res(size, static_cast<T>(0));
				
		const std::size_t end = start + size; 
		
		auto it1 = std::ranges::lower_bound(this->index_, start);
		auto it2 = std::ranges::lower_bound(this->index_, end);

		for (; it1 != it2; ++it1) {

			std::size_t pos = std::distance(this->index_.begin(), it1);

			res[*it1 - start] = this->data_[pos];
		}

		return res;
	}

	SparseVector<T, Index>& operator *=(T val){

		std::ranges::for_each(this->data_, [val](auto&& ele) {
			ele *= val;
		});

		return *this;
	}
};

