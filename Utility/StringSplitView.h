#pragma once
#include <string>
#include <vector>

class StringSplitView {

public:

    class Iterator {

    private:

        const StringSplitView& view_;
        std::size_t loc_;

    public:

        Iterator(const StringSplitView& view, std::size_t loc) : view_(view), loc_(loc) {};

        std::string_view operator*() const {

            return std::string_view(view_.s_.data() + view_.d_[loc_], view_.d_[loc_ + 1] - view_.d_[loc_] - 1);
        }

        Iterator& operator++() {

            ++this->loc_;

            return *this;
        }

        bool operator!=(const Iterator& other) const {

            return loc_ != other.loc_;
        }

    };

private:
    const std::string& s_;
    std::vector<size_t> d_;

public:
    StringSplitView(const std::string& str, char delim) : s_(str) {
        std::size_t start{ 0 };
        std::size_t end{ 0 };
        while ((end = s_.find(delim, start)) != std::string::npos) {
            d_.push_back(start);
            start = end + 1;
        }
        d_.push_back(start);
        d_.push_back(s_.size() + 1);
    }

    std::string_view operator[](std::size_t loc) const {

        return std::string_view(this->s_.data() + this->d_[loc], this->d_[loc + 1] - this->d_[loc] - 1);
    }

    std::size_t size() const {
        return this->d_.size() - 1;
    }

    Iterator begin() const {
        return Iterator(*this, 0);
    }

    Iterator end() const {
        return Iterator(*this, this->size());
    }
};


