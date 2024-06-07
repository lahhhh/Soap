#pragma once

#include "Identifier.h"

#include <vector>
#include <string>

class TranscriptModel
{
public:
    G_CLASS_FUNCTION_DEFAULT(TranscriptModel);

    std::vector<int> loc_;

    std::string gene_name_;

    int start_ = 0;

    int end_ = 0;

    bool operator<(const TranscriptModel& rhs) const{
        return this->start_ == rhs.start_ ? (this->end_ < rhs.end_) : (this->start_ < rhs.start_);
    }

    template<typename String>
    void reset(String&& gene_name) {

        this->loc_.clear();
        this->gene_name_ = std::forward<String>(gene_name);
        this->start_ = 0;
        this->end_ = 0;
    }

    std::size_t size() const{
        return this->loc_.size();
    }

    bool empty() const{
        return this->loc_.empty();
    }

    void append_exon(int exon_start, int exon_end) {
        this->loc_.push_back(exon_start);
        this->loc_.push_back(exon_end);
    }

    void finalize();

    // return value - 0 : no overlap; 1 : part intron; 2 : exon; assume that loc and segments are sorted.
    int match(const std::vector<int>& segments) const;

};

