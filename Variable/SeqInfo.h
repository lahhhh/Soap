#pragma once

#include "Identifier.h"

#include <QVector>
#include <QString>

class SeqInfo
{
public:
	G_CLASS_FUNCTION_DEFAULT(SeqInfo);

	QStringList sequence_names_;

	QVector<int> sequence_lengths_;

	QVector<bool> is_circular_;

	QStringList genome_;

	void append(const SeqInfo& seq_info);

	void clear() {
		this->sequence_names_.clear();
		this->sequence_lengths_.clear();
		this->is_circular_.clear();
		this->genome_.clear();
	}
};

