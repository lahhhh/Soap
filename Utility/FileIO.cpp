#include "FileIO.h"

#include <fstream>

#include "Custom.h"
#include "CustomMatrix.h"

std::map<QString, PatternWeightMatrix> read_motif_database(const QString& file_name) {

	std::map<QString, PatternWeightMatrix> ret;

	QFile file(file_name);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return ret;

	QTextStream in(&file);
	QString line = in.readLine();

	if (line.isNull()) {

		return ret;
	}

	while (!line.isNull()) {
		QStringList line_info = line.split(' ');

		QString motif_name = line_info[0];
		QString transcriptional_factor_name = line_info[1];
		auto& tf = ret[motif_name];
		tf.data_type_ = PatternWeightMatrix::DataType::Frequency;
		tf.motif_name_ = motif_name;
		tf.transcriptional_factor_name_ = transcriptional_factor_name;

		int pattern_length = (line_info.size() - 2) / 4;
		tf.weight_.mat_.resize(4, pattern_length);

		for (int i = 0; i < pattern_length; ++i) {
			tf.weight_.mat_(0, i) = line_info[4 * i + 2].toInt();
			tf.weight_.mat_(1, i) = line_info[4 * i + 3].toInt();
			tf.weight_.mat_(2, i) = line_info[4 * i + 4].toInt();
			tf.weight_.mat_(3, i) = line_info[4 * i + 5].toInt();
		}

		tf.weight_.rownames_ = { "A", "C", "G", "T" };
		tf.weight_.colnames_ = custom::cast<QString>(custom::seq_n(1, pattern_length));

		line = in.readLine();
	}

	return ret;
};

void write_sv(const CustomMatrix& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	bool keep_rownames,
	bool keep_colnames) {

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(mat.colnames_, delimiter, quote) << Qt::endl;
	}

	int nrow = mat.rows();
	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + mat.rownames_[i] + "\"" << delimiter;
			}
			else {
				out << mat.rownames_[i] << delimiter;
			}
		}

		out << custom::merge_to_string(mat.get_row(i), delimiter, quote) << Qt::endl;
	}

	out.flush();
};

void write_sv(
	const QStringList& s,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	bool keep_rownames,
	bool keep_colnames) {

	write_sv(s, out, delimiter, quote);
};

void write_sv(
	const QStringList& s,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = s.size();
	const int ncol{ 1 };

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << colnames[0] << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		if (quote) {
			out << "\"" + s[i] + "\"" << Qt::endl;
		}
		else {
			out << s[i] << Qt::endl;
		}
	}
};

void write_sv(
	const MotifPosition& mp,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	Eigen::MatrixXi mat = mp.get_motif_matrix();

	write_sv(mat, out, delimiter, quote, rownames, colnames);
};

void write_sv(
	const GenomicRange& gr,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = gr.rows();
	const int ncol = gr.cols();

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(colnames, delimiter, quote) << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		if (quote) {
			out << "\"" + gr.get_qstring(i, 0) + "\"";

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + gr.get_qstring(i, j) + "\"");
			}
		}
		else {
			out << gr.get_qstring(i, 0);

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + gr.get_qstring(i, j);
			}
		}

		out << Qt::endl;
	}
};

void write_sv(
	const Eigen::ArrayXXd& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = mat.rows();
	const int ncol = mat.cols();

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(colnames, delimiter, quote) << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		if (quote) {
			if (ncol > 0) {
				out << "\"" + QString::number(mat(i, 0)) + "\"";
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + QString::number(mat(i, j)) + "\"");
			}
		}
		else {

			if (ncol > 0) {
				out << QString::number(mat(i, 0));
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + QString::number(mat(i, j));
			}
		}

		out << Qt::endl;
	}
};

void write_sv(
	const Eigen::MatrixXd& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = mat.rows();
	const int ncol = mat.cols();

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(colnames, delimiter, quote) << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		if (quote) {
			if (ncol > 0) {
				out << "\"" + QString::number(mat(i, 0)) + "\"";
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + QString::number(mat(i, j)) + "\"");
			}
		}
		else {

			if (ncol > 0) {
				out << QString::number(mat(i, 0));
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + QString::number(mat(i, j));
			}
		}

		out << Qt::endl;
	}
};

void write_sv(
	const Eigen::SparseMatrix<double>& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = mat.rows();
	const int ncol = mat.cols();

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(colnames, delimiter, quote) << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		Eigen::ArrayXd row = mat.row(i);

		if (quote) {
			if (ncol > 0) {
				out << "\"" + QString::number(row(0)) + "\"";
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + QString::number(row(j)) + "\"");
			}
		}
		else {
			if (ncol > 0) {
				out << QString::number(row(0));
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + QString::number(row(j));
			}
		}

		out << Qt::endl;
	}
};


void write_sv(
	const Eigen::ArrayXXi& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = mat.rows();
	const int ncol = mat.cols();

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(colnames, delimiter, quote) << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		if (quote) {
			if (ncol > 0) {
				out << "\"" + QString::number(mat(i, 0)) + "\"";
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + QString::number(mat(i, j)) + "\"");
			}
		}
		else {
			if (ncol > 0) {
				out << QString::number(mat(i, 0));
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + QString::number(mat(i, j));
			}
		}

		out << Qt::endl;
	}
};

void write_sv(
	const Eigen::MatrixXi& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = mat.rows();
	const int ncol = mat.cols();

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(colnames, delimiter, quote) << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		if (quote) {
			if (ncol > 0) {
				out << "\"" + QString::number(mat(i, 0)) + "\"";
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + QString::number(mat(i, j)) + "\"");
			}
		}
		else {
			if (ncol > 0) {
				out << QString::number(mat(i, 0));
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + QString::number(mat(i, j));
			}
		}

		out << Qt::endl;
	}
};

void write_sv(
	const Eigen::SparseMatrix<int>& mat,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	const QStringList& rownames,
	const QStringList& colnames)
{
	const int nrow = mat.rows();
	const int ncol = mat.cols();

	bool keep_rownames = (rownames.size() == nrow);
	bool keep_colnames = (colnames.size() == ncol);

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(colnames, delimiter, quote) << Qt::endl;
	}

	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + rownames[i] + "\"" << delimiter;
			}
			else {
				out << rownames[i] << delimiter;
			}
		}

		Eigen::ArrayXi row = mat.row(i);

		if (quote) {
			if (ncol > 0) {
				out << "\"" + QString::number(row(0)) + "\"";
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + QString::number(row(j)) + "\"");
			}
		}
		else {
			if (ncol > 0) {
				out << QString::number(row(0));
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + QString::number(row(j));
			}
		}

		out << Qt::endl;
	}
};

void write_sv(
	const DenseDouble& dd,
	QTextStream& out,
	const QChar delimiter,
	bool quote,
	bool keep_rownames,
	bool keep_colnames) {

	if (keep_colnames) {
		if (keep_rownames) {
			out << delimiter;
		}
		out << custom::merge_to_string(dd.colnames_, delimiter, quote) << Qt::endl;
	}

	int nrow = dd.mat_.rows(), ncol = dd.mat_.cols();
	for (int i = 0; i < nrow; ++i) {

		if (keep_rownames) {
			if (quote) {
				out << "\"" + dd.rownames_[i] + "\"" << delimiter;
			}
			else {
				out << dd.rownames_[i] << delimiter;
			}
		}
		if (quote) {
			if (ncol > 0) {
				out << "\"" + QString::number(dd.mat_(i, 0)) + "\"";
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + ("\"" + QString::number(dd.mat_(i, j)) + "\"");
			}
		}
		else {
			if (ncol > 0) {
				out << QString::number(dd.mat_(i, 0));
			}

			for (int j = 1; j < ncol; ++j) {
				out << delimiter + QString::number(dd.mat_(i, j));
			}
		}

		out << Qt::endl;
	}

	out.flush();
};

CustomMatrix* read_table(const QString& file_name, bool fast) {

	QFile file(file_name);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return nullptr;

	QTextStream in(&file);
	QString line = in.readLine();

	if (line.isNull()) {

		return nullptr;
	}

	QString line2 = in.readLine();
	if (line2.isNull()) {

		return nullptr;
	}

	QStringList data_table;
	auto delimiter = custom::detect_delimiter(line, line2);

	file.close();

	if (fast) {
		return read_sv_fast(file_name, delimiter.toLatin1());
	}
	else {
		return read_sv(file_name, delimiter);
	}
};

CustomMatrix* read_sv_fast(const QString& file_name, const char delimiter) {

	std::ifstream ifs(file_name.toStdString());

	if (!ifs.is_open()) {
		return nullptr;
	}

	std::string content((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

	if (content.empty()) {
		return nullptr;
	}


	std::vector<std::pair<std::size_t, std::size_t>> line_start_end;
	std::size_t size = content.size();

	if (size == 0 || content[0] == '\n') {
		return nullptr;
	}

	std::size_t index{ 0 };
	for (std::size_t i = 0; i < size; ++i) {
		if (content[i] == '\n') {
			std::size_t line_end = i;
			if (content[i - 1] == '\r') {
				--line_end;
			}
			line_start_end.emplace_back(index, line_end);
			index = i + 1;
		}
	}
	if (size > index) {
		line_start_end.emplace_back(index, size);
	}

	int nrow = line_start_end.size();
	if (nrow < 2) {
		return nullptr;
	}

	std::vector<std::vector<std::pair<std::size_t, std::size_t>>> string_loc(nrow);

#pragma omp parallel for
	for (int i = 0; i < nrow; ++i) {
		custom::digest(content, delimiter, line_start_end[i].first, line_start_end[i].second, string_loc[i]);
	}

	if (!custom::is_same(custom::sapply(string_loc, [](auto&& v) {return v.size(); }))) {
		return nullptr;
	}

	int ncol = string_loc[0].size();
	if (ncol == 0) {
		return nullptr;
	}

	bool indent{ false };
	if (string_loc[0][0].first == string_loc[0][0].second) {

		if (ncol == 1) {
			return nullptr;
		}

		indent = true;
	}

	QStringList colnames, rownames;
	CustomMatrix* ret = new CustomMatrix();
	if (indent) {

		for (int i = 1; i < nrow; ++i) {
			rownames << QString::fromStdString(content.substr(string_loc[i][0].first, string_loc[i][0].second - string_loc[i][0].first));
		}

		for (int i = 1; i < ncol; ++i) {
			colnames << QString::fromStdString(content.substr(string_loc[0][i].first, string_loc[0][i].second - string_loc[0][i].first));
		}

		colnames = custom::make_unique(colnames);

		ret->set_rownames(rownames);

	#pragma omp parallel for
		for (int i = 1; i < ncol; ++i) {

			QStringList col(nrow - 1);

			for (int j = 1; j < nrow; ++j) {
				col[j - 1] = QString::fromStdString(content.substr(string_loc[j][i].first, string_loc[j][i].second - string_loc[j][i].first));
			}
		#pragma omp critical
			{
				ret->update(colnames[i - 1], col);
			}
		}
	}
	else {

		rownames = custom::paste("Row ", custom::cast<QString>(custom::integer_linspaced(nrow, 1, nrow)));

		ret->set_rownames(rownames);

		colnames = custom::paste("Column ", custom::cast<QString>(custom::integer_linspaced(ncol, 1, ncol)));

	#pragma omp parallel for
		for (int i = 0; i < ncol; ++i) {

			QStringList col(nrow);

			for (int j = 0; j < nrow; ++j) {
				col[j] = QString::fromStdString(content.substr(string_loc[j][i].first, string_loc[j][i].second - string_loc[j][i].first));
			}

		#pragma omp critical
			{
				ret->update(colnames[i], col);
			}
		}
	}

	return ret;

};

CustomMatrix* read_sv(const QString& file_name, const QChar& delimiter) {

	QFile file(file_name);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		return nullptr;
	}

	QTextStream in(&file);	

	QString content = in.readAll();
	if (content.isNull()) {
		return nullptr;
	}

	QVector<QPair<qsizetype, qsizetype>> line_start_end;
	qsizetype size = content.size();

	if (size == 0 || content[0] == '\n') {
		return nullptr;
	}

	qsizetype index{ 0 };
	for (qsizetype i = 0; i < size; ++i) {
		if (content[i] == '\n') {
			qsizetype line_end = i;
			if (content[i - 1] == '\r') {
				--line_end;
			}
			line_start_end << qMakePair(index, line_end);
			index = i + 1;
		}
	}
	if (size > index) {
		line_start_end << qMakePair(index, size);
	}

	int nrow = line_start_end.size();
	if (nrow < 2) {
		return nullptr;
	}

	std::vector<std::vector<std::pair<qsizetype, qsizetype>>> string_loc(nrow);

#pragma omp parallel for
	for (int i = 0; i < nrow; ++i) {
		custom::digest(content, delimiter, line_start_end[i].first, line_start_end[i].second, string_loc[i]);
	}

	if (!custom::is_same(custom::sapply(string_loc, [](auto&& v) {return v.size(); }))) {
		return nullptr;
	}

	int ncol = string_loc[0].size();
	if (ncol == 0) {
		return nullptr;
	}

	bool indent{ false };
	if (string_loc[0][0].first == string_loc[0][0].second) {

		if (ncol == 1) {
			return nullptr;
		}

		indent = true;
	}

	QStringList colnames, rownames;
	CustomMatrix* ret = new CustomMatrix();
	if (indent) {

		for (int i = 1; i < nrow; ++i) {
			rownames << content.sliced(string_loc[i][0].first, string_loc[i][0].second - string_loc[i][0].first);
		}

		for (int i = 1; i < ncol; ++i) {
			colnames << content.sliced(string_loc[0][i].first, string_loc[0][i].second - string_loc[0][i].first);
		}

		colnames = custom::make_unique(colnames);

		ret->set_rownames(rownames);

	#pragma omp parallel for
		for (int i = 1; i < ncol; ++i) {

			QStringList col(nrow - 1);

			for (int j = 1; j < nrow; ++j) {
				col[j - 1] = content.sliced(string_loc[j][i].first, string_loc[j][i].second - string_loc[j][i].first);
			}
		#pragma omp critical
			{
				ret->update(colnames[i - 1], col);
			}
		}
	}
	else {

		rownames = custom::paste("Row ", custom::cast<QString>(custom::integer_linspaced(nrow, 1, nrow)));

		ret->set_rownames(rownames);

		colnames = custom::paste("Column ", custom::cast<QString>(custom::integer_linspaced(ncol, 1, ncol)));

	#pragma omp parallel for
		for (int i = 0; i < ncol; ++i) {

			QStringList col(nrow);

			for (int j = 0; j < nrow; ++j) {
				col[j] = content.sliced(string_loc[j][i].first, string_loc[j][i].second - string_loc[j][i].first);
			}

		#pragma omp critical
			{
				ret->update(colnames[i], col);
			}
		}
	}

	return ret;
};
