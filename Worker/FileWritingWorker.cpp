#include "FileWritingWorker.h"

#include <QFile>

#include "FileIO.h"
#include "SingleCellRna.h"
#include "DataFrame.h"

FileWritingWorker::FileWritingWorker(
	void* data,
	soap::VariableType data_type, 
	const QString& file_name, 
	soap::FileType file_type,
	const QStringList& settings
) : 
	data_(data), 
	data_type_(data_type), 
	file_name_(file_name), 
	file_type_(file_type),
	settings_(settings)
{}

void FileWritingWorker::run() {

	QFile file(this->file_name_);
	
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		G_TASK_WARN("Can not create file " + this->file_name_ + " .");
		G_TASK_END;
	}
	
	QTextStream out(&file);

	if (this->data_type_ == soap::VariableType::Metadata) {
		Metadata* metadata = (Metadata*)this->data_;
		if (this->file_type_ == soap::FileType::CSV) {
			write_csv(metadata->mat_, out, this->settings_);
		}
		else if (this->file_type_ == soap::FileType::TSV) {
			write_tsv(metadata->mat_, out, this->settings_);
		}
	}
	else if (this->data_type_ == soap::VariableType::DataFrame) {
		DataFrame* df = (DataFrame*)this->data_;
		if (this->file_type_ == soap::FileType::CSV) {
			write_csv(df->mat_, out, this->settings_);
		}
		else if (this->file_type_ == soap::FileType::TSV) {
			write_tsv(df->mat_, out, this->settings_);
		}
	}
	else if (this->data_type_ == soap::VariableType::Enrichment) {
		Enrichment* en = (Enrichment*)this->data_;
		if (this->file_type_ == soap::FileType::CSV) {
			write_csv(en->mat_, out, this->settings_);
		}
		else if (this->file_type_ == soap::FileType::TSV) {
			write_tsv(en->mat_, out, this->settings_);
		}
	}
	else if (this->data_type_ == soap::VariableType::DifferentialAnalysis) {
		DifferentialAnalysis* da = (DifferentialAnalysis*)this->data_;
		if (this->file_type_ == soap::FileType::CSV) {
			write_csv(da->mat_, out, this->settings_);
		}
		else if (this->file_type_ == soap::FileType::TSV) {
			write_tsv(da->mat_, out, this->settings_);
		}
	}
	else if (this->data_type_ == soap::VariableType::GSEA) {
		GSEA* gsea = (GSEA*)this->data_;
		if (this->file_type_ == soap::FileType::CSV) {
			write_csv(gsea->mat_, out, this->settings_);
		}
		else if (this->file_type_ == soap::FileType::TSV) {
			write_tsv(gsea->mat_, out, this->settings_);
		}
	}
	else if (this->data_type_ == soap::VariableType::Embedding) {
		Embedding* embedding = (Embedding*)this->data_;
		if (this->file_type_ == soap::FileType::CSV) {
			write_csv(embedding->data_, out, this->settings_);
		}
		else if (this->file_type_ == soap::FileType::TSV) {
			write_tsv(embedding->data_, out, this->settings_);
		}
	}
	else {
		G_TASK_WARN("Unsupported File Type.");
		G_TASK_END;
	}

	G_TASK_LOG("File writing succeeds.");
	G_TASK_END;
}
