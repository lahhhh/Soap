#include "ItemIOWorker.h"

#include <fstream>

#include "Serialization.h"

const QStringList valid_editions = {SOAP_EDITION_1_0};

ItemIOWorker::ItemIOWorker(
	const QString& file_path, 
	soap::VariableType data_type, 
	void* data, 
	const QString& item_name
) :
	mode_(WorkMode::Write),
	file_path_(file_path),
	item_name_(item_name),
	data_type_(data_type),
	data_(data)
{}

ItemIOWorker::ItemIOWorker(const QString& file_path) :
	mode_(WorkMode::Read), 
	file_path_(file_path)
{}

void ItemIOWorker::write() {

	std::ofstream ofs(this->file_path_.toStdString(), std::ios::binary);

	if (!ofs.is_open()) {
		G_TASK_WARN("File open failed.");
		return;
	}

	int id1 = ITEM_IDENTIFIER_1;
	int id2 = ITEM_IDENTIFIER_2;
	QString edition = SOAP_EDITION_NOW;

	swrite(ofs, id1);
	swrite(ofs, id2);
	swrite(ofs, edition);
	swrite(ofs, this->item_name_);

	if (this->data_type_ == soap::VariableType::SingleCellRna) {
		swrite(ofs, SingleCellRna::g_identifier());
		swrite(ofs, *static_cast<SingleCellRna*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::SingleCellAtac) {
		swrite(ofs, SingleCellAtac::g_identifier());
		swrite(ofs, *static_cast<SingleCellAtac*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::CellChat) {
		swrite(ofs, CellChat::g_identifier());
		swrite(ofs, *static_cast<CellChat*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::CNV) {
		swrite(ofs, CNV::g_identifier());
		swrite(ofs, *static_cast<CNV*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::DataFrame) {
		swrite(ofs, DataFrame::g_identifier());
		swrite(ofs, *static_cast<DataFrame*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::DifferentialAnalysis) {
		swrite(ofs, DifferentialAnalysis::g_identifier());
		swrite(ofs, *static_cast<DifferentialAnalysis*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::Pando) {
		swrite(ofs, Pando::g_identifier());
		swrite(ofs, *static_cast<Pando*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::DenseDouble) {
		swrite(ofs, DenseDouble::g_identifier());
		swrite(ofs, *static_cast<DenseDouble*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::Embedding) {
		swrite(ofs, Embedding::g_identifier());
		swrite(ofs, *static_cast<Embedding*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::Enrichment) {
		swrite(ofs, Enrichment::g_identifier());
		swrite(ofs, *static_cast<Enrichment*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::GSEA) {
		swrite(ofs, GSEA::g_identifier());
		swrite(ofs, *static_cast<GSEA*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::Metadata) {
		swrite(ofs, Metadata::g_identifier());
		swrite(ofs, *static_cast<Metadata*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::SparseDouble) {
		swrite(ofs, SparseDouble::g_identifier());
		swrite(ofs, *static_cast<SparseDouble*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::SparseInt) {
		swrite(ofs, SparseInt::g_identifier());
		swrite(ofs, *static_cast<SparseInt*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::SingleCellMultiome) {
		swrite(ofs, SingleCellMultiome::g_identifier());
		swrite(ofs, *static_cast<SingleCellMultiome*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::GenomicRange) {
		swrite(ofs, GenomicRange::g_identifier());
		swrite(ofs, *static_cast<GenomicRange*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::MotifPosition) {
		swrite(ofs, MotifPosition::g_identifier());
		swrite(ofs, *static_cast<MotifPosition*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::CoverageTrack) {
		swrite(ofs, CoverageTrack::g_identifier());
		swrite(ofs, *static_cast<CoverageTrack*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::Footprint) {
		swrite(ofs, Footprint::g_identifier());
		swrite(ofs, *static_cast<Footprint*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::StringVector) {
		swrite(ofs, StringVector::g_identifier());
		swrite(ofs, *static_cast<StringVector*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::GeneName) {
		swrite(ofs, GeneName::g_identifier());
		swrite(ofs, *static_cast<GeneName*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::NumericMatrix) {
		swrite(ofs, NumericMatrix::g_identifier());
		swrite(ofs, *static_cast<NumericMatrix*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::DenseInt) {
		swrite(ofs, DenseInt::g_identifier());
		swrite(ofs, *static_cast<DenseInt*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::DenseDouble) {
		swrite(ofs, DenseDouble::g_identifier());
		swrite(ofs, *static_cast<DenseDouble*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::ChromVAR) {
		swrite(ofs, ChromVAR::g_identifier());
		swrite(ofs, *static_cast<ChromVAR*>(this->data_), edition);
	}
	else if (this->data_type_ == soap::VariableType::BulkRna) {
		swrite(ofs, BulkRna::g_identifier());
		swrite(ofs, *static_cast<BulkRna*>(this->data_), edition);
	}
	else {
		G_TASK_WARN("Unknown type");
	}
	G_TASK_LOG("Item saved.");
};

void ItemIOWorker::read() {

	std::ifstream ifs(this->file_path_.toStdString(), std::ios::binary);

	if (!ifs.is_open()) {
		G_TASK_WARN("File open failed.");
		return;
	}

	int id1{ 0 };
	int id2{ 0 };
	QString edition;

	sread(ifs, id1);
	sread(ifs, id2);
	sread(ifs, edition);
	
	if (id1 != ITEM_IDENTIFIER_1 || id2 != ITEM_IDENTIFIER_2) {
		G_TASK_WARN("File is broken.");
		return;
	}

	if (!valid_editions.contains(edition)) {
		G_TASK_WARN("Unrecognized edition! Item Edition : [" + edition + "]");
		return;
	}
	
	sread(ifs, this->item_name_);

	QString type;

	sread(ifs, type);

	if (type == SingleCellRna::g_identifier()) {
		SingleCellRna* data = new SingleCellRna();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::SingleCellRna, this->item_name_);
	}
	else if (type == SingleCellAtac::g_identifier()) {
		SingleCellAtac* data = new SingleCellAtac();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::SingleCellAtac, this->item_name_);
	}
	else if (type == SingleCellMultiome::g_identifier()) {
		SingleCellMultiome* data = new SingleCellMultiome();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::SingleCellMultiome, this->item_name_);
	}
	else if (type == CellChat::g_identifier()) {
		CellChat* data = new CellChat();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::CellChat, this->item_name_);
	}
	else if (type == CNV::g_identifier()) {
		CNV* data = new CNV();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::CNV, this->item_name_);
	}
	else if (type == DataFrame::g_identifier()) {
		DataFrame* data = new DataFrame();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::DataFrame, this->item_name_);
	}
	else if (type == DifferentialAnalysis::g_identifier()) {
		DifferentialAnalysis* data = new DifferentialAnalysis();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::DifferentialAnalysis, this->item_name_);
	}
	else if (type == DenseDouble::g_identifier()) {
		DenseDouble* data = new DenseDouble();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::DenseDouble, this->item_name_);
	}
	else if (type == Embedding::g_identifier()) {
		Embedding* data = new Embedding();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::Embedding, this->item_name_);
	}
	else if (type == Enrichment::g_identifier()) {
		Enrichment* data = new Enrichment();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::Enrichment, this->item_name_);
	}
	else if (type == GSEA::g_identifier()) {
		GSEA* data = new GSEA();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::GSEA, this->item_name_);
	}
	else if (type == Metadata::g_identifier()) {
		Metadata* data = new Metadata();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::Metadata, this->item_name_);
	}
	else if (type == SparseDouble::g_identifier()) {
		SparseDouble* data = new SparseDouble();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::SparseDouble, this->item_name_);
	}
	else if (type == SparseInt::g_identifier()) {
		SparseInt* data = new SparseInt();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::SparseInt, this->item_name_);
	}
	else if (type == GenomicRange::g_identifier()) {
		GenomicRange* data = new GenomicRange();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::GenomicRange, this->item_name_);
	}
	else if (type == MotifPosition::g_identifier()) {
		MotifPosition* data = new MotifPosition();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::MotifPosition, this->item_name_);
	}
	else if (type == CoverageTrack::g_identifier()) {
		CoverageTrack* data = new CoverageTrack();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::CoverageTrack, this->item_name_);
	}
	else if (type == Footprint::g_identifier()) {
		Footprint* data = new Footprint();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::Footprint, this->item_name_);
	}
	else if (type == StringVector::g_identifier()) {
		StringVector* data = new StringVector();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::StringVector, this->item_name_);
	}
	else if (type == GeneName::g_identifier()) {
		GeneName* data = new GeneName();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::GeneName, this->item_name_);
	}
	else if (type == Pando::g_identifier()) {
		Pando* data = new Pando();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::Pando, this->item_name_);
	}
	else if (type == NumericMatrix::g_identifier()) {
		NumericMatrix* data = new NumericMatrix();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::NumericMatrix, this->item_name_);
	}
	else if (type == DenseInt::g_identifier()) {
		DenseInt* data = new DenseInt();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::DenseInt, this->item_name_);
	}
	else if (type == DenseDouble::g_identifier()) {
		DenseDouble* data = new DenseDouble();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::DenseDouble, this->item_name_);
	}
	else if (type == ChromVAR::g_identifier()) {
		ChromVAR* data = new ChromVAR();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::ChromVAR, this->item_name_);
	}
	else if (type == BulkRna::g_identifier()) {
		BulkRna* data = new BulkRna();
		sread(ifs, *data, edition);
		emit x_data_create_soon(data, soap::VariableType::BulkRna, this->item_name_);
	}
	else {
		G_TASK_WARN("Unrecognized item type : " + type);
	}

	G_TASK_LOG("Item Reading finished.");
};

bool ItemIOWorker::work() {

	if (this->mode_ == WorkMode::Read) {
		this->read();
	}
	else {
		this->write();
	}

	return true;
};

void ItemIOWorker::run() {

	if (!this->work()) {
		G_TASK_END;
	}

	G_TASK_END;
}
