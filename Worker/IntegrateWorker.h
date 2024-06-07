#pragma once

#include "Identifier.h"
#include <QTreeWidgetItem>

#include "BulkRna.h"
#include "SingleCellRna.h"
#include "SingleCellAtac.h"
#include "SingleCellMultiome.h"

class IntegrateWorker 
	: public QObject
{
	Q_OBJECT
public:

	IntegrateWorker(
		const QList<const BulkRna*>& data,
		soap::Species species
	) :
		bulkrna_data_(data),
		species_(species),
		mode_(WorkMode::BulkRna)
	{};

	IntegrateWorker(
		const QList<const SingleCellRna*>& data,
		soap::Species species,
		bool distinguish_barcode
	) :
		scrna_data_(data),
		species_(species),
		distinguish_barcode_(distinguish_barcode),
		mode_(WorkMode::SingleCellRna)
	{}

	IntegrateWorker(
		const QList<const SingleCellMultiome* >& data,
		bool distinguish_barcode,
		soap::Species species) :
		scmultiome_data_(data),
		species_(species),
		distinguish_barcode_(distinguish_barcode),
		mode_(WorkMode::SingleCellMultiome)
	{}


	IntegrateWorker(
		const QList<const SingleCellAtac* >& data,
		bool distinguish_barcode,
		soap::Species species) :
		scatac_data_(data),
		species_(species),
		distinguish_barcode_(distinguish_barcode),
		mode_(WorkMode::SingleCellAtac)
	{}

	enum class WorkMode { BulkRna, SingleCellRna, SingleCellAtac, SingleCellMultiome };

	QList<const BulkRna*> bulkrna_data_;

	QList<const SingleCellRna*> scrna_data_;

	QList<const SingleCellMultiome*> scmultiome_data_;

	QList<const SingleCellAtac*> scatac_data_;

	soap::Species species_;

	bool distinguish_barcode_ = true;

	WorkMode mode_ = WorkMode::SingleCellRna;

	std::pair<QStringList, QList<QVector<int>>> get_index(const QList<QStringList>& components, bool distinguish);

	void scrna_mode();

	void scmultiome_mode();

	void scatac_mode();

	void bulkrna_mode();

private:

	void rna_quality_control(SingleCellRna& single_cell_rna);

	void atac_quality_control(SingleCellAtac& single_cell_atac);

	void multiome_quality_control(SingleCellMultiome& single_cell_multiome);

	void integrate_metadata(Metadata& to, QList<const Metadata*> froms);

	bool integrate_fragments_object(
		Fragments& to,
		const QList< QVector<int> >& col_map,
		const QList<const Fragments* >& fragmentss
	);

	void integrate_velocyto(VelocytoBase& to, QList<const VelocytoBase*> froms, bool distinguish_barcode);

	void map_sparse_int_value(
		SparseInt& to,
		QList<SparseInt const*> froms,
		const QList<QVector<int>>& row_maps,
		const QList<QVector<int>>& col_maps,
		const QStringList& rownames,
		const QStringList& colnames);

	void map_dense_int_value(
		DenseInt& to,
		QList<DenseInt const*> froms,
		const QList<QVector<int>>& row_maps,
		const QList<QVector<int>>& col_maps,
		const QStringList& rownames,
		const QStringList& colnames);


public slots:

	void run();

signals:

	void x_message(QString, int);

	void x_results_ready();

	void x_scrna_ready(SingleCellRna*, QList<const SingleCellRna*>);

	void x_scmultiome_ready(SingleCellMultiome*, QList<const SingleCellMultiome*>);

	void x_scatac_ready(SingleCellAtac*, QList<const SingleCellAtac*>);

	void x_bulkrna_ready(BulkRna*, QList<const BulkRna*>);
};
