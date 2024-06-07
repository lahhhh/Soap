#pragma once

class SingleCellRna;
class SingleCellAtac;
class SingleCellMultiome;
class BulkRna;
class DataField;
class DenseDouble;
class DenseInt;
class SparseDouble;
class SparseInt;
class Metadata;
class Embedding;
class Enrichment;
class DataFrame;
class GSEA;
class CellChat;
class CNV;
class GenomicRange;
class MotifPosition;
class Footprint;
class CoverageTrack;
class Fragments;
class VelocytoBase;
class VelocityEstimate;
class ScveloEstimate;
class StringVector;
class IntegerVector;
class NumericVector;
class NumericMatrix;
class GeneName;
class Pando;
class DifferentialAnalysis;
class Monocle3;
class ChromVAR;
class Cicero;

class VariableItem;
class SingleCellRnaItem;
class BulkRnaItem;
class SingleCellAtacItem;
class SingleCellMultiomeItem;
class DataFieldItem;
class DenseDoubleItem;
class DenseIntItem;
class SparseDoubleItem;
class SparseIntItem;
class MetadataItem;
class EmbeddingItem;
class EnrichmentItem;
class DataFrameItem;
class GSEAItem;
class CellChatItem;
class CNVItem;
class GenomicRangeItem;
class MotifPositionItem;
class FootprintItem;
class CoverageTrackItem;
class FragmentsItem;
class VelocytoBaseItem;
class VelocityEstimateItem;
class ScveloEstimateItem;
class StringVectorItem;
class NumericMatrixItem;
class GeneNameItem;
class PandoItem;
class DifferentialAnalysisItem;
class Monocle3Item;
class ChromVARItem;
class CiceroItem;

namespace soap {

	enum class VariableType {
		AnyVariable, 
		BulkRna,
		SingleCellRna, 
		SingleCellMultiome,
		DataField,
		DenseDouble, 
		DenseInt, 
		SparseDouble, 
		SparseInt, 
		Metadata, 
		Embedding, 
		Enrichment, 
		DataFrame, 
		GSEA, 
		CellChat, 
		CNV,
		GenomicRange, 
		MotifPosition, 
		Footprint,
		CoverageTrack, 
		Fragments,
		VelocytoBase,
		VelocityEstimate,
		ScveloEstimate,
		StringVector,
		IntegerVector,
		NumericVector,
		GeneName,
		Pando,
		DifferentialAnalysis,
		NumericMatrix,
		Monocle3,
		ChromVAR,
		Cicero,
		SingleCellAtac
	};

	template <typename T> inline constexpr VariableType type() { return VariableType::AnyVariable; }
	template <> inline constexpr VariableType type<BulkRna>() { return VariableType::BulkRna; }
	template <> inline constexpr VariableType type<SingleCellRna>() { return VariableType::SingleCellRna; }
	template <> inline constexpr VariableType type<SingleCellAtac>() { return VariableType::SingleCellAtac; }
	template <> inline constexpr VariableType type<SingleCellMultiome>() { return VariableType::SingleCellMultiome; }
	template <> inline constexpr VariableType type<DataField>() { return VariableType::DataField; }
	template <> inline constexpr VariableType type<DenseDouble>() { return VariableType::DenseDouble; }
	template <> inline constexpr VariableType type<DenseInt>() { return VariableType::DenseInt; }
	template <> inline constexpr VariableType type<SparseDouble>() { return VariableType::SparseDouble; }
	template <> inline constexpr VariableType type<SparseInt>() { return VariableType::SparseInt; }
	template <> inline constexpr VariableType type<Metadata>() { return VariableType::Metadata; }
	template <> inline constexpr VariableType type<Embedding>() { return VariableType::Embedding; }
	template <> inline constexpr VariableType type<Enrichment>() { return VariableType::Enrichment; }
	template <> inline constexpr VariableType type<DataFrame>() { return VariableType::DataFrame; }
	template <> inline constexpr VariableType type<GSEA>() { return VariableType::GSEA; }
	template <> inline constexpr VariableType type<CellChat>() { return VariableType::CellChat; }
	template <> inline constexpr VariableType type<CNV>() { return VariableType::CNV; }
	template <> inline constexpr VariableType type<GenomicRange>() { return VariableType::GenomicRange; }
	template <> inline constexpr VariableType type<MotifPosition>() { return VariableType::MotifPosition; }
	template <> inline constexpr VariableType type<Footprint>() { return VariableType::Footprint; }
	template <> inline constexpr VariableType type<CoverageTrack>() { return VariableType::CoverageTrack; }
	template <> inline constexpr VariableType type<Fragments>() { return VariableType::Fragments; }
	template <> inline constexpr VariableType type<VelocytoBase>() { return VariableType::VelocytoBase; }
	template <> inline constexpr VariableType type<VelocityEstimate>() { return VariableType::VelocityEstimate; }
	template <> inline constexpr VariableType type<ScveloEstimate>() { return VariableType::ScveloEstimate; }
	template <> inline constexpr VariableType type<StringVector>() { return VariableType::StringVector; }
	template <> inline constexpr VariableType type<IntegerVector>() { return VariableType::IntegerVector; }
	template <> inline constexpr VariableType type<NumericVector>() { return VariableType::NumericVector; }
	template <> inline constexpr VariableType type<GeneName>() { return VariableType::GeneName; }
	template <> inline constexpr VariableType type<Pando>() { return VariableType::Pando; }
	template <> inline constexpr VariableType type<DifferentialAnalysis>() { return VariableType::DifferentialAnalysis; }
	template <> inline constexpr VariableType type<NumericMatrix>() { return VariableType::NumericMatrix; }
	template <> inline constexpr VariableType type<Monocle3>() { return VariableType::Monocle3; }
	template <> inline constexpr VariableType type<ChromVAR>() { return VariableType::ChromVAR; }
	template <> inline constexpr VariableType type<Cicero>() { return VariableType::Cicero; }

	template <typename T> class get_item_type { public:	using type = VariableItem; };

	template <> class get_item_type<SingleCellRna> { public: using type = SingleCellRnaItem; };
	template <> class get_item_type<BulkRna> { public: using type = BulkRnaItem; };
	template <> class get_item_type<SingleCellAtac> { public: using type = SingleCellAtacItem; };
	template <> class get_item_type<SingleCellMultiome> { public: using type = SingleCellMultiomeItem; };
	template <> class get_item_type<DataField> { public: using type = DataFieldItem; };
	template <> class get_item_type<DenseDouble> { public: using type = DenseDoubleItem; };
	template <> class get_item_type<DenseInt> { public: using type = DenseIntItem; };
	template <> class get_item_type<SparseDouble> { public: using type = SparseDoubleItem; };
	template <> class get_item_type<SparseInt> { public: using type = SparseIntItem; };
	template <> class get_item_type<Metadata> { public: using type = MetadataItem; };
	template <> class get_item_type<Embedding> { public: using type = EmbeddingItem; };
	template <> class get_item_type<Enrichment> { public: using type = EnrichmentItem; };
	template <> class get_item_type<DataFrame> { public: using type = DataFrameItem; };
	template <> class get_item_type<GSEA> { public: using type = GSEAItem; };
	template <> class get_item_type<CellChat> { public: using type = CellChatItem; };
	template <> class get_item_type<CNV> { public: using type = CNVItem; };
	template <> class get_item_type<GenomicRange> { public: using type = GenomicRangeItem; };
	template <> class get_item_type<MotifPosition> { public: using type = MotifPositionItem; };
	template <> class get_item_type<Footprint> { public: using type = FootprintItem; };
	template <> class get_item_type<CoverageTrack> { public: using type = CoverageTrackItem; };
	template <> class get_item_type<Fragments> { public: using type = FragmentsItem; };
	template <> class get_item_type<VelocytoBase> { public: using type = VelocytoBaseItem; };
	template <> class get_item_type<VelocityEstimate> { public: using type = VelocityEstimateItem; };
	template <> class get_item_type<ScveloEstimate> { public: using type = ScveloEstimateItem; };
	template <> class get_item_type<StringVector> { public: using type = StringVectorItem; };
	template <> class get_item_type<NumericMatrix> { public: using type = NumericMatrixItem; };
	template <> class get_item_type<GeneName> { public: using type = GeneNameItem; };
	template <> class get_item_type<Pando> { public: using type = PandoItem; };
	template <> class get_item_type<DifferentialAnalysis> { public: using type = DifferentialAnalysisItem; };
	template <> class get_item_type<Monocle3> { public: using type = Monocle3Item; };
	template <> class get_item_type<ChromVAR> { public: using type = ChromVARItem; };
	template <> class get_item_type<Cicero> { public: using type = CiceroItem; };
													 
	template <typename T>
	using item_type = get_item_type<T>::type;

	template <typename T> class get_data_type { public:	using type = void; };

	template <> class get_data_type<SingleCellRnaItem> { public: using type = SingleCellRna; };
	template <> class get_data_type<BulkRnaItem> { public: using type = BulkRna; };
	template <> class get_data_type<SingleCellAtacItem> { public: using type = SingleCellAtac; };
	template <> class get_data_type<SingleCellMultiomeItem> { public: using type = SingleCellMultiome; };
	template <> class get_data_type<DataFieldItem> { public: using type = DataField; };
	template <> class get_data_type<DenseDoubleItem> { public: using type = DenseDouble; };
	template <> class get_data_type<DenseIntItem> { public: using type = DenseInt; };
	template <> class get_data_type<SparseDoubleItem> { public: using type = SparseDouble; };
	template <> class get_data_type<SparseIntItem> { public: using type = SparseInt; };
	template <> class get_data_type<MetadataItem> { public: using type = Metadata; };
	template <> class get_data_type<EmbeddingItem> { public: using type = Embedding; };
	template <> class get_data_type<EnrichmentItem> { public: using type = Enrichment; };
	template <> class get_data_type<DataFrameItem> { public: using type = DataFrame; };
	template <> class get_data_type<GSEAItem> { public: using type = GSEA; };
	template <> class get_data_type<CellChatItem> { public: using type = CellChat; };
	template <> class get_data_type<CNVItem> { public: using type = CNV; };
	template <> class get_data_type<GenomicRangeItem> { public: using type = GenomicRange; };
	template <> class get_data_type<MotifPositionItem> { public: using type = MotifPosition; };
	template <> class get_data_type<FootprintItem> { public: using type = Footprint; };
	template <> class get_data_type<CoverageTrackItem> { public: using type = CoverageTrack; };
	template <> class get_data_type<FragmentsItem> { public: using type = Fragments; };
	template <> class get_data_type<VelocytoBaseItem> { public: using type = VelocytoBase; };
	template <> class get_data_type<VelocityEstimateItem> { public: using type = VelocityEstimate; };
	template <> class get_data_type<ScveloEstimateItem> { public: using type = ScveloEstimate; };
	template <> class get_data_type<StringVectorItem> { public: using type = StringVector; };
	template <> class get_data_type<NumericMatrixItem> { public: using type = NumericMatrix; };
	template <> class get_data_type<GeneNameItem> { public: using type = GeneName; };
	template <> class get_data_type<PandoItem> { public: using type = Pando; };
	template <> class get_data_type<DifferentialAnalysisItem> { public: using type = DifferentialAnalysis; };
	template <> class get_data_type<Monocle3Item> { public: using type = Monocle3; };
	template <> class get_data_type<ChromVARItem> { public: using type = ChromVAR; };
	template <> class get_data_type<CiceroItem> { public: using type = Cicero; };

	template <typename T>
	using data_type = get_data_type<T>::type;

	enum class TableType { 
		EigenArrayXXd, 
		EigenArrayXXi, 
		EigenMatrixXd, 
		EigenMatrixXi, 
		CustomMatrix, 
		EigenSparseMatrixInt, 
		EigenSparseMatrixDouble, 
		GenomicRange, 
		MotifPosition,
		StringList
	};

	enum class Species { 
		Undefined = 0, 
		Human = 1, 
		Mouse = 2 
	};

	enum class InputStyle {
		IntegerLineEdit,
		NumericLineEdit,
		StringLineEdit,
		ComboBox,
		SwitchButton, 
		SimpleChoice, 
		FactorChoice, 
		MultipleLineEdit, 
		MultipleLineEditWithCompleter,
		LineEditWithCompleter,
		LineEditWithCompleter2,
		SimpleChoiceWithLineEdit, 
		FactorChoiceWithLineEdit,
		MultipleDoubleLineEditWithCompleterLayout,
		SimpleChoiceWithLineEditAndColorChoiceLayout, 
		FactorDoubleLineEditWithCompleterLayout,
		MultiFile,
		MultiLineEditWithColorChoiceLayout,
		MultiCheckBox,
		TextEdit,
		CompareLayout,
		ColorChoice,
		ChooseOpenFile,
		ChooseSaveFile,
		LogicLayout
	};

	enum class FileType { 
		Table,
		CSV, 
		TSV 
	};
}