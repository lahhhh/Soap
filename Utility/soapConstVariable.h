#pragma once

#include <QSize>
#include <QStringList>

namespace soap {

	constexpr int LargeWidth = 250;
	constexpr int MiddleWidth = 150;
	constexpr int SmallWidth = 100;

	constexpr int LargeHeight = 50;
	constexpr int MiddleHeight = 30;
	constexpr int SmallHeight = 20;

	constexpr QSize LargeSize = QSize(LargeWidth, LargeHeight);
	constexpr QSize MiddleSize = QSize(MiddleWidth, MiddleHeight);
	constexpr QSize SmallSize = QSize(SmallWidth, SmallHeight);

	const QStringList HumanChromosomeOrder{ 
		"1", 
		"2", 
		"3",
		"4", 
		"5", 
		"6", 
		"7",
		"8", 
		"9", 
		"10", 
		"11", 
		"12", 
		"13", 
		"14",
		"15",
		"16",
		"17", 
		"18",
		"19", 
		"20", 
		"21", 
		"22" 
	};

	const QStringList MouseChromosomeOrder{ 
		"1",
		"2",
		"3",
		"4",
		"5",
		"6",
		"7",
		"8",
		"9",
		"10",
		"11",
		"12",
		"13",
		"14",
		"15",
		"16",
		"17",
		"18",
		"19" 
	};
};