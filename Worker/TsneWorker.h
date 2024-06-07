#pragma once

#include "Identifier.h"

class TsneWorker
    : public QObject
{
    Q_OBJECT
public:

    TsneWorker(
        const Eigen::MatrixXd& mat,
        unsigned int random_state = 1997,
        int max_iter = 1000,
        double perplexity = 30.0,
        double theta = 0.5,
        double momentum = 0.5,
        double final_momentum = 0.8,
        int stop_lying_iter = 250,
        int momentum_switch_iter = 250,
        double eta = 200.0,
        double exaggeration_factor = 12.0
    );

	Eigen::MatrixXd mat_;

    unsigned int random_state_;

	int max_iter_;
    double perplexity_;
    double theta_;
	double momentum_;
	double final_momentum_;
	int stop_lying_iter_;
	int momentum_switch_iter_;
	double eta_;
	double exaggeration_factor_;

public slots:

    void run();

signals:

    void x_message(QString, int);

    void x_results_ready();

    void x_tsne_ready(Eigen::MatrixXd);
};

