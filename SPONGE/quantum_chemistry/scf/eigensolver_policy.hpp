#pragma once

enum QC_SCF_Eigensolver_Stage
{
    QC_SCF_EIGENSOLVER_OVERLAP = 0,
    QC_SCF_EIGENSOLVER_FOCK = 1,
    QC_SCF_EIGENSOLVER_RI_LOEWNER = 2,
    QC_SCF_EIGENSOLVER_WORKSPACE = 3,
    QC_SCF_EIGENSOLVER_DENSITY_FACTOR = 4,
    QC_SCF_EIGENSOLVER_ENSEMBLE_OCCUPATION = 5,
};

enum QC_SCF_Eigensolver_Channel
{
    QC_SCF_EIGENSOLVER_CHANNEL_OVERLAP = 0,
    QC_SCF_EIGENSOLVER_CHANNEL_ALPHA = 1,
    QC_SCF_EIGENSOLVER_CHANNEL_BETA = 2,
    QC_SCF_EIGENSOLVER_CHANNEL_AUXILIARY = 3,
    QC_SCF_EIGENSOLVER_CHANNEL_SHARED = 4,
};

struct QC_SCF_Eigensolver_Failure
{
    QC_SCF_Eigensolver_Stage stage;
    QC_SCF_Eigensolver_Channel channel;
    const char* stage_name;
    const char* channel_name;
    int dimension;
    int api_status;
    int info;
};

struct QC_SCF_Eigensolver_Workspace_Failure
{
    QC_SCF_Eigensolver_Stage stage;
    QC_SCF_Eigensolver_Channel channel;
    const char* stage_name;
    const char* channel_name;
    int dimension;
    int api_status;
    int workspace_size;
    bool workspace_available;
};

static inline const char* QC_SCF_Eigensolver_Stage_Name(
    QC_SCF_Eigensolver_Stage stage)
{
    switch (stage)
    {
        case QC_SCF_EIGENSOLVER_OVERLAP:
            return "overlap orthogonalization";
        case QC_SCF_EIGENSOLVER_FOCK:
            return "SCF Fock diagonalization";
        case QC_SCF_EIGENSOLVER_RI_LOEWNER:
            return "RI Coulomb-metric Loewner factorization";
        case QC_SCF_EIGENSOLVER_WORKSPACE:
            return "shared double-precision eigensolver workspace query";
        case QC_SCF_EIGENSOLVER_DENSITY_FACTOR:
            return "RI-K spin-density factorization";
        case QC_SCF_EIGENSOLVER_ENSEMBLE_OCCUPATION:
            return "ensemble natural-occupation diagonalization";
    }
    return "unknown eigensolver stage";
}

static inline const char* QC_SCF_Eigensolver_Channel_Name(
    QC_SCF_Eigensolver_Channel channel)
{
    switch (channel)
    {
        case QC_SCF_EIGENSOLVER_CHANNEL_OVERLAP:
            return "overlap";
        case QC_SCF_EIGENSOLVER_CHANNEL_ALPHA:
            return "alpha";
        case QC_SCF_EIGENSOLVER_CHANNEL_BETA:
            return "beta";
        case QC_SCF_EIGENSOLVER_CHANNEL_AUXILIARY:
            return "auxiliary";
        case QC_SCF_EIGENSOLVER_CHANNEL_SHARED:
            return "shared";
    }
    return "unknown";
}

// Keep the policy independent of CONTROLLER so the failure path can be fault
// injected.  Production passes a handler which raises a formatted fatal error;
// tests pass a recorder and verify that no solver output would be consumed.
template <typename FailureHandler>
static inline bool QC_SCF_Require_Eigensolver_Success(
    QC_SCF_Eigensolver_Stage stage, QC_SCF_Eigensolver_Channel channel,
    int dimension, int api_status, int info, FailureHandler failure_handler)
{
    if (api_status == 0 && info == 0) return true;
    const QC_SCF_Eigensolver_Failure failure = {
        stage,
        channel,
        QC_SCF_Eigensolver_Stage_Name(stage),
        QC_SCF_Eigensolver_Channel_Name(channel),
        dimension,
        api_status,
        info,
    };
    failure_handler(failure);
    return false;
}

template <typename FailureHandler>
static inline bool QC_SCF_Require_Eigensolver_Workspace(
    int dimension, int api_status, int workspace_size,
    bool workspace_available, FailureHandler failure_handler)
{
    if (api_status == 0 && workspace_size > 0 && workspace_available)
        return true;
    const QC_SCF_Eigensolver_Workspace_Failure failure = {
        QC_SCF_EIGENSOLVER_WORKSPACE,
        QC_SCF_EIGENSOLVER_CHANNEL_SHARED,
        QC_SCF_Eigensolver_Stage_Name(QC_SCF_EIGENSOLVER_WORKSPACE),
        QC_SCF_Eigensolver_Channel_Name(QC_SCF_EIGENSOLVER_CHANNEL_SHARED),
        dimension,
        api_status,
        workspace_size,
        workspace_available,
    };
    failure_handler(failure);
    return false;
}
