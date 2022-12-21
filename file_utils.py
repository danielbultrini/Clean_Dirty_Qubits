def exp_filename(exp_params):
    """Returns filename for dataframe

    Args:
        exp_params (dict of experimental parameters): list of experimental parameters from simulator function

    Returns:
        str: filename to be used
    """
    return "results/S_{}-{}_L_{}-{}_Q_{}_{}_{}".format(
        exp_params["SimSeed"][0],
        exp_params["SimSeed"][-1],
        exp_params["layers"][0],
        exp_params["layers"][-1],
        exp_params["qubits"][0],
        exp_params["noise_type"][0],
        exp_params["Err"][0],
    )
