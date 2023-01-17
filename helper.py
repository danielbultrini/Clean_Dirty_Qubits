def total_error_realistic(
    num_qubits: int, clean_qubits: int, layers: int, error_rates: dict = None
):
    """
    calculates total error in realistic model

    Args:
        num_qubits (int): number of qubits
        layers (int): layers in circuit.
        error_rates (dict, optional): Dictionary with 'single' and 'double'
            entries, which have a dict of noise probabilities. Defaults to None
            which defaults to the experiment's model.
    """
    if error_rates is None:
        error_rates = {
            "single": {
                "p_d": 1.5e-4,
                "p_dep": 8e-4,
                "p_d1": 7.5e-4,
                "p_a": 1e-4,
                "p_idle": 8e-4,
            },
            "double": {"p_xx": 1e-3, "p_h": 1.25e-3},
        }
    # total number of gates in circuit
    two_qgate_no = layers * num_qubits
    dirty_singles = layers * (num_qubits - clean_qubits)
    # calculate number of dirty gates for noise
    if clean_qubits == num_qubits:
        dirty_doubles = 0
    elif clean_qubits == 0:
        dirty_doubles = two_qgate_no
    else:
        dirty_doubles = two_qgate_no - (clean_qubits + 1) * layers
    single_error = dirty_singles * sum(error_rates["single"].values())
    double_error = dirty_doubles * sum(error_rates["double"].values())
    return single_error + double_error
