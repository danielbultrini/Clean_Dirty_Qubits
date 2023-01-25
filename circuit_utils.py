from qiskit import QuantumCircuit, QuantumRegister


def hva_circuit(n_qubits, n_layers, params, measure=False):
    """Generated the HVA ansatz

    Args:
        n_qubits (int): number of qubits
        n_layers (int): layers of ansatz
        params (List[floats]): variational parameters
        measure (bool, optional): Perform measurement?. Defaults to False.

    Returns:
        QuantumCircuit: circuit to run for HVA instance
    """
    qr = QuantumRegister(n_qubits)
    qc = QuantumCircuit(qr)

    for i in range(0, 2 * n_layers, 2):
        for j in range(n_qubits - 1):
            qc.rxx(params[i], j, j + 1)
            qc.rx(0,j)
            qc.rx(0,j+1)

            # qc.rx(0,j)
            # qc.rx(0,j+1)

        qc.rxx(params[i], -1, 0)  # periodic boundary conditions
        # qc.rx(0,-1)
        # qc.rx(0,0)
        qc.rx(0,0)
        qc.rx(0,-1)

        qc.rz(params[i + 1], range(n_qubits))

    if measure:
        return qc.measure_all()
    else:
        return qc


def param_shift_hva_circuit(
    n_qubits, n_layers, params, num_dummy_params, target_param, measure=False
):
    """HVA circuit for parameter shift rule

    Args:
        n_qubits (int): number of qubits in ansatz
        n_layers (int): number of layers
        params (list[float]): list of parameters
        num_dummy_params (int): additional dummy parameters
        target_param (int): which parameter to target?
        measure (bool, optional): Measure?. Defaults to False.

    Returns:
        QuantumCircuit: circuit of HVA for parameter shift rule.
    """
    qr = QuantumRegister(n_qubits)
    qc = QuantumCircuit(qr)

    # construct Hamiltonian variational ansatz for TFQIM in 1D with dummy parameters in sub layer with target parameter
    for j in range(target_param):
        if (j + 1) % 2 == 0:
            qc.rz(params[j], range(n_qubits))
        else:
            for k in range(
                n_qubits - 1
            ):  # add Molmer-Sorenson gates https://en.wikipedia.org/wiki/M%C3%B8lmer%E2%80%93S%C3%B8rensen_gate
                qc.rxx(params[j], k, k + 1)
                qc.rx(0,k)
                qc.rx(0,k+1)

            qc.rxx(params[j], -1, 0)  # periodic boundary conditions
            qc.rx(0,-1)
            qc.rx(0,0)

    if (target_param + 1) % 2 == 0:  # Z layer
        for j in range(target_param, target_param + num_dummy_params):
            qc.rz(params[j], j - target_param)
    else:  # XX layer
        for j in range(target_param, target_param + num_dummy_params - 1):
            qc.rxx(params[j], j - target_param, j - target_param + 1)
            qc.rx(0,j - target_param)
            qc.rx(0,j - target_param+1)

        qc.rxx(
            params[target_param + num_dummy_params - 1], -1, 0
        )
        qc.rx(0,-1)
        qc.rx(0,0)

    for j in range(
        target_param + num_dummy_params, 2 * n_layers - 1 + num_dummy_params
    ):
        if num_dummy_params % 2 == 0:
            if j % 2 == 0:
                qc.rz(params[j], range(n_qubits))
            else:
                for k in range(n_qubits - 1):
                    qc.rxx(params[j], k, k + 1)
                    qc.rx(0,k)
                    qc.rx(0,k+1)

                qc.rxx(params[j], -1, 0)  # periodic boundary conditions
                qc.rx(0,-1)
                qc.rx(0,0)


        else:
            if (j + 1) % 2 == 0:
                qc.rz(params[j], range(n_qubits))
            else:
                for k in range(n_qubits - 1):
                    qc.rxx(params[j], k, k + 1)
                    qc.rx(0,k)
                    qc.rx(0,k+1)

                qc.rxx(params[j], -1, 0)  # periodic boundary conditions
                qc.rx(0,-1)
                qc.rx(0,0)

    if measure:
        return qc.measure_all()
    else:
        return qc


def opstr_to_meas_circ(op_str):
    """Takes a list of operator strings and makes circuit with the correct post-rotations for measurements.

    Parameters:
        op_str (list): List of strings representing the operators needed for measurements.

    Returns:
        list: List of circuits for measurement post-rotations
    """
    num_qubits = len(op_str[0])
    circs = []
    for op in op_str:
        qc = QuantumCircuit(num_qubits)
        for idx, item in enumerate(op):
            if item == "X":
                qc.h(num_qubits - idx - 1)
            elif item == "Y":
                qc.sdg(num_qubits - idx - 1)
                qc.h(num_qubits - idx - 1)
        circs.append(qc)
    return circs
