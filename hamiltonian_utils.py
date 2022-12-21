from qiskit.opflow import I, X, Z, SummedOp, TensoredOp


def tfim_1d_ops(n_qubits):
    """Returns list of measurements in operator format"""
    ops = []

    for n in range(n_qubits - 1):
        ops.append("".join(["Z" if i == n else "I" for i in range(n_qubits)]))
        ops.append(
            "".join(["X" if (i == n or i == n + 1) else "I" for i in range(n_qubits)])
        )

    ops.append("".join(["Z" if i == n_qubits - 1 else "I" for i in range(n_qubits)]))
    ops.append(
        "".join(
            ["X" if (i == 0 or i == n_qubits - 1) else "I" for i in range(n_qubits)]
        )
    )

    return ops


def tfim_1d_H(n_qubits, g):
    """Returns hamiltonian for HVA

    Args:
        n_qubits (int): num of qubits
        g (float): scaling factor for z terms

    Returns:
        Hamiltonian: Returns qiskit form hamiltonian for the system.
    """
    Z_terms = g * SummedOp(
        [
            TensoredOp([Z if i == n else I for i in range(n_qubits)])
            for n in range(n_qubits)
        ]
    )
    XX_terms = SummedOp(
        [
            TensoredOp([X if (i == n or i == n + 1) else I for i in range(n_qubits)])
            for n in range(n_qubits - 1)
        ]
    )
    if n_qubits > 2:
        XX_terms += TensoredOp(
            [X if (i == 0 or i == n_qubits - 1) else I for i in range(n_qubits)]
        )  # periodic boundary conditions
    H = -XX_terms - Z_terms
    # return Operator(H.reduce())
    return H.reduce()
