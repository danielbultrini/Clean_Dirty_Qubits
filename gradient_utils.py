from qiskit import transpile
from qiskit.circuit import Parameter
from qiskit.opflow import StateFn, CircuitStateFn
from qiskit.opflow.gradients import Gradient
from qiskit.quantum_info import Statevector

import mthree
import numpy as np

import hamiltonian_utils as hu
import circuit_utils as cu


def qiskit_builtin_gradients(H, n_qubits, n_layers, values):
    """Returns parameter shift gradients from qiskit internal functions.

    Args:
        H (qiskit hamiltonian): Hamiltonian of the system
        n_qubits (int): number of qubits
        n_layers (int): layers of circuit
        values (parameter values): values of parameters for HAV circuit

    Returns:
        state_grad_results: list of gradients of result
    """
    params = [Parameter("\\phi" + str(i)) for i in range(2 * n_layers)]
    param_values = {param: value for (param, value) in zip(params, values)}
    qc = cu.hva_circuit(n_qubits, n_layers, params, measure=False)
    op = ~StateFn(H) @ CircuitStateFn(primitive=qc, coeff=1.0)
    state_grad = Gradient(grad_method="param_shift").convert(operator=op, params=params)
    # Assign the parameters and evaluate the gradient
    state_grad_result = state_grad.assign_parameters(param_values).eval()
    return state_grad_result


def parameter_shift_gradients_hva(
    backend, n_qubits, n_layers, initial_param_values, shift, prefactor
):
    """Returns gradients for HVA circuit

    Args:
        backend (qskit.provideri): backend used
        n_qubits (int): num qubits
        n_layers (int): number of layers
        initial_param_values (list): parameters for the circuit
        shift (list): shift values for parameter shift rule
        prefactor (float): prefactor multiple

    Returns:
        list: list of gradients for each parameter
    """
    # This should go somewhere else
    basis_gates = ["rx", "ry", "rz", "rxx"]

    # Construct operator strings to compute expectation values
    ops = hu.tfim_1d_ops(n_qubits)
    coeffs = -np.ones(2 * n_qubits)
    meas_qcs = cu.opstr_to_meas_circ(ops)
    meas_ops = [string.replace("X", "Z").replace("Y", "Z") for string in ops]

    H = hu.tfim_1d_H(n_qubits, 1).to_matrix_op().to_matrix()
    gradients_for_each_param = []
    num_dummy_params = n_qubits

    for i in range(2 * n_layers):  # i.e. for each circuit parameter
        # Create list of Parameter objects with distinct parameters in the target sub-layer
        params = (
            [Parameter("\\theta" + str(p)) for p in range(i)]
            + [
                Parameter("\\theta" + str(i) + "_" + str(p))
                for p in range(i, i + num_dummy_params)
            ]
            + [Parameter("\\theta" + str(p)) for p in range(i + 1, 2 * n_layers)]
        )
        # Create list of parameter values with the target sub-layer parameters duplicated
        param_values = (
            [v for v in initial_param_values[:i]]
            + [initial_param_values[i] for _ in range(num_dummy_params)]
            + [v for v in initial_param_values[i + 1 :]]
        )

        # the gradient for the ith parameter is defined as the sum over the gradients for the 'num_dummy_params' dummy parameters associated with the ith parameter (product rule)
        dummy_parameter_gradient_sum = 0
        # construct the HVA circuit with distinct parameters in the target sub-layer
        qc = cu.param_shift_hva_circuit(n_qubits, n_layers, params, num_dummy_params, i)
        plus_qcs = []
        minus_qcs = []

        # compute the gradient of the cost function wrt each dummy parameter in the target sub-layer
        if backend.name() == "statevector_simulator":
            for j in range(num_dummy_params):
                # shift the target dummy parameter
                plus_params = {
                    param: value + shift if k == i + j else value
                    for param, value, k in zip(params, param_values, range(len(params)))
                }
                minus_params = {
                    param: value - shift if k == i + j else value
                    for param, value, k in zip(params, param_values, range(len(params)))
                }

                # assign the shifted parameters, producing copies of the circuit
                param_shift_plus_qc = qc.assign_parameters(plus_params)
                param_shift_minus_qc = qc.assign_parameters(minus_params)

                # construct measurement circuits
                full_plus_qcs = [
                    param_shift_plus_qc.compose(meas_qc) for meas_qc in meas_qcs
                ]
                full_minus_qcs = [
                    param_shift_minus_qc.compose(meas_qc) for meas_qc in meas_qcs
                ]

                exact_svs_plus = [
                    Statevector(backend.run(full_plus_qc).result().get_statevector())
                    for full_plus_qc in full_plus_qcs
                ]
                exact_svs_minus = [
                    Statevector(backend.run(full_minus_qc).result().get_statevector())
                    for full_minus_qc in full_minus_qcs
                ]

                exact_probs_plus = [
                    mthree.classes.ProbDistribution(
                        {
                            f"{index:0{n_qubits}b}": prob
                            for (index, prob) in zip(
                                range(2**n_qubits), sv.probabilities()
                            )
                        }
                    )
                    for sv in exact_svs_plus
                ]
                exact_probs_minus = [
                    mthree.classes.ProbDistribution(
                        {
                            f"{index:0{n_qubits}b}": prob
                            for (index, prob) in zip(
                                range(2**n_qubits), sv.probabilities()
                            )
                        }
                    )
                    for sv in exact_svs_minus
                ]

                # convert counts to expectation values
                exp_plus = np.array(
                    [
                        exact_prob_plus.expval(meas_op)
                        for (exact_prob_plus, meas_op) in zip(
                            exact_probs_plus, meas_ops
                        )
                    ]
                )
                exp_minus = np.array(
                    [
                        exact_prob_minus.expval(meas_op)
                        for (exact_prob_minus, meas_op) in zip(
                            exact_probs_minus, meas_ops
                        )
                    ]
                )

                # sum to compute cost
                cost_plus = np.sum(coeffs * exp_plus)
                cost_minus = np.sum(coeffs * exp_minus)

                dummy_parameter_gradient_sum += prefactor * (cost_plus - cost_minus)
        elif backend.name() == "aer_simulator":
            # compute the gradient of the cost function wrt each dummy parameter in the target sub-layer
            for j in range(num_dummy_params):
                # shift the target dummy parameter
                plus_params = {
                    param: value + shift if k == i + j else value
                    for param, value, k in zip(params, param_values, range(len(params)))
                }
                minus_params = {
                    param: value - shift if k == i + j else value
                    for param, value, k in zip(params, param_values, range(len(params)))
                }

                # assign the shifted parameters, producing copies of the circuit
                param_shift_plus_qc = transpile(
                    qc.assign_parameters(plus_params),
                    backend=backend,
                    basis_gates=basis_gates,
                    optimization_level=0
                )
                param_shift_minus_qc = transpile(
                    qc.assign_parameters(minus_params),
                    backend=backend,
                    basis_gates=basis_gates,
                    optimization_level=0
                )

                # construct measurement circuits
                full_plus_qcs = [
                    param_shift_plus_qc.compose(meas_qc).measure_all(inplace=False)
                    for meas_qc in meas_qcs
                ]
                full_minus_qcs = [
                    param_shift_minus_qc.compose(meas_qc).measure_all(inplace=False)
                    for meas_qc in meas_qcs
                ]
                plus_qcs += full_plus_qcs
                minus_qcs += full_minus_qcs

            result = backend.run(plus_qcs + minus_qcs, shots).result()

            for j in range(num_dummy_params):
                counts_plus = result.get_counts()[
                    j * len(meas_qcs) : j * len(meas_qcs) + len(meas_qcs)
                ]
                counts_minus = result.get_counts()[
                    j * len(meas_qcs)
                    + num_dummy_params * len(meas_qcs) : j * len(meas_qcs)
                    + num_dummy_params * len(meas_qcs)
                    + len(meas_qcs)
                ]

                # convert counts to expectation values
                exp_plus = mthree.utils.expval(counts_plus, meas_ops)
                exp_minus = mthree.utils.expval(counts_minus, meas_ops)

                # sum to compute cost
                cost_plus = np.sum(coeffs * exp_plus).real
                cost_minus = np.sum(coeffs * exp_minus).real

                dummy_parameter_gradient_sum += prefactor * (cost_plus - cost_minus)
        elif backend.name() == "aer_simulator_density_matrix":
            for j in range(num_dummy_params):
                # shift the target dummy parameter
                plus_params = {
                    param: value + shift if k == i + j else value
                    for param, value, k in zip(params, param_values, range(len(params)))
                }
                minus_params = {
                    param: value - shift if k == i + j else value
                    for param, value, k in zip(params, param_values, range(len(params)))
                }

                # assign the shifted parameters, producing copies of the circuit
                param_shift_plus_qc = qc.assign_parameters(plus_params)
                param_shift_minus_qc = qc.assign_parameters(minus_params)

                param_shift_plus_qc.save_density_matrix()
                param_shift_minus_qc.save_density_matrix()

                full_exact_dm_plus = (
                    backend.run(param_shift_plus_qc).result().data()["density_matrix"]
                )
                full_exact_dm_minus = (
                    backend.run(param_shift_minus_qc).result().data()["density_matrix"]
                )
                cost_plus = np.trace(np.matmul(H, full_exact_dm_plus)).real
                cost_minus = np.trace(np.matmul(H, full_exact_dm_minus)).real

                dummy_parameter_gradient_sum += prefactor * (cost_plus - cost_minus)

        gradients_for_each_param.append(dummy_parameter_gradient_sum)

    return gradients_for_each_param
