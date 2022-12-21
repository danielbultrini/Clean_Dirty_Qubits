# %%
import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
from typing import Callable

# from matlab_caller import run_parameters
# Code to read all available dataframes and concatinating them in order
# to look for missing steps and then plot the resutls

readable_names = {
    "x_s0i0": "Pauli X",
    "y_s0i0": "Pauli Y",
    "z_s0i0": "Pauli Z",
    "depol_s0i0": "Depolarizing Noise",
    "random_s0i0": "Random Pauli Noise 1",
    "random_s1i0": "Random Pauli Noise 2",
    "random_s2i0": "Random Pauli Noise 3",
    "random_s3i0": "Random Pauli Noise 4",
    "random_s4i0": "Random Pauli Noise 5",
    "realistic": "Realistic Noise",
    "NumCleanQubits": "Clean Qubits",
    "Err": "Error Rate",
    "AbsGrad": "Absolute Gradient",
    "noise_type_ID": "Noise Type",
    "layerXdirtyXerr": "Total Error",
    "NumDirtyQubits": "Noisy Qubits",
    "layers": "Layers",
}


def load_data(folder: str, filter: str = ""):
    df = pd.DataFrame()

    # for file in os.listdir(folder):
    #     if not file.endswith(".txt"):
    #         df = pd.concat([df, pd.read_pickle(os.path.join(folder, file))])
    files = glob.iglob(os.path.join(folder, f"*{filter}*.pkl"))
    df = pd.concat((pd.read_pickle(f) for f in files), ignore_index=True)
    df = df.reset_index(drop=True)
    df["AbsGrad"] = pd.to_numeric(np.divide(np.abs(df["Grad"]), df["qubits"]))
    df["NumDirtyQubits"] = pd.to_numeric(df["qubits"] - df["NumCleanQubits"])
    df["Qdirty*p"] = pd.to_numeric(np.multiply(df["NumDirtyQubits"], df["layers"]))
    df["drityXerr"] = pd.to_numeric(np.multiply(df["Err"], df["NumDirtyQubits"]))
    df["noise_type_ID"] = (
        df.noise_type.map(str)
        + "_s"
        + df.seed_error.map(str)
        + "i"
        + df.seed_iter.map(str)
    )
    if "ealistic" in df["noise_type"][0]:
        df["layerXdirtyXerr"] = df.apply(
            lambda row: total_error_realistic(
                row["qubits"], row["NumCleanQubits"], row["layers"]
            )
            * row["Err"],
            axis=1,
        )
        df["Scaled Error"] = df["layerXdirtyXerr"] / df["layerXdirtyXerr"].max()
    else:
        df["layerXdirtyXerr"] = pd.to_numeric(
            3 * np.multiply(np.multiply(df["Err"], df["NumDirtyQubits"]), df["layers"])
        )
        temp_group = df.groupby(["qubits"])
        df["Scaled Error"] = df["layerXdirtyXerr"] / temp_group[
            "layerXdirtyXerr"
        ].transform("max")

    return df


def facet_plot(
    dataframe,
    x: str,
    y: str,
    hue: str = None,
    style: str = None,
    col: int = None,
    row: int = None,
    kind: str = "line",
    col_wrap: int = None,
    estimator = np.median,
    ci: str = None,
    height: float = 5,
    aspect: float = 1.5,
):
    """Returns a facet plot with taylored presets

    Args:
        dataframe (pandas dataframe): Project data
        x (str): x axis
        y (str): y axis
        hue (str, optional): color key from dataframe. Defaults to None.
        style (str, optional): style key from dataframe. Defaults to None.
        col (int, optional): number of columns. Defaults to None.
        row (int, optional): number of rows. Defaults to None.
        kind (str, optional): type of graph. Defaults to "line".
        col_wrap (int, optional): wrapping, not usable with col/row. Defaults to None.
        estimator (function or string, optional): type of computation to perform on the
                                                different seeds. Defaults to np.median.

    Returns:
        axis object: The plot from sns as an axis object.
    """
    ax = sns.relplot(
        hue=hue,
        data=dataframe,
        x=x,
        y=y,
        col=col,
        row=row,
        kind=kind,
        style=style,
        col_wrap=col_wrap,
        estimator=estimator,
        ci=ci,
        height=height,
        aspect=aspect,
        palette='flare'
    )

    return ax


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
                "p_a": 1e-4,
            },
            "double": {
                "p_xx": 1e-3,
                "p_h": 1.25e-3,
                "p_d12": 2 * (7.5e-4),
                "p_dep": 2 * (8e-4),
            },
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


def prepare_dataframe(folder_path: str, filter: str = ""):
    """Prepares all files in a folder into a friendly dataframes

    Args:
        folder_path (str): folder with generated pickle files.

    Returns:
        dataframe: dataframe for further processing/plotting
    """
    df = load_data(folder_path, filter)
    df = df.rename(columns=readable_names)
    for old, new in readable_names.items():
        df["Noise Type"] = df["Noise Type"].str.replace(old, new, regex=False)
    return df

