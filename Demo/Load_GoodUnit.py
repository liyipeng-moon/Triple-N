# Load one GoodUnit .mat file and extract image response data.
# Compute per-unit split-half noise ceiling from trial raster data.
# Developer: 王柯盛(Kesheng Wang, https://github.com/JustMuteAll)

import os
from typing import Optional
import h5py
import numpy as np
from scipy.stats import pearsonr

def compute_noise_ceiling(
    resp_mat: np.ndarray,
    stim_idx: np.ndarray,
    n_splits: int = 100,
    random_state: Optional[int] = None,
) -> np.ndarray:
    """
    Compute split-half noise ceiling for each unit and apply
    Spearman-Brown correction.

    Parameters
    ----------
    resp_mat : np.ndarray
        Shape (n_trials, n_units). Trial-by-unit response matrix.
    stim_idx : np.ndarray
        Shape (n_trials,). Stimulus identity for each trial.
    n_splits : int, default=100
        Number of random split-half repetitions.
    random_state : int or None, default=None
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Shape (n_units,). Spearman-Brown corrected noise ceiling
        for each unit.
    """
    resp_mat = np.asarray(resp_mat, dtype=float)
    stim_idx = np.asarray(stim_idx)

    if resp_mat.ndim != 2:
        raise ValueError(f"resp_mat must be 2D, got shape {resp_mat.shape}")
    if stim_idx.ndim != 1:
        raise ValueError(f"stim_idx must be 1D, got shape {stim_idx.shape}")
    if resp_mat.shape[0] != stim_idx.shape[0]:
        raise ValueError(
            f"Number of trials mismatch: resp_mat has {resp_mat.shape[0]} trials, "
            f"but stim_idx has {stim_idx.shape[0]} entries."
        )

    rng = np.random.default_rng(random_state)

    n_trials, n_units = resp_mat.shape
    unique_stim = np.unique(stim_idx)
    n_stim = len(unique_stim)

    corr_all = np.full((n_splits, n_units), np.nan, dtype=float)

    for split in range(n_splits):
        resp1 = np.full((n_stim, n_units), np.nan, dtype=float)
        resp2 = np.full((n_stim, n_units), np.nan, dtype=float)

        for i, stim in enumerate(unique_stim):
            trials = np.where(stim_idx == stim)[0]
            if len(trials) < 2:
                continue

            trials = rng.permutation(trials)

            # ensure equal-sized halves
            n_use = (len(trials) // 2) * 2
            trials = trials[:n_use]
            half = n_use // 2

            trials1 = trials[:half]
            trials2 = trials[half:]

            resp1[i] = resp_mat[trials1].mean(axis=0)
            resp2[i] = resp_mat[trials2].mean(axis=0)

        for u in range(n_units):
            x = resp1[:, u]
            y = resp2[:, u]

            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 2:
                continue

            x_valid = x[valid]
            y_valid = y[valid]

            if np.std(x_valid) == 0 or np.std(y_valid) == 0:
                continue

            corr_all[split, u] = pearsonr(x_valid, y_valid)[0]

    noise_ceiling = np.nanmean(corr_all, axis=0)
    noise_ceiling_sb = (2 * noise_ceiling) / (1 + noise_ceiling)

    return noise_ceiling_sb


def load_data_from_good_unit(
    file_path: str,
    start_time: int = 60,
    end_time: int = 220,
    n_splits: int = 100,
    random_state: Optional[int] = None,
) -> dict:
    """
    Load response data and compute noise ceiling from one .mat file.

    Parameters
    ----------
    file_path : str
        Path to a single .mat file.
    start_time : int, default=60
        Start time relative to stimulus onset, in ms.
    end_time : int, default=220
        End time relative to stimulus onset, in ms.
    n_splits : int, default=100
        Number of split-half repetitions for noise ceiling.
    random_state : int or None, default=None
        Random seed for reproducibility.

    Returns
    -------
    dict
        A dictionary with keys:
        - "data": array of unit responses, shape (n_units, n_images)
        - "nc":   array of noise ceilings, shape (n_units,)
    """
    file_path = str(file_path)

    if not file_path.endswith(".mat"):
        raise ValueError(f"Expected a .mat file, got: {file_path}")

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    print(f"Loading: {file_path}")

    with h5py.File(file_path, "r") as f:
        n_units = len(f["GoodUnitStrc"]["response_matrix_img"])

        global_params_key = next(
            k for k in f.keys() if k.lower() == "global_params"
        )
        pre_onset_key = next(
            k for k in f[global_params_key].keys() if k.lower() == "pre_onset"
        )
        pre_onset = int(np.ravel(np.array(f[global_params_key][pre_onset_key]))[0])

        cur_start_time = pre_onset + start_time
        cur_end_time = pre_onset + end_time

        trial_valid_idx = np.array(f["meta_data"]["trial_valid_idx"]).reshape(-1).astype(np.int32)
        dataset_valid_idx = np.array(f["meta_data"]["dataset_valid_idx"]).reshape(-1).astype(bool)
        trial_idx = trial_valid_idx[dataset_valid_idx]

        response_list = []
        raster_list = []

        for unit_idx in range(n_units):
            resp_mat = f[f["GoodUnitStrc"]["response_matrix_img"][unit_idx][0]]
            raster = f[f["GoodUnitStrc"]["Raster"][unit_idx][0]]

            cur_resp = np.mean(resp_mat[cur_start_time:cur_end_time], axis=0)
            cur_raster = np.mean(raster[cur_start_time:cur_end_time], axis=0) * 1000

            response_list.append(cur_resp)
            raster_list.append(cur_raster)

        response_arr = np.array(response_list)   # shape: (n_units, n_images)
        raster_arr = np.array(raster_list)       # shape: (n_units, n_trials)

        noise_ceiling = compute_noise_ceiling(
            raster_arr.T,
            trial_idx,
            n_splits=n_splits,
            random_state=random_state,
        )

    return {"data": response_arr, "nc": noise_ceiling}


def main():
    input_path = r"/path/to/your/data_file.mat"
    result = load_data_from_good_unit(
        file_path=input_path,
        start_time=60,
        end_time=220,
        n_splits=100,
        random_state=42,
    )

    print("Response data shape:", result["data"].shape)
    print("Noise ceiling shape:", result["nc"].shape)
    print("Mean noise ceiling:", np.nanmean(result["nc"]))


if __name__ == "__main__":
    main()
