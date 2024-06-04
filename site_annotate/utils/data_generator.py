import pandas as pd
import numpy as np


def generate_test_data(n_rows=16, seed=43):
    """
    this is tmt test data
    for labelfree test data there simply is no "sample_xx" additional tmt measurement,
    the `Intensity` column is all that is used/needed
    """
    # Initialize random seed for reproducibility
    np.random.seed(seed)

    # Generate DataFrame with specified columns and random data
    data = {
        "Peptide": [
            str.join(
                "",
                np.random.choice(
                    list("ACDEFGHIKLMPRSTVWY"),
                    replace=True,
                    size=np.random.randint(7, 30),
                ),
            )
            for _ in np.arange(0, n_rows)
        ],
        "Peptide Length": np.random.randint(
            7, 13, size=n_rows
        ),  # Integer between 7 and 12
        "Charge": np.random.choice(
            [2, 3], size=n_rows, p=[0.2, 0.8]
        ),  # Mostly 3, some 2
        "Retention": np.random.normal(
            1549, 50, size=n_rows
        ),  # Normal distribution around 1549
        "Observed Mass": np.random.uniform(
            1248, 1793, size=n_rows
        ),  # Uniform distribution
        "Calibrated Observed Mass": np.random.uniform(1248, 1793, size=n_rows),
        "Observed M/Z": np.random.uniform(417, 743, size=n_rows),
        "Calibrated Observed M/Z": np.random.uniform(417, 743, size=n_rows),
        "Calculated Peptide Mass": np.random.uniform(1248, 1792, size=n_rows),
        "Calculated M/Z": np.random.uniform(417, 743, size=n_rows),
        "Delta Mass": np.random.uniform(
            -0.001, 1.005, size=n_rows
        ),  # Covering specified range
        "SpectralSim": np.random.uniform(0.375, 1.0, size=n_rows),
        "RTScore": np.random.uniform(0.004, 23, size=n_rows),
        "Expectation": np.random.uniform(0.004, 23, size=n_rows),
        "Hyperscore": 30 * np.random.beta(5, 2, size=n_rows),
        "Nextscore": 30 * np.random.beta(2, 2, size=n_rows), #np.random.uniform(0, 25, size=n_rows),
        "PeptideProphet Probability": np.random.uniform(0.85, 1, size=n_rows),
        "Number of Enzymatic Termini": np.random.choice([2], size=n_rows),  # Always 2
        "Number of Missed Cleavages": np.random.randint(
            0, 3, size=n_rows
        ),  # Integer between 0 and 2
        "Protein Start": np.random.randint(0, 890, size=n_rows),  # Random integer
        #"Protein End": np.random.randint(112, 899, size=n_rows),  # Random integer  not guaranteed to be greater than start
        "Intensity": np.random.uniform(
            2e6, 5.6e9, size=n_rows
        ),  # Wide range of intensities
        "STY:79.9663 Best Localization": np.random.beta(3, 1, size=n_rows),
        "Purity": 0.6 * np.random.beta(3, 1, size=n_rows) + .4,
        # Add sample columns
        **{
            f"sample-0{i+1}": np.random.uniform(0, 300000, size=n_rows)
            for i in range(16)
        },
    }

    # Create DataFrame
    df = pd.DataFrame(data)
    df["Protein End"] = df["Protein Start"] + df['Peptide'].apply(len)

    return df
