import pandas as pd
import numpy as np



def generate_test_data(seed=43, n_samples=16):
    """
    this is tmt test data
    for labelfree test data there simply is no "sample_xx" additional tmt measurement,
    the `Intensity` column is all that is used/needed
    """
    # Initialize random seed for reproducibility
    np.random.seed(seed)

    # Generate DataFrame with specified columns and random data
    data = {
        'Peptide Length': np.random.randint(7, 13, size=n_samples),  # Integer between 7 and 12
        'Charge': np.random.choice([2, 3], size=n_samples, p=[0.2, 0.8]),  # Mostly 3, some 2
        'Retention': np.random.normal(1549, 50, size=n_samples),  # Normal distribution around 1549
        'Observed Mass': np.random.uniform(1248, 1793, size=n_samples),  # Uniform distribution
        'Calibrated Observed Mass': np.random.uniform(1248, 1793, size=n_samples),
        'Observed M/Z': np.random.uniform(417, 743, size=n_samples),
        'Calibrated Observed M/Z': np.random.uniform(417, 743, size=n_samples),
        'Calculated Peptide Mass': np.random.uniform(1248, 1792, size=n_samples),
        'Calculated M/Z': np.random.uniform(417, 743, size=n_samples),
        'Delta Mass': np.random.uniform(-0.001, 1.005, size=n_samples),  # Covering specified range
        'SpectralSim': np.random.uniform(0.375, 0.867, size=n_samples),
        'RTScore': np.random.uniform(0.004, 23, size=n_samples),
        'Expectation': np.random.uniform(0.004, 23, size=n_samples),
        'Hyperscore': np.random.uniform(5, 25, size=n_samples),
        'Nextscore': np.random.uniform(0, 25, size=n_samples),
        'PeptideProphet Probability': np.random.uniform(0.85, 1, size=n_samples),
        'Number of Enzymatic Termini': np.random.choice([2], size=n_samples),  # Always 2
        'Number of Missed Cleavages': np.random.randint(0, 3, size=n_samples),  # Integer between 0 and 2
        'Protein Start': np.random.randint(106, 890, size=n_samples),  # Random integer
        'Protein End': np.random.randint(112, 899, size=n_samples),  # Random integer
        'Intensity': np.random.uniform(2e6, 5.6e9, size=n_samples),  # Wide range of intensities
        'STY:79.9663 Best Localization': np.random.uniform(0.41, 1, size=n_samples),
        'Purity': np.random.uniform(0.55, 1, size=n_samples),
        # Add sample columns
        **{f'sample-0{i+1}': np.random.uniform(0, 300000, size=n_samples) for i in range(16)}
    }

    # Create DataFrame
    df = pd.DataFrame(data)

    return df

