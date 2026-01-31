import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import norm
import argparse

# Function to generate a random correlation matrix
def generate_random_correlation_matrix(n_vars):
    """
    Generates a random, symmetric, positive semi-definite correlation matrix
    of size n_vars x n_vars.

    Args:
        n_vars (int): Number of variables (dimensions of the correlation matrix).

    Returns:
        np.ndarray: A random positive semi-definite correlation matrix.
    """
    random_matrix = np.random.uniform(-1, 1, (n_vars, n_vars))
    correlation_matrix = (random_matrix + random_matrix.T) / 2
    np.fill_diagonal(correlation_matrix, 1.0)

    eigvals, eigvecs = np.linalg.eigh(correlation_matrix)
    eigvals = np.maximum(eigvals, 0)
    correlation_matrix = eigvecs @ np.diag(eigvals) @ eigvecs.T

    return correlation_matrix

# Function to handle z-score calculation with zero standard deviation handling
def safe_zscore(x, axis=0):
    """
    Computes the z-score while handling zero standard deviation cases.

    Args:
        x (np.ndarray or pd.DataFrame): The input data.
        axis (int): Axis along which to compute the z-score. Default is 0 (columns).

    Returns:
        np.ndarray: The z-scored data, with zero standard deviation cases set to zero.
    """
    mean = np.mean(x, axis=axis)
    std = np.std(x, axis=axis)
    std = np.where(std == 0, 1, std)  # Replace zero standard deviations with 1
    return (x - mean) / std


def replace_neg9_with_sampled(row):
    if (row == -9).sum() > 0:
        # print("Filling in missing.")
        valid_values = row[row != -9]  # Extract valid values (continuous between 0 and 2)
        
        if valid_values.empty:  # If all values are -9, return row as-is
            return row
        
        # Sample missing values directly from the observed non -9 values
        replace_values = np.random.choice(valid_values, size=(row == -9).sum(), replace=True)
        
        row[row == -9] = replace_values  # Replace -9 with sampled values
    
    return row

# Load input data and preprocess
def load_data(fname):
    """
    Loads and preprocesses the input data from a tab-separated file.

    Args:
        fname (str): Path to the input file.

    Returns:
        pd.DataFrame: Preprocessed DataFrame.
        np.ndarray: Transposed genotype matrix G.
    """
    df = pd.read_csv(fname, sep='\t')
    df.set_index('SNP', inplace=True)
    G = df.iloc[:, 5:] #setting SNP as index removed it from the data - originally we want the columns after the first six
    # assert G.shape[1]==32017
    # replace missing
    G = G.apply(replace_neg9_with_sampled, axis=1)
    G = G.T
    return df, G

# Main function to generate and save the data
def generate_and_save_data(fname, n_vars=5, h2=0.5, theta=0.1):
    """
    Generates the random variables and saves them to a CSV file.

    Args:
        fname (str): Path to the input file.
        n_vars (int): Number of variables for the random correlation matrix.
        h2 (float): Heritability value for sigma_a.
        theta (float): Threshold for generating ybool.
    """
    print("Loading genotypes.")
    DF, G = load_data(fname)
    
    print("Generating correlations for covars.")
    random_correlation_matrix = generate_random_correlation_matrix(n_vars)

    print("Covariates Correlation Matrix:")
    print(np.round(random_correlation_matrix, 2))

    eigvals, _ = np.linalg.eigh(random_correlation_matrix)
    mean = np.zeros(n_vars)
    COV = np.random.multivariate_normal(mean, random_correlation_matrix, size=G.shape[0])

    sigma_a = np.sqrt(h2)
    sigma_cov = np.sqrt(0.1)
    sigma_e = np.sqrt(1 - h2 - 0.1)

    x_ = np.random.binomial(1, 0.01, size=G.shape[1])
    
    ## Save x_ to a CSV file
    x_values_df = pd.DataFrame(x_, columns=["x_values"])
    x_values_df['variant IDs'] = DF.index.values
    output_filename = fname.split('.')[0]
    output_file_path = f"{output_filename}_known_g_effects.csv"
    x_values_df.to_csv(output_file_path, index=False)
    print(f"x_ values saved to {output_file_path}")
    ######
    
    print("Running GLM model.")
    
    z_ = np.random.binomial(1, 1, size=COV.shape[1])

    Z = safe_zscore(G, axis=0)
    
    assert np.isnan(Z.values).sum() == 0
        
    E = np.random.normal(0, 1, G.shape[0])[:, None]

    y = sigma_a * safe_zscore(Z @ x_[:, None]) + \
        sigma_cov * (COV @ z_[:, None]) + \
        sigma_e * E

    threshold = norm.ppf(1 - theta)
    ybool = (y > threshold).astype(int) + 1

    count_2 = np.sum(ybool == 2)
    frequency_2 = count_2 / G.shape[0]
    print(f"\nybool case frequency: {frequency_2}")

    IID = np.array([i.split("0_")[1] for i in G.index.values])

    output_df = pd.DataFrame(
        np.hstack([np.zeros([len(y), 1]).astype(int), IID[:, None], y, ybool, COV]),
        columns=["#FID"] +
                ["IID"] +
                ["y"] +
                ['ybool'] +
                [f'COV_{i+1}' for i in range(COV.shape[1])]
    )

    print("\nFirst few rows of the output data:")
    print(output_df.head())

    print("\nSummary statistics of the output data:")
    print(output_df.describe())

    output_filename = fname.split('.')[0]
    output_file_path = f"{output_filename}_phenocov.csv"

    output_df.to_csv(output_file_path, index=False, sep="\t")
    print(f"Data saved to {output_file_path}")

# Command-line interface using argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate random variables and save data.")
    parser.add_argument("fname", type=str, help="Path to the input file.")
    parser.add_argument("--n_vars", type=int, default=5, help="Number of variables for the random correlation matrix.")
    parser.add_argument("--h2", type=float, default=0.5, help="Heritability value for sigma_a.")
    parser.add_argument("--theta", type=float, default=0.1, help="Threshold for generating ybool.")

    args = parser.parse_args()
    generate_and_save_data(args.fname, args.n_vars, args.h2, args.theta)
