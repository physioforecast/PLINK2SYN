import pandas as pd
import numpy as np
import random
import string
import os
import argparse
import tqdm
import gc

def get_var_protect(genotype_data):
    # Identify rows without missing values (-9)
    rows_without_missing = 1 - np.any(np.isnan(genotype_data.iloc[:, 6:].values.astype(float)), axis=1)

    # Compute differences to find transitions (start and end of blocks)
    diff = np.diff(rows_without_missing.astype(int))

    # Find start and end indices of contiguous blocks
    start_indices = np.where(diff == 1)[0] + 1  # Add 1 because diff shifts by 1
    end_indices = np.where(diff == -1)[0] + 1  # End indices are where -1 occurs

    # Handle edge cases: if contiguous block starts or ends at the boundary
    if rows_without_missing[0]:
        start_indices = np.insert(start_indices, 0, 0)
    if rows_without_missing[-1]:
        end_indices = np.append(end_indices, len(rows_without_missing))

    # Pair start and end indices to form intervals
    intervals = list(zip(start_indices, end_indices))[1:]

    # Create a protect mask for all contiguous non-missing
    protect_mask = np.zeros(genotype_data.shape[0])
    for ix0, ix1 in intervals:
        if ix1 - ix0 > 1:
            protect_mask[ix0:ix1] = 1

    return protect_mask

def augment_traw(fname, fname2, nsamples=32017, mean=0, std_dev=0.1, miss=False, p=0.1, outputflag=1):
    """
    Augments genetic data to match the specified number of samples. Generates new traw and psam files 
    for use with PLINK. This function adds Gaussian noise to the genotypic data, applies bounds between 0 and 2,
    and optionally introduces missing values.

    Parameters:
    - fname (str): Path to the input genotype file (tab-delimited).
    - fname2 (str): Path to the input sample file (psam, tab-delimited).
    - nsamples (int): The number of samples to generate (default is 32017).
    - mean (float): Mean for the Gaussian noise distribution (default is 0).
    - std_dev (float): Standard deviation for the Gaussian noise distribution (default is 0.1).
    - miss (bool): Whether to randomly introduce missing values (default is False).
    - p (float): Probability of setting a genotype value to missing (if `miss=True`).

    Returns:
    - modified_sample_data (DataFrame): The modified sample data (psam).
    - modified_genotype_data (DataFrame): The modified genotype data (traw).
    """
    

    
    # Step 1: Read and filter the sample data (psam file)
    sample_data = pd.read_csv(fname2, sep="\t", low_memory=False)
    sample_data = sample_data.query("(PAT=='0') & (MAT=='0')")  # Filter for samples with unknown parents
    sample_data = sample_data.sample(n=nsamples, replace=True, random_state=42)  # Randomly sample nsamples rows
    sample_data.reset_index(drop=True, inplace=True)  # Reset index for consistency
    sample_data.index = range(len(sample_data))  # Reindex to start from 0
    modified_sample_data = sample_data.copy()  # Create a copy for modifying sample IDs

    # Step 2: Read and prepare the genotype data (traw file)
    with open(fname, 'r') as f:
        header = f.readline().strip().split('\t')

    # Define dtype for columns (set columns 6 and onward as float)
    dtype_dict = {header[i]: float for i in range(6, len(header))}

    # Load the genotype data
    genotype_data = pd.read_csv(fname, sep='\t', header=0, low_memory=False, dtype=dtype_dict)

    # Step 3: Generate new sample IDs based on the psam file
    sample_ids = ["0_" + i for i in sample_data['#IID']]  # Generate new sample IDs based on the psam file

    # Step 4: Get the protect mask for variants
    if miss:
        protect_mask = get_var_protect(genotype_data)
        # print(np.sum(protect_mask)/len(protect_mask))

    # Step 5: Prepare new traw and psam file names
    base_name = os.path.splitext(fname)[0]  # Get the base name (without extension) from the original file
    new_fname_traw = f"{base_name}_nonintdose_{outputflag}.traw"  # New traw file name
    new_fname_psam = f"{base_name}_nonintdose_{outputflag}.psam"  # New psam file name
    
    
    # Step 6: Modify genotype data for each sample
    g_columns = []
    g_columns_names = []
    for ix, sid in tqdm.tqdm(enumerate(sample_ids), total=len(sample_ids), desc="Processing samples"):
        # print(f"SID:{sid}")
        # Generate a new random sample ID
        new_sid = "0_" + ''.join(random.choices(string.ascii_letters + string.digits, k=6))  # Random ID
        modified_sample_data.loc[ix, '#IID'] = new_sid.split("_")[1]  # Update sample ID
        
        # Modify genotype data for the sample
        g_ = genotype_data[sid].astype(float)#np.float16)  # Get the genotypic data for the sample
        
        file_path = "np_isnan_debug.txt"
        with open(file_path, "w") as file:
            file.write("\n".join(list(np.isnan(g_).astype(str))))
       
        # x = bool(np.sum(np.isnan(g_))>0)
        # print(f"g_ NAN COUNT:{x}")
        
        gaussian_noise = np.random.normal(mean, std_dev, len(g_))  # Generate Gaussian noise
        noisy_genotype_data = g_ + gaussian_noise  # Add noise
        noisy_genotype_data = np.clip(noisy_genotype_data, a_min=0, a_max=2)  # Clip values to the range [0, 2]
        noisy_genotype_data[np.isnan(noisy_genotype_data)] = -9  # Replace NaN with -9

        # Introduce missingness if `miss` is True
        if miss:
            # print(f"MISS:{miss}")
            # print(f"MISS P:{p}")
            add_miss_mask = (np.random.rand(len(noisy_genotype_data)) > p).astype(int)
            # mask =  #add_miss_mask * 
            mask = ((1 - protect_mask)*add_miss_mask).astype(bool)  
            noisy_genotype_data[mask] = -9
            
            
        
        
        g_columns.append(noisy_genotype_data.astype(np.float32))
        g_columns_names.append(new_sid)
    
    G_ = pd.DataFrame(g_columns).T  # Create a DataFrame from the list and transpose it
    del(g_columns)
    gc.collect()
    G_.columns = g_columns_names
    G_ = pd.concat([genotype_data.iloc[:, :6], G_], axis=1)

        # noisy_genotype_data = noisy_genotype_data.reshape(-1, 1)  # Convert to 2D array

        # # Save modified genotype data for this sample
        # df = pd.DataFrame(noisy_genotype_data.values, columns=[new_sid])
        
        # df.to_csv("tmp.geno", index=False, sep="\t")

        # # Append the modified genotype data to the traw file
        # os.system(f"paste {new_fname_traw} tmp.geno > temp_output.txt && mv temp_output.txt {new_fname_traw}")

    # Save the modified sample data (psam) to file
    G_.to_csv(new_fname_traw, sep='\t', index=False)  # Save the modified genotype data
    modified_sample_data.to_csv(new_fname_psam, sep='\t', index=False)

    return 

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="""Data augmentation for genotype data. Augments genetic data to match the specified number of samples. Generates new traw and psam files 
    for use with PLINK. This function adds Gaussian noise to the genotypic data, applies bounds between 0 and 2,
    and optionally introduces missing values.""")
    
    # Required positional arguments
    parser.add_argument('--fname', type=str, help="Path to the genotype data file (traw format)")
    parser.add_argument('--fname2', type=str, help="Path to the sample data file (psam format)")
    
    # Optional arguments with default values
    parser.add_argument('--nsamples', type=int, default=32017, help="Number of samples to generate (default: 32017)")
    parser.add_argument('--mean', type=float, default=0, help="Mean for noise (default: 0)")
    parser.add_argument('--std_dev', type=float, default=0.1, help="Standard deviation for noise (default: 0.1)")
    parser.add_argument('--miss', action='store_true', help="Whether to introduce missing values (default: False)")
    parser.add_argument('--p', type=float, default=0.1, help="Probability of missing values on top of original missingness (default: 0.1)")
    parser.add_argument('--outputflag', type=str, default=1, help="Flag for output filename.")

    # Parse the arguments
    args = parser.parse_args()

    # Call the augment_traw function with the parsed arguments
    augment_traw(args.fname, args.fname2, args.nsamples, args.mean, args.std_dev, args.miss, args.p,args.outputflag)

if __name__ == "__main__":
    main()
