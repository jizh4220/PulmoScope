import math
import numpy as np
import scanpy as sc

def bitpacking_compression(adata):
    # Convert the matrix to a dense array
    data = adata.X.toarray().flatten()
    
    # Find the minimum and maximum values in the data
    min_val = min(data)
    max_val = max(data)

    # Compute the range and number of bits needed for each value
    data_range = max_val - min_val
    num_bits = math.ceil(math.log2(data_range + 1))

    # Initialize the compressed data as an empty list
    compressed_data = []

    # Encode the input data
    for value in data:
        # Subtract the minimum value and encode using the computed number of bits
        encoded_value = value - min_val
        compressed_data.append(encoded_value)

    # Return the compressed data, minimum value, number of bits and shape of original matrix
    return compressed_data, min_val, num_bits, adata.X.shape

# Example usage
adata = sc.read_10x_mtx(
    './filtered_feature_bc_matrix/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)                                # write a cache file for faster subsequent reading

compressed_data, min_val, num_bits, shape = bitpacking_compression(adata)
print(f"Compressed data: {compressed_data}")
print(f"Minimum value: {min_val}")
print(f"Number of bits: {num_bits}")
print(f"Original shape: {shape}")