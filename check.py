import os
import subprocess
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_row(row, enzyme_tuple, index):
    """
    Function to process a single row in parallel.
    """
    try:
        combined_smiles = f"{row.predicted_substrate_smiles}>>{row.predicted_product_smiles}"
        command = f'java -jar ReactionDecoder.jar -Q SMI -q "{combined_smiles}" -g -c -j ANNOTATE -f XML -p "run{index}"'

        # Run the command, suppress stderr
        result = subprocess.run(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=30)

        if result.returncode != 0:
            return {'index': index, 'error': f"Command failed for row {index}: suppressed stderr output."}

        xml_path = f'run{index}_ECBLAST_smiles_ANNONATE.xml'
        if not os.path.exists(xml_path):
            return {'index': index, 'error': f"XML file not found for row {index}"}

        # Extract SMILES and calculate likelihoods
        unique_matches, substrate_smiles, product_smiles = reaction_likelihoods.extract_and_analyze_SMILES(xml_path)
        substrate_som_labels, product_som_labels, max_substrate_likelihood, max_product_likelihood = reaction_likelihoods.compare_molecules(
            substrate_smiles, product_smiles, enzyme_tuple, unique_matches
        )

        # Return the calculated likelihoods
        return {'index': index, 'substrate_likelihood': max_substrate_likelihood, 'product_likelihood': max_product_likelihood}

    except Exception as e:
        return {'index': index, 'error': str(e)}

# Parallelized main script
results = []
max_workers = 8  # Adjust this based on the number of available threads/cores
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = []
    for row in tqdm(reactions_df.itertuples(), total=reactions_df.shape[0], desc='Calculating reaction likelihoods'):
        enzyme_tuple = tuple(map(int, row.enzyme_ec.split('.')[:2]))
        futures.append(executor.submit(process_row, row, enzyme_tuple, row.Index))

    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing Results"):
        result = future.result()
        if 'error' in result:
            print(f"Error for row {result['index']}: {result['error']}")
        else:
            results.append(result)

# Update DataFrame
for result in results:
    reactions_df.loc[result['index'], 'substrate_likelihood'] = result.get('substrate_likelihood')
    reactions_df.loc[result['index'], 'product_likelihood'] = result.get('product_likelihood')
