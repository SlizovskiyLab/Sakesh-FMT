import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def process_disease_data(df_full, aitchison_file_path, disease_type_name):
    """
    Filters data for a specific disease, calculates Aitchison distances, and bins them.

    Args:
        df_full (pd.DataFrame): The complete metadata DataFrame.
        aitchison_file_path (str): Path to the Aitchison distance matrix for the disease.
        disease_type_name (str): The name of the disease to filter by (e.g., 'rCDI').

    Returns:
        pd.DataFrame: A DataFrame with binned distances for plotting.
    """
    print(f"\n--- Processing: {disease_type_name} ---")
    try:
        df_aitchison = pd.read_csv(aitchison_file_path, index_col=0)
        print(f"Successfully loaded Aitchison matrix: {aitchison_file_path}")
    except FileNotFoundError:
        print(f"Error: Aitchison file not found at '{aitchison_file_path}'. Skipping this disease.")
        return pd.DataFrame() # Return empty DataFrame if file is missing

    # --- Step 1: Filter and clean the metadata DataFrame ---
    df_filtered = df_full[df_full['Disease_type'] == disease_type_name].copy()
    print(f"Filtered for '{disease_type_name}'. Shape of data is now: {df_filtered.shape}")

    if len(df_filtered) <= 5:
        print("Warning: DataFrame has 5 or fewer rows.")

    # Standardize 'donor_pre_post' column
    replacements = {'Pre-FMT': 'PreFMT', 'Post-FMT': 'PostFMT'}
    df_filtered['donor_pre_post'] = df_filtered['donor_pre_post'].replace(replacements)

    # Clean 'timepoint' column
    df_filtered['timepoint'] = pd.to_numeric(df_filtered['timepoint'], errors='coerce')
    df_filtered.dropna(subset=['timepoint'], inplace=True)
    df_filtered['timepoint'] = df_filtered['timepoint'].astype(int)

    # --- Step 2: Map donors to patients ---
    donors_df = df_filtered[df_filtered['donor_pre_post'] == 'Donor']
    donor_map = donors_df.set_index('Patient')['run_accession'].to_dict()
    print(f"Created a map for {len(donor_map)} donors.")

    # --- Step 3: Extract distances and bin them ---
    bins = {'0': [], '1-30': [], '31-60': [], '>60': []}
    patient_samples_df = df_filtered[df_filtered['donor_pre_post'].isin(['PreFMT', 'PostFMT'])]

    for _, sample_row in patient_samples_df.iterrows():
        patient_id = sample_row['Patient']
        patient_run_accession = sample_row['run_accession']
        timepoint = sample_row['timepoint']
        donor_run_accession = donor_map.get(patient_id)

        if donor_run_accession and patient_run_accession in df_aitchison.index and donor_run_accession in df_aitchison.index:
            distance = df_aitchison.loc[patient_run_accession, donor_run_accession]
            if timepoint == 0:
                bins['0'].append(distance)
            elif 1 <= timepoint <= 30:
                bins['1-30'].append(distance)
            elif 31 <= timepoint <= 60:
                bins['31-60'].append(distance)
            elif timepoint > 60:
                bins['>60'].append(distance)

    # --- Step 4: Convert binned data to a long-form DataFrame for plotting ---
    plot_data = []
    for bin_name, distances in bins.items():
        for distance in distances:
            plot_data.append({
                'Timepoint Bin': bin_name,
                'Aitchison Distance': distance,
                'Disease Type': disease_type_name
            })
    return pd.DataFrame(plot_data)


# --- Main Script ---

# --- Step 0: Define File Paths ---
metadata_file = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_paired.csv"
# A dictionary mapping disease types to their respective Aitchison distance files
disease_files = {
    'rCDI': "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/aitchison_rcdi.csv",
    'Melanoma': "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/aitchison_melanoma.csv",
    'MDRB': "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/aitchison_mdrb.csv"
}

try:
    df_full_metadata = pd.read_csv(metadata_file, low_memory=False)
    print("Successfully loaded the main metadata file.")
except FileNotFoundError:
    print(f"Error: Main metadata file not found at '{metadata_file}'.")
    exit()

# --- Process each disease and collect the results ---
all_plot_data = []
for disease, aitchison_file in disease_files.items():
    disease_df = process_disease_data(df_full_metadata, aitchison_file, disease)
    if not disease_df.empty:
        all_plot_data.append(disease_df)

# Combine all data into a single DataFrame
if not all_plot_data:
    print("\nNo data was processed successfully. Exiting.")
    exit()

final_plot_df = pd.concat(all_plot_data, ignore_index=True)

# --- Step 5: Create Stacked Violin Plots ---
print("\nGenerating stacked violin plots...")

# Get the list of diseases that were successfully processed
processed_diseases = final_plot_df['Disease Type'].unique()
num_diseases = len(processed_diseases)

fig, axes = plt.subplots(num_diseases, 1, figsize=(12, 6 * num_diseases), sharex=True, sharey=True)
if num_diseases == 1: # If only one plot, axes is not a list
    axes = [axes]

bin_order = ['0', '1-30', '31-60', '>60']

for i, disease_name in enumerate(processed_diseases):
    ax = axes[i]
    disease_data = final_plot_df[final_plot_df['Disease Type'] == disease_name]
    
    sns.violinplot(x='Timepoint Bin', y='Aitchison Distance', data=disease_data, order=bin_order, ax=ax, palette='viridis')
    
    ax.set_title(f'Distribution of Aitchison Distances for {disease_name}', fontsize=16)
    ax.set_ylabel('Aitchison Distance to Donor', fontsize=12)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

# Set common X-axis label
axes[-1].set_xlabel('Timepoint Bin (Days Post-FMT)', fontsize=12)
plt.xticks(rotation=45)

plt.tight_layout()
plt.show()

