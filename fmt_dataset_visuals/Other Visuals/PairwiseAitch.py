import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from skbio.stats.composition import clr
from sklearn.impute import KNNImputer
from scipy.stats import sem
from scipy.stats import chi2
from matplotlib.patches import Ellipse

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
        return pd.DataFrame()

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

fig, axes = plt.subplots(num_diseases, 1, figsize=(12, 6 * num_diseases), sharex=True)
fig.suptitle('Resistome-Donor Distance by Disease', fontsize=18, y=0.95)
if num_diseases == 1:
    axes = [axes]

bin_order = ['0', '1-30', '31-60', '>60']
palette = sns.color_palette("viridis", n_colors=len(bin_order))

# Calculate medians for plotting points
medians_df = final_plot_df.groupby(['Disease Type', 'Timepoint Bin'])['Aitchison Distance'].median().reset_index()

for i, disease_name in enumerate(processed_diseases):
    ax = axes[i]
    disease_data = final_plot_df[final_plot_df['Disease Type'] == disease_name]
    
    # --- Aesthetic Changes ---
    sns.violinplot(x='Timepoint Bin', y='Aitchison Distance', hue='Timepoint Bin', data=disease_data, order=bin_order, ax=ax, 
                   palette=palette, inner=None, linewidth=2, legend=False)
    
    sns.boxplot(x='Timepoint Bin', y='Aitchison Distance', data=disease_data, order=bin_order, ax=ax, 
                color="white", width=0.2, boxprops={'zorder': 2},
                medianprops={'visible': False})
                
    sns.stripplot(x='Timepoint Bin', y='Aitchison Distance', data=disease_data, order=bin_order, ax=ax, 
                  color="black", alpha=0.3, jitter=True)
    
    # Plot a horizontal line for the median
    current_medians = medians_df[medians_df['Disease Type'] == disease_name]
    sns.scatterplot(x='Timepoint Bin', y='Aitchison Distance', data=current_medians, ax=ax,
                    color='black', s=600, marker='_', zorder=3, legend=False)

    for violin in ax.collections[:len(bin_order)]:
        for path in violin.get_paths():
            vertices = path.vertices
            vertices[:, 1] = np.maximum(0, vertices[:, 1])
    
    ax.set_title(f'Distribution of Aitchison Distances for {disease_name}', fontsize=11)
    
    if i == num_diseases // 2:
        ax.set_ylabel('Aitchison Distance to Donor', fontsize=12)
    else:
        ax.set_ylabel('')
    
    ax.grid(False)
    ax.set_ylim(bottom=0) 

# Set common X-axis label
axes[-1].set_xlabel('Timepoint Bin (Days Post-FMT)', fontsize=12)
plt.xticks(rotation=45)

# Adjust vertical spacing between plots and remove tight_layout
# plt.subplots_adjust(hspace=0.4)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/resistome_pairwise_aichison.png", format='png', dpi=600, bbox_inches='tight', transparent=True)

plt.show()

# --- START: Code for Mobilome Sister Figure ---

print("\n\n" + "="*50)
print("Generating Sister Figure for Mobilome Data")
print("="*50)


# --- Step 1: Define New File Paths for Mobilome ---
mobilome_disease_files = {
    'rCDI': "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/aitchison_rcdi.csv",
    'Melanoma': "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/aitchison_melanoma.csv",
    'MDRB': "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/aitchison_mdrb.csv"
}

# --- Step 2: Process each disease for Mobilome data ---
all_mobilome_plot_data = []
for disease, aitchison_file in mobilome_disease_files.items():
    disease_df = process_disease_data(df_full_metadata, aitchison_file, disease)
    if not disease_df.empty:
        all_mobilome_plot_data.append(disease_df)

# Combine all mobilome data into a single DataFrame
if not all_mobilome_plot_data:
    print("\nNo Mobilome data was processed successfully. Exiting.")
    exit()

final_mobilome_plot_df = pd.concat(all_mobilome_plot_data, ignore_index=True)

# --- Step 3: Create Stacked Violin Plots for Mobilome ---
print("\nGenerating stacked violin plots for Mobilome...")

# Get the list of diseases that were successfully processed
processed_mobilome_diseases = final_mobilome_plot_df['Disease Type'].unique()
num_mobilome_diseases = len(processed_mobilome_diseases)

fig_mobilome, axes_mobilome = plt.subplots(num_mobilome_diseases, 1, figsize=(12, 6 * num_mobilome_diseases), sharex=True)
fig_mobilome.suptitle('Mobilome-Donor Distance by Disease', fontsize=18, y=0.95)

if num_mobilome_diseases == 1:
    axes_mobilome = [axes_mobilome]

# Re-calculate medians for the new dataset
medians_mobilome_df = final_mobilome_plot_df.groupby(['Disease Type', 'Timepoint Bin'])['Aitchison Distance'].median().reset_index()

for i, disease_name in enumerate(processed_mobilome_diseases):
    ax = axes_mobilome[i]
    disease_data = final_mobilome_plot_df[final_mobilome_plot_df['Disease Type'] == disease_name]
    
    # Apply aesthetics
    sns.violinplot(x='Timepoint Bin', y='Aitchison Distance', hue='Timepoint Bin', data=disease_data, order=bin_order, ax=ax, 
                   palette=palette, inner=None, linewidth=2, legend=False)
    
    sns.boxplot(x='Timepoint Bin', y='Aitchison Distance', data=disease_data, order=bin_order, ax=ax, 
                color="white", width=0.2, boxprops={'zorder': 2},
                medianprops={'visible': False})
                
    sns.stripplot(x='Timepoint Bin', y='Aitchison Distance', data=disease_data, order=bin_order, ax=ax, 
                  color="black", alpha=0.3, jitter=True)
    
    # Plot a horizontal line for the median
    current_medians = medians_mobilome_df[medians_mobilome_df['Disease Type'] == disease_name]
    sns.scatterplot(x='Timepoint Bin', y='Aitchison Distance', data=current_medians, ax=ax,
                    color='black', s=600, marker='_', zorder=3, legend=False)

    for violin in ax.collections[:len(bin_order)]:
        for path in violin.get_paths():
            vertices = path.vertices
            vertices[:, 1] = np.maximum(0, vertices[:, 1])
    
    ax.set_title(f'Distribution of Aitchison Distances for {disease_name} (Mobilome)', fontsize=11)
    
    if i == num_mobilome_diseases // 2:
        ax.set_ylabel('Aitchison Distance to Donor', fontsize=12)
    else:
        ax.set_ylabel('')
    
    ax.grid(False)
    ax.set_ylim(bottom=0) 

# Set common X-axis label
axes_mobilome[-1].set_xlabel('Timepoint Bin (Days Post-FMT)', fontsize=12)
plt.xticks(rotation=45)

# Adjust vertical spacing between plots
# plt.subplots_adjust(hspace=0.4)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/mobilome_pairwise_aichison.png", format='png', dpi=600, bbox_inches='tight', transparent=True)
plt.show()

# --- END: Code for Mobilome Sister Figure ---