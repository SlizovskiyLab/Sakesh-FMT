import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from skbio.stats.composition import clr
from sklearn.impute import KNNImputer
from scipy.stats import chi2
from matplotlib.patches import Ellipse
from matplotlib.patches import Patch

# File paths
mge_matrix_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\Telcomb_MGE_analytical_matrix.xlsx - Sheet1.csv"
mge_annotations_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\MGE_total_classification.xlsx - Sheet1.csv"
# --- Use the new binned dataset ---
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_paired_binned.csv"

# Loading datasets
mge_matrix = pd.read_csv(mge_matrix_path)
mge_annotations = pd.read_csv(mge_annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# Filter for disease and clean data
fmt_dataset = fmt_dataset.iloc[:-5]
fmt_dataset = fmt_dataset.dropna(subset=['Patient', 'donor_pre_post_binned'])
fmt_dataset = fmt_dataset[fmt_dataset['Patient'].astype(str).str.strip() != '']
fmt_dataset = fmt_dataset[fmt_dataset['Disease_type'] == 'rCDI']

# Renaming column for merging
mge_annotations = mge_annotations.rename(columns={'IDs': 'gene_accession'})

# Merging MGE matrix with annotations
mge_merged = mge_matrix.merge(
    mge_annotations[['gene_accession', 'final_classification']],
    on='gene_accession',
    how='left'
)

# Aggregating by 'final_classification'
mobilome_features = mge_merged.drop(columns=['gene_accession']).groupby('final_classification').sum().T

# Resetting index and renaming
mobilome_features.reset_index(inplace=True)
mobilome_features.rename(columns={'index': 'ID'}, inplace=True)
fmt_dataset['run_accession'] = fmt_dataset['run_accession'].astype(str).str.strip()
mobilome_features['ID'] = mobilome_features['ID'].astype(str).str.strip()

# --- SIMPLIFIED MERGE ---
# Merge directly with the new 'donor_pre_post_binned' column
merged_mobilome_df = mobilome_features.merge(fmt_dataset[['run_accession', 'donor_pre_post_binned', 'Patient']], 
                                             left_on='ID', right_on='run_accession', how='left').drop(columns=['run_accession'])
merged_mobilome_df = merged_mobilome_df.dropna(subset=['donor_pre_post_binned'])

# Applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_mobilome_df.drop(columns=['ID', 'donor_pre_post_binned']))

# Applying CLR transformation
pseudocount = 1e-6
clr_input = imputed_data
clr_input[clr_input <= 0] = pseudocount
clr_transformed = clr(clr_input)
clr_transformed[np.isinf(clr_transformed)] = 0
clr_transformed = np.nan_to_num(clr_transformed)

# Computing Aitchison distances
aitchison_distances = squareform(pdist(clr_transformed, metric='euclidean'))

# Performing PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(aitchison_distances)
merged_mobilome_df['PC1'] = pca_result[:, 0]
merged_mobilome_df['PC2'] = pca_result[:, 1]

# Function to compute and draw 95% confidence ellipses
def confidence_ellipse(x, y, ax, color, n_std=1.96):
    if len(x) < 2:
        return
    mean_x, mean_y = np.mean(x), np.mean(y)
    cov = np.cov(x, y)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    chi2_val = np.sqrt(chi2.ppf(0.95, df=2))
    width, height = 2 * chi2_val * np.sqrt(eigvals)
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    ellipse = Ellipse(xy=(mean_x, mean_y), width=width, height=height, angle=angle, edgecolor=color, facecolor=color, alpha=0.2)
    ax.add_patch(ellipse)

# Creating scatter plot
plt.figure(figsize=(10, 6))
unique_groups = merged_mobilome_df['donor_pre_post_binned'].unique()
group_colors = {
    'Donor': '#34301f',
    'PreFMT': '#726732',
    'PostFMT (1-30d)': '#b3cde0',
    'PostFMT (31-60d)': '#6497b1',
    'PostFMT (60d+)': '#005b96',
}

ax = sns.scatterplot(
    x='PC1', 
    y='PC2', 
    hue='donor_pre_post_binned', 
    data=merged_mobilome_df, 
    palette=group_colors, 
    alpha=0.7, 
    edgecolor='k',
    legend=False
)

# Computing confidence ellipses
for group in unique_groups:
    subset = merged_mobilome_df[merged_mobilome_df['donor_pre_post_binned'] == group]
    confidence_ellipse(subset['PC1'], subset['PC2'], ax, group_colors.get(group, 'gray'))

plt.xlim(merged_mobilome_df['PC1'].min() - 150, merged_mobilome_df['PC1'].max() + 70)
plt.ylim(merged_mobilome_df['PC2'].min() - 100, merged_mobilome_df['PC2'].max() + 70)

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.grid(False)

# Create and add the custom legend
legend_handles = [Patch(facecolor=color, edgecolor='k', label=label) for label, color in group_colors.items()]
lgd = plt.legend(
    handles=legend_handles, 
    title='Sample Type', 
    bbox_to_anchor=(1.7, 0.5), 
    loc='right', 
    markerscale=2, 
    fontsize=20
)

# Adjust layout to make room for the legend
plt.subplots_adjust(right=0.7)

# Save plot and show
plt.savefig(
    "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/pca_rcdi.svg", 
    format='svg', 
    dpi=600, 
    bbox_extra_artists=(lgd,),
    bbox_inches='tight', 
    transparent=True
)
plt.savefig(
    "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/pca_rcdi.png", 
    format='png', 
    dpi=600, 
    bbox_extra_artists=(lgd,),
    bbox_inches='tight', 
    transparent=True
)
plt.show()

# Save metadata with the correct column
metadata_df = merged_mobilome_df[['ID', 'donor_pre_post_binned', 'Patient']]
metadata_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/metadata_rcdi.csv", index=False)

# Save Aitchison distance matrix
aitchison_df = pd.DataFrame(
    aitchison_distances,
    index=merged_mobilome_df['ID'],
    columns=merged_mobilome_df['ID']
)
aitchison_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/aitchison_rcdi.csv")
