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

# File paths
mge_matrix_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\Telcomb_MGE_analytical_matrix.xlsx - Sheet1.csv"
mge_annotations_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\MGE_total_classification.xlsx - Sheet1.csv"
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_paired.csv"

# Loading datasets
mge_matrix = pd.read_csv(mge_matrix_path)
mge_annotations = pd.read_csv(mge_annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

fmt_dataset = fmt_dataset[fmt_dataset['Disease_type'] == 'MDRB']

fmt_dataset = fmt_dataset.iloc[:-5]

# Standardize values in donor_pre_post
fmt_dataset['donor_pre_post'] = fmt_dataset['donor_pre_post'].replace({
    'Post-FMT': 'PostFMT',
    'Pre-FMT': 'PreFMT',
    'Pre-Abx/FMT': 'PreFMT'
})
fmt_dataset = fmt_dataset[fmt_dataset['donor_pre_post'].isin(['PostFMT', 'PreFMT', 'Donor'])]

# Remove rows where 'Patient' is missing or blank
fmt_dataset = fmt_dataset.dropna(subset=['Patient'])
fmt_dataset = fmt_dataset[fmt_dataset['Patient'].astype(str).str.strip() != '']


# Renaming column for merging
mge_annotations = mge_annotations.rename(columns={'IDs': 'gene_accession'})

# Merging MGE matrix with annotations, keeping only the 'final_classification' column
mge_merged = mge_matrix.merge(
    mge_annotations[['gene_accession', 'final_classification']],
    on='gene_accession',
    how='left'
)

# Aggregating by 'final_classification' (summation across samples)
mobilome_features = mge_merged.drop(columns=['gene_accession']).groupby('final_classification').sum().T

# Resetting index and renaming for merging with FMT dataset
mobilome_features.reset_index(inplace=True)
mobilome_features.rename(columns={'index': 'ID'}, inplace=True)
fmt_dataset['run_accession'] = fmt_dataset['run_accession'].astype(str).str.strip()
mobilome_features['ID'] = mobilome_features['ID'].astype(str).str.strip()
# Merging mobilome features with donor_pre_post column
merged_mobilome_df = mobilome_features.merge(fmt_dataset[['run_accession', 'donor_pre_post', 'Patient']], 
                                             left_on='ID', right_on='run_accession', how='left').drop(columns=['run_accession'])

merged_mobilome_df = merged_mobilome_df[merged_mobilome_df['donor_pre_post'].notna()]

# Applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_mobilome_df.drop(columns=['ID', 'donor_pre_post']))

# Applying CLR transformation with zero replacement
pseudocount = 1e-6
clr_input = imputed_data
clr_input[clr_input <= 0] = pseudocount
clr_transformed = clr(clr_input)

# Handling NaN and infinite values
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
        return  # Skip groups with insufficient points

    mean_x, mean_y = np.mean(x), np.mean(y)
    cov = np.cov(x, y)  # Compute covariance matrix

    # Eigen decomposition to get ellipse parameters
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]  # Sort eigenvalues (largest first)
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # Compute width and height of the ellipse (scaled by chi2 for 95% CI)
    chi2_val = np.sqrt(chi2.ppf(0.95, df=2))  # Scaling factor for 95% confidence
    width, height = 2 * chi2_val * np.sqrt(eigvals)

    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))  # Compute rotation angle

    ellipse = Ellipse(
        xy=(mean_x, mean_y),
        width=width,
        height=height,
        angle=angle,
        edgecolor=color,
        facecolor=color,
        alpha=0.2
    )
    ax.add_patch(ellipse)

# Creating scatter plot and explicitly setting color palette
plt.figure(figsize=(10, 6))
unique_groups = merged_mobilome_df['donor_pre_post'].unique()
group_colors = {
    'Donor': 'blue',
    'PreFMT': 'red',
    'PostFMT': 'green'
}

ax = sns.scatterplot(x='PC1', y='PC2', hue='donor_pre_post', data=merged_mobilome_df, palette=group_colors, alpha=0.7, edgecolor='k')

# Computing confidence ellipses
for group in unique_groups:
    subset = merged_mobilome_df[merged_mobilome_df['donor_pre_post'] == group]
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
legend = plt.legend()
legend.set_visible(False)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/pca_mdrb.svg", format='svg', dpi=600, bbox_inches='tight', transparent=True)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/pca_mdrb.png", format='png', dpi=600, bbox_inches='tight', transparent=True)
plt.show()

# Save metadata: ID, fmt_prep, Patient
metadata_df = merged_mobilome_df[['ID', 'donor_pre_post', 'Patient']]
metadata_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/metadata_mdrb.csv", index=False)

# Save Aitchison distance matrix
aitchison_df = pd.DataFrame(
    aitchison_distances,
    index=merged_mobilome_df['ID'],
    columns=merged_mobilome_df['ID']
)
aitchison_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/aitchison_mdrb.csv")
