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
from skbio.stats.distance import DistanceMatrix, permanova
from itertools import combinations
from statsmodels.stats.multitest import multipletests

# File paths
mge_matrix_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\Telcomb_MGE_analytical_matrix.xlsx - Sheet1.csv"
mge_annotations_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\MGE_total_classification.xlsx - Sheet1.csv"
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset.csv"

# Loading datasets
mge_matrix = pd.read_csv(mge_matrix_path)
mge_annotations = pd.read_csv(mge_annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# Remove blank/missing entries from 'fmt_prep'
fmt_dataset = fmt_dataset.dropna(subset=['fmt_prep'])
fmt_dataset = fmt_dataset[fmt_dataset['fmt_prep'].str.strip() != '']

# Filter MGE matrix columns to include only IDs present in fmt_dataset
valid_ids = set(fmt_dataset['run_accession'])
columns_to_keep = ['gene_accession'] + [col for col in mge_matrix.columns if col in valid_ids]
mge_matrix = mge_matrix[columns_to_keep]

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

# Merging mobilome features with fmt_prep column
merged_mobilome_df = mobilome_features.merge(fmt_dataset[['run_accession', 'fmt_prep']], 
                                             left_on='ID', right_on='run_accession', how='left').drop(columns=['run_accession'])

# Standardizing 'fmt_prep' values
merged_mobilome_df['fmt_prep'] = merged_mobilome_df['fmt_prep'].replace({'Frozen': 'frozen'})

# Applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_mobilome_df.drop(columns=['ID', 'fmt_prep']))

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

# Ensuring no NaN values in 'fmt_prep'
merged_mobilome_df['fmt_prep'] = merged_mobilome_df['fmt_prep'].fillna('Unknown')


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
unique_groups = merged_mobilome_df['fmt_prep'].unique()
palette = sns.color_palette('tab10', len(unique_groups))
group_colors = {group: palette[i] for i, group in enumerate(unique_groups)}

ax = sns.scatterplot(x='PC1', y='PC2', hue='fmt_prep', data=merged_mobilome_df, palette=group_colors, alpha=0.7, edgecolor='k')

# Computing confidence ellipses
for group in unique_groups:
    subset = merged_mobilome_df[merged_mobilome_df['fmt_prep'] == group]
    confidence_ellipse(subset['PC1'], subset['PC2'], ax, group_colors.get(group, 'gray'))

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of Aitchison Distances for Mobilome Samples (FMT Preparation)')
plt.legend(title='FMT Preparation', bbox_to_anchor=(1, 1))
plt.grid(True)
plt.show()

#PERMANOVA

pairwise_results = []
unique_groups = merged_mobilome_df['fmt_prep'].dropna().unique()

for group1, group2 in combinations(unique_groups, 2):
    subset_df = merged_mobilome_df[merged_mobilome_df['fmt_prep'].isin([group1, group2])]
    subset_ids = subset_df['ID'].values
    subset_grouping = subset_df['fmt_prep'].values

    # Get corresponding rows in the CLR-transformed matrix
    subset_indices = merged_mobilome_df.index[merged_mobilome_df['fmt_prep'].isin([group1, group2])]
    subset_dist_matrix = squareform(pdist(clr_transformed[subset_indices, :], metric='euclidean'))
    dm = DistanceMatrix(subset_dist_matrix, ids=subset_ids)

    # Run PERMANOVA
    result = permanova(dm, grouping=subset_grouping, permutations=999, seed=42)
    pairwise_results.append((f"{group1} vs {group2}", result['p-value']))

# FDR correction
labels, raw_pvals = zip(*pairwise_results)
fdr_corrected = multipletests(raw_pvals, method='fdr_bh')[1]

# Print results
print("\nPairwise PERMANOVA with FDR Correction (Benjamini-Hochberg):")
print("{:<30} {:<15} {:<15}".format("Comparison", "Raw p-value", "FDR-adjusted p"))
for label, raw, fdr in zip(labels, raw_pvals, fdr_corrected):
    print("{:<30} {:<15.4f} {:<15.4f}".format(label, raw, fdr))