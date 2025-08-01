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

fmt_dataset = fmt_dataset[fmt_dataset['Disease_type'] == 'Melanoma']

fmt_dataset['DNA_extraction_kit'] = fmt_dataset['DNA_extraction_kit'].replace({
    'DNeasy_PowerSoil_Pro_Kit': 'DNeasy Pro',
    'MagCore_Genomic_DNA_Tissue_Kit_MODIFIED': 'MagCore',
    'PureLink_Genomic_DNA_Mini_Kit': 'PureLink mini',
    'DNeasy_PowerSoil_Kit': 'DNeasy',
    'Purelink_Microbiome_DNA_purification_Kit': 'PureLink'
})

# Clean FMT dataset
fmt_dataset = fmt_dataset.dropna(subset=['Patient'])
fmt_dataset = fmt_dataset[fmt_dataset['Patient'].astype(str).str.strip() != '']
fmt_dataset = fmt_dataset.dropna(subset=['DNA_extraction_kit'])
fmt_dataset = fmt_dataset[fmt_dataset['DNA_extraction_kit'].str.strip() != '']

# Rename for merging
mge_annotations = mge_annotations.rename(columns={'IDs': 'gene_accession'})

# Merge and aggregate
mge_merged = mge_matrix.merge(
    mge_annotations[['gene_accession', 'final_classification']],
    on='gene_accession',
    how='left'
)
mobilome_features = mge_merged.drop(columns=['gene_accession']).groupby('final_classification').sum().T
mobilome_features.reset_index(inplace=True)
mobilome_features.rename(columns={'index': 'ID'}, inplace=True)

# Merge with metadata
merged_mobilome_df = mobilome_features.merge(
    fmt_dataset[['run_accession', 'DNA_extraction_kit', 'Patient']],
    left_on='ID', right_on='run_accession', how='left'
).drop(columns=['run_accession'])

merged_mobilome_df = merged_mobilome_df[merged_mobilome_df['DNA_extraction_kit'].notna()]


# KNN imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_mobilome_df.drop(columns=['ID', 'DNA_extraction_kit']))

# CLR transformation
pseudocount = 1e-6
clr_input = imputed_data
clr_input[clr_input <= 0] = pseudocount
clr_transformed = clr(clr_input)
clr_transformed[np.isinf(clr_transformed)] = 0
clr_transformed = np.nan_to_num(clr_transformed)

# Aitchison distance
aitchison_distances = squareform(pdist(clr_transformed, metric='euclidean'))

# PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(aitchison_distances)
merged_mobilome_df['PC1'] = pca_result[:, 0]
merged_mobilome_df['PC2'] = pca_result[:, 1]

# Unmatched IDs
unmatched_ids = merged_mobilome_df[merged_mobilome_df['DNA_extraction_kit'].isna()]['ID']
print("Unmatched IDs:\n", unmatched_ids)

# Confidence ellipse
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
    ellipse = Ellipse(
        xy=(mean_x, mean_y),
        width=width, height=height, angle=angle,
        edgecolor=color, facecolor=color, alpha=0.2
    )
    ax.add_patch(ellipse)

# Plotting
plt.figure(figsize=(10, 6))
unique_kits = merged_mobilome_df['DNA_extraction_kit'].unique()
kit_colors = {
    'DNeasy Pro': 'red',
    'MagCore': 'blue',
    'PureLink mini': 'green',
    'DNeasy': 'orange',
    'QIAamp' : 'Purple',
    'MoBio PS': 'brown',
    'NucleoSpin': 'yellow',
    'MoBio MB': 'Pink',
    'Phenol': 'lightblue',
    'Zymo': 'maroon',
    'PureLink' : 'coral'
}

ax = sns.scatterplot(
    x='PC1', y='PC2',
    hue='DNA_extraction_kit',
    data=merged_mobilome_df,
    palette=kit_colors,
    alpha=0.7,
    edgecolor='k'
)

# Ellipses
for kit in unique_kits:
    subset = merged_mobilome_df[merged_mobilome_df['DNA_extraction_kit'] == kit]
    confidence_ellipse(subset['PC1'], subset['PC2'], ax, kit_colors.get(kit, 'gray'))

# Expand axis limits
plt.xlim(merged_mobilome_df['PC1'].min() - 100, merged_mobilome_df['PC1'].max() + 100)
plt.ylim(merged_mobilome_df['PC2'].min() - 100, merged_mobilome_df['PC2'].max() + 100)

# Adjust legend
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

# Save plot as SVG
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/pca_melanoma.svg", format='svg', dpi=600, bbox_inches='tight', transparent=True)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/pca_melanoma.png", format='png', dpi=600, bbox_inches='tight', transparent=True)
plt.show()

# Save metadata: ID, fmt_prep, Patient
metadata_df = merged_mobilome_df[['ID', 'DNA_extraction_kit', 'Patient']]
metadata_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/metadata_melanoma.csv", index=False)

# Save Aitchison distance matrix
aitchison_df = pd.DataFrame(
    aitchison_distances,
    index=merged_mobilome_df['ID'],
    columns=merged_mobilome_df['ID']
)
aitchison_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/aitchison_melanoma.csv")
