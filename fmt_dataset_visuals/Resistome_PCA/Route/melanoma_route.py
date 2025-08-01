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

# File paths
amr_matrix_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\dedup_AMR_analytic_matrix.csv"
annotations_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\megares_annotations_v3.00.csv"
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_paired.csv"

# Loading datasets
amr_matrix = pd.read_csv(amr_matrix_path)
annotations = pd.read_csv(annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# Remove rows where 'Patient' is missing or blank
fmt_dataset = fmt_dataset.dropna(subset=['Patient'])
fmt_dataset = fmt_dataset[fmt_dataset['Patient'].astype(str).str.strip() != '']

fmt_dataset = fmt_dataset[fmt_dataset['Disease_type'] == 'Melanoma']


# Removing rows where 'gene_accession' contains "RequiresSNPConfirmation"
amr_matrix_filtered = amr_matrix[~amr_matrix['gene_accession'].str.contains("RequiresSNPConfirmation", na=False)]

# Merging with megares annotations to categorize resistance genes
amr_merged = amr_matrix_filtered.merge(annotations[['header', 'group']], left_on='gene_accession', right_on='header', how='left').drop(columns=['header'])

# Summing values across genes for each ID, preserving individual resistance features
resistance_features = amr_merged.drop(columns=['gene_accession', 'group']).groupby(amr_merged['group']).sum().T
resistance_features.reset_index(inplace=True)
resistance_features.rename(columns={'index': 'ID'}, inplace=True)

# Standardizing fmt_route values
fmt_dataset['fmt_route'] = fmt_dataset['fmt_route'].str.strip()
fmt_dataset = fmt_dataset.dropna(subset=['fmt_route'])

# Merging with fmt_route data from FMT dataset
merged_df = resistance_features.merge(fmt_dataset[['run_accession', 'fmt_route', 'Patient']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)

# Removing rows with missing fmt_route values
merged_df = merged_df[merged_df['fmt_route'].notna()]
merged_df = merged_df[merged_df['fmt_route'].astype(str).str.strip() != '']

# Applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_df.drop(columns=['ID', 'fmt_route']))

# Applying CLR transformation using scikit-bio with zero replacement
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
merged_df['PC1'] = pca_result[:, 0]
merged_df['PC2'] = pca_result[:, 1]

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

# creating scatter plot and explicitly setting color palette
plt.figure(figsize=(10, 6))
unique_diseases = merged_df['fmt_route'].unique()
palette = sns.color_palette('tab10', len(unique_diseases))
disease_colors = {disease: palette[i] for i, disease in enumerate(unique_diseases)}

ax = sns.scatterplot(x='PC1', y='PC2', hue='fmt_route', data=merged_df, palette=disease_colors, alpha=0.7, edgecolor='k')

# computing confidence ellipses
for disease in unique_diseases:
    subset = merged_df[merged_df['fmt_route'] == disease]
    confidence_ellipse(subset['PC1'], subset['PC2'], ax, disease_colors.get(disease, 'gray'))

plt.xlim(merged_df['PC1'].min() - 600, merged_df['PC1'].max() + 550)
plt.ylim(merged_df['PC2'].min() - 100, merged_df['PC2'].max() + 200)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='FMT Route', bbox_to_anchor=(1, 1), loc='upper right')
plt.grid(True)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/pca_melanoma.svg", format='svg', dpi=600, bbox_inches='tight')

plt.show()

# Save metadata: ID, fmt_prep, Patient
metadata_df = merged_df[['ID', 'fmt_route', 'Patient']]
metadata_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/metadata_melanoma.csv", index=False)

# Save Aitchison distance matrix
aitchison_df = pd.DataFrame(
    aitchison_distances,
    index=merged_df['ID'],
    columns=merged_df['ID']
)
aitchison_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/aitchison_melanoma.csv")