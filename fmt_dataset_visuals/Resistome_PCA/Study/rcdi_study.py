import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from skbio.stats.composition import clr
from sklearn.impute import KNNImputer
from scipy.stats import sem
from matplotlib.patches import Ellipse
from scipy.stats import chi2

# File paths
amr_matrix_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\dedup_AMR_analytic_matrix.csv"
annotations_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\megares_annotations_v3.00.csv"
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_paired.csv"

# loading datasets
amr_matrix = pd.read_csv(amr_matrix_path)
annotations = pd.read_csv(annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# Remove rows where 'Patient' is missing or blank
fmt_dataset = fmt_dataset.dropna(subset=['Patient'])
fmt_dataset = fmt_dataset[fmt_dataset['Patient'].astype(str).str.strip() != '']

fmt_dataset = fmt_dataset[fmt_dataset['Disease_type'] == 'rCDI']

# removing rows where 'gene_accession' contains "RequiresSNPConfirmation"
amr_matrix_filtered = amr_matrix[~amr_matrix['gene_accession'].str.contains("RequiresSNPConfirmation", na=False)]

# merging with megares annotations to categorize resistance genes
amr_merged = amr_matrix_filtered.merge(annotations[['header', 'group']], left_on='gene_accession', right_on='header', how='left').drop(columns=['header'])

# summing values across genes for each ID, preserving individual resistance features
resistance_features = amr_merged.drop(columns=['gene_accession', 'group']).groupby(amr_merged['group']).sum().T
resistance_features.reset_index(inplace=True)
resistance_features.rename(columns={'index': 'ID'}, inplace=True)

# merging with study labels from FMT_dataset
merged_df = resistance_features.merge(fmt_dataset[['run_accession', 'study_data', 'Patient']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)

merged_df = merged_df[merged_df['study_data'].notna()]
merged_df = merged_df[merged_df['study_data'].astype(str).str.strip() != '']

# applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_df.drop(columns=['ID', 'study_data']))

# applying CLR transformation using scikit-bio with zero replacement
pseudocount = 1e-6
clr_input = imputed_data
clr_input[clr_input <= 0] = pseudocount 

clr_transformed = clr(clr_input)
if np.any(np.isnan(clr_transformed)) or np.any(np.isinf(clr_transformed)):
    raise ValueError("CLR transformation produced NaN or infinite values.")

# handling NaN and infinite values
clr_transformed[np.isinf(clr_transformed)] = 0
clr_transformed = np.nan_to_num(clr_transformed)

# computing Aitchison distances
aitchison_distances = squareform(pdist(clr_transformed, metric='euclidean'))

# performing PCA
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
unique_diseases = merged_df['study_data'].unique()
palette = sns.color_palette('tab10', len(unique_diseases))
disease_colors = {disease: palette[i] for i, disease in enumerate(unique_diseases)}

ax = sns.scatterplot(x='PC1', y='PC2', hue='study_data', data=merged_df, palette=disease_colors, alpha=0.7, edgecolor='k')

# computing confidence ellipses
for disease in unique_diseases:
    subset = merged_df[merged_df['study_data'] == disease]
    confidence_ellipse(subset['PC1'], subset['PC2'], ax, disease_colors.get(disease, 'gray'))

plt.xlim(merged_df['PC1'].min() - 1300, merged_df['PC1'].max() + 1200)
plt.ylim(merged_df['PC2'].min() - 700, merged_df['PC2'].max() + 700)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Study Data', bbox_to_anchor=(1, 1), loc='upper right')
plt.grid(True)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/pca_rcdi.svg", format='svg', dpi=600, bbox_inches='tight')

plt.show()

# Save metadata: ID, fmt_prep, Patient
metadata_df = merged_df[['ID', 'study_data', 'Patient']]
metadata_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/metadata_rcdi.csv", index=False)

# Save Aitchison distance matrix
aitchison_df = pd.DataFrame(
    aitchison_distances,
    index=merged_df['ID'],
    columns=merged_df['ID']
)
aitchison_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/aitchison_rcdi.csv")