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
# --- Use the new binned dataset ---
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_paired_binned.csv"

# Loading datasets
amr_matrix = pd.read_csv(amr_matrix_path)
annotations = pd.read_csv(annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# Filter for disease and clean data
fmt_dataset = fmt_dataset.dropna(subset=['Patient'])
fmt_dataset = fmt_dataset[fmt_dataset['Patient'].astype(str).str.strip() != '']
fmt_dataset = fmt_dataset[fmt_dataset['Disease_type'] == 'MDRB']
fmt_dataset = fmt_dataset.iloc[:-5]
fmt_dataset = fmt_dataset.dropna(subset=['donor_pre_post_binned'])

# Removing rows where 'gene_accession' contains "RequiresSNPConfirmation"
amr_matrix_filtered = amr_matrix[~amr_matrix['gene_accession'].str.contains("RequiresSNPConfirmation", na=False)]

# Merging with megares annotations
amr_merged = amr_matrix_filtered.merge(annotations[['header', 'group']], left_on='gene_accession', right_on='header', how='left').drop(columns=['header'])

# Summing values across genes for each ID
resistance_features = amr_merged.drop(columns=['gene_accession', 'group']).groupby(amr_merged['group']).sum().T
resistance_features.reset_index(inplace=True)
resistance_features.rename(columns={'index': 'ID'}, inplace=True)

# --- SIMPLIFIED MERGE ---
# Merge directly with the new 'donor_pre_post_binned' column
merged_df = resistance_features.merge(fmt_dataset[['run_accession', 'donor_pre_post_binned', 'Patient']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)
merged_df = merged_df.dropna(subset=['donor_pre_post_binned'])

# Applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
# Use the new binned column for imputation
imputed_data = imputer.fit_transform(merged_df.drop(columns=['ID', 'donor_pre_post_binned']))

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
unique_groups = merged_df['donor_pre_post_binned'].unique()
group_colors = {
    'Donor': '#003771',
    'PreFMT': '#726732',
    'PostFMT (1-30d)': '#b3cde0', # Light Blue
    'PostFMT (31-60d)': '#6497b1', # Medium Blue
    'PostFMT (60d+)': '#005b96', # Dark Blue
}

ax = sns.scatterplot(x='PC1', y='PC2', hue='donor_pre_post_binned', data=merged_df, palette=group_colors, alpha=0.7, edgecolor='k')

# computing confidence ellipses
for group in unique_groups:
    subset = merged_df[merged_df['donor_pre_post_binned'] == group]
    confidence_ellipse(subset['PC1'], subset['PC2'], ax, group_colors.get(group, 'gray'))

plt.xlim(merged_df['PC1'].min() - 500, merged_df['PC1'].max() + 550)
plt.ylim(merged_df['PC2'].min() - 100, merged_df['PC2'].max() + 500)

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
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/pca_mdrb.svg", format='svg', dpi=600, bbox_inches='tight', transparent = True)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/pca_mdrb.png", format='png', dpi=600, bbox_inches='tight', transparent = True)

plt.show()

# Save metadata with the correct column
metadata_df = merged_df[['ID', 'donor_pre_post_binned', 'Patient']]
metadata_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/metadata_mdrb.csv", index=False)

# Save Aitchison distance matrix
aitchison_df = pd.DataFrame(
    aitchison_distances,
    index=merged_df['ID'],
    columns=merged_df['ID']
)
aitchison_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/aitchison_mdrb.csv")
