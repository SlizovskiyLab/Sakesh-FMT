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
from matplotlib.patches import Patch

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

# Ignoring bottom five rows
fmt_dataset = fmt_dataset.iloc[:-5]

# Cleaning donor_pre_post column
fmt_dataset['donor_pre_post'] = fmt_dataset['donor_pre_post'].replace({
    'Pre-FMT': 'PreFMT', 
    'Post-FMT': 'PostFMT',
    'Pre-Abx/FMT': 'PreFMT'
})
fmt_dataset = fmt_dataset[fmt_dataset['donor_pre_post'].isin(['PostFMT', 'PreFMT', 'Donor'])]
fmt_dataset = fmt_dataset.dropna(subset=['donor_pre_post'])

# Removing rows where 'gene_accession' contains "RequiresSNPConfirmation"
amr_matrix_filtered = amr_matrix[~amr_matrix['gene_accession'].str.contains("RequiresSNPConfirmation", na=False)]

# Merging with megares annotations
amr_merged = amr_matrix_filtered.merge(annotations[['header', 'group']], left_on='gene_accession', right_on='header', how='left').drop(columns=['header'])

# Summing values across genes for each ID
resistance_features = amr_merged.drop(columns=['gene_accession', 'group']).groupby(amr_merged['group']).sum().T
resistance_features.reset_index(inplace=True)
resistance_features.rename(columns={'index': 'ID'}, inplace=True)

# --- START: MODIFIED SECTION ---

# Merging with metadata, INCLUDING timepoint
merged_df = resistance_features.merge(fmt_dataset[['run_accession', 'donor_pre_post', 'Patient', 'timepoint']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)
merged_df = merged_df.dropna(subset=['donor_pre_post'])

# Clean the timepoint column
merged_df['timepoint'] = pd.to_numeric(merged_df['timepoint'], errors='coerce')
merged_df.dropna(subset=['timepoint'], inplace=True)
merged_df['timepoint'] = merged_df['timepoint'].astype(int)

# Split 'PostFMT' into time-based bins
merged_df.loc[(merged_df['donor_pre_post'] == 'PostFMT') & (merged_df['timepoint'].between(1, 30)), 'donor_pre_post'] = 'PostFMT (1-30d)'
merged_df.loc[(merged_df['donor_pre_post'] == 'PostFMT') & (merged_df['timepoint'].between(31, 60)), 'donor_pre_post'] = 'PostFMT (31-60d)'
merged_df.loc[(merged_df['donor_pre_post'] == 'PostFMT') & (merged_df['timepoint'] > 60), 'donor_pre_post'] = 'PostFMT (60d+)'

# Applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
# Drop timepoint before imputation
imputed_data = imputer.fit_transform(merged_df.drop(columns=['ID', 'donor_pre_post', 'timepoint']))

# --- END: MODIFIED SECTION ---

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
unique_groups = merged_df['donor_pre_post'].unique()
group_colors = {
    'Donor': '#003771',
    'PreFMT': '#726732',
    'PostFMT (1-30d)': '#b3cde0', # Light Blue
    'PostFMT (31-60d)': '#6497b1', # Medium Blue
    'PostFMT (60d+)': '#005b96', # Dark Blue
}

ax = sns.scatterplot(x='PC1', y='PC2', hue='donor_pre_post', data=merged_df, palette=group_colors, alpha=0.7, edgecolor='k', legend=False)

# computing confidence ellipses
for group in unique_groups:
    subset = merged_df[merged_df['donor_pre_post'] == group]
    confidence_ellipse(subset['PC1'], subset['PC2'], ax, group_colors.get(group, 'gray'))

plt.xlim(merged_df['PC1'].min() - 600, merged_df['PC1'].max() + 550)
plt.ylim(merged_df['PC2'].min() - 300, merged_df['PC2'].max() + 350)

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
legend_handles = [Patch(facecolor=color, edgecolor='k', label=label)
                  for label, color in group_colors.items()]

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

plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/pca_melanoma.svg", format='svg', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight', transparent = True)
plt.savefig("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/pca_melanoma.png", format='png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight', transparent = True)

plt.show()

# Save metadata
metadata_df = merged_df[['ID', 'donor_pre_post', 'Patient', 'timepoint']]
metadata_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/metadata_melanoma.csv", index=False)

# Save Aitchison distance matrix
aitchison_df = pd.DataFrame(
    aitchison_distances,
    index=merged_df['ID'],
    columns=merged_df['ID']
)
aitchison_df.to_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/aitchison_melanoma.csv")