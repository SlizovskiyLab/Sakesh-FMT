import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from skbio.stats.composition import clr
from sklearn.impute import KNNImputer
from scipy.stats import sem

# File paths
amr_matrix_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\dedup_AMR_analytic_matrix.csv"
annotations_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\megares_annotations_v3.00.csv"
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"

# Loading datasets
amr_matrix = pd.read_csv(amr_matrix_path)
annotations = pd.read_csv(annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# Ignoring bottom five rows
fmt_dataset = fmt_dataset.iloc[:-5]

# Cleaning donor_pre_post column
fmt_dataset['donor_pre_post'] = fmt_dataset['donor_pre_post'].replace({
    'Pre-FMT': 'PreFMT', 
    'Post-FMT': 'PostFMT'
})

# Filtering out 'Pre-Abx/FMT'
fmt_dataset = fmt_dataset[fmt_dataset['donor_pre_post'] != 'Pre-Abx/FMT']

# Drop any NaN values in donor_pre_post before merging
fmt_dataset = fmt_dataset.dropna(subset=['donor_pre_post'])

# Removing rows where 'gene_accession' contains "RequiresSNPConfirmation"
amr_matrix_filtered = amr_matrix[~amr_matrix['gene_accession'].str.contains("RequiresSNPConfirmation", na=False)]

# Merging with megares annotations
amr_merged = amr_matrix_filtered.merge(annotations[['header', 'group']], left_on='gene_accession', right_on='header', how='left').drop(columns=['header'])

# Summing values across genes for each ID
resistance_features = amr_merged.drop(columns=['gene_accession', 'group']).groupby(amr_merged['group']).sum().T
resistance_features.reset_index(inplace=True)
resistance_features.rename(columns={'index': 'ID'}, inplace=True)

# Merging with donor_pre_post labels
merged_df = resistance_features.merge(fmt_dataset[['run_accession', 'donor_pre_post']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)

# Drop remaining NaN values in donor_pre_post
merged_df = merged_df.dropna(subset=['donor_pre_post'])

# Keep only valid categories
valid_categories = {'PreFMT', 'PostFMT', 'Donor'}
merged_df = merged_df[merged_df['donor_pre_post'].isin(valid_categories)]

# Applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_df.drop(columns=['ID', 'donor_pre_post']))

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

# Creating scatter plot
plt.figure(figsize=(10, 6))
ax = sns.scatterplot(x='PC1', y='PC2', hue='donor_pre_post', data=merged_df, palette='tab10', alpha=0.7, edgecolor='k')

# Extracting unique categories for coloring
handles, labels = ax.get_legend_handles_labels()
category_colors = {label: handle.get_color() for label, handle in zip(labels, handles)}

# Adding confidence ellipses
for category in merged_df['donor_pre_post'].unique():
    subset = merged_df[merged_df['donor_pre_post'] == category]
    if len(subset) < 2:
        continue
    mean_x, mean_y = subset['PC1'].mean(), subset['PC2'].mean()
    ci_x, ci_y = sem(subset['PC1']) * 1.96, sem(subset['PC2']) * 1.96
    color = category_colors.get(category, 'gray')
    circle = plt.Circle((mean_x, mean_y), max(ci_x, ci_y), color=color, alpha=0.2)
    plt.gca().add_patch(circle)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of Aitchison Distances for Resistome Samples (Donor Pre/Post)')
plt.legend(title='Donor Pre/Post', bbox_to_anchor=(1, 1))
plt.grid(True)
plt.show()