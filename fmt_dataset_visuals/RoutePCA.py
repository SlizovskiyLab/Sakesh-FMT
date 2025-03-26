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
fmt_dataset_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset.csv"

# Loading datasets
amr_matrix = pd.read_csv(amr_matrix_path)
annotations = pd.read_csv(annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

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
merged_df = resistance_features.merge(fmt_dataset[['run_accession', 'fmt_route']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)

# Removing rows with missing fmt_route values
merged_df = merged_df.dropna(subset=['fmt_route'])

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

# Creating scatter plot
plt.figure(figsize=(10, 6))
ax = sns.scatterplot(x='PC1', y='PC2', hue='fmt_route', data=merged_df, palette='tab10', alpha=0.7, edgecolor='k')

# Extract colors properly from the scatterplot
handles, labels = ax.get_legend_handles_labels()
legend_colors = {label: handle.get_color() for label, handle in zip(labels, handles)}

# Computing confidence intervals and adding correctly colored bubbles
for route in merged_df['fmt_route'].unique():
    subset = merged_df[merged_df['fmt_route'] == route]
    if len(subset) < 2:
        continue
    
    mean_x, mean_y = subset['PC1'].mean(), subset['PC2'].mean()
    ci_x, ci_y = sem(subset['PC1']) * 1.96, sem(subset['PC2']) * 1.96 
    color = legend_colors.get(route, 'gray') 
    circle = plt.Circle((mean_x, mean_y), max(ci_x, ci_y), color=color, alpha=0.2)
    plt.gca().add_patch(circle)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of Aitchison Distances for Resistome Samples (FMT Route-based)')
plt.legend(title='FMT Route', bbox_to_anchor=(1, 1))
plt.grid(True)
plt.show()
