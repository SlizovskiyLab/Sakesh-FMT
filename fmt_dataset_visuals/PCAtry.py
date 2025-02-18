import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform

# File paths
amr_matrix_path = "/mnt/data/dedup_AMR_analytic_matrix.csv"
annotations_path = "/mnt/data/megares_annotations_v3.00.csv"
fmt_dataset_path = "/mnt/data/FMT_dataset.csv"

# Load datasets
amr_matrix = pd.read_csv(amr_matrix_path)
annotations = pd.read_csv(annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# Step 1: Remove rows where 'gene_accession' contains "RequiresSNPConfirmation"
amr_matrix_filtered = amr_matrix[~amr_matrix['gene_accession'].str.contains("RequiresSNPConfirmation", na=False)]

# Step 2: Sum values across all sample ID columns per ID
sample_id_cols = amr_matrix_filtered.columns[1:]
summed_data = amr_matrix_filtered[sample_id_cols].sum()

# Convert summed data to a dataframe
summed_df = summed_data.reset_index()
summed_df.columns = ['ID', 'Total_Resistance']

# Step 3: Merge with study labels from FMT_dataset
merged_df = summed_df.merge(fmt_dataset[['run_accession', 'study_data']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)

# Step 4: Compute Aitchison distances
clr_transformed = np.log1p(merged_df['Total_Resistance'])
aitchison_distances = squareform(pdist(clr_transformed.values.reshape(-1, 1), metric='euclidean'))

# Step 5: Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(aitchison_distances)
merged_df['PC1'] = pca_result[:, 0]
merged_df['PC2'] = pca_result[:, 1]

# Step 6: Plot PCA results
plt.figure(figsize=(10, 6))
sns.scatterplot(x='PC1', y='PC2', hue='study_data', data=merged_df, palette='tab10', alpha=0.7, edgecolor='k')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of Aitchison Distances for Resistome Samples')
plt.legend(title='Study', bbox_to_anchor=(1, 1))
plt.grid(True)
plt.show()
