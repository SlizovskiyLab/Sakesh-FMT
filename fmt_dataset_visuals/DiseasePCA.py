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

# loading datasets
amr_matrix = pd.read_csv(amr_matrix_path)
annotations = pd.read_csv(annotations_path)
fmt_dataset = pd.read_csv(fmt_dataset_path)

# removing rows where 'gene_accession' contains "RequiresSNPConfirmation"
amr_matrix_filtered = amr_matrix[~amr_matrix['gene_accession'].str.contains("RequiresSNPConfirmation", na=False)]

# merging with megares annotations to categorize resistance genes
amr_merged = amr_matrix_filtered.merge(annotations[['header', 'group']], left_on='gene_accession', right_on='header', how='left').drop(columns=['header'])

# summing values across genes for each ID, preserving individual resistance features
resistance_features = amr_merged.drop(columns=['gene_accession', 'group']).groupby(amr_merged['group']).sum().T
resistance_features.reset_index(inplace=True)
resistance_features.rename(columns={'index': 'ID'}, inplace=True)

# merging with disease type labels from FMT_dataset
merged_df = resistance_features.merge(fmt_dataset[['run_accession', 'Disease_type']], left_on='ID', right_on='run_accession', how='left')
merged_df.drop(columns=['run_accession'], inplace=True)

# applying Bayesian Missing Data Imputation
imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(merged_df.drop(columns=['ID', 'Disease_type']))

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

# ensuring no NaN values in 'Disease_type'
merged_df['Disease_type'] = merged_df['Disease_type'].fillna('Unknown')

# creating scatter plot and extract color mapping from Seaborn
plt.figure(figsize=(10, 6))
ax = sns.scatterplot(x='PC1', y='PC2', hue='Disease_type', data=merged_df, palette='tab10', alpha=0.7, edgecolor='k')

# extracting colors properly from scatterplot
unique_diseases = merged_df['Disease_type'].unique()
facecolors = ax.collections[0].get_facecolors()
disease_colors = {disease: facecolors[i] for i, disease in enumerate(unique_diseases)}

# computing confidence intervals and add matching colored bubbles
for disease in unique_diseases:
    subset = merged_df[merged_df['Disease_type'] == disease]
    if len(subset) < 2:
        continue

    mean_x, mean_y = subset['PC1'].mean(), subset['PC2'].mean()
    ci_x, ci_y = sem(subset['PC1']) * 1.96, sem(subset['PC2']) * 1.96 
    color = disease_colors.get(disease, 'gray') 
    circle = plt.Circle((mean_x, mean_y), max(ci_x, ci_y), color=color, alpha=0.2)
    plt.gca().add_patch(circle)

plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of Aitchison Distances for Resistome Samples (Disease Type)')
plt.legend(title='Disease Type', bbox_to_anchor=(1, 1))
plt.grid(True)
plt.show()
