import pandas as pd
import matplotlib.pyplot as plt
# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
df = pd.read_csv(file_path)
# Define the sequencing depth using the number of bases in gigabases
df_filtered = df.dropna(subset=["number_bases_gigabases", "Disease_type"])
# Grouping data by disease
disease_types = df_filtered["Disease_type"].unique()
groups = [df_filtered[df_filtered["Disease_type"] == condition]["number_bases_gigabases"] 
          for condition in disease_types]
# Set up for the plot
plt.figure(figsize=(8, 6))
plt.boxplot(groups, labels=disease_types, patch_artist=True)
# Labels and title
plt.xlabel("Condition (Disease Type)")
plt.ylabel("Sequencing Depth (Gbp)")
plt.title("Sequencing Depth by Condition (Disease Type)")
plt.xticks(rotation=45)  # Rotate x-axis labels if needed
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.show()