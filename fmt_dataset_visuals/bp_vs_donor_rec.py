import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
df = pd.read_csv(file_path)

# Define the sequencing depth using the number of bases in gigabases
df_filtered = df.dropna(subset=["number_bases_gigabases", "donor_pre_post"])

# Replace 'Pre-FMT' with 'PreFMT' and 'Post-FMT' with 'PostFMT'
df_filtered['donor_pre_post'] = df_filtered['donor_pre_post'].replace({'Pre-FMT': 'PreFMT', 'Post-FMT': 'PostFMT'})

# Grouping data by disease
donor_patient_types = df_filtered["donor_pre_post"].unique()
groups = [df_filtered[df_filtered["donor_pre_post"] == condition]["number_bases_gigabases"] 
          for condition in donor_patient_types]

# Set up for the plot
plt.figure(figsize=(8, 6))
plt.boxplot(groups, tick_labels=donor_patient_types, patch_artist=True)

# Labels and title
plt.xlabel("Donor, PreFMT, PostFMT")
plt.ylabel("Sequencing Depth (Gbp)")
plt.title("Sequencing Depth by Condition (Disease Type)")
plt.xticks(rotation=45)  # Rotate x-axis labels if needed
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.show()
