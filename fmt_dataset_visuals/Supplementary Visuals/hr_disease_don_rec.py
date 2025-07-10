import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset.csv"
df = pd.read_csv(file_path)
df_filtered = df.iloc[:264]

# clean donor_pre_post column
df_filtered = df_filtered[df_filtered["donor_pre_post"] != "Pre-Abx/FMT"]
df_filtered["donor_pre_post"] = df_filtered["donor_pre_post"].replace(
    {"Pre-FMT": "PreFMT", "Post-FMT": "PostFMT"}
)

# new column combining Disease_type and donor_pre_post
df_filtered["X_Label"] = df_filtered["Disease_type"] + ", " + df_filtered["donor_pre_post"]

# Convert host_removal_input_reads to numeric and divide by a million
df_filtered["host_removal_input_reads"] = pd.to_numeric(df_filtered["host_removal_input_reads"], errors="coerce")
df_filtered["host_removal_input_reads_millions"] = df_filtered["host_removal_input_reads"] / 1e6

# violin plot
plt.figure(figsize=(12, 6))
ax = sns.violinplot(x="X_Label", y="host_removal_input_reads_millions", data=df_filtered, palette="muted", inner=None, linewidth =3)
sns.boxplot(x="X_Label", y = "host_removal_input_reads_millions", data=df_filtered, color = "white", width=0.2)
sns.stripplot(x="X_Label", y="host_removal_input_reads_millions", data=df_filtered, color="black", alpha=0.3, jitter=True)

# This loop goes through the violin shapes and adjusts their coordinates
for violin in ax.collections:
    for path in violin.get_paths():
        vertices = path.vertices
        vertices[:, 1] = np.maximum(0, vertices[:, 1])

# plot readability
plt.xticks(rotation=20, ha="right")
plt.ylim(-10, 210)
plt.xlabel("Disease Type & Donor/Recipient Status")
plt.ylabel("Host Removal Input Reads (Millions)")
plt.title("Violin Plot of Host Removal Input Reads by Disease Type and Donor/Recipient Status")
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.show()
