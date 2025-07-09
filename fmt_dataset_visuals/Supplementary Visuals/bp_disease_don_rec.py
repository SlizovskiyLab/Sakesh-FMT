import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Load the dataset
# You will need to update this path to your local file
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_with_throughput.csv"
df = pd.read_csv(file_path)
df_filtered = df.iloc[:264]

# Clean donor_pre_post column
df_filtered = df_filtered[df_filtered["donor_pre_post"] != "Pre-Abx/FMT"]
df_filtered["donor_pre_post"] = df_filtered["donor_pre_post"].replace(
    {"Pre-FMT": "PreFMT", "Post-FMT": "PostFMT"}
)

# Create a new column combining Disease_type and donor_pre_post
df_filtered["X_Label"] = df_filtered["Disease_type"] + ", " + df_filtered["donor_pre_post"]

# Convert the throughput column to a numeric type
throughput_col = "Throughput"
df_filtered[throughput_col] = df_filtered[throughput_col].astype(str).str.replace(',', '')
df_filtered[throughput_col] = pd.to_numeric(df_filtered[throughput_col], errors='coerce')
df_filtered.dropna(subset=[throughput_col], inplace=True)

# Create a new column for throughput in Gigabases
df_filtered["Throughput_GB"] = df_filtered[throughput_col] / 1e9

# --- Plotting ---
plt.figure(figsize=(12, 8))

# capture the axes object returned by the plot for later manipulation
ax = sns.violinplot(x="X_Label", y="Throughput_GB", data=df_filtered, hue="X_Label", palette="muted", inner=None, linewidth=2, legend=False)
sns.boxplot(x="X_Label", y="Throughput_GB", data=df_filtered, color="white", width=0.2, boxprops=dict(alpha=0.7), ax=ax)
sns.stripplot(x="X_Label", y="Throughput_GB", data=df_filtered, color="black", alpha=0.4, jitter=True, ax=ax, marker='.', size=8)

# This loop goes through the violin shapes and adjusts their coordinates
for violin in ax.collections:
    for path in violin.get_paths():
        vertices = path.vertices
        vertices[:, 1] = np.maximum(0, vertices[:, 1])

# plot readability
plt.yticks(np.arange(0, 26, 5))
plt.ylim(-1, 25)
plt.xticks(rotation=45, ha="right")
plt.xlabel("Disease Type & Donor/Recipient Status")
plt.ylabel("Throughput (Gigabases)")
plt.title("Throughput by Disease Type and Donor/Recipient Status")
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.savefig("C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\bp_disease_donor_rec.png", format='png', dpi=600, bbox_inches='tight')

plt.show()