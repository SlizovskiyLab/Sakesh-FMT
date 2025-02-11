import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
df = pd.read_csv(file_path)

# Filter out PreFMT samples
df_filtered = df[df["donor_pre_post"] != "PreFMT"].copy()

# Ensure Post-FMT is labeled consistently
df_filtered["donor_pre_post"] = df_filtered["donor_pre_post"].replace({"Post-FMT": "PostFMT"})

# Remove data points above 200 days and at 0 days
df_filtered = df_filtered[(df_filtered["timepoint"] > 0) & (df_filtered["timepoint"] <= 200)]

# Convert host removal input reads to numeric and scale to millions
df_filtered["host_removal_input_reads"] = pd.to_numeric(df_filtered["host_removal_input_reads"], errors='coerce')
df_filtered["host_removal_input_reads_millions"] = df_filtered["host_removal_input_reads"] / 1e6

# Initialize a wider figure
plt.figure(figsize=(16, 6))

# Scatter plot using seaborn scatterplot
sns.scatterplot(x="timepoint", y="host_removal_input_reads_millions", hue="Disease_type", data=df_filtered, alpha=0.6)

# Fit LOESS curve for smoothing
unique_diseases = df_filtered["Disease_type"].unique()
colors = sns.color_palette("tab10", len(unique_diseases))

for disease, color in zip(unique_diseases, colors):
    subset = df_filtered[df_filtered["Disease_type"] == disease]
    subset = subset.dropna(subset=["timepoint", "host_removal_input_reads_millions"])
    
    if not subset.empty:
        x_sorted = np.sort(subset["timepoint"])
        y_sorted = subset["host_removal_input_reads_millions"].iloc[np.argsort(subset["timepoint"])]
        lowess_result = sm.nonparametric.lowess(y_sorted, x_sorted, frac=0.5)  # Increased smoothing
        plt.plot(lowess_result[:, 0], lowess_result[:, 1], color=color, linewidth=2, label=f"{disease} LOESS")

# Customize x-axis for better spacing
plt.xticks(rotation=45, ha="right")
plt.xlabel("Timepoint (Days)")
plt.ylabel("Host Removal Input Reads (Millions)")
plt.title("Host Removal Input Reads Over Time by Disease Type")
plt.legend()
plt.grid(alpha=0.3)

# Show the updated plot
plt.show()
