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

# Convert sequencing depth to log scale
df_filtered["log_number_bases_gigabases"] = np.log10(df_filtered["number_bases_gigabases"])

plt.figure(figsize=(16, 6))

sns.scatterplot(x="timepoint", y="log_number_bases_gigabases", hue="Disease_type", data=df_filtered, alpha=0.6)

# Fit LOESS curve for smoothing and compute confidence intervals using bootstrapping
unique_diseases = df_filtered["Disease_type"].unique()
colors = sns.color_palette("tab10", len(unique_diseases))

n_boot = 1000

for disease, color in zip(unique_diseases, colors):
    subset = df_filtered[df_filtered["Disease_type"] == disease].dropna(subset=["timepoint", "log_number_bases_gigabases"])
    
    if not subset.empty:
        x_sorted = np.sort(subset["timepoint"])
        y_sorted = subset["log_number_bases_gigabases"].iloc[np.argsort(subset["timepoint"])]

        # Compute LOESS fit
        lowess_result = sm.nonparametric.lowess(y_sorted, x_sorted, frac=0.5)
        loess_x, loess_y = lowess_result[:, 0], lowess_result[:, 1]

        # Bootstrapping for confidence interval
        y_boot_samples = []
        for _ in range(n_boot):
            sample_indices = np.random.choice(len(subset), len(subset), replace=True)
            x_sampled = subset["timepoint"].iloc[sample_indices]
            y_sampled = subset["log_number_bases_gigabases"].iloc[sample_indices]

            lowess_boot = sm.nonparametric.lowess(y_sampled, x_sampled, frac=0.5)
            y_boot_samples.append(np.interp(loess_x, lowess_boot[:, 0], lowess_boot[:, 1]))

        # Compute percentiles for 95% CI
        y_boot_samples = np.array(y_boot_samples)
        lower_bound = np.percentile(y_boot_samples, 2.5, axis=0)
        upper_bound = np.percentile(y_boot_samples, 97.5, axis=0)

        # Plot LOESS curve
        plt.plot(loess_x, loess_y, color=color, linewidth=2, label=f"{disease} LOESS")
        
        # Plot confidence interval
        plt.fill_between(loess_x, lower_bound, upper_bound, color=color, alpha=0.2)

# Customize x-axis for better spacing
plt.xticks(rotation=45, ha="right")
plt.xlabel("Timepoint (Days)")
plt.ylabel("Log Sequencing Depth (Gigabases)")
plt.title("Log Sequencing Depth Over Time by Disease Type with 95% CI")
plt.legend()
plt.grid(alpha=0.3)
plt.show()
