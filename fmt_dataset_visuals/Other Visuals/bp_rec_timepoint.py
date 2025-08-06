import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_with_throughput.csv"
df = pd.read_csv(file_path)

df_filtered = df[df["donor_pre_post"] != "Pre-Abx/FMT"]
df_filtered["donor_pre_post"] = df_filtered["donor_pre_post"].replace(
    {"Pre-FMT": "PreFMT", "Post-FMT": "PostFMT"}
)

# Filter out PreFMT samples
df_filtered = df_filtered[df_filtered["donor_pre_post"] != "PreFMT"].copy()

# Remove data points above 200 days and at 0 days
df_filtered = df_filtered[(df_filtered["timepoint"] > 0) & (df_filtered["timepoint"] <= 200)]

# Cap MDRB timepoints at 175 days
df_filtered = df_filtered[~((df_filtered['Disease_type'] == 'MDRB') & (df_filtered['timepoint'] > 175))]

# Convert the throughput column to a numeric type
throughput_col = "Throughput"
df_filtered[throughput_col] = df_filtered[throughput_col].astype(str).str.replace(',', '')
df_filtered[throughput_col] = pd.to_numeric(df_filtered[throughput_col], errors='coerce')
df_filtered.dropna(subset=[throughput_col], inplace=True)


# --- Direct Logging Data Transformation ---
# Filter out non-positive values before logging to prevent errors
df_filtered = df_filtered[df_filtered[throughput_col] > 0].copy()

# Create a new log-transformed column for throughput ONLY
df_filtered["log_throughput"] = np.log10(df_filtered[throughput_col])


# --- Plotting Section ---
plt.figure(figsize=(16, 6))

# Use linear timepoint for x-axis and log_throughput for y-axis
sns.scatterplot(x="timepoint", y="log_throughput", hue="Disease_type", data=df_filtered, alpha=0.6)

# Fit LOESS curve and compute confidence intervals
unique_diseases = df_filtered["Disease_type"].unique()
colors = sns.color_palette("tab10", len(unique_diseases))

n_boot = 1000

for disease, color in zip(unique_diseases, colors):
    subset = df_filtered[df_filtered["Disease_type"] == disease].dropna(subset=["timepoint", "log_throughput"])
    
    if not subset.empty:
        # Sort by the linear timepoint
        x_sorted = np.sort(subset["timepoint"])
        y_sorted = subset["log_throughput"].iloc[np.argsort(subset["timepoint"])]

        # Compute LOESS fit on semi-log data
        lowess_result = sm.nonparametric.lowess(y_sorted, x_sorted, frac=0.8)
        loess_x, loess_y = lowess_result[:, 0], lowess_result[:, 1]

        # Bootstrapping for confidence interval
        y_boot_samples = []
        for _ in range(n_boot):
            sample_indices = np.random.choice(len(subset), len(subset), replace=True)
            x_sampled = subset["timepoint"].iloc[sample_indices]
            y_sampled = subset["log_throughput"].iloc[sample_indices]

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

# Customize plot labels and title for the semi-log scale
plt.xticks(rotation=45, ha="right")
plt.xlabel("Timepoint (Days)")
plt.ylabel("Log10(Throughput (Bases))")
plt.title("Throughput Over Time by Disease Type with 95% CI")
plt.legend()
plt.grid(alpha=0.3)
plt.show()