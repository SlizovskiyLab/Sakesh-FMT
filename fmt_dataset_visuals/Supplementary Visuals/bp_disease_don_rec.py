import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm

file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_with_throughput.csv"
df = pd.read_csv(file_path)

df_filtered = df[df["donor_pre_post"] != "Pre-Abx/FMT"]
df_filtered["donor_pre_post"] = df_filtered["donor_pre_post"].replace(
    {"Pre-FMT": "PreFMT", "Post-FMT": "PostFMT"}
)
df_filtered = df_filtered[df_filtered["donor_pre_post"] != "PreFMT"].copy()
df_filtered = df_filtered[(df_filtered["timepoint"] > 0) & (df_filtered["timepoint"] <= 200)]

throughput_col = "Throughput"
df_filtered[throughput_col] = df_filtered[throughput_col].astype(str).str.replace(',', '')
df_filtered[throughput_col] = pd.to_numeric(df_filtered[throughput_col], errors='coerce')
df_filtered.dropna(subset=[throughput_col], inplace=True)

df_filtered["log_number_bases"] = np.log10(df_filtered[throughput_col])

plt.figure(figsize=(16, 6))
sns.scatterplot(x="timepoint", y="log_number_bases", hue="Disease_type", data=df_filtered, alpha=0.6)

unique_diseases = df_filtered["Disease_type"].unique()
colors = sns.color_palette("tab10", len(unique_diseases))
n_boot = 1000

for disease, color in zip(unique_diseases, colors):
    subset = df_filtered[df_filtered["Disease_type"] == disease].dropna(subset=["timepoint", "log_number_bases"])
    
    if not subset.empty:
        x_sorted = np.sort(subset["timepoint"])
        y_sorted = subset["log_number_bases"].iloc[np.argsort(subset["timepoint"])]

        lowess_result = sm.nonparametric.lowess(y_sorted, x_sorted, frac=0.7)
        loess_x, loess_y = lowess_result[:, 0], lowess_result[:, 1]

        y_boot_samples = []
        for _ in range(n_boot):
            sample_indices = np.random.choice(len(subset), len(subset), replace=True)
            x_sampled = subset["timepoint"].iloc[sample_indices]
            y_sampled = subset["log_number_bases"].iloc[sample_indices]

            lowess_boot = sm.nonparametric.lowess(y_sampled, x_sampled, frac=0.7)
            y_boot_samples.append(np.interp(loess_x, lowess_boot[:, 0], lowess_boot[:, 1]))

        y_boot_samples = np.array(y_boot_samples)
        lower_bound = np.percentile(y_boot_samples, 2.5, axis=0)
        upper_bound = np.percentile(y_boot_samples, 97.5, axis=0)

        plt.plot(loess_x, loess_y, color=color, linewidth=2, label=f"{disease} LOESS")
        plt.fill_between(loess_x, lower_bound, upper_bound, color=color, alpha=0.2)

plt.xticks(rotation=45, ha="right")
plt.xlabel("Timepoint (Days)")
plt.ylabel("Log Sequencing Depth (Bases)")
plt.title("Log Sequencing Depth Over Time by Disease Type with 95% CI")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig(
    "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\increase_loess.png",
    format='png', dpi=600, bbox_inches='tight'
)
plt.show()
