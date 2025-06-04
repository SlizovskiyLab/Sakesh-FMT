import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load datasets
fmt_file = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
mge_file = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\Telcomb_MGE_analytical_matrix.xlsx - Sheet1.csv"
amr_file = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\dedup_AMR_analytic_matrix.csv"

fmt_df = pd.read_csv(fmt_file)
mge_df = pd.read_csv(mge_file)
amr_df = pd.read_csv(amr_file)

# data cleaning
fmt_df = fmt_df.iloc[:264]
fmt_df["donor_pre_post"] = fmt_df["donor_pre_post"].replace({"Pre-FMT": "PreFMT", "Post-FMT": "PostFMT"})
fmt_df = fmt_df[~fmt_df["donor_pre_post"].isin(["Pre-Abx/FMT", "Max:"])]

# Sum MGE and AMR values based on column IDs
mge_sums = mge_df.set_index(mge_df.columns[0]).sum(axis=0)
amr_sums = amr_df.set_index(amr_df.columns[0]).sum(axis=0)

# Combine MGE and AMR data into a single dataframe
combined_df = pd.DataFrame({"ID": mge_sums.index, "MGE": mge_sums.values, "AMR": amr_sums.values})

# Merge with FMT dataset based on matching IDs
fmt_df = fmt_df.set_index(fmt_df.columns[0])
merged_df = fmt_df.merge(combined_df, left_index=True, right_on="ID")

# Sort by number_bases_gigabases and assign percentile bins
merged_df = merged_df.sort_values("number_bases_gigabases").reset_index()
merged_df["size_bin"] = pd.qcut(merged_df.index, 10, labels=False)  # 10 equal bins

# Scale sizes for better visualization
merged_df["size"] = (merged_df["size_bin"] + 1) * 10 

fig, ax = plt.subplots(figsize=(14, 8))  # Increased figure size

# Define color and shape mapping
palette = sns.color_palette("tab10", len(merged_df["study_data"].unique()))
marker_styles = {"PreFMT": "o", "PostFMT": "s", "Donor": "^"}

# Scatter plot with different colors for study_data and shapes for donor_pre_post
for (study, donor_status), subset in merged_df.groupby(["study_data", "donor_pre_post"]):
    ax.scatter(
        subset["AMR"], subset["MGE"],
        s=subset["size"],  # Set sizes based on percentile bins
        label=f"{study}, {donor_status}",
        color=palette[list(merged_df["study_data"].unique()).index(study)],
        marker=marker_styles.get(donor_status, "o"), alpha=0.7
    )

ax.set_xlabel("AMR (Summed)")
ax.set_ylabel("MGE (Summed)")
ax.set_title("MGE vs. AMR by Study Data and Donor Status (Sized by Percentile Bins)")
plt.subplots_adjust(left=0.1, right=0.75) 
ax.legend(title="Study & Donor Status", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small', frameon=True)
ax.grid(alpha=0.3)
plt.show()