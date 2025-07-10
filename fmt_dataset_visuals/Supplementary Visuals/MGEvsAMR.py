import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D  # For custom legend handles

# Load datasets
fmt_file = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset.csv"
mge_file = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\Telcomb_MGE_analytical_matrix.xlsx - Sheet1.csv"
amr_file = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\dedup_AMR_analytic_matrix.csv"

fmt_df = pd.read_csv(fmt_file)
mge_df = pd.read_csv(mge_file)
amr_df = pd.read_csv(amr_file)

# Data cleaning
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
study_labels = merged_df["study_data"].unique()
palette = sns.color_palette("tab10", len(study_labels))
marker_styles = {"PreFMT": "o", "PostFMT": "s", "Donor": "^"}

# Scatter plot
for (study, donor_status), subset in merged_df.groupby(["study_data", "donor_pre_post"]):
    ax.scatter(
        subset["AMR"], subset["MGE"],
        s=subset["size"],
        label=f"{study}, {donor_status}",
        color=palette[list(study_labels).index(study)],
        marker=marker_styles.get(donor_status, "o"),
        alpha=0.7
    )

# Set axis labels, log scale, title
ax.set_xlabel("Log AMR (Summed)")
ax.set_ylabel("Log MGE (Summed)")
ax.set_title("Log-scaled MGE vs. AMR by Study Data and Donor Status (Sized by Percentile Bins)")
ax.set_yscale("log")
ax.set_xscale("log")
plt.subplots_adjust(left=0.1, right=0.75) 

# Create custom legends
color_handles = [
    Line2D([0], [0], marker='o', color='w',
           markerfacecolor=palette[i], markersize=10, label=label)
    for i, label in enumerate(study_labels)
]

shape_handles = [
    Line2D([0], [0], marker=marker_styles[status], color='k',
           linestyle='None', markersize=10, label=status)
    for status in marker_styles
]

# Add legends
legend1 = ax.legend(handles=color_handles, title="Study", bbox_to_anchor=(1.05, 1), loc='upper left')
legend2 = ax.legend(handles=shape_handles, title="Donor Status", bbox_to_anchor=(1.05, 0.4), loc='upper left')
ax.add_artist(legend1)

ax.grid(alpha=0.3)
ax.grid(False)
plt.show()
