import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
df = pd.read_csv(file_path)

# Step 1: Filter the first 264 rows
df_filtered = df.iloc[:264]

# Step 2: Clean 'donor_pre_post' column
df_filtered = df_filtered[df_filtered["donor_pre_post"] != "Pre-Abx/FMT"]
df_filtered["donor_pre_post"] = df_filtered["donor_pre_post"].replace(
    {"Pre-FMT": "PreFMT", "Post-FMT": "PostFMT"}
)

# Step 3: Create a new column combining Disease_type and donor_pre_post
df_filtered["X_Label"] = df_filtered["Disease_type"] + ", " + df_filtered["donor_pre_post"]

# Convert host_removal_input_reads to numeric and divide by a million
df_filtered["host_removal_input_reads"] = pd.to_numeric(df_filtered["host_removal_input_reads"], errors="coerce")
df_filtered["host_removal_input_reads_millions"] = df_filtered["host_removal_input_reads"] / 1e6

# Step 4: Plot the violin plot
plt.figure(figsize=(12, 6))
sns.violinplot(x="X_Label", y="host_removal_input_reads_millions", data=df_filtered, palette="muted", inner=None)
sns.stripplot(x="X_Label", y="host_removal_input_reads_millions", data=df_filtered, color="black", alpha=0.6, jitter=True)

# Enhance plot readability
plt.xticks(rotation=45, ha="right")
plt.xlabel("Disease Type & Donor/Recipient Status")
plt.ylabel("Host Removal Input Reads (Millions)")
plt.title("Violin Plot of Host Removal Input Reads by Disease Type and Donor/Recipient Status")
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show the plot
plt.show()
