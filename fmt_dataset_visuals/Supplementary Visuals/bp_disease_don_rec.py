import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset.csv"
df = pd.read_csv(file_path)
df_filtered = df.iloc[:264]

# Clean donor_pre_post column
df_filtered = df_filtered[df_filtered["donor_pre_post"] != "Pre-Abx/FMT"]
df_filtered["donor_pre_post"] = df_filtered["donor_pre_post"].replace(
    {"Pre-FMT": "PreFMT", "Post-FMT": "PostFMT"}
)

# Create a new column combining Disease_type and donor_pre_post
df_filtered["X_Label"] = df_filtered["Disease_type"] + ", " + df_filtered["donor_pre_post"]

# Plot the violin plot
plt.figure(figsize=(12, 6))
sns.violinplot(x="X_Label", y="number_bases_gigabases", data=df_filtered, palette="muted", inner=None, linewidth = 3)
sns.boxplot(x="X_Label", y = "number_bases_gigabases", data=df_filtered, color = "white", width=0.2)
sns.stripplot(x="X_Label", y="number_bases_gigabases", data=df_filtered, color="black", alpha=0.3, jitter=True)

# plot readability
plt.xticks(rotation=45, ha="right")
plt.xlabel("Disease Type & Donor/Recipient Status")
plt.ylabel("Number Bases (Gigabases)")
plt.title("Violin Plot of Number Bases by Disease Type and Donor/Recipient Status")
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.show()
