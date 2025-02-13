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

# Step 4: Plot the violin plot
plt.figure(figsize=(12, 6))
sns.violinplot(x="X_Label", y="number_bases_gigabases", data=df_filtered, palette="muted", inner=None, cut = 0)
sns.boxplot(x="X_Label", y = "number_bases_gigabases", data=df_filtered, color = "gray", width=0.2)
sns.stripplot(x="X_Label", y="number_bases_gigabases", data=df_filtered, color="black", alpha=0.6, jitter=True)

# Enhance plot readability
plt.xticks(rotation=45, ha="right")
plt.xlabel("Disease Type & Donor/Recipient Status")
plt.ylabel("Number Bases (Gigabases)")
plt.title("Violin Plot of Number Bases by Disease Type and Donor/Recipient Status")
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show the plot
plt.show()
