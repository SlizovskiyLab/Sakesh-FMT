import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
df = pd.read_csv(file_path)

# Group by Disease_type and sum the host reads
df_grouped = df.groupby("Disease_type")["trimmomatic_input_reads"].sum() / 1e6  # Convert to millions

# Plot the bar chart
plt.figure(figsize=(10, 6))
plt.bar(df_grouped.index, df_grouped, color="skyblue", edgecolor="black")

# Customize the plot
plt.xlabel("Disease Type")
plt.ylabel("Host Reads (Millions)")
plt.title("Number of Host Reads by Disease Type")
plt.xticks(rotation=45)
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show the plot
plt.show()
