import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
df = pd.read_csv(file_path)
df = df[:-5]
# Replace 'Pre-FMT' with 'PreFMT' and 'Post-FMT' with 'PostFMT'
df['donor_pre_post'] = df['donor_pre_post'].replace({'Pre-FMT': 'PreFMT', 'Post-FMT': 'PostFMT'})

# Group by donor_pre_post and sum the host reads
df_grouped = df.groupby("donor_pre_post")["trimmomatic_input_reads"].sum() / 1e6  # Convert to millions

# Plot the bar chart
plt.figure(figsize=(10, 6))
plt.bar(df_grouped.index, df_grouped, color="skyblue", edgecolor="black")

# Customize the plot
plt.xlabel("Donor/PreFMT/PostFMT")
plt.ylabel("Host Reads (Millions)")
plt.title("Number of Host Reads by Donor/Recipient Types")
plt.xticks(rotation=45)
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show the plot
plt.show()
