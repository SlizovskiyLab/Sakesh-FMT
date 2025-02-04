import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Load the dataset
file_path = "C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_dataset.csv"
df = pd.read_csv(file_path)

# Replace values for clarity
df['donor_pre_post'] = df['donor_pre_post'].replace({'Pre-FMT': 'PreFMT', 'Post-FMT': 'PostFMT'})

# Define the sequencing depth using the number of bases in gigabases
# Filter out rows with missing values in these columns
df_filtered = df.dropna(subset=["number_bases_gigabases", "donor_pre_post"])

# Filter only preFMT and postFMT samples
df_filtered = df_filtered[df_filtered["donor_pre_post"].isin(["PreFMT", "PostFMT"])]

# Convert timepoint (days) to numeric
df_filtered["timepoint"] = pd.to_numeric(df_filtered["timepoint"], errors="coerce")

# Sort values by time for better visualization
df_filtered = df_filtered.sort_values("timepoint")

# Initialize plot
plt.figure(figsize=(10, 6))

# Define colors for each group
colors = {"PreFMT": "blue", "PostFMT": "red"}

# Plot scatter points
for condition, color in colors.items():
    subset = df_filtered[df_filtered["donor_pre_post"] == condition]
    plt.scatter(subset["timepoint"], subset["number_bases_gigabases"], 
                color=color, alpha=0.7, label=f"{condition} Samples")

# Approximate LOESS manually using a rolling average
window_size = 5  # Adjust for smoothness
for condition, color in colors.items():
    subset = df_filtered[df_filtered["donor_pre_post"] == condition].copy()  # Avoid warnings
    subset["Smoothed"] = subset["number_bases_gigabases"].rolling(window=window_size, center=True).mean()
    
    # Plot LOESS approximation
    plt.plot(subset["timepoint"], subset["Smoothed"], color=color, linewidth=2, linestyle="-", 
             label=f"{condition} LOESS Approx.")

# Customize x-axis ticks
plt.xticks(rotation=45)  # Rotate for readability
plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))  # Ensure enough ticks

plt.xlim(0, 400)  
plt.ylim(0, 35)   
plt.xticks(range(0, 401, 50))  
plt.yticks(range(0, 36, 5))

# Customize plot
plt.xlabel("Days Until Collection (timepoint)")
plt.ylabel("Sequencing Depth (Gigabases)")
plt.title("Sequencing Depth Over Time for preFMT and postFMT Samples")
plt.legend()
plt.grid(alpha=0.3)

# Show plot
plt.show()
