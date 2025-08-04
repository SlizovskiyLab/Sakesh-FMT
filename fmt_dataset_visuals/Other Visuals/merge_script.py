import pandas as pd


# The main file to add the column to
df_main = pd.read_csv('C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_paired.csv')

# The lookup file containing the data
df_lookup = pd.read_csv('C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\seqkit_stats_results.csv')

# we'll keep only the first entry for each 'run_accession'.
df_lookup_unique = df_lookup.drop_duplicates(subset=['run_accession'], keep='first')

throughput_map = df_lookup_unique.set_index('run_accession')['total_sum_len']

# to the 'sum_len' values in our map.
df_main['Throughput(base_level,paired_end)'] = df_main['run_accession'].map(throughput_map)

# --- 4. Save the Result ---
output_filename = 'C:\\Users\\asake\\OneDrive\\Desktop\\Homework\\FMT\\FMT_full_dataset_with_throughput.csv'
df_main.to_csv(output_filename, index=False)

print(f"The new file is saved as: {output_filename}")

