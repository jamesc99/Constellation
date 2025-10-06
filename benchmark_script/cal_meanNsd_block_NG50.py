import pandas as pd
import sys

def analyze_block_lengths(file_path):
    output_filename = "meanNmedian_sd.stat.txt"
    try:
        with open(file_path, 'r') as f:
            header_line = f.readline().strip()
        column_names = header_line.lstrip('#').split()
        df = pd.read_csv(file_path, sep=r'\s+', skiprows=1, names=column_names)
        df['length'] = (df['to'] - df['from']) + 1
        block_lengths = df['length']
        mean_length = block_lengths.mean()
        std_dev_length = block_lengths.std()
        median_length = block_lengths.median()
        output_content = (
            f"Descriptive statistics for individual phase block lengths from '{file_path}':\n\n"
            f"Total number of blocks: {len(block_lengths):,}\n"
            f"Mean block length: {mean_length:,.2f} bp\n"
            f"Median block length: {median_length:,.0f} bp\n"
            f"Standard Deviation of block length: {std_dev_length:,.2f} bp\n"
            f"Longest block: {block_lengths.max():,} bp\n"
            f"Shortest block: {block_lengths.min():,} bp\n"
        )
        with open(output_filename, 'w') as f_out:
            f_out.write(output_content)
        print(f"Successfully wrote statistics to '{output_filename}'.")
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except KeyError as e:
        print(f"Error: A required column was not found: {e}.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <path_to_block_list_file>")
        sys.exit(1)
    analyze_block_lengths(sys.argv[1])

