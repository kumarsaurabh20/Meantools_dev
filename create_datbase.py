import pandas as pd
import sqlite3

# Load the CSV file into a Pandas DataFrame
csv_file = 'Linoleic_acid_pathway_lotus.csv'  # Replace with your CSV file path
df = pd.read_csv(csv_file)

# Connect to SQLite database (creates a new database if it doesn't exist)
database_file = 'roots_7dpi_mgh_met_up.sqlite'  # Replace with your database file path
conn = sqlite3.connect(database_file)

# Write the DataFrame to SQLite database
table_name = 'linoleic_acid_pathway'  # Replace with your desired table name
df.to_sql(table_name, conn, if_exists='replace', index=False)

# Close the database connection
conn.close()

print(f"Data from {csv_file} has been stored in {database_file} under the table '{table_name}'.")
