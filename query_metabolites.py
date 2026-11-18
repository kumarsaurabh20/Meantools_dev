import sqlite3
import pandas as pd

metabolites_df = pd.read_csv("jeon_metabolites.csv")

# Extract metabolites from the CSV
metabolites = metabolites_df['metabolites'].iloc[0]  # Assuming a single row in CSV
metabolite_list = metabolites.split(',')

# Connect to the SQLite database
conn = sqlite3.connect("Jeon_tomato/Jeon_bac.sqlite")
cursor = conn.cursor()

query = f"""
SELECT metabolite, gene, correlation, P
FROM Jeon_bac_correlations
WHERE gene = 'Solyc12g100250.1' AND metabolite IN ({', '.join('?' for _ in metabolite_list)})
"""
# Execute the query
cursor.execute(query, metabolite_list)
results = cursor.fetchall()

# Convert the results to a DataFrame for better readability
columns = ['metabolite', 'gene', 'correlation', 'P']
results_df = pd.DataFrame(results, columns=columns)

# Display the results
print(results_df)

# Close the database connection
conn.close()