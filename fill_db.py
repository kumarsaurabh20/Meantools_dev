import sqlite3
import pandas as pd

# Connect to SQLite database
conn = sqlite3.connect('/Users/singh018/Documents/Meantools_v1/PAL.sqlite')
cursor = conn.cursor()

# Read data from CSV file
csv_file = '/Users/singh018/Documents/Meantools_v1/PAL_related_metabolite_lotus.tsv'
df = pd.read_csv(csv_file, delimiter=',', quotechar="'")
print(df)

data_to_insert = df[[
    'lotus_id', 'wikidata_id', 'smiles', 'inchi', 'inchikey', 'iupac_name',
    'molecular_formula', 'molecular_weight', 'NPclassifierPathway', 'NPclassifierSuperclass',
    'NPclassifierClass', 'ClassyfireKingdom', 'ClassyfireSuperclass', 'ClassyfireClass',
    'ClassyfireDirectParent', 'allTaxa'
]]

# Insert data into the SQLite database
data_to_insert.to_sql('Jeon_dummy_PAL', conn, if_exists='replace', index=False)

# Commit the changes and close the connection
conn.commit()
conn.close()

