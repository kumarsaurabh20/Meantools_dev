import sqlite3

import pandas as pd

con = sqlite3.connect("jeon_tomato.sqlite")
cur = con.cursor()

sql1="select * from fungal_correlations"
query = cur.execute(sql1)
#s_mnx_names = []
#for each in query:
#    s_mnx_names.append(''.join(list(each)))


#kri = pd.DataFrame(s_mnx_names, columns=["metabolite", "gene", "correlation", "p_value"])
results= pd.DataFrame.from_records(data = query.fetchall(), columns=["metabolite", "gene", "correlation", "p_value"])
print(results)
results.to_csv("./Jeon_tomato/fungal_correlations.csv", index=False)

con.close()
