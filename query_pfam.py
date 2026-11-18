import requests

# Base URL for Pfam API
base_url = "https://www.ebi.ac.uk/interpro/api/pfam/"

# List of Pfam IDs
pfam_ids = ['PF00067','PF00097','PF00106','PF00107','PF00155','PF00171','PF00173','PF00176','PF00240','PF00248','PF00266','PF00270','PF00271','PF00441','PF00487','PF00636','PF00646','PF00724','PF00743','PF00903','PF01370','PF01565','PF01596','PF02770','PF02771','PF02852','PF03171','PF04030','PF04116','PF04408','PF05050','PF07717','PF07992','PF08031','PF08240','PF08241','PF10607','PF14226','PF21010','PF21791']

# Function to query Pfam API
def query_pfam(pfam_id):
    url = f"{base_url}/{pfam_id}?output=json"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error: Could not fetch data for {pfam_id}. HTTP Status: {response.status_code}")
        return None

# Fetch information for each Pfam ID
for pfam_id in pfam_ids:
    pfam_data = query_pfam(pfam_id)
    if pfam_data:
        print(f"Pfam ID: {pfam_id}")
        print(f"Description: {pfam_data.get('description', 'N/A')}")
        print(f"EC Numbers: {pfam_data.get('ec_numbers', 'N/A')}")
        print("---")
