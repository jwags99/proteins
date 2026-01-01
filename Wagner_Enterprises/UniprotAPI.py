import requests, sys, json

# Documentation: https://www.uniprot.org/help/api
WEBSITE_API = "https://rest.uniprot.org/"

# Documentation: https://www.ebi.ac.uk/proteins/api/doc/
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

# Helper function to download data
def get_url(url, **kwargs):
  response = requests.get(url, **kwargs)

  if not response.ok:
    print(response.text)
    response.raise_for_status()
    sys.exit()

  return response

r = get_url(f"{WEBSITE_API}/uniprotkb/search?query=*")

data = r.json()

# print the number of results in the payload
n_results = len(data["results"])
print(f"Number of results: {n_results}\n")

# print all headers in the server response
for (key, value) in r.headers.items():
  print(f"{key}: {value}")