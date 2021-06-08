from pygbif import occurrences
months = "7,8"
gbif_id = "2496287"
years = "2010, 2012"
records = occurrences.search(gbif_id, year=years, month=months)
months = []
years = []
for x in records["results"]:
    months.append(x["month"])
    years.append(x["year"])
print(set(months))
print(set(years))
