###################################################  INSERT INTO DB (big)
########################################################################
biggin = datetime.now()
'''  # This is an alternate way to insert records
sql1 = """INSERT INTO occurrences ('occ_id', 'taxon_id', 'source',
                                   'latitude', 'longitude',
                                   'coordinateUncertaintyInMeters',
                                   'occurrenceDate', 'request_id',
                                   'filter_id', 'generalizations',
                                   'remarks')
          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"""
for x in df9.index:
    insert2 = [df9.loc[x,"id"], taxon_id, df9.loc[x,"source"],
               df9.loc[x,"decimalLatitude"], df9.loc[x,"decimalLongitude"],
               df9.loc[x,"coordinateUncertaintyInMeters"],
               df9.loc[x,"eventDate"], request_id, filter_id,
               df9.loc[x,"dataGeneralizations"], df9.loc[x,"remarks"]]
    cursor.execute(sql1, [(insert2)])
conn.commit()
'''
df9.to_sql(name='occurrences', con = conn, if_exists='replace',
           chunksize=2000)
sql_toad = '''SELECT AddGeometryColumn('occurrences', 'geom_xy4326', 4326,
                                       'POINT', 'XY');'''
cursor.execute(sql_toad)
print("Inserted records into table: " + str(datetime.now() - biggin))
