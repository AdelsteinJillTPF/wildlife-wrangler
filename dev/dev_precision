import pandas as pd

df1 = pd.read_csv("T:/Temp/dev_acc_prec.csv")
df2 = df1[["record_id", "eventDate", "gps_accuracy_m", "source"]]
print(df2.head())


# Set GPS accuracy to 100 m if record is pre-2000.
df2["gps_accuracy_m"] = np.where(df2["eventDate"].apply(lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S").year) < 2000, 100, df2["gps_accuracy_m"])
