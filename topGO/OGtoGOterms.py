import pandas as pd

i = 0
m = 0
df = pd.read_csv("ANN.txt")
df["Orthogroups GO terms"] = df["Orthogroups GO terms"].str.split()


print(
#print(df["Orthogroups GO terms"][1])


#for line in df.itertuples():
#	print(line.str.strip())
	#df = df.str.split("\t")
	#print(df)
#print(len(df)-1)
"""
hako = []
tyuukei = pd.DataFrame([])
wanted = pd.DataFrame([], index = ["1","2","3","4","5","6","7"])	#OGは何個ある？


for i in range(len(df)-1):
	hako.append(df[i][2])
	if df[i][1] != df[i+1][1]:
#		hako = pd.Series(hako)
		wanted.loc[m] = hako
		hako.clear()
		print(wanted)
		m += 1
#		tyuukei = pd.concat([tyuukei,hako.T])
#		hako = hako.values.tolist()
#		hako.clear()
	i += 1
hako.append(df[i][2])

hako = pd.Series(hako)
#tyuukei = pd.concat([tyuukei,hako.T])


#print(tyuukei)


#print(type(tyuukei))
i = 0
for line in tyuukei.itertuples():
	if line[0]
#	i += 1		
"""
