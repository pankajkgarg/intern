import json

f = open("geneCoord.txt", "rb")
rawCoord = f.read().split("\n")
f.close()

f = open("tempLabels.txt", "rb")
rawLabels = f.read().split("\n")
f.close()

f = open("tempScores.txt", "rb")
rawScores = f.read().split("\n")
f.close()

coordinates = []
for i, row in enumerate(rawCoord):
	temp = row.split(",")
	if len(temp) != 2:
		continue
	x,y = temp[0].strip(), temp[1].strip()
	coordinates.append({
		"i": i,
		"x": x,
		"y": y,
		"label": rawLabels[i].strip(),
		"score": rawScores[i].strip()
	})
	

f = open("jsonOutput2.js", "wb")
text = "var data = " + json.dumps(coordinates)
f.write(text)
f.close()
