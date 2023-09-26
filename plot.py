#!/usr/env/bin python3

# Imports
import matplotlib.pyplot as plt
import os

# Load all files in the data directory
dataDir = "./data/"
dataFiles = os.listdir(dataDir)
info = []
dimsUsed = []
for i in range(len(dataFiles)):

    # Open the file 
    with open(dataDir + dataFiles[i], "r") as f:

        # Parse the file
        data = f.read()
        data = data.replace("\r", "\n")
        data = data.split("\n")

        # Get the last non-empty line
        lastLine = ""
        for j in range(1, len(data)):
            if len(data[-j]) > 0:
                lastLine = data[-j]
                break

        # Parse the line
        lastLine = lastLine.strip().replace("  ", " ")
        lineSplit = lastLine.split(" ")
        info.append({
                "filename": dataFiles[i],
                "numGlobal": int(lineSplit[3]),
                "numTotal": int(lineSplit[5]),
                "numLocal": int(lineSplit[5]) - int(lineSplit[3]),
                "numChains": int(lineSplit[6]),
                "numFreeChains": int(lineSplit[7]),
                "numFixedChains": int(lineSplit[6]) - int(lineSplit[7]),
                "numFreeVars": int(lineSplit[8]),
                "dimension": int(int(lineSplit[8]) / int(lineSplit[7])),
                "setSizes": dataFiles[i][:-4].split("-")
            })
        if info[-1]["dimension"] not in dimsUsed:
            dimsUsed.append(info[-1]["dimension"])

        # Check for sanity
        print("Info: ", info[-1])

# Plot the data
# pointTypes = {
        # 2: (".", "red"),
        # 3: (".", "blue"),
        # 4: (".", "green"),
        # 5: (".", "orange"),

        # 6: ("x", "red"),
        # 7: ("x", "green"),
        # 8: ("x", "blue"),
        # 9: ("x", "orange"),

        # 10: ("+", "red"),
        # 16: ("+", "green"),
        # 18: ("+", "blue"),

        # 20: ("*", "red"),
        # 50: ("*", "green"),
        # 100: ("*", "blue"),
    # }
# plt.figure(figsize=(10, 6))
# for dim in dimsUsed:
    # xData = [d["numChains"] for d in info if d["dimension"] == dim]
    # yData = [d["numLocal"]  / float(d["numTotal"]) for d in info if d["dimension"] == dim]
    # plt.scatter(xData, yData, marker=pointTypes[dim][0], color=pointTypes[dim][1], label="d="+str(dim))
# plt.legend()
# plt.show()
import plotly.express as px
xData = [d["numChains"] for d in info]
yData = [d["numLocal"]  / float(d["numTotal"]) for d in info]
labels = [d["filename"] for d in info]
fig = px.scatter(x=xData, y=yData, hover_name=labels)
fig.update_traces(textposition='top left')
fig.update_layout(
    xaxis_title="Number of chains",
    yaxis_title="Fraction of local variables",
)
fig.show()


