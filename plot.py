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
    if dataFiles[i] == "temp.log":
        continue

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
        infoObj = {
                "filename": dataFiles[i],
                "numGlobal": int(lineSplit[3]),
                "numTotal": int(lineSplit[5]),
                "numLocal": int(lineSplit[5]) - int(lineSplit[3]),
                "numChains": int(lineSplit[6]),
                "numFreeChains": int(lineSplit[7]),
                "numFixedChains": int(lineSplit[6]) - int(lineSplit[7]),
                "numFreeVars": int(lineSplit[8]),
                "averageObjective": float(lineSplit[9]),
                "averageObjectiveLog": float(lineSplit[10]),
                "dimension": int(int(lineSplit[8]) / int(lineSplit[7])),
                "setSizes": [int(a) for a in dataFiles[i][:-4].split("-")],
                "knownInfeasible": False,
            }
        if len(infoObj["setSizes"]) > infoObj["dimension"]+1:
            infoObj["knownInfeasible"] = True
        for size in infoObj["setSizes"]:
            if size > infoObj["dimension"]:
                infoObj["knownInfeasible"] = True
        info.append(infoObj)
        if infoObj["dimension"] not in dimsUsed and not infoObj["knownInfeasible"]:
            dimsUsed.append(infoObj["dimension"])

import plotly.express as px
import plotly.graph_objs as go
import numpy as np

x = [d["numChains"] for d in info]
y = [d["averageObjective"] for d in info]
labels = [d["filename"] for d in info]

# Creating the dataset, and generating the plot
trace1 = go.Scatter(
                  x=x,
                  y=y,
                  mode='markers',
                  hovertext=labels
                  )

gObjs = [trace1]
for dim in dimsUsed:
    xNoInfeasible = [d["numChains"] for d in info if not d["knownInfeasible"] and d["dimension"] == dim]
    yNoInfeasible = [d["averageObjective"] for d in info if not d["knownInfeasible"] and d["dimension"] == dim]
    z = np.polyfit(xNoInfeasible, yNoInfeasible, 2)
    f = np.poly1d(z)
    print("for d=", dim)
    print(f)

    # calculate new x's and y's
    x_new = np.linspace(min(xNoInfeasible), max(xNoInfeasible), 50)
    y_new = f(x_new)

    trace2 = go.Scatter(
                      x=x_new,
                      y=y_new,
                      mode='lines',
                      name='fit for ' + str(dim)
                      )
    gObjs.append(trace2)

layout = go.Layout(
                )

fig = go.Figure(data=gObjs, layout=layout)
fig.show()

# fig = px.scatter(x=xData, y=yData, hover_name=labels)
# fig.update_traces(textposition='top left')
# fig.update_layout(
    # xaxis_title="Number of chains",
    # yaxis_title="Average Objective",
# )
# fig.show()


