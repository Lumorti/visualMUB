#!/usr/env/bin python3

# Imports
import matplotlib.pyplot as plt
import os
import plotly.express as px
import plotly.graph_objs as go
import numpy as np

# Load all files in the data directory
dataDir = "./data/"
dataFiles = os.listdir(dataDir)
info = []
dimsUsed = []
for i in range(len(dataFiles)):
    if "temp" in dataFiles[i]:
        continue

    # Open the file 
    with open(dataDir + dataFiles[i], "r") as f:

        # Parse the file
        data = f.read()
        data = data.replace("\r", "#DELETE\n")
        data = data.split("\n")

        # Get the last non-empty line
        lastLine = ""
        for j in range(1, len(data)):
            if len(data[-j]) > 0 and "DELETE" not in data[-j]:
                lastLine = data[-j]
                break

        # Get each unique run
        finalValues = []
        for j in range(len(data)):
            if "inf" in data[j]:
                if "sqr=" in data[j] and "DELETE" not in data[j] and "nan" not in data[j]:
                    finalValues.append(float(data[j].split(" ")[1]))
            else:
                if "s=" in data[j] and "DELETE" not in data[j] and "nan" not in data[j]:
                    finalValues.append(float(data[j].split(" ")[0][2:]))

        # Calculate statistics
        avg = sum(finalValues) / len(finalValues)
        std = (sum([(a-avg)**2 for a in finalValues]) / len(finalValues))**0.5
        err = std / (len(finalValues)**0.5)
        maxVal = max(finalValues)
        minVal = min(finalValues)

        # Parse the line
        lastLine = lastLine.strip().replace("  ", " ")
        lineSplit = lastLine.split(" ")
        if "inf" in lastLine:
            numLocal = sum([1 for a in finalValues if abs(minVal-a) < 1e-5])
            infoObj = {
                    "filename": dataFiles[i],
                    "numTotal": len(finalValues),
                    "numGlobal": numLocal,
                    "numLocal": numLocal,
                    "numChains": int(lineSplit[8]),
                    "numFreeChains": int(lineSplit[9]),
                    "numFixedChains": int(lineSplit[8]) - int(lineSplit[9]),
                    "numFreeVars": int(lineSplit[10]),
                    "averageObjective": avg,
                    "standardDeviation": std,
                    "error": err,
                    "max": maxVal,
                    "min": minVal,
                    "dimension": int(int(lineSplit[10]) / int(lineSplit[9])),
                    "setSizes": [int(a) for a in dataFiles[i][:-4].split("-")],
                    "knownInfeasible": False,
                }
        else:
            infoObj = {
                    "filename": dataFiles[i],
                    "numGlobal": int(lineSplit[3]),
                    "numTotal": int(lineSplit[5]),
                    "numLocal": int(lineSplit[5]) - int(lineSplit[3]),
                    "numChains": int(lineSplit[6]),
                    "numFreeChains": int(lineSplit[7]),
                    "numFixedChains": int(lineSplit[6]) - int(lineSplit[7]),
                    "numFreeVars": int(lineSplit[8]),
                    "averageObjective": avg,
                    "standardDeviation": std,
                    "error": err,
                    "max": maxVal,
                    "min": minVal,
                    "dimension": int(int(lineSplit[8]) / int(lineSplit[7])),
                    "setSizes": [int(a) for a in dataFiles[i][:-4].split("-")],
                    "knownInfeasible": False,
                }
        if len(infoObj["setSizes"]) > infoObj["dimension"]+1:
            infoObj["knownInfeasible"] = True
        for size in infoObj["setSizes"]:
            if size > infoObj["dimension"]:
                infoObj["knownInfeasible"] = True
        if minVal > 1e-3:
            infoObj["knownInfeasible"] = True
        info.append(infoObj)
        if infoObj["dimension"] not in dimsUsed and not infoObj["knownInfeasible"]:
            dimsUsed.append(infoObj["dimension"])


# Creating the dataset, and generating the plot
gObjs = []
for dim in dimsUsed:
    xAll = [d["numChains"] for d in info if d["dimension"] == dim]
    yAll = [d["averageObjective"] for d in info if d["dimension"] == dim]
    xInfeasible = [d["numChains"] for d in info if d["knownInfeasible"] and d["dimension"] == dim]
    yInfeasible = [d["averageObjective"] for d in info if d["knownInfeasible"] and d["dimension"] == dim]
    xFeasible = [d["numChains"] for d in info if not d["knownInfeasible"] and d["dimension"] == dim]
    yFeasible = [d["averageObjective"] for d in info if not d["knownInfeasible"] and d["dimension"] == dim]
    labelsFeasible = [d["filename"] for d in info if not d["knownInfeasible"] and d["dimension"] == dim]
    labelsInfeasible = [d["filename"] for d in info if d["knownInfeasible"] and d["dimension"] == dim]
    deltaPlusFeasible = [abs(d["max"]-d["averageObjective"]) for d in info if not d["knownInfeasible"] and d["dimension"] == dim]
    deltaPlusInfeasible = [abs(d["max"]-d["averageObjective"]) for d in info if d["knownInfeasible"] and d["dimension"] == dim]
    deltaMinusFeasible = [abs(d["min"]-d["averageObjective"]) for d in info if not d["knownInfeasible"] and d["dimension"] == dim]
    deltaMinusInfeasible = [abs(d["min"]-d["averageObjective"]) for d in info if d["knownInfeasible"] and d["dimension"] == dim]
    trace1 = go.Scatter(
                      x=xFeasible,
                      y=yFeasible,
                      error_y={"array": deltaPlusFeasible, "arrayminus": deltaMinusFeasible, "symmetric": False},
                      mode='markers',
                      hovertext=labelsFeasible,
                      fillcolor="blue",
                      name='feasible dim=' + str(dim)
                      )
    trace2 = go.Scatter(
                      x=xInfeasible,
                      y=yInfeasible,
                      mode='markers',
                      hovertext=labelsInfeasible,
                      error_y={"array": deltaPlusInfeasible, "arrayminus": deltaMinusInfeasible, "symmetric": False},
                      fillcolor="red",
                      name='infeasible dim=' + str(dim)
                      )
    z = np.polyfit(xAll, yAll, 2)
    f = np.poly1d(z)
    x_new = np.linspace(min(xAll), max(xAll), 50)
    y_new = f(x_new)
    trace3 = go.Scatter(
                      x=x_new,
                      y=y_new,
                      mode='lines',
                      fillcolor="green",
                      name='fit for ' + str(dim)
                      )
    gObjs.append(trace1)
    gObjs.append(trace2)
    gObjs.append(trace3)

layout = go.Layout(
                )

fig = go.Figure(data=gObjs, layout=layout)
fig.show()


