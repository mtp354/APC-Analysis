# APC-Analysis
Some Template Code for analyzing APC effects.

# Rough User Guide
Run a SEERstat Rate session such as the included example. Copy session matrix to excel and save as .csv file.
Run R code, selecting the file path that points to the new SEER data.
Export plot and Mutual Information Matrix via clipbard to excel.

# Understanding Mutual Information Matrix
The elements in the matrix correspond amount of information shared between the row variable and the column variable. Information here is a measure of the amount of variance reduced in one Random Variable by observations of another Random Variable. Diagonal entries should therefore be the largest elements in the array.
