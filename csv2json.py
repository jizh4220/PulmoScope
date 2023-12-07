import csv
import json
import glob
import os
 
# Function to convert a CSV to JSON
# Takes the file paths as arguments
def make_json(csvFilePath, jsonFilePath):
     
    data = []
     
    # Open a csv reader called DictReader
    with open(csvFilePath, encoding='utf-8') as csvf:
        csvReader = csv.DictReader(csvf)
        # head_row = next(csvReader)
        # Convert each row into a dictionary
        # and add it to data
        for rows in csvReader:             
            rows.pop('')
            data.append(rows)
 
    # Open a json writer, and use the json.dumps()
    # function to dump data
    with open(jsonFilePath, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))


def bulk_csv2json(path):
    all_csv = glob.glob(path, recursive=True)
    for i in all_csv:
        csvFilePath = i
        fp = os.path.split(csvFilePath)
        fname = fp[1].split('.', 1)[0]
        jsonFilePath = os.path.join(fp[0], (fname + '.json'))
        jsonFilePath

        make_json(csvFilePath, jsonFilePath)

