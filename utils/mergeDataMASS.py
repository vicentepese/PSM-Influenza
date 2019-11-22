import csv
import os 
import json


def main():

    # Read options 
    with open("options.json", 'r') as inFile:
        options = json.load(inFile)
    
    # Merge (colapse) files
    colData = list()
    for file in os.listdir(options['folders']['FiltDataMASS']):
        with open(options['folders']['FiltDataMASS'] + file, 'r') as inFile:
            reader = csv.reader(inFile)
            for row in reader:
                colData.append(row)

    # Write files
    with open(options['files']['mergedFiltMASS'], 'w') as outFile:
        writer = csv.writer(outFile)
        writer.writerows(colData)


if __name__ == "__main__":
    main()