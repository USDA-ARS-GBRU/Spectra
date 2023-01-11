#!/usr/bin/env python3

# Take a given image file and define a spectra-compatible color palette from evenly spaced pixels in the image

from PIL import Image
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Load matrix of pixel colors from image')
parser.add_argument('-i', '--input', dest='input_file', type=str, help='Input image file', required=True)
parser.add_argument('-o', '--output', dest='output_file', type=str, help='Output csv file', required=True)
parser.add_argument('-x', '--width', dest='x', type=int, help='Sample X times across width of image', default=8)
parser.add_argument('-y', '--height', dest='y', type=int, help='Sample Y times across height of image', default=8)
args = parser.parse_args()

with Image.open(args.input_file) as image:
    imagePixels = image.load()

    sampleSpacingX = image.size[0] / args.x
    sampleSpacingY = image.size[1] / args.y

    outputMatrix = []
    # sample the image at X intervals of Y rows
    for y in range(args.y):
        currentRow = []
        for x in range(args.x):
            currentSamplePixel = (
                round((x * sampleSpacingX) + sampleSpacingX / 2),
                round((y * sampleSpacingY) + sampleSpacingY / 2)
            )
            try:
                r, g, b, a = imagePixels[currentSamplePixel[0], currentSamplePixel[1]]
                currentRow.append(f"#{r:02x}{g:02x}{b:02x}")
            except ValueError:
                try:
                    r, g, b = imagePixels[currentSamplePixel[0], currentSamplePixel[1]]
                    currentRow.append(f"#{r:02x}{g:02x}{b:02x}")
                except TypeError:
                    currentRow.append(f"#FFFFFF")

        outputMatrix.append(currentRow)
    outputDF = pd.DataFrame(outputMatrix, columns=[f'c{a}' for a in range(len(outputMatrix[1]))])
    outputDF.to_csv(args.output_file, index=False)
