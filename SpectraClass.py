import csv
import scipy.stats as sp

def combineWindows(windows):
    if len(windows) > 1:
        combinedWindow = Window(library=windows[0].library, sequence=windows[0].sequence, coordinates=(windows[0].coords[0], windows[-1].coords[1]))
        countNames = list(windows[0].counts.keys())
        combinedWindow.setCounts(countNames, [sum([window.counts[name] for window in windows]) for name in countNames])
        combinedWindow.setFrequenciesFromCounts()
        return combinedWindow
    else:
        return windows

def writeToFile(spectra, outputFilePath):
    with open(outputFilePath, 'w', newline='') as outputFile:
        csvWriter = csv.writer(outputFile, delimiter='\t')
        csvWriter.writerow(["Library", "Sequence", "Start", "End"] + list(spectra.windows[0].counts.keys()))
        for window in spectra.windows:
            if isinstance(window, list):
                print(window)
                csvWriter.writerow(windowToList(window[0]))
            else:
                csvWriter.writerow(windowToList(window))

def windowToList(window):
    return [window.library, window.sequence, window.coords[0], window.coords[1]] + [window.counts[a] for a in window.counts]

def frequencyFromValues(values):
    total = sum(values)
    if total == 0:
        return [0] * len(values)
    else:
        return [a/total for a in values]

class Spectra:
    def __init__(self, spectraTsv, frequencies=False):
        self.spectraTsv = spectraTsv
        self.isFrequencies = frequencies
        self.windows = []
        self.globalFrequencies = {}
        self.globalCounts = {}
        self.mers = []
        self.loadSequences()

    def loadSequences(self):
        with open(self.spectraTsv, 'r') as fileInput:
            tsvReader = csv.reader(fileInput, delimiter='\t')
            headers = next(tsvReader)
            self.mers = headers[-64:]

            for spectraWindow in tsvReader:
                newWindow = Window(spectraWindow[0], spectraWindow[1], (int(spectraWindow[2]), int(spectraWindow[3])))
                if self.isFrequencies:
                    newWindow.setFrequencies(self.mers, [float(a) for a in spectraWindow[-64:]])
                    newWindow.setCountsFromFrequencies()
                else:
                    newWindow.setCounts(self.mers, [int(a) for a in spectraWindow[-64:]])
                    newWindow.setFrequenciesFromCounts()
                self.windows.append(newWindow)

    def getMerSubset(self, library=None, sequence=None, coordinates=0):
        if not sequence:
            sequence = []
        if not library:
            library = []
        toReturn = [a for a in self.windows]
        if library:
            toReturn = [a for a in toReturn if a.library in library]
        if sequence:
            toReturn = [a for a in toReturn if a.sequence in sequence]
        if isinstance(coordinates, tuple) and len(coordinates) > 1:
            toReturn = [a for a in toReturn if a.coords[0] >= coordinates[0] and a.coords[1] <= coordinates[1]]
        return toReturn

    def getLibraryAndSequenceList(self):
        librariesDict = {}
        for window in self.windows:
            if window.library not in librariesDict:
                librariesDict[window.library] = [window.sequence]
            elif window.sequence not in librariesDict[window.library]:
                librariesDict[window.library].append(window.sequence)
        return librariesDict

    def setGlobalCounts(self):
        for mer in self.mers:
            self.globalCounts[mer] = sum([a.counts[mer] for a in self.windows])
        pass

    def getGlobalFrequencies(self):
        return self.globalFrequencies

    def setGlobalFrequencies(self):
        totalLength = sum([a.coords[1] - a.coords[0] - 1 for a in self.windows])
        for mer in self.mers:
            self.globalFrequencies[mer] = sum([a.frequencies[mer] * (a.coords[1] - a.coords[0] - 1) for a in self.windows]) / totalLength

    def resizeWindows(self, targetFactor=1):
        targetFactor = int(targetFactor)
        replacementSpectraWindows = []
        librarySequences = self.getLibraryAndSequenceList()
        for library in librarySequences:
            for sequence in librarySequences[library]:
                windowSubset = self.getMerSubset(library=library, sequence=sequence)
                for index in range(0, len(windowSubset), targetFactor):
                    replacementSpectraWindows.append(combineWindows(windowSubset[index:index+targetFactor]))
        self.windows = replacementSpectraWindows

    def reduceByGlobalFrequencies(self):
        if not self.globalFrequencies:
            self.setGlobalFrequencies()
        for window in self.windows:
            window.modifyFrequencies([-self.globalFrequencies[a] for a in self.globalFrequencies])

    def removeOutliers(self):
        if not self.globalFrequencies:
            self.setGlobalFrequencies()
        globalFreq = [self.globalFrequencies[a] for a in self.globalFrequencies]
        if sum(globalFreq) < 1:
            globalFreq = frequencyFromValues(globalFreq)
        toRemove = []
        for i in range(len(self.windows)):

            window = self.windows[i]
            windowFreq = [window.frequencies[a] for a in window.frequencies]
            if sum(windowFreq) == 0:
                toRemove.append(i)
            elif sum(windowFreq) < 1:
                windowFreq = frequencyFromValues(windowFreq)
                chiResult = sp.chisquare(windowFreq, f_exp=globalFreq)
                if chiResult[1] >= .999:
                    toRemove.append(i)
            else:
                chiResult = sp.chisquare(windowFreq, f_exp=globalFreq)
                if chiResult[1] >= .999:
                    toRemove.append(i)
        for remove in toRemove[::-1]:
            self.windows.pop(remove)

class Window:
    def __init__(self, library, sequence, coordinates):
        self.counts = {}
        self.coordinates = coordinates
        self.library = library
        self.sequence = sequence
        self.frequencies = {}

    def modifyFrequencies(self, modifiers):
        mers = list(self.counts.keys())
        for i in range(len(modifiers)):
            self.frequencies[mers[i]] = self.frequencies[mers[i]] + modifiers[i] if self.frequencies[mers[i]] + modifiers[i] > 0 else 0
        self.setCountsFromFrequencies()

    def setFrequenciesFromCounts(self):
        maxCount = self.coordinates[1] - self.coordinates[0] - 1
        for count in self.counts:
            self.frequencies[count] = self.counts[count] / maxCount

    def setCountsFromFrequencies(self):
        maxCount = self.coordinates[1] - self.coordinates[0] - 1
        for frequency in self.frequencies:
            self.counts[frequency] = int(self.frequencies[frequency] * maxCount)

    def setCounts(self, headers, counts):
        for i in range(len(counts)):
            self.counts[headers[i]] = counts[i]

    def setFrequencies(self, headers, frequencies):
        for i in range(len(frequencies)):
            self.frequencies[headers[i]] = frequencies[i]

    def setLibrary(self, library):
        self.library = library

    def setSequence(self, sequence):
        self.sequence = sequence

    def setCoordinates(self, coordinates):
        self.coordinates = (coordinates[0], coordinates[1])
