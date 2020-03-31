from dataSets import dataSets, dataCombos

COLOURS = {'white':0, 'black':1, 'gray':920, 'red':632, 'green':416,
	   'blue':600, 'yellow':400, 'magenta':616, 'cyan':432, 'orange':800,
	   'spring':820, 'teal':840, 'azure':860, 'violet':880, 'pink':900}

def getInput():
    """
    Function to get input list of decay chains from the user and verify them
    Returns an array of keys to analyse and a boolean with their validity
    """

    # get list of decay chains from user
    userInput = input() # get input

    # check whether to read list from file or user input
    if "text" in userInput: 
        # read from file, remove newlines if present
        with open("input.txt","r") as infile:
            userInput = infile.read().replace("\n","")

    src = list(map(lambda x: x.strip(), userInput.split(","))) # get list of series of chains from string
    series, colourcodes = [], []
    for ch in src:
        s,c = ch.split(" ")
        series.append(s)
        c = c.lower()
        cc = COLOURS[c] if c in COLOURS else 632
        colourcodes.append(cc)

    # new 2d list to hold chains
    # first index is series, second is chain in the series
    chains = [[]]*len(series)

    for i in range(len(series)):
       chains[i] = series[i].split("+") 

    # check that all the decay chains are valid
    chainsValid = True

    for i in range(len(chains)):
        for j in range(len(chains[i])):
            if not ((chains[i][j] in dataSets.keys()) or (chains[i][j] in dataCombos.keys())):
                print("Sorry I don't recognise " + chains[i][j] + " as a valid data set.")
                chainsValid = False

    # return the chains and their validity
    return chains, colourcodes, chainsValid 
