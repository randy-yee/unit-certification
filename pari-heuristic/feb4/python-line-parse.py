FILENAME = "Q-experiment-0-4"
with open(FILENAME,"r") as file:
    with open("results"+FILENAME, "w") as outfile:
        for line in file:
            if(line[0:3] == "Max"):
                outfile.write(line)
