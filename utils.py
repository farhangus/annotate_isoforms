# utils.py


def detectDelimiter(csvFile):
    with open(csvFile, "r") as myCsvfile:
        header = myCsvfile.readline()
        if header.find(",") != -1:
            return ","
        elif header.find("\t") != -1:
            return "\t"
        elif header.find(" ") != -1:
            return " "
        else:
            return None
