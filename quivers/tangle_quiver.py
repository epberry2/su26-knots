from helpers.findloops import findpathwithloops_sector

def writhe_of_path(num, denom, j, i):
    "assume i < j and num >= denom"
    if i == j: return 0
    path, loops = findpathwithloops_sector(num, denom, j, i)
    i_left = i <= num
    j_left = j <= num
    writhe = 0
    on_left = j_left
    

    if i_left and j_left:
        if loops[0] == "L":
            writhe += -1 if j > i else 0
        elif loops[0] == "C":
            writhe += 0 if j > i else -1
        elif loops[0] == "T":
            writhe += 0 if j > i else -1
    elif (not i_left) and (not j_left):
        if loops[0] == "R":
            writhe += -1 if i > j else 0
        elif loops[0] == "T":
            writhe += 0 if i > j else -1
            
    for next_point, path_type in zip(path[1:-1], loops[0:-1]):
        if i_left:
            match path_type:
                case "L":
                    writhe += 1 if next_point > i else -1
                case "C":
                    writhe += 1 if next_point < i else -1
                case "T":
                    if not on_left:
                        writhe += 1 if next_point < i else -1
                    on_left = not on_left
        else:
            match path_type:
                case "R":
                    writhe += 1 if next_point < i else -1
                case "T":
                    if on_left:
                        writhe += -1 if next_point < i else 1
                    on_left = not on_left
    return writhe

print(writhe_of_path(5, 2, 5, 6))
            