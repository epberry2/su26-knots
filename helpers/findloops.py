from collections import defaultdict
from helpers.tanglestate import get_state

def findpath(num, denom):
    tupleslist = []
    for i in range(num // 2):
        tupleslist.append((1 + i, num - i))
    for i in range(denom // 2):
        tupleslist.append((num + 1 + i, num + denom - i))
    if num >= denom:
        for i in range(denom):
            tupleslist.append((num - i, num + denom - i))
        for i in range((num - denom) // 2):
            tupleslist.append((i + 1, num - denom - i))
    else:
        for i in range(num):
            tupleslist.append((num - i, num + denom - i))
        for i in range((denom - num) // 2):
            tupleslist.append((num + 1 + i, denom - i))
    adj = defaultdict(list)
    for a, b in tupleslist:
        adj[a].append(b)
        adj[b].append(a)    
    start = 0
    if num >= denom:
        if (num % 2 == 1):
            start = (num + 1) // 2
        else:
            start = num + ((denom + 1) // 2)
    else:
        if (denom % 2 == 1):
            start = num + ((denom + 1) // 2)
        else:
            start = (num + 1) // 2
    prev = None
    path = [start]
    while len(path) < num + denom:
        cur = path[-1]
        nxt = next(x for x in adj[cur] if x != prev)
        path.append(nxt)
        prev = cur
    points = get_state(num, denom)[1]
    x_min_pos = next(index for index, point in enumerate(points) if point == "X-")
    match x_min_pos:
        case 0:
            if not (path[0] == (num + 1) // 2):
                path.reverse()
        case 1:
            if not (path[0] == (num - denom + 1) // 2):
                path.reverse()
        case 2:
            if not (path[0] == num + (denom + 1) // 2):
                path.reverse()
    return path
    

def findpathwithloops(num, denom):
    path = findpath(num, denom)
    pathwithloops = []
    iter = 1
    end = len(path) - 1
    index = 0
    while path[index] != path[end]:
        nextindex = index + iter
        cur = path[index]
        nxt = path[nextindex]
        index += iter
        if cur > num and nxt > num:
            if num > denom or cur + nxt == denom + (2 * num) + 1:
                pathwithloops.append("R")
                continue
            pathwithloops.append("C")
            continue
        if (cur > num and nxt <= num) or (cur <= num and nxt > num):
            pathwithloops.append("T")
            continue
        if  denom > num or cur + nxt == num + 1:
            pathwithloops.append("L")
            continue
        pathwithloops.append("C")
    return path, pathwithloops


def findpathwithloops_sector(num, denom, start, end):
    path = findpath(num, denom)
    index = path.index(start)
    end_idx = path.index(end)
    startbeforeend = index < end_idx
    loops = []
    path_block = [path[index]]
    if startbeforeend:
        iter = 1
    else:
        iter = -1
    print(f"end: {end}")
    while path[index] != end:
        nextindex = index + iter
        cur = path[index]
        print(nextindex)
        nxt = path[nextindex]
        path_block.append(nxt)
        index += iter
        if cur > num and nxt > num:
            if num > denom or cur + nxt == denom + (2 * num) + 1:
                loops.append("R")
                continue
            loops.append("C")
            continue
        if (cur > num and nxt <= num) or (cur <= num and nxt > num):
            loops.append("T")
            continue
        if  denom > num or cur + nxt == num + 1:
            loops.append("L")
            continue
        loops.append("C")
    return path_block, loops

# print(findpathwithloops(3, 1))
