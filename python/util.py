from nzmath import vector

def ffvector(compo, field):
    ffcmp = map(field.createElement, compo)
    return vector.Vector(ffcmp)

def range1(start, stop=None, step=None):
    if stop is None:
        return range(1, start+1)
    elif step is None:
        return range(start, stop+1)
    else:
        return range(start, stop+1, step)
        
if __name__ == "__main__":
    print range1(10)
    print range1(2,10)
    print range1(2,10,2)
