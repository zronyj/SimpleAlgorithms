import random

def validate(iters, typ=0):
    ctrl = 0
    try:
        if typ == 0:
            ctrl = int(iters)
        else:
            ctrl = float(iters)
        return False
    except ValueError:
        print "No ingresaste bien el numero."
        return True

def get_data(get_message,retry_message,typ=0):
    data = raw_input(get_message)
    while validate(data, typ):
        data = raw_input(retry_message)
    if typ == 0:
        return int(data)
    else:
        return float(data)

def get_range():
    a = get_data("Por favor ingresa el limite inferior de integracion: ",\
        "Por favor reingresa el numero: ", 1)
    b = get_data("Por favor ingresa el limite superior de integracion: ",\
        "Por favor reingresa el numero: ", 1)
    return [a,b]

def get_iters():
    n = get_data("Por favor ingresa el numero de iteraciones: ",\
    "Por favor reingresa el numero: ", 0)
    return n

def find_max(f, rang):
    lim = 10.**(-7)
    pp = rang[0]
    np = rang[1]
    while abs(np - pp) > lim:
        mid = (np - pp) / 2 + pp
        if f(np) > f(pp) and f(mid) > f(pp):
            if f(np) > f(mid):
                pp = mid
            elif f(np) < f(mid):
                pp = np
                np = mid
            else:
                return np
        elif f(pp) > f(np) and f(mid) > f(np):
            if f(pp) > f(mid):
                np = mid
            elif f(pp) < f(mid):
                np = pp
                pp = mid
            else:
                return pp
        else:
            if f(pp) == f(np) and f(mid) > f(pp):
                left = find_max(f, [pp, mid])
                right = find_max(f, [mid,np])
                if left > right:
                    return left
                elif left < right:
                    return right
                elif left == right:
                    return left
                else:
                    return mid
            else:
                return pp
    return np


def monte_carlo(f, iters, rang, f_max, box_x):
    inner_points = 0
    for i in range(iters):
        x = random.random() * box_x + rang[0]
        y = random.random() * f_max
        v = f(x)
        if ( v > 0 and y < v ) or ( v < 0 and y > v ):
            inner_points += 1
    box_area = abs(box_x) * abs(f_max)
    point_relation = inner_points * 1./iters
    area = box_area * point_relation
    return area

def mmc(f):
    area = 0
    iters = get_iters()
    rang = get_range()
    f_max = f(find_max(f, rang))
    box_x = rang[1] - rang[0]
    for i in range(100):
        area += monte_carlo(f, iters, rang, f_max, box_x)
    area /= 100
    return area
