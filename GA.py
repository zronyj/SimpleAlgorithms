import random, string

# Creates the first generation (type 1 = integers, type 0 = floats)
def first_gen(population=100,dimension=10,interval=1,tipo=0):
    unit = []
    for j in range(population):
        chromosome = []
        if tipo == 0:
            for i in range(dimension): chromosome.append(random.random()*interval)
        elif tipo == 1:
            for i in range(dimension): chromosome.append(int(random.random()*interval))
        unit.append(chromosome)
    return unit

# Makes crossovers between chromosomes (0 = single point, 1 = double point, 2 = uniform)
def crossover(unit, operator=0):
    largo = len(unit)
    cromo = len(unit[0])
    if operator == 0:
        segment = 0
        while segment >= cromo or segment == 0:
            print "Crossover point out of range.\n\n Recalculating ...\n\n"
            segment = int(random.random()*cromo)
        offsprings = []
        for i in range(largo):
            for j in range(largo):
                if i != j:
                    offsprings.append( unit[i][:segment] + unit[j][segment:] )
        return unit + offsprings
    elif operator == 1:
        segments = []
        for i in range(2):
            segment = 0
            while segment >= cromo or segment == 0:
                print "Crossover point out of range.\n\n Recalculating ...\n\n"
                segment = int(random.random()*cromo)
            segments.append(segment)
        segments.sort()
        offsprings = []
        for i in range(largo):
            for j in range(largo):
                for k in range(largo):
                    if i != j and i != k and j != k:
                        offsprings.append( unit[i][:segments[0]] + unit[j][segments[0]:segments[1]] + unit[k][segments[1]:])
        return unit + offsprings
    elif operator == 2:
        for i in range(largo):
            for j in range(largo):
                if i != j:
                    new = []
                    for k in range(cromo):
                        new.append( random.choice([ unit[i][k], unit[j][k] ]) )
                    offsprings.append(new[:])
        return unit + offsprings

# Introduces mutations in the chromosomes
def mutation(unit, mut_factor=0.1, interval=1):
    xmen = unit[:]
    for mutant in range(len(xmen)):
        for fac in range( int( len(xmen[0]) * mut_factor ) ):
            gene = int(random.random()*len(xmen[mutant]))
            if type(xmen[mutant][gene]) == float:
                xmen[mutant][gene] = random.random()
            elif type(xmen[mutant][gene]) == int:
                xmen[mutant][gene] = int(random.random()*interval)
    return unit + xmen

# Deletes chromosomes which are the same (everyone is unique!)
def multi_del(unit):
    unique = [unit[0]]
    for i in unit:
        ctrl = 0
        for j in unique:
            if i != j:
                ctrl += 1
        if ctrl == len(unique):
            unique.append(i)
    return unique

# Function defining which chromosome is the best (CHANGES FOR DIFFERENT APPLICATIONS)
def scoring(chromosome):
    ctrl = 0
    for i in chromosome:
        ctrl += i
    return ctrl

# Kills the chromosomes that are too weak
def natural_sel(unit, limit=10):
    scores = []
    for i in unit:
        scores.append(scoring(i))
    top = []
    for j in range(limit):
        high = 0
        posi = 0
        for k in list(set(range(len(scores))) - set(top)):
            if scores[k] > high:
                high = scores[k]
                posi = k
        top.append(posi)
    best = []
    for l in range(len(top)):
        best.append(unit[top[l]])
    return best

# Just something to visualize the best chromosome
def results(unit):
    print "*"*60
    print "\t\t\tResults"
    print "*"*60
    print "\n\tScore\tChromosome\n"
    for i in range(len(unit)):
        print "\t" + str(scoring(unit[i])) + "\t" + str(unit[i])
    print "\n" + "*"*60

# Commented genetic algorithm
def GAv(initial_pop=100, generations=5, chrom_len=10, interval=1, kind=1, mut=0.1, expected=50, final_pop=10):
    print "*"*60 + "\n\t\tGenetic Algorithm\n" + "*"*60 + "\n\n"
    print "- Creating first population:", initial_pop
    unit = first_gen(initial_pop, chrom_len, interval, kind)
    for i in range(generations):
        print "\n\n- Making generation:", i+1, "\n"
        unit = crossover(unit, 0)
        print "-- Recombined\tPopulation:", len(unit)
        unit = mutation(unit, mut, interval)
        print "-- Mutated!\tMutants:",len(unit)/2
        unit = multi_del(unit)
        print "-- Filtered!\tPopulation:", len(unit)
        unit = natural_sel(unit, expected)
        print "-- Selected!\tPopulation:", len(unit)
    unit = natural_sel(unit, final_pop)
    print "\n\n- Final population:", len(unit)
    print "\n-Process finished!\n\n" + "*"*60
    results(unit)
    return unit

# Uncommented genetic algorithm
def GA(initial_pop=100, generations=5, chrom_len=10, interval=1, kind=1, mut=0.1, expected=50, final_pop=10):
    unit = first_gen(initial_pop, chrom_len, interval, kind)
    for i in range(generations):
        unit = crossover(unit, 0)
        unit = mutation(unit, mut)
        unit = multi_del(unit)
        unit = natural_sel(unit, expected)
    unit = natural_sel(unit, final_pop)
    return unit
