import numpy
import sklearn
import sklearn.metrics
import sklearn.metrics.cluster
import logging
import sys
import gc
import random
import itertools
import pandas
import pulp
from threading import Thread
import time
import math
try:
    import cplex
    import cplex.exceptions
    from cplex.exceptions import CplexError
except:
    pass

logger = logging.getLogger(__name__)
stdout = logging.StreamHandler(sys.stdout)
logger.handlers = []
logger.addHandler(stdout)
#logger.setLevel(logging.INFO)

def mixed_triplets(d):
    for i, j, k in itertools.combinations(range(d.shape[0]), 3):
        if d[i,j] <= 0 and d[i,k] <= 0 and d[j,k] <= 0:
            pass
        elif d[i,j] >= 0 and d[i,k] >= 0 and d[j,k] >= 0:
            pass
        else:
            yield i,j,k

def same_sign_triplets(d):
    for i, j, k in itertools.combinations(range(d.shape[0]), 3):
        if (d[i,j] <= 0 and d[i,k] <= 0 and d[j,k] <= 0) or (d[i,j] >= 0 and d[i,k] >= 0 and d[j,k] >= 0):
            yield i,j,k


# def processProblemWithPuLP(weights, all_constraints):
#     logger.debug("Creating problem instance ... ")
#     prob = pulp.LpProblem("problem", pulp.LpMinimize)
#     N = weights.shape[0]
#     numVariables = N * (N - 1) // 2
#     #variables = [pulp.LpVariable("x" + str(i), cat='Binary') for i in range(numVariables)]
#     variables = [pulp.LpVariable("x" + str(i), 0, 1) for i in range(numVariables)]
#     allPairs = list(itertools.combinations(range(N), 2))
#     mapDict = {pair : i for i, pair in enumerate(allPairs)}
#     #for i, j in allPairs:
#         #prob.solverModel.getVars()[mapDict[i, j]].start = start_solution[i, j]
#         #variables[mapDict[i, j]].setInitialValue(start_solution[i, j])
#     allWeights = [weights[i][j] for i, j in allPairs]
#     triplets = itertools.combinations(range(N), 3) if all_constraints else mixed_triplets(weights)
#     for i, j, k in triplets:
#         x1 = variables[mapDict[i, j]]
#         x2 = variables[mapDict[i, k]]
#         x3 = variables[mapDict[j, k]]
#         prob += x1 <= x2 + x3, "C_%d,%d,%d_1" % (i, j, k)
#         prob += x2 <= x1 + x3, "C_%d,%d,%d_2" % (i, j, k)
#         prob += x3 <= x1 + x2, "C_%d,%d,%d_3" % (i, j, k)
#     prob += pulp.lpDot(variables, allWeights)
#     gc.collect()
#     logger.debug("Solving ... ")
#     solver = pulp.solvers.COIN()
#     prob.setSolver(solver)
#     #prob.solver.buildSolverModel(prob)
#     #for i, j in allPairs:
#     #    prob.solverModel.getVars()[mapDict[i, j]].start = start_solution[i, j]
#     while True:
#         status = prob.solve(pulp.COIN())
#         logger.debug("Solution status: %s" % pulp.LpStatus[status])
#         solMatrix = numpy.zeros((N, N))
#         for i, pair in enumerate(allPairs):
#             solMatrix[pair[0]][pair[1]] = pulp.value(variables[i])
#             solMatrix[pair[1]][pair[0]] = solMatrix[pair[0]][pair[1]]
#         solMatrix[solMatrix < 0] = 0
#         solMatrix[solMatrix > 1] = 1
#         if all_constraints:
#             break
#         logger.debug("Processing solution ... ")
#         # check if any same sign triangle is violated:
#         violated = False
#         for i, j, k in same_sign_triplets(weights):
#             a, b, c = sorted([solMatrix[i][j], solMatrix[i][k], solMatrix[j][k]])
#             EPSILON = 10**-8
#             if c - a - b > EPSILON:
#                 logger.debug("Constraint violated for triplet %s, %s and %s." % (i, j, k))
#                 violated = True
#                 x1 = variables[mapDict[i, j]]
#                 x2 = variables[mapDict[i, k]]
#                 x3 = variables[mapDict[j, k]]
#
#                 prob += x1 <= x2 + x3, "C_%d,%d,%d_1" % (i, j, k)
#                 prob += x2 <= x1 + x3, "C_%d,%d,%d_2" % (i, j, k)
#                 prob += x3 <= x1 + x2, "C_%d,%d,%d_3" % (i, j, k)
#         if not violated:
#             break
#         logger.debug("Re-optimizing with all violated constraints added ...")
#     logger.debug("OBJ value: %.f" % prob.objective.value())
#     logger.debug("Finished PuLP solving.")
#     return solMatrix

def processProblem(Distances, all_constraints, start_solution=None):
    logger.debug("Creating problem instance ... ")
    my_prob = cplex.Cplex()
    N = Distances.shape[0]
    numConstraints = populateByNonZero(my_prob, Distances) if all_constraints else populateByNonZero_only_mixed(my_prob, Distances)
    gc.collect()
    #if not presolve:
    if start_solution is not None:
        start_vector = [start_solution[i, j] for i, j in itertools.combinations(range(N), 2)]
        my_prob.MIP_starts.add([range(len(start_vector)), start_vector], my_prob.MIP_starts.effort_level.solve_MIP) # CHANGE HERE!!

    my_prob.parameters.preprocessing.presolve.set(0) # try without this also.
    my_prob.parameters.emphasis.memory.set(1)  # try without this also.
    my_prob.parameters.timelimit.set(3600 * 2)
    #my_prob.parameters.workmem.set(1024)
    #my_prob.parameters.mip.strategy.file.set(3)
    my_prob.parameters.simplex.display.set(2)
    # set optimality gap to 1 over sum of all the negative weights
    #sum_neg = sum(Distances[Distances < 0])
    #my_prob.parameters.mip.tolerances.mipgap.set(1/abs(sum_neg))
    num_iterations = 0
    logger.debug("Solving ... ")
    while True:
        num_iterations += 1
        try:
            my_prob.solve()
            #print("iterations:", my_prob.solution.progress.get_num_iterations())
            
        except CplexError as exc:
            if exc.args[2] == cplex.exceptions.error_codes.CPXERR_NO_MEMORY:
                return
        logger.debug(" Processing solution ... ")
        X = my_prob.solution.get_values()
        solMatrix = [[0 for i in range(N)] for j in range(N)]
        for ind, pair in enumerate(itertools.combinations(range(N),2)):
            solMatrix[pair[0]][pair[1]] = X[ind]
            solMatrix[pair[1]][pair[0]] = X[ind]
        # check if any same sign triangle is violated:
        violated = False
        # map from x_(i,j) pair to a column:
        mapDict = {pair : i for i, pair in enumerate(itertools.combinations(range(N), 2))}
        for i,j,k in same_sign_triplets(Distances):
            a, b, c = sorted([solMatrix[i][j], solMatrix[i][k], solMatrix[j][k]])
            if c > a+b: # test triangle ineq.
                logger.debug(" Constraint violated for triplet %s, %s and %s." % (i,j,k))
                violated = True
                my_prob.linear_constraints.add(rhs = [0,0,0], senses = ["G","G","G"])
                my_prob.linear_constraints.set_coefficients(zip([numConstraints]*3, [mapDict[i,j], mapDict[i,k], mapDict[j,k]], [1,1,-1]))
                my_prob.linear_constraints.set_coefficients(zip([numConstraints+1]*3, [mapDict[i,j], mapDict[i,k], mapDict[j,k]], [1,-1,1]))
                my_prob.linear_constraints.set_coefficients(zip([numConstraints+2]*3, [mapDict[i,j], mapDict[i,k], mapDict[j,k]], [-1,1,1]))
                numConstraints += 3
        if not violated:
            break
        logger.debug("Re-optimizing with all violated constraints added ...")
    logger.debug("OBJ value: %.f" % my_prob.solution.get_objective_value())
    print(my_prob.solution.MIP.get_mip_relative_gap())
    print(my_prob.solution.MIP.get_best_objective())
    print("objective", my_prob.solution.get_objective_value())
    print("constraints", numConstraints)
    print("iterations", num_iterations)
    print("status", my_prob.solution.status[my_prob.solution.get_status()])
    logger.debug("Finished CPLEX solving.")
    return solMatrix

def populateByNonZero(prob, Distances):
    logger.debug(" Creating problem instance with all constraints")
    Distances = Distances.astype(float)
    N = Distances.shape[0]
    numVariables = N * (N - 1) // 2
    numConstraints = N * (N - 1) * (N - 2) // 2
    rowIndices = range(N)
    allPairs = itertools.combinations(rowIndices, 2)
    mapDict = {pair : i for i, pair in enumerate(allPairs)}
    allPairs = itertools.combinations(rowIndices, 2)
    allValues = [Distances[pair[0],pair[1]] for pair in allPairs]
    lowerBounds = [0] * numVariables
    upperBounds = [1] * numVariables
    my_rhs = [0] * numConstraints
    my_sense = "G" * numConstraints
    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.linear_constraints.add(rhs = my_rhs, senses = my_sense) #, names = my_rownames)
    # prob.variables.add(obj = allValues, ub = upperBounds, lb = lowerBounds) #, names = my_colnames)
    t = prob.variables.type
    prob.variables.add(obj = allValues, types = [t.binary] * numVariables)
    numBlocks = int(numConstraints/3)
    myRange = range(numBlocks)
    rows1 = itertools.chain.from_iterable(itertools.repeat(3 * x, 3) for x in myRange)
    rows2 = itertools.chain.from_iterable(itertools.repeat(3 * x + 1, 3) for x in myRange)
    rows3 = itertools.chain.from_iterable(itertools.repeat(3 * x + 2, 3) for x in myRange)
    vals1 = [1, 1, -1] * numBlocks
    vals2 = [1, -1, 1] * numBlocks
    vals3 = [-1, 1, 1] * numBlocks
    cols1 = itertools.chain.from_iterable((mapDict[x[0], x[1]], mapDict[x[0], x[2]], mapDict[x[1], x[2]]) for x in itertools.combinations(rowIndices, 3))
    prob.linear_constraints.set_coefficients(zip(rows1, cols1, vals1))
    cols2 = itertools.chain.from_iterable((mapDict[x[0], x[1]], mapDict[x[0], x[2]], mapDict[x[1], x[2]]) for x in itertools.combinations(rowIndices, 3))
    prob.linear_constraints.set_coefficients(zip(rows2, cols2, vals2))
    cols3 = itertools.chain.from_iterable((mapDict[x[0], x[1]], mapDict[x[0], x[2]], mapDict[x[1], x[2]]) for x in itertools.combinations(rowIndices, 3))
    prob.linear_constraints.set_coefficients(zip(rows3, cols3, vals3))
    return numConstraints

def populateByNonZero_only_mixed(prob, Distances):
    logger.debug(" Creating problem instance only with mixed constraints")
    N = Distances.shape[0]
    numVariables = int(N * (N - 1) / 2)
    mixed_trips = list(mixed_triplets(Distances))
    numConstraints = len(mixed_trips) * 3
    if numConstraints == 0:
        raise ValueError("There are zero mixed triplet constraints. Recommend using all constraints or adjusting threshold value.")
        sys.exit(1)
    rowIndices = range(N)
    allPairs = itertools.combinations(rowIndices, 2)
    mapDict = {pair : i for i, pair in enumerate(allPairs)}
    allPairs = itertools.combinations(rowIndices, 2)
    allValues = [Distances[pair[0],pair[1]] for pair in allPairs]
    lowerBounds = [0] * numVariables
    upperBounds = [1] * numVariables
    my_rhs = [0] * numConstraints
    my_sense = "G" * numConstraints
    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.linear_constraints.add(rhs = my_rhs, senses = my_sense) #, names = my_rownames)
    #prob.variables.add(obj = allValues, ub = upperBounds, lb = lowerBounds) #, names = my_colnames)
    t = prob.variables.type
    prob.variables.add(obj = allValues, types = [t.binary] * numVariables)
    numBlocks = int(numConstraints/3)
    myRange = range(numBlocks)
    rows1 = itertools.chain.from_iterable(itertools.repeat(3 * x, 3) for x in myRange)
    rows2 = itertools.chain.from_iterable(itertools.repeat(3 * x + 1, 3) for x in myRange)
    rows3 = itertools.chain.from_iterable(itertools.repeat(3 * x + 2, 3) for x in myRange)
    vals1 = [1, 1, -1] * numBlocks
    vals2 = [1, -1, 1] * numBlocks
    vals3 = [-1, 1, 1] * numBlocks
    cols1 = itertools.chain.from_iterable((mapDict[x[0], x[1]], mapDict[x[0], x[2]], mapDict[x[1], x[2]]) for x in mixed_trips)
    prob.linear_constraints.set_coefficients(zip(rows1, cols1, vals1))
    cols2 = itertools.chain.from_iterable((mapDict[x[0], x[1]], mapDict[x[0], x[2]], mapDict[x[1], x[2]]) for x in mixed_trips)
    prob.linear_constraints.set_coefficients(zip(rows2, cols2, vals2))
    cols3 = itertools.chain.from_iterable((mapDict[x[0], x[1]], mapDict[x[0], x[2]], mapDict[x[1], x[2]]) for x in mixed_trips)
    prob.linear_constraints.set_coefficients(zip(rows3, cols3, vals3))
    return numConstraints

def f_plus(x,a=0.19, b=0.5095):
    if x < a:
        return 0
    if x >= b:
        return 1
    return ((x-a)/(b-a)) ** 2

def f_minus(x):
    return x

def prob(u, v, sol_matrix, weight_matrix):
    return f_plus(sol_matrix[u][v]) if weight_matrix[u][v] >= 0 else f_minus(sol_matrix[u][v])

def e_cost(u, v, w, sol_matrix, weight_matrix, probs):
    return probs[u][w] * (1 - probs[v][w]) + (1 - probs[u][w]) * probs[v][w] if weight_matrix[u][v] >= 0 else (1 - probs[u][w]) * (1 - probs[v][w])

def e_lp(u, v, w, sol_matrix, weight_matrix, probs):
    return (1 - probs[u][w] * probs[v][w]) * sol_matrix[u][v] if weight_matrix[u][v] >= 0 else (1 - probs[u][w] * probs[v][w]) * (1 - sol_matrix[u][v])

def alg(u, v, w, sol_matrix, weight_matrix, probs):
    return e_cost(u, v, w, sol_matrix, weight_matrix, probs) + e_cost(w, u, v, sol_matrix, weight_matrix, probs) + e_cost(v, w, u, sol_matrix, weight_matrix, probs)

def lp(u, v, w, sol_matrix, weight_matrix, probs):
    return e_lp(u, v, w, sol_matrix, weight_matrix, probs) + e_lp(w, u, v, sol_matrix, weight_matrix, probs) + e_lp(v, w, u, sol_matrix, weight_matrix, probs)

def cost_func(u, v, w, sol_matrix, weight_matrix, probs, alpha):
    return alpha * lp(u, v, w, sol_matrix, weight_matrix, probs) - alg(u, v, w, sol_matrix, weight_matrix, probs)

def sum_cost_func(u, v, v_set, sol_matrix, weight_matrix, probs, alpha):
    c = 0
    for w in v_set:
        c += cost_func(u, v, w, sol_matrix, weight_matrix, probs, alpha)
    return c

def fixed_pivot_cost_func(u, v, pivot, sol_matrix, weight_matrix, probs, alpha):
    return alpha * e_lp(u, v, pivot, sol_matrix, weight_matrix, probs) - e_cost(u, v, pivot, sol_matrix, weight_matrix, probs)

def sum_fixed_pivot_cost_func(v_set, pivot, sol_matrix, weight_matrix, probs, alpha):
    c = 0
    for u, v in itertools.combinations(v_set, 2):
        c += fixed_pivot_cost_func(u, v, pivot, sol_matrix, weight_matrix, probs, alpha)
    return c

def best_prob(u, v, v_set, sol_matrix, weight_matrix, probs, alpha=2.06):
    uv_prob = probs[u][v]
    probs[u][v] = 0
    probs[v][u] = 0
    c0 = sum_cost_func(u, v, v_set, sol_matrix, weight_matrix, probs, alpha)
    probs[u][v] = 1
    probs[v][u] = 1
    c1 = sum_cost_func(u, v, v_set, sol_matrix, weight_matrix, probs, alpha)
    probs[u][v] = uv_prob
    probs[v][u] = uv_prob
    return 0 if c0 > c1 else 1

def best_pivot(v_set, sol_matrix, weight_matrix, probs, alpha=2.06):
    max_score = -1
    max_pivot = -1
    for pivot in v_set:
        score = sum_fixed_pivot_cost_func(v_set, pivot, sol_matrix, weight_matrix, probs, alpha)
        if score > max_score:
            max_score = score
            max_pivot = pivot
    return max_pivot

def add_element_test(u, v, sol_matrix, weight_matrix):
    # depending on the weight_matrix value for u,v, use f_plus or f_minus to find the probability of adding this element,
    # then sample on [0,1] and return true if larger than the probability
    return random.random() < 1 - prob(u, v, sol_matrix, weight_matrix)

def chawla_rounding(sol_matrix, weight_matrix):
    logger.debug(" Running Chawla rounding ... ")
    v_set = set(range(len(sol_matrix)))
    clustering = []
    while len(v_set) > 0:
        pivot = random.sample(v_set, 1)[0]
        current_cluster = [v for v in v_set if add_element_test(pivot, v, sol_matrix, weight_matrix)]
        for v in current_cluster:
            v_set.remove(v)
        clustering.append(sorted(current_cluster))
    return clustering

def derandomized_chawla_rounding(sol_matrix, weight_matrix):
    logger.debug(" Running derandomized Chawla rounding ... ")
    v_set = set(range(len(sol_matrix)))

    # setting all the probablities p_uv using function f
    probs = [[prob(u, v, sol_matrix, weight_matrix) for v in v_set] for u in v_set]

    # changing the value of each p_uv to either 0 or 1 to maximize the cost function
    for u, v in itertools.combinations(v_set, 2):
        probs[u][v] = best_prob(u, v, v_set, sol_matrix, weight_matrix, probs)
        probs[v][u] = probs[u][v]

    # picking the pivot to maximize the cost function and clustering like the randomized version
    clustering = []

    while len(v_set) > 0:
        pivot = best_pivot(v_set, sol_matrix, weight_matrix, probs)
        current_cluster = {v for v in v_set if probs[pivot][v] == 0}
        v_set = v_set.difference(current_cluster)
        clustering.append(sorted(current_cluster))

    return clustering

def clustering_to_pandas(list_of_clusters,samples):
    '''
    Returns representation of a clustering (represented as a list of lists) as a cluster assignment
    vector, represented as a Pandas Dataframe
    '''
    logger.debug(" Organizing clustering results ... ")
    clustering_temp = []
    for cluster_num,cluster in enumerate(list_of_clusters):
        for sample_num in cluster:
            clustering_temp.append( [samples[sample_num],cluster_num+1] )
    clustering = pandas.DataFrame(clustering_temp,columns=['Sample','Cluster'])
    clustering.set_index('Sample',inplace=True)
    return clustering

# based on https://stackoverflow.com/a/2785908/1056345
def wait_until(somepredicate, timeout, period=0.25, *args, **kwargs):
    must_end = time.time() + timeout
    while time.time() < must_end:
        if somepredicate(*args, **kwargs):
            return True
        time.sleep(period)
        logger.debug('waiting for ', must_end - time.time())
    return False

def createCluster(v, m, pi, pi_dict, G, clusterIDs):
    clusterIDs[v] = pi_dict[v]
    for u in pi[m:]:
        if G[u, v] == 1:
            clusterIDs[u] = min(clusterIDs[u], pi_dict[v])

def isCenter(v, pi, pi_dict, G, clusterIDs, is_center_dict):
    if v in is_center_dict:
        return is_center_dict[v]
    for u in pi:
        if pi_dict[u] < pi_dict[v]:
            if G[u, v] == 1:
                if not wait_until(lambda x, idx: x[idx] != math.inf, 5, 0.1, clusterIDs, u):
                    logger.debug('Timeout!', u, clusterIDs[u], clusterIDs[u] != math.inf)
                if isCenter(u, pi, pi_dict, G, clusterIDs, is_center_dict):
                    is_center_dict[v] = 0
                    return 0
        else:
            break
    is_center_dict[v] = 1
    return 1

def attemptCluster(v, m, pi, pi_dict, G, clusterIDs, is_center_dict):
    if clusterIDs[v] == math.inf and isCenter(v, pi, pi_dict, G, clusterIDs, is_center_dict):
        createCluster(v, m, pi, pi_dict, G, clusterIDs)

def c4(G, epsilon):
    n = G.shape[0]
    clusterIDs = math.inf * numpy.ones(n)
    pi = numpy.random.permutation(n)
    pi_dict = {pi[i]:i for i in range(n)}
    max_deg = G.sum(axis=1).max()
    min_deg = G.sum(axis=1).min()
    delta = max_deg
    round = 0
    while len(pi) > 0:
        min_deg = 1 if min_deg == 0 else min_deg
        ratio = max_deg / min_deg
        ratio = 1.1 if ratio <= 1 else ratio
        if round > (2 / epsilon) * math.log(n * math.log2(ratio)):
            round = 0
            delta /= 2
        else:
            round += 1
        # delta = G[pi][:,pi].sum(axis=1).max()
        m = math.ceil(epsilon * n / delta)
        A = set(pi[:m])
        jobs = []
        is_center_dict = dict()
        while A: # parallel
            v = A.pop()
            j = Thread(target = attemptCluster, args=(v, m, pi, pi_dict, G, clusterIDs, is_center_dict))
            j.start()
            jobs += [j]
        for j in jobs:
            j.join()
        pi = pi[m:]
    return clusterIDs

def c4_correlation(distance_matrix, threshold):
    threshold = float(threshold)
    samples = distance_matrix.columns.values
    n = distance_matrix.shape[0]
    weight_matrix = threshold - distance_matrix
    G = numpy.where(weight_matrix > 0, numpy.ones((n, n)), numpy.zeros((n, n)))
    clusterIDs = c4(G, 0.5)
    clusters_dict = dict()
    for idx, c in enumerate(clusterIDs):
        if c in clusters_dict:
            clusters_dict[c].append(idx)
        else:
            clusters_dict[c] = [idx]
    list_of_clusters = sorted(clusters_dict.values(), key=lambda x:x[0])
    clustering = clustering_to_pandas(list_of_clusters,samples)
    return clustering

def obj_func(distance_matrix, threshold, clustering):
    samples = distance_matrix.columns.values
    sim_matrix = threshold - distance_matrix
    obj_value = 0
    for i, j in itertools.combinations(samples, 2):
        s_ij = sim_matrix[i][j]
        x_ij = 0 if clustering[i] == clustering[j] else 1
        obj_value += x_ij * s_ij
    return obj_value

def multiple_c4(distance_matrix, threshold, repeat=20):
    best_clustering = None
    for i in range(repeat):
        clustering = c4_correlation(distance_matrix, threshold)
        obj_val = obj_func(distance_matrix, threshold, clustering['Cluster'])
        if best_clustering is None or obj_val < min_obj_val:
            min_obj_val = obj_val
            best_clustering = clustering
    return best_clustering

def dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            neighbors = set(numpy.where(graph[vertex] == 0)[0])
            
            stack.extend(neighbors - visited)
    return visited

def make_clustering(sol_matrix):
    v_set = set(range(len(sol_matrix)))
    clustering = []
    while len(v_set) > 0:
        start = v_set.pop()
        current_cluster = dfs(numpy.array(sol_matrix), start)
        v_set -= current_cluster
        clustering.append(sorted(current_cluster))

    return clustering

def correlation(distance_matrix, threshold, all_constraints=False, method='ILP'):
    '''
    Given a distance matrix as a Pandas DataFrame and a distance threshold, solve a correlation
    clustering problem instance LP problem and then apply the Chawla et al. 2015 rounding algorithm,
    guaranteeing a 2-epsilon approximation ratio.
    @param distance_matrix: distance matrix represented as a Pandas DataFrame object, doubly indexed
                            by sample names
    @param threshold: a threshold value, in form of an int or a float
    @param all_constraints: boolean indicating whether all triangle inequality constraints should be
                            used in the CPLEX problem
    @param method: the method that is used to solve the correlation clustering, which is one of these: 'C4', 'ILP', 'C4+ILP'
    @rvalue clustering: the approximate optimal clustering represented as a Pandas DataFrame
    '''
    threshold = float(threshold)
    samples = distance_matrix.columns.values
    weight_matrix = threshold - distance_matrix    
    if method == 'C4':
        #clustering = c4_correlation(distance_matrix, threshold)
        clustering = multiple_c4(distance_matrix, threshold)
    else:
        #logger.info("Solving instance for threshold value " + str(threshold) + " ...")
        #if solver == 'cplex':
        #    if not sol_matrix:
        #        raise CplexError
        #elif solver == 'pulp':
        #    sol_matrix = processProblemWithPuLP(weight_matrix.values, all_constraints)
        #else:
        #    print("Error: unsupported solver %s" % (solver))
        #    sys.exit(1) 
        #    if not sol_matrix:
        print("threshold", threshold)
        print("method", method)
        print("all_con", all_constraints)
        
        c4_clustering = c4_correlation(distance_matrix, threshold)
        indexes = distance_matrix.index
        N = distance_matrix.shape[0]
        start_solution = numpy.ones((N, N))
        numpy.fill_diagonal(start_solution, 0)
        for i, j in itertools.combinations(range(N), 2):
                if c4_clustering['Cluster'].loc[indexes[i]] == c4_clustering['Cluster'].loc[indexes[j]]:
                    start_solution[i, j] = 0
                    start_solution[j, i] = 0

        sol_matrix = processProblem(weight_matrix.values, all_constraints, start_solution)
        #sol_matrix = processProblem(weight_matrix.values, all_constraints)
        if not sol_matrix:
            raise CplexError
        list_of_clusters = sorted(make_clustering(sol_matrix), key=lambda x:x[0])
        clustering = clustering_to_pandas(list_of_clusters,samples)
    
    #logger.info("Applying Chawla rounding ...")
    #list_of_clusters = sorted(derandomized_chawla_rounding(sol_matrix,weight_matrix.values),
    #                          key=lambda x:x[0])
    #logger.info("Done! %d clusters found" % clustering['Cluster'].values.max())
    return clustering

def multiple_correlation(distance_matrix, thresholds, all_constraints=False,solver='pulp', method='C4+ILP'):
    '''
    Perform correlation clustering on a list of thresholds
    @param distance_matrix: distance matrix represented as a Pandas DataFrame object, doubly indexed
                            by sample names
    @param thresholds: a list of threshold values, where the values can be ints or floats
    @param all_constraints: boolean indicating whether all triangle inequality constraints should be
                            used in the CPLEX problem
    @param solver: the solver to use to solve the correlation clustering instance
    @rvalue clusterings: a dictionary of clusterings (represented by Pandas DataFrames), indexed by
                         threshold values in thresholds
    '''
    clusterings = {}
    for threshold in thresholds:
        clustering = correlation(distance_matrix,threshold,all_constraints,solver, method)
        clusterings[threshold] = clustering
    return clusterings

def create_contingency_table(A_vect,B_vect):
    '''
    Construct the contingency table as defined in the paper
    "A Confidence Interval for the Wallace Coefficient of Concordance
    and Its Applications to Microbial Typing Methods" by Pinto, Melo-Cristino, and Ramirez
    @param A: left-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @param B: right-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @rvalue ct: the contingency table, definition as in the paper.
    '''
    assert(set(A_vect.index.values) == set(B_vect.index.values)),\
        "Clusterings A and B do not describe the same set of samples"

    A_copy = A_vect.copy().iloc[:,0]
    B_copy = B_vect.copy().iloc[:,0]
    samples = list(A_vect.index.values)
    A_clusters = set(A_vect.iloc[:,0].tolist())
    B_clusters = set(B_vect.iloc[:,0].tolist())

    # Index clusterings by cluster
    A = {i: set([sample for sample in samples if A_copy[sample] == i])
         for i in A_clusters}
    A_num_samples = sum([len(A[cluster]) for cluster in A.keys()])

    B = {i: set([sample for sample in samples if B_copy[sample] == i])
         for i in B_clusters}
    B_num_samples = sum([len(B[cluster]) for cluster in B.keys()])

    assert(A_num_samples == B_num_samples == len(samples))
    # contingency table
    ct = pandas.DataFrame(numpy.zeros(shape=(len(A_clusters),len(B_clusters)),dtype=int),
                          index = list(A_clusters), columns = list(B_clusters))
    for i in A_clusters:
        A_cluster = A[i]
        for j in B_clusters:
            B_cluster = B[j]
            intersection = (A_cluster & B_cluster)
            ct.loc[i,j] = len(intersection)
    assert( ct.values.sum() == len(samples) )
    return ct

def compute_wallace_coefficient(A,B,ct=None):
    '''
    Calculate the Wallace coefficient from clustering A to clustering B.
    Definition taken from "A Confidence Interval for the Wallace Coefficient of Concordance
    and Its Applications to Microbial Typing Methods" by Pinto, Melo-Cristino, and Ramirez
    @param A: left-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @param B: right-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @param ct: the contingency table, definition as in the paper.
    @rvalue W: the Wallace coefficient of the two clusterings
    '''
    if ct is None:
        ct = create_contingency_table(A,B)
    numerator = (ct*(ct-1)).values.sum()
    row_sums = ct.sum(axis=1)
    denominator = (row_sums*(row_sums-1)).values.sum()
    assert(numerator <= denominator)
    W = numerator/denominator
    return W


def compute_expected_wallace_coefficient(A,B,ct=None):
    '''
    Calculate the expected value of the Wallace coefficient of clusterings A,B if A,B are
    independent.
    Definition taken from "A Confidence Interval for the Wallace Coefficient of Concordance
    and Its Applications to Microbial Typing Methods" by Pinto, Melo-Cristino, and Ramirez
    @param A: left-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @param B: right-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @param ct: the contingency table, definition as in the paper.
    @rvalue exp_W: the expected Wallace coefficient of the two clusterings
    '''
    if ct is None:
        ct = create_contingency_table(A,B)
    numerator = (ct*(ct-1)).values.sum()
    num_samples = len(list(A.index.values))
    denominator = num_samples * (num_samples-1)
    assert(numerator <= denominator)
    simpson = 1 - (numerator/denominator)
    exp_W = 1 - simpson
    return exp_W

def compute_adjusted_wallace_coefficient(A,B):
    '''
    Compute the adjusted wallace coefficient, as defined in the paper
    "Adjusted Wallace Coefficient as a Measure of Congruence between Typing Methods" by
    Severiano et al 2011.
    @param A: left-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @param B: right-hand input clustering, in the form of a Pandas DataFrame where the lone
                       column is labelled 'Cluster', and the index is sample names
    @rvalue awc: the adjusted wallace coefficient of the two clusterings
    '''
    ct = create_contingency_table(A,B)
    wallace = compute_wallace_coeffient(A,B,ct)
    expected_wallace = compute_expected_wallace_coeffient(A,B,ct)
    numerator = wallace - expected_wallace
    denominator = 1 - expected_wallace
    if denominator == 0:
        awc = 1
    else:
        awc = numerator/denominator
    return awc

def compute_neighbor_adjusted_wallace_coefficients(clusterings,descending=True):
    '''
    Compute the neighbor adjusted wallace coefficients of a dictionary of clusterings, where the keys
    are threshold values.
    Definition taken from paper "Rapid Identification of Stable Clusters in Bacterial Populations Using
    the Adjusted Wallace Coefficient" by Barker et al, 2018.
    @param clusterings: a dictionary of clusterings (represented by Pandas DataFrames), indexed by
                         threshold values
    @param descending: boolean indicating whether to compute nAWC in descending order (if True), or
                       descending order otherwise.
    @rvalue nawc_values: a Pandas DataFrame object containing the nAWC values. Indexed by threshold.
    '''
    thresholds = sorted(clusterings.keys(),key=float)
    if descending:
        thresholds = thresholds[::-1]
    num_thresholds = len(thresholds)
    nawc_values = pandas.DataFrame(numpy.zeros(shape=(len(thresholds),1), columns=['nAWC'],
                                               index=thresholds))
    for i in range(1,len(thresholds)):
        from_threshold = thresholds[i-1]
        to_threshold = thresholds[i]
        from_clust = clusterings[from_threshold]
        to_clust = clusterings[to_threshold]
        nawc = compute_adjusted_wallace_coefficient(from_clust,to_clust)
        nawc_values.loc[to_threshold,'nAWC'] = nawc
    return nawc_values

def cluster_vector_to_matrix(cluster_vector):
    '''
    Returns a matrix-representation of the input vector-representation of a clustering.
    Inputs:
    - cluster_vector: a dictionary where the keys are sample names and the output is an integer
                      representing the assigned cluster.
    Outputs:
    - cluster_matrix: a matrix represented as a pandas dataframe
                      Indexed by sample names, and the output is 0 if the samples are in the same
                      cluster, 0 otherwise.
                      e.g.
                      cluster_matrix[sample1][sample2] = 0 # they are in the same cluster
                      cluster_matrix[sample1][sample3] = 1 # they are not in the same cluster
    '''
    samples = cluster_vector.index.values
    cluster_matrix = pandas.DataFrame(data=1,index=samples,columns=samples,dtype=float)
    inverted_cluster_vector = {}
    sample_cluster_pairs = [(sample,cluster_vector.loc[sample,'Cluster']) for sample in samples]

    for sample,cluster in sample_cluster_pairs:
        if cluster not in inverted_cluster_vector:
            inverted_cluster_vector[cluster] = []
        inverted_cluster_vector[cluster].append(sample)

    # reflexivity of cluster matrix
    for sample in samples:
        cluster_matrix[sample][sample] = 0

    for cluster in inverted_cluster_vector:
        for sample1, sample2 in itertools.combinations(inverted_cluster_vector[cluster],2):
            # symmetry of cluster matrix
            cluster_matrix[sample1][sample2] = 0
            cluster_matrix[sample2][sample1] = 0

    return cluster_matrix

def normalize_distances(distances):
    '''
    Normalize input distance matrices
    '''
    normal_distances = {}
    for clustering in distances:
        max_value = numpy.amax(distances[clustering])
        normal_distances[clustering] = distances[clustering]/max_value
    return normal_distances

def construct_consensus_D(clusterings,normal_distances,fine_clusterings):
    '''
    Construct the D matrix as described in the document describing consensus clustering.
    '''
    samples = clusterings[list(clusterings.keys())[0]].columns.values
    num_samples = len(samples)
    D = pandas.DataFrame(numpy.zeros(shape=(num_samples,num_samples),dtype=float),
                         index=samples,columns=samples)
    for sample1,sample2 in itertools.combinations(samples,2):
        for clustering in normal_distances:
            if clusterings[clustering][sample1][sample2] == 0 or clustering in fine_clusterings:
                D[sample1][sample2] += normal_distances[clustering][sample1][sample2]
                D[sample2][sample1] = D[sample1][sample2]
    return D

def construct_consensus_weights(clustering_vectors,distances,fine_clusterings):
    '''
    @parameter clustering_vectors: dictionary of pandas dataframe representing clusterings as vectors
    @parameter distances: dictionary of pandas dataframes representing distance matrices of
                          different data types
    @parameter fine_clusterings: list of key values for clustering/distances corresponding to "finest"
                                 clusterings
    @rvalue S: the consensus clustering weight matrix represented as a Pandas Dataframe object
    '''
    clusterings = {key: cluster_vector_to_matrix(clustering_vectors[key]) 
                   for key in clustering_vectors.keys()}
    normal_distances = normalize_distances(distances)
    # now construct the Pi and D matrices
    samples = clusterings[list(clusterings.keys())[0]].columns.values
    num_samples = len(samples)
    unlabeled_pi = numpy.full(shape=(num_samples,num_samples),fill_value=len(clusterings.keys()),
                              dtype=float)
    Pi = pandas.DataFrame(unlabeled_pi,index=samples,columns=samples)
    for clustering in clusterings:
        Pi = Pi.subtract(clusterings[clustering])
    D = construct_consensus_D(clusterings,normal_distances,fine_clusterings)
    S = Pi.subtract(D)
    return S

def consensus(distances,clusterings,fine_clusterings,weight_matrix=None, all_constraints=False, method="C4"):
    '''
    Solve an instane of consensus clustering.
    @param clusterings: dictionary of pandas dataframe representing multiple clusterings as vectors
    @param distances: dictionary of pandas dataframes representing distance matrices of
                      different data types
    @param fine_clusterings: list of key values for clustering/distances corresponding to "finest"
                             clusterings
    @param weight_matrix (optional): precomputed consensus clustering weight matrix
    @param all_constraints: boolean variable indicating whether to use all_constraints, or only those
                            involving mixed triplets
    @rvalue clustering: a Pandas DataFrame
    '''
    if weight_matrix is None:
        #clustering_matrices = {key: cluster_vector_to_matrix(clusterings[key]) for key in clusterings.keys()}
        weight_matrix = construct_consensus_weights(clusterings,distances,fine_clusterings)
    clustering = correlation(-weight_matrix, 0, all_constraints, method)
    '''
    samples = weight_matrix.columns.values
    
    if solver == 'cplex':
        sol_matrix = processProblem(weight_matrix.values,all_constraints)
        if not sol_matrix:
            raise CplexError
    elif solver == 'pulp':
        if solver == 'pulp':
        sol_matrix = processProblemWithPuLP(weight_matrix.values,all_constraints)
    else:
        print("Error: unsupported solver %s" % (solver))
        sys.exit(1)
    list_of_clusters = sorted(derandomized_chawla_rounding(sol_matrix,weight_matrix.values)
                               ,key=lambda x:x[0])
    # Turn the list of clusters into pandas data frame
    clustering = clustering_to_pandas(list_of_clusters,samples)
    logger.info(" Done! %d clusters found" % clustering['Cluster'].values.max())
    '''    
    clustering.columns = ['Consensus']
    return clustering

def summarize_clusterings(consensus_clustering,clusterings):
    '''
    Returns a single dataframe containing all clusterings
    '''
    clusterings_copy = {key: clusterings[key].copy(deep=True) for key in clusterings.keys()}
    for cluster in clusterings_copy:
        clusterings_copy[cluster].columns = [cluster]
    all_clusterings = [consensus_clustering] \
                    + [clusterings_copy[name] for name in clusterings_copy.keys()]
    summary_clustering = pandas.concat(all_clusterings,axis=1,sort=True)
    return summary_clustering

def adjusted_rand_index(clustering1,clustering2):
    '''
    Compares the similarity of two clusterings based on the adjusted rand index.
    '''
    return sklearn.metrics.cluster.adjusted_rand_score(clustering1.sort_index().T.values.flatten(),
                                                       clustering2.sort_index().T.values.flatten())

def cluster_purity(ground_truth, predicted):
    y_true = numpy.array(ground_truth.loc[:,'Cluster'])
    y_pred = numpy.array(predicted.loc[:,'Cluster'])
    contingency_matrix = sklearn.metrics.cluster.contingency_matrix(y_true, y_pred)
    return numpy.sum(numpy.amax(contingency_matrix, axis=0)) / numpy.sum(contingency_matrix)


