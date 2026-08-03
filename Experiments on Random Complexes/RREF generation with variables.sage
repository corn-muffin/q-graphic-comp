###########This seeks to more efficiently generate the possible rrefs by storing entries as variables to be evaluated later
###########Unfortunately, many more variables end up being introduced than what is necessary.  To improve we could look for any variables which only appear once and make the corresponding entry into a.  We'd then want to reindex so that the variables are x0,...,x(k-1) again.  This would speed up the second part but not all_RREFS().
###########There may still be errors.  As the improvements don't seem to bring computing [[1,2,7],[1,6,7],[2,7,8],[2,8,3],[3,4,8],[4,8,0],[4,5,0],[5,6,0],[6,7,0],[7,8,0]] within computational reach, however, it may not be worth the time to fully correct these.

import numpy

def normal_vectors(delta): #vertices of delta must be 0,1,...,n-1
    n = len(delta.vertices())
    vecs = []
    
    for i in range(n):
        faces = [F for F in delta.face_iterator()
             if i in F and all(v >= i for v in F)]#faces with i but not 1,...,i-1
        delta2 = SimplicialComplex(faces)
        if(len(delta2.vertices())==0):
            continue
        for F in delta2.facets():
            v = ["0"] * n
            v[F[0]]="1"
            for i in F[1:]:
                v[i]="a"
            vecs.append(numpy.matrix([v],dtype=object))
    return vecs

def row_reduce(M,variable_name = 0, assuming = [],pivot_row = 0):
    R = M.copy()
    m, n = R.shape

    for pivot_col in range(0,n):
        # Find pivot
        pivot = None
        for r in range(pivot_row, m):
            if R[r, pivot_col] != "0":
                if R[r, pivot_col] == "1":
                    pivot = r
                    break
        if pivot is None: #no pivot with leading 1
            for r in range(pivot_row, m):
                if R[r, pivot_col] != "0":
                    #expression must be in x0,... or an "a"
                    #we break into 2 cases
                    leading_term = R[r, pivot_col]
                    
                    R0 = R.copy()
                    R0[r, pivot_col] = "0"
                    if R[r, pivot_col] == "a":
                        R0_output = row_reduce(R0,variable_name,assuming,pivot_row)
                    else:
                        R0_output = row_reduce(R0,variable_name,assuming + [[leading_term,"0"]],pivot_row)
                    
                    
                    Rneq0 = R.copy()
                    if R[r, pivot_col] == "a":
                        leading_term = "x" + str(variable_name)
                        variable_name = variable_name + 1
                    for entry in Rneq0[r].flat:
                        if entry != "0" and entry != "a":
                            entry = "(" + entry + ")/(" + leading_term + ")"
                    Rneq0[r,pivot_col] = "1"
                    Rneq0_output = row_reduce(Rneq0,variable_name,assuming + [[leading_term,"1"]],pivot_row)#1 signifies neq 0
                    
                    return R0_output +  Rneq0_output
                continue
        
        if pivot is None:
            continue

        # Swap rows
        if pivot != pivot_row:
            R[[pivot_row, pivot]] = R[[pivot, pivot_row]]


        # Eliminate all other rows
        for r in range(m):
            if r != pivot_row and R[r,pivot_col] != "0":
                for c in range(pivot_col+1,n):
                    if R[pivot_row,c] != "0":
                        if R[r,c] == "a":
                            continue #if entry is free it remains free
                        else:
                            #we need to give variable names if they don't yet exist
                            if R[r, pivot_col] == "a":#we know its not 0
                                R[r, pivot_col] = "-x" + str(variable_name)#the negative sign reduces the number of negative signs we will get by the end
                                variable_name = variable_name + 1
                            if R[pivot_row, c] == "a":#we know its not 0
                                R[pivot_row, c] = "x" + str(variable_name)
                                variable_name = variable_name + 1
                            R[r,c] = "(" + R[r,c] + "-" + R[r, pivot_col] + "*" + R[pivot_row,c] +")"
                R[r,pivot_col] = "0"
                    
        pivot_row += 1
        if pivot_row == m:
            break

    return [tuple([R,variable_name,assuming])]

def all_RREFs(delta):
    rrefs = []
    vecs = normal_vectors(delta)
    
    vecs_immutable = [tuple(map(tuple, M.tolist())) for M in vecs]
    
    for A in Subsets(vecs_immutable):
        if len(A)==0:
            continue
        print("Number of hyperplane classes being intersected",len(A)," / ",len(vecs_immutable))
        rrefs = rrefs + row_reduce(numpy.vstack([numpy.matrix(M,dtype=object) for M in A]))
    return rrefs
            
def print_matrix_aligned(A): #takes numpy matrix and prints
    m, n = A.shape
    # Convert numpy array to ordinary Python lists
    A = A.tolist()

    # Find width of each column
    widths = [max(len(str(row[j])) for row in A) for j in range(n)]

    # Print rows
    for row in A:
        print("  ".join(f"{str(SR(x)):>{widths[j]}}" for j, x in enumerate(row)))
        
def get_variable_terms(R):
    output = []
    for entry in R.flat:
        if entry != "0" and entry != "1" and entry != "a":
            
            output.append(entry)
    return output
        
def get_num_free_entries(R):
    output = 0
    for entry in R.flat:
        if entry == "a":
            output = output + 1
    return output
def num_of_evals(var_terms_SR,num_vars,assumptions,q):
    K = GF(q)
    R = PolynomialRing(K, [f'x{i}' for i in range(num_vars)])
    var_terms = [R(f) for f in var_terms_SR]
    
    
    V = VectorSpace(K, num_vars)
    
    S = set()
    
    for v in V:
        subs_dict = {
            SR.var(f'x{i}'): v[i]
            for i in range(num_vars)
        }
        #verify assumptions are satisfied
        assumptions_satisfied = true
        for f,val in assumptions:
            eval_f = R(f)(*v)
            if eval_f == 0 and val == "1":
                assumptions_satisfied = false
                break
            if eval_f != 0 and val == "0":
                assumptions_satisfied = false
                break
        if assumptions_satisfied:
            S.add(tuple(entry(*v) for entry in var_terms))
    return len(S)





#The "percent of predicted" refers to the assumption that if one rref with a given distribution of 0 entries exists then all such rrefs must appear
#The program will likely not terminate

#D=SimplicialComplex([[0,1,2],[0,2,3],[1,2,3]])
D=SimplicialComplex([[1,2,7],[1,6,7],[2,7,8],[2,8,3],[3,4,8],[4,8,0],[4,5,0],[5,6,0],[6,7,0],[7,8,0]])

all_the_rrefs = all_RREFs(D)

print("RREFs shapes Calculated")

for qfixed in range(2,100):
    if is_prime_power(qfixed):
        print("q=",qfixed)
        existent = 0
        total = 0
        for rref,num_vars,assumptions in all_the_rrefs:
            #print_matrix_aligned(rref)
            #print()
            free_entries = get_num_free_entries(rref)
            var_terms = get_variable_terms(rref)
            
            existent = existent + num_of_evals(var_terms,num_vars,assumptions,qfixed)*qfixed^free_entries
            total = total + qfixed^len(var_terms)*qfixed^free_entries

        print("\tPercent of predicted found: ",float((total-existent)/total))#should only be used to check if things change.  To count true RREF the sum would need to include factors for the a
