import numpy as np
from quivers.knot_quiver import colored_homfly_vectors_and_quiver

def create_quiver_table(n, path="table"):
    # writes to text file all quivers of knots up to numerator n
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    S, A, Q = colored_homfly_vectors_and_quiver(i,j)
                    Q_numpy = np.array(Q.tolist(), dtype=int) 
                    f.write(f"Quiver for K_{i}/{j}\n")
                    np.savetxt(f, Q_numpy, fmt="%3d", delimiter=" ")
                    f.write("\n")
                except:
                    pass

def create_A_table(n, path="a_table"):
    # writes to text file all quivers of knots up to numerator n
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    S, A, Q = colored_homfly_vectors_and_quiver(i, j)
                    A = np.array(A)
                    f.write(f"A vector for K_{i}/{j}\n")
                    np.savetxt(f, [A], fmt="%3d", delimiter=" ")
                    f.write("\n")
                
                except:
                    pass    

def create_S_table(n, path="s_table"):
    # writes to text file all quivers of knots up to numerator n
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    S, A, Q = colored_homfly_vectors_and_quiver(i, j)
                    S = np.array(S)
                    f.write(f"S vector for K_{i}/{j}\n")
                    np.savetxt(f, [S], fmt="%3d", delimiter=" ")
                    f.write("\n")
                
                except:
                    pass  

def create_diag_table(n, path="diag_table"):
    # writes to text file all quivers of knots up to numerator n
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    S, A, Q = colored_homfly_vectors_and_quiver(i, j)
                    Q_numpy = np.diag(np.array(Q.tolist(), dtype=int))  
                    f.write(f"Q diagonal vector for K_{i}/{j}\n")
                    np.savetxt(f, [Q_numpy], fmt="%3d", delimiter=" ")
                    f.write("\n")
                    
                except:
                    pass   