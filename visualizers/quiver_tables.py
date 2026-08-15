import numpy as np
from quivers.knot_quiver import colored_homfly_vectors_and_quiver
from quivers.tail_quivs import colored_homfly_head_vectors_and_quiver, colored_homfly_tail_vectors_and_quiver

def create_quiver_table(n, path="table"):
    """Writes to text file all quivers of rational knots up to numerator n."""
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
    """Writes to text file A vectors of rational knots up to numerator n."""
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
    """Writes to text file S vectors of rational knots up to numerator n."""
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
    """Writes to text file all quiver diagonals of rational knots up to numerator n."""
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

def diff_quiv(i, j):
    _, _, Q = colored_homfly_vectors_and_quiver(i, j)
    D = Q[1:, 1:] - Q[1:, :-1] - Q[:-1, 1:] + Q[:-1, :-1]
    return D

def create_diff_table(n, path="diff_table"):
    """Writes to text file all difference quivers of rational knots up to numerator n."""
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    Q = diff_quiv(i,j)
                    Q_numpy = np.array(Q.tolist(), dtype=int)
                    f.write(f"D matrix for K_{i}/{j}\n")
                    np.savetxt(f, Q_numpy, fmt="%3d", delimiter=" ")
                    f.write("\n")
                    
                except:
                    pass   

def create_head_table(n, path="head_table"):
    """Writes to text file all head subquivers of rational knots up to numerator n."""
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    S, A, Q = colored_homfly_head_vectors_and_quiver(i,j)
                    S = np.array(S)
                    A = np.array(A)
                    Q_numpy = np.array(Q.tolist(), dtype=int)
                    f.write(f"Head subquiver data for K_{i}/{j}\n")
                    f.write("S:\n")
                    np.savetxt(f, [S], fmt="%3d", delimiter=" ")
                    f.write("\nA:\n")
                    np.savetxt(f, [A], fmt="%3d", delimiter=" ")
                    f.write("\nQ:\n")
                    np.savetxt(f, Q_numpy, fmt="%3d", delimiter=" ")
                    f.write("\n")
                    
                except:
                    pass   

def create_tail_table(n, path="tail_table"):
    """Writes to text file all tail subquivers of rational knots up to numerator n."""
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    S, A, Q = colored_homfly_tail_vectors_and_quiver(i,j)
                    S = np.array(S)
                    A = np.array(A)
                    Q_numpy = np.array(Q.tolist(), dtype=int)
                    f.write(f"Tail subquiver data for K_{i}/{j}\n")
                    f.write("S:\n")
                    np.savetxt(f, [S], fmt="%3d", delimiter=" ")
                    f.write("\nA:\n")
                    np.savetxt(f, [A], fmt="%3d", delimiter=" ")
                    f.write("\nQ:\n")
                    np.savetxt(f, Q_numpy, fmt="%3d", delimiter=" ")
                    f.write("\n")
                    
                except:
                    pass   
create_S_table(50)