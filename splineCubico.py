import numpy as np

def spline_cubico_natural(x_points, y_points, x_eval):

    N = len(x_points) - 1
    h = np.diff(x_points) # h_i = x_i - x_{i-1}
    
    # 1. Montagem do Sistema Linear A * M = b para encontrar M (segundas derivadas)
    # Dimensão do sistema interno é (N-1) pois M_0=0 e M_N=0
    A = np.zeros((N-1, N-1))
    b = np.zeros(N-1)
    
    for i in range(1, N):
        idx = i - 1 
        
        # Lado direito da equação (diferenças divididas)
        term1 = (y_points[i+1] - y_points[i]) / h[i]
        term2 = (y_points[i] - y_points[i-1]) / h[i-1]
        b[idx] = term1 - term2
        
        # Diagonal Principal: (h_i + h_{i+1}) / 3
        A[idx, idx] = (h[i-1] + h[i]) / 3
        
        # Diagonais secundárias
        if idx > 0:
            A[idx, idx-1] = h[i-1] / 6
        if idx < N-2:
            A[idx, idx+1] = h[i] / 6
            
    # Resolver para M (nos pontos interiores)
    M_interior = np.linalg.solve(A, b)
    
    # Adicionar condições de fronteira natural M0=0, Mn=0 (Slide 113)
    M = np.concatenate(([0], M_interior, [0]))
    
    # 2. Avaliação do Spline no ponto x_eval
    # Encontrar em qual intervalo [x_{i-1}, x_i] o ponto x_eval está
    for i in range(1, N + 1):
        if x_points[i-1] <= x_eval <= x_points[i]:
            xi_1 = x_points[i-1]
            xi = x_points[i]
            fi_1 = y_points[i-1]
            fi = y_points[i]
            hi = h[i-1]
            Mi_1 = M[i-1]
            Mi = M[i]
            
            # Fórmula do polinómio S_i(x)
            term_M1 = Mi_1 * (xi - x_eval)**3 / (6 * hi)
            term_M2 = Mi * (x_eval - xi_1)**3 / (6 * hi)
            term_F1 = (fi_1 - Mi_1 * hi**2 / 6) * (xi - x_eval) / hi
            term_F2 = (fi - Mi * hi**2 / 6) * (x_eval - xi_1) / hi
            
            return term_M1 + term_M2 + term_F1 + term_F2
            
    return None # Fora do intervalo

# --- Exemplo Slide 114 ---
x_spline = [1, 2, 3, 4]
y_spline = [1, 1/2, 1/3, 1/4]

# Avaliar em 1.5
val_spline = spline_cubico_natural(x_spline, y_spline, 1.5)
print(f"Spline Cúbico Natural em 1.5 = {val_spline:.5f}")