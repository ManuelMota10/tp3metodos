import numpy as np
def diferencas_divididas(x_points, y_points):
    n = len(x_points)
    coef = np.zeros([n, n])
    # A primeira coluna são os valores de y
    coef[:, 0] = y_points
    
    for j in range(1, n):
        for i in range(n - j):
            coef[i][j] = (coef[i + 1][j - 1] - coef[i][j - 1]) / (x_points[i + j] - x_points[i])
            
    return coef[0, :] # Retorna a primeira linha (coeficientes a0, a1, ...)

def newton_slides(coef, x_points, x_eval):
    n = len(coef)
    p = coef[0]

    for k in range(1, n):
        termo = coef[k]
        for j in range(k):
            termo *= (x_eval - x_points[j])
        p += termo
    return p
    
def newton_poly(coef, x_points, x_eval):
    #Avalia o polinómio de Newton num ponto x_eval.
    n = len(x_points) - 1 
    p = coef[n]
    for k in range(1, n + 1):
        p = coef[n - k] + (x_eval - x_points[n - k]) * p
    return p

x_ex = [2.1, 2.2, 2.3, 2.4, 2.5]
y_ex = [0.32222, 0.34242, 0.36173, 0.38021, 0.39794]

coeficientes = diferencas_divididas(x_ex, y_ex)
print("Coeficientes (Diferenças Divididas):", coeficientes)

val = newton_poly(coeficientes, x_ex, 2.15)
print(f"Newton: Valor estimado em 2.15 = {val:.5f}")
# Esperado próximo de 0.33243 (Slide 102)