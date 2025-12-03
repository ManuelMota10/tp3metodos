def lagrange_interpolation(x_points, y_points, x_eval):
    n = len(x_points)
    result = 0.0
    
    for k in range(n):
        # Calcula o termo l_k(x)
        l_k = 1.0
        for i in range(n):
            if i != k:
                # Fórmula do produto (Slide 90)
                l_k *= (x_eval - x_points[i]) / (x_points[k] - x_points[i])
        
        # Soma ponderada
        result += l_k * y_points[k]
        
    return result

# --- Exemplo do Slide 86 (Estimar e^0.826) ---
x_dados = [0.80, 0.81, 0.82, 0.83, 0.84] 
y_dados = [2.225541, 2.247908, 2.270500, 2.293319, 2.316367]
x_alvo = 0.826

valor_estimado = lagrange_interpolation(x_dados, y_dados, x_alvo)
print(f"Lagrange (n=1): Valor estimado em {x_alvo} = {valor_estimado:.10f}")
# Valor esperado slide 86: ~2.284164