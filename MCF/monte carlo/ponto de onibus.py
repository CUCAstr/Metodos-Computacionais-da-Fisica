import numpy as np
import numpy.random as rd

tot_esp = 0
n = 0
onibus = 1000
pessoa = 1000


for i in range (onibus):

    i += 1
    x = rd.rand()
    
    for j in range (pessoa):
        
        j += 1
        n += 1
        y = rd.rand()
        
        if y == 0:
            tot_esp += 0

        elif y <= x:
            tot_esp += x-y

        else:
            tot_esp += 1-y
        
        
print(f"\nTempo médio de espera para {onibus*pessoa} iterações: {tot_esp/n}\n")