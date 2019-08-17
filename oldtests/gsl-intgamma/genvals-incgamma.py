import random
import math

jmu = 2
mu = 5

for i in range(10000):
    x = abs(random.expovariate(1/mu))
    j = math.floor(abs(random.expovariate(1/jmu)))
    print(j, x)
    
            
