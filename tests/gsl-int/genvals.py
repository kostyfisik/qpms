import random
import math

nmu = 2
mu = 1e-5
ksigma = 1e7

for i in range(10000):
    R = abs(random.expovariate(1/mu))
    k = random.normalvariate(0,ksigma)
    n = math.ceil(abs(random.expovariate(1/nmu)))
    eta = abs(k/2/6**0.5)# abs(random.expovariate(8/k))
    print(n, k, R, eta)
    
            
