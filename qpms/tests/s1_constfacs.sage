def s1_constfacs(m, n):
    if (m+n) % 2 == 1:
        return []
    s1consts = list()
    for j in range((n-abs(m))+1):
        s1consts.append(
            -I**(n+1)/2 * (-1)**((n+m)/2) * sqrt((2*n+1)*factorial(n-m)*factorial(n+m))
            * (-1)**j / 2**(n-2*j)
            / (factorial(j) * factorial((n-m)/2-j) * factorial((n+m)/2-j)) / 2**(2*j-1)
            )
    return s1consts

for l in range(1, 11):
    for m in range (-l, l+1):
        print(l, m, [N(x, prec=66) for x in s1_constfacs(m,l)])

