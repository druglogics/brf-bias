# from a contingency table, calculate cohen's kappa
# see https://en.wikipedia.org/wiki/Cohen%27s_kappa
# a b
# c d

a = 6
b = 1
c = 1
d = 4
total = a + b + c + d

# observed, percent agreement
p_obs = (a+d)/total

p_yes = ((a+b)/total) * ((a+c)/total)
p_no  = ((c+d)/total) * ((b+d)/total)
p_exp = p_yes + p_no # expected

k = (p_obs - p_exp)/(1 - p_exp)
print(k)
