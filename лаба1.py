# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 10:52:20 2022

@author: RIPLECS
"""


import numpy as np
import random
from scipy.stats import chi2
from scipy import integrate 

gamma = 0.99
z_gamma = 2.575
epsilon = 0.01

def factorial(x):
    return 1 if x< 2 else x*factorial(x - 1)
'''
A.	Побудувати довірчий інтервал для математичного сподівання   у припущенні, 
що спостерігаються в.в.  , які мають нормальний розподіл, але дисперсія   невідома.

B.	Побудувати довірчий інтервал для математичного сподівання   у припущенні, 
що спостерігаються в.в.  , розподіл яких невідомий.

C.	Побудувати довірчий інтервал для дисперсії   у припущенні, що спостерігаються 
в.в.  , які мають нормальний розподіл.

'''
'''

for n in (10**2, 10**3, 10**4):
    X = np.random.normal(0, 1, n)
    print(f'\nn = {n}')
    print('\nПрипущення: x_i мають нормальний розподіл, дисперсія невідома.')
    a = sum(X)/n
    sigma = sum((x - a)**2 for x in X)/n
    print('a^ = ', a)
    interval = (a - z_gamma*sigma/np.sqrt(n-1), a + z_gamma*sigma/np.sqrt(n-1))
    print(f'a ∈ ({interval[0]}, {interval[1]})')
    print('Довжина інтервалу:', interval[1]-interval[0])
    print('\nПрипущення: розподіл x_i невідомий.')
    print('a^ = ', a)
    interval = (a - z_gamma*sigma/np.sqrt(n), a + z_gamma*sigma/np.sqrt(n))
    print(f'a ∈ ({interval[0]}, {interval[1]})')
    print('Довжина інтервалу:', interval[1]-interval[0])
    print('\nПрипущення: розподіл x_i нормальний.')
    print('sigma^ = ', sigma)
    z_1 = chi2.ppf(1 - gamma, n - 1)
    z_2 = chi2.ppf(gamma, n - 1)
    interval = (n*sigma/z_2, n*sigma/z_1)
    print(f'sigma ∈ ({interval[0]}, {interval[1]})')
    print('Довжина інтервалу:', interval[1]-interval[0])
    
'''
m = 10

def prob(u):
    return np.exp(-u)*u**(m - 1)/((1 + u)*factorial(m-1))

print('A. Точне значення імовірності:')
print('Q_m = ', integrate.quad(prob, 0, np.infty)[0])
print(f'\nB. Метод Монте-Карло:')
X = [random.expovariate(1) for i in range(m)]
eta = 1/random.uniform(0, 1) - 1
n = 2
q = 1 if eta > sum(X) else 0
sum_q = q
sum_sq_q = q
while True:
    X = [random.expovariate(1) for i in range(m)]
    eta = 1/random.uniform(0, 1) - 1
    q = 1 if eta > sum(X) else 0
    sum_q += q
    sum_sq_q += q**2
    Q = sum_q/n
    s = (sum_sq_q - n*Q**2)/(n - 1)
    n += 1
    if Q != 0:
        if n >= (z_gamma**2*s)/(epsilon**2*Q**2):
            break
print('Q_m = ', Q)
print('n* = ', n)

print(f'\nC. Метод 2:')
X = [random.expovariate(1) for i in range(m)]
n = 2
q = 1/(sum(X) + 1)
sum_q = q
sum_sq_q = q
while True:
    X = [random.expovariate(1) for i in range(m)]
    q = 1/(sum(X) + 1)
    sum_q += q
    sum_sq_q += q**2
    Q = sum_q/n
    s = (sum_sq_q - n*Q**2)/(n - 1)
    n += 1
    if Q != 0:
        if n >= (z_gamma**2*s)/(epsilon**2*Q**2):
            break
print('Q_m = ', Q)
print('n* = ', n)


def summ(x):
    res = 0
    for i in range(len(x)):
        res += x[i]**i*np.exp(-x[i])/factorial(i)
    return res

print(f'\nC. Метод 3:')
X = [1/random.uniform(0, 1) - 1 for i in range(m)]
n = 2
q = 1 - summ(X)
sum_q = q
sum_sq_q = q
while True:
    X = [random.expovariate(1) for i in range(m)]
    q = 1 - summ(X)
    sum_q += q
    sum_sq_q += q**2
    Q = sum_q/n
    s = (sum_sq_q - n*Q**2)/(n - 1)
    n += 1
    if Q != 0:
        if n >= (z_gamma**2*s)/(epsilon**2*Q**2):
            break
print('Q_m = ', Q)
print('n* = ', n)