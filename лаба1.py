# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 10:52:20 2022

@author: RIPLECS
"""


import numpy as np
import random
from scipy.stats import chi2
from scipy import integrate 
from sympy import factorial
import math

gamma = 0.99
z_gamma = 2.575
epsilon = 0.01
epsilon2 = epsilon**2
z_gamma2 = z_gamma**2 



print('\n____________ЗАВДАННЯ 1____________')

for n in (10**2, 10**3, 10**4):
    X = np.random.normal(0, 1, n)
    print(f'\nn = {n}')
    print('\nПрипущення: x_i мають нормальний розподіл, дисперсія невідома.')
    a = sum(X)/n
    sigma = sum((x - a)**2 for x in X)/n
    print('a^ = ', a)
    sub = z_gamma*sigma/np.sqrt(n - 1)
    interval = (a - sub, a + sub)
    print(f'a ∈ ({interval[0]}, {interval[1]})')
    print('Довжина інтервалу:', interval[1] - interval[0])
    print('\nПрипущення: розподіл x_i невідомий.')
    print('a^ = ', a)
    sub = z_gamma*sigma/np.sqrt(n)
    interval = (a - sub, a + sub)
    print(f'a ∈ ({interval[0]}, {interval[1]})')
    print('Довжина інтервалу:', interval[1] - interval[0])
    print('\nПрипущення: розподіл x_i нормальний.')
    print('sigma^ = ', sigma)
    z_1 = chi2.ppf(1 - gamma, n - 1)
    z_2 = chi2.ppf(gamma, n - 1)
    interval = (n*sigma/z_2, n*sigma/z_1)
    print(f'sigma ∈ ({interval[0]}, {interval[1]})')
    print('Довжина інтервалу:', interval[1] - interval[0])
    


def true_prob(u):
    return math.exp(-u)*u**(m_ - 1)/((1 + u)*factorial(m_ - 1))


def culc_prob(f, m, n = 2):
    if f == f3 and m == 1:
        print('Метод не застосовний!')
        return '-', '-'
    ksi = [random.expovariate(1) for i in range(m)]
    eta = 1/random.uniform(0, 1) - 1
    q = f(eta, ksi)
    sum_q = q
    sum_sq_q = q**2
    while True:
        ksi = [random.expovariate(1) for i in range(m)]
        eta = 1/random.uniform(0, 1) - 1
        q = f(eta, ksi)
        sum_q += q
        sum_sq_q += q**2
        Q = sum_q/n
        sigma = (sum_sq_q - sum_q**2/n)/(n - 1)
        if Q == 0:
            continue
        if n >= (z_gamma2*sigma)/(epsilon2*Q**2) and n > 3:
            break
        n += 1
    return Q, n

def f1(eta, ksi):
    return int(eta > sum(ksi))

def f2(eta, ksi):
    return 1/(sum(ksi) + 1)

def f3(eta, ksi):
    s = sum(ksi[:-1])
    return s/((1 + s)*(len(ksi) - 1))
 

print('\n____________ЗАВДАННЯ 2____________')

for m_ in(1, 50, 100):
    print('\nm = ', m_)
    if m_ < 60:
        print('\nA. Точне значення імовірності:')
        print('Q_m = ', integrate.quad(true_prob, 0, np.infty)[0])
    print('\nB. Метод 1:')
    start_time = datetime.now()
    Q_m, n_, disp = culc_prob(f1, m_)
    print('Q_m = ', Q_m)
    print('n* = ', n_)
    print('sigma = ', disp)
    print(datetime.now() - start_time)
    print('\nC. Метод 2:')
    Q_m, n_, disp  = culc_prob(f2, m_)
    print('Q_m = ', Q_m)
    print('n* = ', n_)
    print('sigma = ', disp)
    print('\nC. Метод 3:')
    Q_m, n_, disp  = culc_prob(f3, m_)
    print('Q_m = ', Q_m)
    print('n* = ', n_)
    print('sigma = ', disp)
