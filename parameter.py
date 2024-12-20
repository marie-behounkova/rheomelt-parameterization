#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import math
import numpy as np
import scipy.special as spec #erf
from scipy import optimize as opt
import matplotlib as mpl
import matplotlib.pyplot as plt


#=======================================
# DEFAULT VALUES
B_einstein=2.5 #einstein coefficient
poiss = 0.25 #poisson ratio
gamma = 5. #thickness of transition
phi_c = 0.3 #critical melt fraction
#numerical parameters
epsilon = 1e-8


#order of variables:
#0 liquidus value
#1 delta
#2 xi
#3 gamma
#4 phi_star


#parameters from Kervazo et al., 2021, phistar in paper meant as 1-phistar due to a typo
#for comparison
plot_kervazo    = False
par_eta_kervazo = [1.,25.7,1.17e-9,5.,1-0.569]
par_mu_kervazo  = [10.,2.10,7.08e-7,5.,1-0.597]
par_k_kervazo   = [1e9,2.62,.102,5.,1-0.712]


#parameters of eta
eta_sol     = 1e19 # solidus value in Pas
eta_liq     = 1. #par_eta_kervazo[0] # liquidus value in Pas, as in Kervazo
eta_gamma   = gamma #par_eta_kervazo[3]
eta_phicrit = phi_c
eta_exp     = 30 #exponential decrease, Kervazo A.1

#parameters of mu
mu_sol     = 60e9 # solidus value in Pa
mu_liq     = 10. #par_mu_kervazo[0] # liquidus value in Pa, as in Kervazo
mu_gamma   = gamma #par_mu_kervazo[3]
mu_phicrit = phi_c
mu_B       = (40.-24*poiss)/15. #kervazo A.2

#parameters of K
k_sol      = 200e9
# k_sol      = 2*mu_sol*(1+poiss)/3./(1-2*poiss) # solidus value in Pa, #computed from combination of shear modulus and poisson ratio
k_liq      = 1e9 #par_k_kervazo[0] # liquidus value in Pa, as in Kervazo
k_gamma    = gamma #par_k_kervazo[3]
k_phicrit  = phi_c
k_B        = (5-4*poiss)/3./(1-2*poiss) #kervazo A.3


#================================================================================================
#-------------------------------------------------------
def mavko(phi,A,B): #see also Kervazo et al. (2021), eqs A.2 and A.3
    return A/(1+phi*B)

#-------------------------------------------------------
def kervazo_eq1(phi,par):
    varl=par[0]
    delta=par[1]
    xi=par[2]
    gamma=par[3]
    phistar=par[4]

    theta=(1-phi)/(1.-phistar)
    f=kervazo_eq3(theta,xi,gamma)
    return varl*(1+theta**delta)/(1-f)**(B_einstein*(1-phistar))

#-------------------------------------------------------

def kervazo_eq3(theta,xi,gamma):
    arg = np.sqrt(math.pi)/(2.*(1-xi))*theta*(1+theta**gamma)
    return (1-xi)*spec.erf(arg)


#separation of trends
#-------------------------------------------------------
def kervazo_eq1_trend1(phi,par):
    varl=par[0]
    delta=par[1]
    xi=par[2]
    gamma=par[3]
    phistar=par[4]

    theta=(1-phi)/(1.-phistar)
    f=kervazo_eq3(theta,xi,gamma)
    #return varl*(1+theta**delta)/(1-f)**(B*(1-phistar))
    return varl*(1)/(1-f)**(B_einstein*(1-phistar))

#-------------------------------------------------------
def kervazo_eq1_trend2(phi,par):
    varl=par[0]
    phistar=par[4]
    delta=par[1]
    xi=par[2]
    gamma=par[3]
    theta=(1-phi)/(1-phistar)
    f=kervazo_eq3(theta,xi,gamma)
    #return varl*(1+theta**delta)/(1-f)**(B*(1-phistar))
    return varl*(1+theta**delta)/xi**(B_einstein*(1-phistar))

#-------------------------------------------------------
def kervazo_diff(phicrit,phi,par,abs_val,f=0.1):
    val = kervazo_eq1(phicrit,par)
    trend0 = kervazo_eq1(phi,par)
    trend2 = kervazo_eq1_trend2(phi,par)
    return trend2-trend0-f*val


#=========================================================================
#computation phi_star from phi_critical, gamma and epsilon
#iter scheme
def find_approx(v_liq,v_sol,phicrit,gamma,for_delta,mod='eta',nmax=200,phi_star_ini=0.5,eps=epsilon):
    phi_star_old = phi_star_ini #initial quess
    phi_star = phi_star_old
    if (mod=='eta'):
        delta = delta_from_exp(for_delta,phicrit,phi_star)
    else:
        delta = delta_from_mavko(for_delta,phicrit,phi_star)
    xi = xi_contrast(v_liq/v_sol,delta,phi_star)

    for i in range(nmax):
        sol=opt.root_scalar(lambda x: kervazo_diff(phicrit,phicrit,[v_liq,delta,xi,gamma,x],v_sol), method='bisect', bracket=[.0*phi_c, 0.99], x0=phi_star_old)
        phi_star = sol.root

        if (mod=='eta'):
            delta = delta_from_exp(for_delta,phicrit,phi_star)
        else:
            delta = delta_from_mavko(for_delta,phicrit,phi_star)
        xi = xi_contrast(v_liq/v_sol,delta,phi_star)

        diff = abs((phi_star_old-phi_star)/phi_star_old)
        par = [v_liq,delta,xi,gamma,phi_star]
        if (diff<eps): break

        phi_star_old=phi_star

    return par

#================================================================================================
def delta_from_exp(a,phi_crit,phistar,n=200,eps=1e-8):
    delta_old = -np.log(2*np.exp(a*phistar)-1)/np.log(1-phistar) #estimate from phistar
    #delta_old=20
    #iterative solution
    for i in range(n):
        # limit
        delta = np.log(np.exp(-a*phi_crit)-(1-phistar)**delta_old)/np.log(1-phi_crit)
        #full
        # delta = np.log(np.exp(-a*phi_crit)*((1-phistar)**delta_old+1)-(1-phistar)**delta_old)/np.log(1-phi_crit)
        if abs((delta-delta_old)/delta_old)<eps: break
        delta_old = delta
    return delta


#================================================================================================
def delta_from_mavko_ver1(B,phistar):
    return -(np.log(1.+2*B*phistar))/(np.log(1-phistar))


#================================================================================================
def delta_from_mavko(B,phi_crit,phistar,n=200,eps=1e-4):
    delta_old = -(np.log(1.+2*B*phistar))/(np.log(1-phistar))
    for i in range(n):
        # full
        delta = np.log(((1-phistar)**delta_old+1)/(B*phi_crit+1)-(1-phistar)**delta_old)/np.log(1-phi_crit)
        if abs((delta-delta_old)/delta_old)<eps: break
        delta_old = delta
    return delta

#================================================================================================
def xi_contrast(contrast,delta,phistar):
    theta = 1./(1-phistar)
    return np.power(contrast*(1+theta**delta),theta/B_einstein)

#================================================================================================
def plot(par,*args,scale='log',name=None,phicrit=None,compare = False):
    fig, (ax) = plt.subplots(nrows=1)
    #im = ax.plot(phi,(kervazo1(phi,par_eta)),c='k')
    if (scale=='lin'):
        im = ax.plot(phi,(kervazo_eq1(phi,par)),c='r',label='automatic assassement')
        if compare:
            im = ax.plot(phi,(kervazo_eq1(phi,args[0])),c='b',label='Kervazo et al., 2021')
    elif (scale=='log'):
        im = ax.plot(phi,np.log10(kervazo_eq1(phi,par)),c='r',label='automatic assassement')
        if compare:
            im = ax.plot(phi,np.log10(kervazo_eq1(phi,args[0])),c='b',label='Kervazo et al., 2021')
    if (phicrit): plt.axvline(x=phicrit)
    plt.legend()
    #plt.gca().invert_yaxis()
    if (name):
        plt.savefig(name+'.png')
    else:
        plt.show()


def format(number):
    return ("%s" % (number))
#=====================================================================================
#for plotting
phi  = np.linspace(0.,1,200,endpoint=True)

#order of variables:
#0 liquidus value
#1 delta
#2 xi
#3 gamma
#4 phi_star


try:
    mode = sys.argv[1]
except:
    print('trying default mode, i.e., one-run')
    mode = 'one-run'


saveplot = True

if saveplot:
    dir_plots = 'plots/'

    if not os.path.exists(dir_plots):
        os.makedirs(dir_plots)


if mode=='one-run':
        par_eta = find_approx(eta_liq,eta_sol,eta_phicrit,eta_gamma,eta_exp,mod='eta')
        par_mu  = find_approx(mu_liq,mu_sol,mu_phicrit,mu_gamma,mu_B,mod='mu')
        par_k   = find_approx(k_liq,k_sol,k_phicrit,k_gamma,k_B,mod='k')

        if saveplot:
            plot(par_eta,par_eta_kervazo,scale= 'log',name=dir_plots+'eta',phicrit=eta_phicrit,compare=True)
            plot(par_mu,par_mu_kervazo,scale='log',name=dir_plots+'mu',phicrit=mu_phicrit,compare=True)
            plot(par_k,par_k_kervazo,scale='log',name=dir_plots+'k',phicrit=k_phicrit,compare=True)

        print('eta',par_eta)
        print('mu',par_mu)
        print('k',par_k)

elif mode=='multi-run':

    #gamma and phi_c the same for all variables
    mu_B       = (40.-24*poiss)/15. #kervazo A.2
    k_B        = (5-4*poiss)/3./(1-2*poiss) #kervazo A.3
    if (len(sys.argv)!=3):
        print('please call python3 parameter.py multi-run data.in')
        exit()

    name_in = sys.argv[2]
    name_out_eta = 'eta-'+sys.argv[2]
    name_out_mu = 'mu-'+sys.argv[2]
    name_out_k = 'k-'+sys.argv[2]

    data_in = np.genfromtxt(name_in,skip_header=1)

    ndata = data_in.shape[0]
    eta_out = np.zeros((ndata,5))
    mu_out = np.zeros((ndata,5))
    k_out = np.zeros((ndata,5))

    print('processing ',ndata, 'datasets')
    for i in range(ndata):
        print ('dataset ',i)
        eta_liq = data_in[i,0]
        eta_sol = data_in[i,1]
        eta_exp = data_in[i,2]
        mu_liq = data_in[i,3]
        mu_sol = data_in[i,4]
        k_liq = data_in[i,5]
        k_sol = data_in[i,6]
        phicrit = data_in[i,7]
        gamma = data_in[i,8]

        par_eta = find_approx(eta_liq,eta_sol,phicrit,gamma,eta_exp,mod='eta')
        par_mu  = find_approx(mu_liq,mu_sol,phicrit,gamma,mu_B,mod='mu')
        par_k   = find_approx(k_liq,k_sol,phicrit,gamma,k_B,mod='k')

        eta_out[i,:] = par_eta
        mu_out[i,:] = par_mu
        k_out[i,:] = par_k

        iden_eta = '_eta_sol'+format(eta_sol)+'-eta_liq'+format(eta_liq)+\
            '-exp'+format(eta_exp)+\
            '-phicrit'+format(phicrit)+'-gamma'+format(gamma)
        iden_k = '_k_sol'+format(k_sol)+'-k_liq'+format(k_liq)+\
            '-phicrit'+format(phicrit)+'-gamma'+format(gamma)
        iden_mu = '_mu_sol'+format(mu_sol)+'-mu_liq'+format(mu_liq)+\
            '-phicrit'+format(phicrit)+'-gamma'+format(gamma)


        if saveplot:
            plot(par_eta,scale='log',name=dir_plots+'eta'+iden_eta,phicrit=phicrit)
            plot(par_mu,scale='log',name=dir_plots+'mu'+iden_mu,phicrit=phicrit)
            plot(par_k,scale='log',name=dir_plots+'k'+iden_k,phicrit=phicrit)

    np.savetxt(name_out_eta,eta_out,header = ' par = [v_liq,delta,xi,gamma,phi_star]')
    np.savetxt(name_out_mu,mu_out,header = ' par = [v_liq,delta,xi,gamma,phi_star]')
    np.savetxt(name_out_k,k_out,header = ' par = [v_liq,delta,xi,gamma,phi_star]')


else:
    print('unknown mode')


