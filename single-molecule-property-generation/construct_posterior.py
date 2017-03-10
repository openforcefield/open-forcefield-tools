# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:09:09 2017

@author: Bryce Manubay
"""

"""
This script represents prototyping for constructing posterior distributions
across parameter space given some initial set of data.

Author: Bryce Manubay
"""

import matplotlib as mpl

mpl.use('Agg')


from smarty.forcefield import *
import openeye
from openeye import oechem
import smarty
from smarty.utils import get_data_filename
from simtk import openmm
from simtk import unit
import numpy as np
import netCDF4 as netcdf
import collections as cl
import pandas as pd
import pymbar
from pymbar import timeseries
import glob
import sys
from smarty.forcefield import generateTopologyFromOEMol
import pdb
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#np.set_printoptions(threshold=np.inf)

import scipy as sp
import seaborn as sns

from scipy.stats import norm
from scipt.stats import multivariate_normal
import sys



sns.set_style('white')
sns.set_context('talk')

np.random.seed(123)

data = np.random.randn(20)

#ax = plt.subplot()
#sns.distplot(data, kde=False, ax=ax)
#_ = ax.set(title='Histogram of observed data', xlabel='x', ylabel='# observations');


def calc_posterior_analytical(data, x, mu_0, sigma_0):
    sigma = 1.
    n = len(data)
    mu_post = (mu_0 / sigma_0**2 + data.sum() / sigma**2) / (1. / sigma_0**2 + n / sigma**2)
    sigma_post = (1. / sigma_0**2 + n / sigma**2)**-1
    return norm(mu_post, np.sqrt(sigma_post)).pdf(x)

ax1 = plt.subplot()
x = np.linspace(-1, 1, 500)
posterior_analytical = calc_posterior_analytical(data, x, 0., 1.)
ax1.plot(x, posterior_analytical)
ax1.set(xlabel='mu', ylabel='belief', title='Analytical posterior');
sns.despine()

def sampler(data, samples=4, mu_init=.5, proposal_width=.5, plot=False, mu_prior_mu=0, mu_prior_sd=1.):
    """
    Outline:
    1)Take data and calculate observable
    2)Reweight observable to different state and calculate observable based on new state
	- smarty move in many parameters
        - will have to start without torsion moves
        - Safe moves in equilibrium bond length and angle is ~3%. For force constants ~5%.
    3)Will have to make decision here:
        a)Continue to sample in order to gather more data? -or-
        b)Attempt to create surrogate models from the data we have? What does that entail?
            i)We want a surrogate model for every observable we have, $O\left(\theta\right)$
            ii)Thus for bonds and angles; we have 4 observables as a function of however many parameters we're working with at the time
            iii)Choice of surrogate model becomes important. Start with splining though
            iv)What is the best surrogate modeling technique to use when we have very sparse data?
    4)
   
            
          

    Other things to consider:
    1)Choice of surrogate models:
        a)Splining
        b)Rich's ideas
        c)Other ideas from Michael he got at conference last week
    2)Choice of likelihood:
        a)Gaussian likelihood
        b)More general based on mean squared error
    3)Prior
        a)Start with uniforms with physically relevant bounds for given parameter
        b)Informationless priors
    """
    mu_current = mu_init
    posterior = [mu_current]
    for i in range(samples):
        # suggest new position
        mu_proposal = norm(mu_current, proposal_width).rvs()

        # Compute likelihood by multiplying probabilities of each data point
        likelihood_current = norm(mu_current, 1).pdf(data).prod()
        likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
        
        # Compute prior probability of current and proposed mu        
        prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
        prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
        
        p_current = likelihood_current * prior_current
        p_proposal = likelihood_proposal * prior_proposal
        
        # Accept proposal?
        p_accept = p_proposal / p_current
        
        # Usually would include prior probability, which we neglect here for simplicity
        accept = np.random.rand() < p_accept
        
        if plot:
            plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accept, posterior, i)
        
        if accept:
            # Update position
            mu_current = mu_proposal
        
        posterior.append(mu_current)
        
    return posterior

# Function to display
def plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accepted, trace, i):
    from copy import copy
    trace = copy(trace)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(16, 4))
    fig.suptitle('Iteration %i' % (i + 1))
    x = np.linspace(-3, 3, 5000)
    color = 'g' if accepted else 'r'
        
    # Plot prior
    prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
    prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
    prior = norm(mu_prior_mu, mu_prior_sd).pdf(x)
    ax1.plot(x, prior)
    ax1.plot([mu_current] * 2, [0, prior_current], marker='o', color='b')
    ax1.plot([mu_proposal] * 2, [0, prior_proposal], marker='o', color=color)
    ax1.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    ax1.set(ylabel='Probability Density', title='current: prior(mu=%.2f) = %.2f\nproposal: prior(mu=%.2f) = %.2f' % (mu_current, prior_current, mu_proposal, prior_proposal))
    
    # Likelihood
    likelihood_current = norm(mu_current, 1).pdf(data).prod()
    likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
    y = norm(loc=mu_proposal, scale=1).pdf(x)
    sns.distplot(data, kde=False, norm_hist=True, ax=ax2)
    ax2.plot(x, y, color=color)
    ax2.axvline(mu_current, color='b', linestyle='--', label='mu_current')
    ax2.axvline(mu_proposal, color=color, linestyle='--', label='mu_proposal')
    #ax2.title('Proposal {}'.format('accepted' if accepted else 'rejected'))
    ax2.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    ax2.set(title='likelihood(mu=%.2f) = %.2f\nlikelihood(mu=%.2f) = %.2f' % (mu_current, 1e14*likelihood_current, mu_proposal, 1e14*likelihood_proposal))
    
    # Posterior
    posterior_analytical = calc_posterior_analytical(data, x, mu_prior_mu, mu_prior_sd)
    ax3.plot(x, posterior_analytical)
    posterior_current = calc_posterior_analytical(data, mu_current, mu_prior_mu, mu_prior_sd)
    posterior_proposal = calc_posterior_analytical(data, mu_proposal, mu_prior_mu, mu_prior_sd)
    ax3.plot([mu_current] * 2, [0, posterior_current], marker='o', color='b')
    ax3.plot([mu_proposal] * 2, [0, posterior_proposal], marker='o', color=color)
    ax3.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    #x3.set(title=r'prior x likelihood $\propto$ posterior')
    ax3.set(title='posterior(mu=%.2f) = %.5f\nposterior(mu=%.2f) = %.5f' % (mu_current, posterior_current, mu_proposal, posterior_proposal))
    
    if accepted:
        trace.append(mu_proposal)
    else:
        trace.append(mu_current)
    ax4.plot(trace)
    ax4.set(xlabel='iteration', ylabel='mu', title='trace')
    plt.tight_layout()
    #plt.legend()
    
np.random.seed(123)

posterior = sampler(data, samples=8000, mu_init=2.)
fig, ax = plt.subplots()
ax.plot(posterior)
_ = ax.set(xlabel='sample', ylabel='mu');

sys.exit()
ax1 = plt.subplot()

sns.distplot(posterior[500:], ax=ax1, label='estimated posterior')
x = np.linspace(-.5, .5, 500)
post = calc_posterior_analytical(data, x, 0, 1)
ax1.plot(x, post, 'g', label='analytic posterior')
_ = ax1.set(xlabel='mu', ylabel='belief');
ax1.legend();


#Each piece of experimental data that we use as evidence (i.e. bond length for 
#one smirks vs any other or variance of a specific bond angle, etc.) has a 
#separate likelihood to evaluate and the product of all of those likelihoods
#the overall likelihood for a particular proposal. So we propose a change in 
#all parameters in one go and then evaluate that full change with all of the 
#evidence we have available.
