import numpy as np
import numpy.random as npr
import matplotlib.pyplot as mp
import pickle
import pystan as ps
import corner as co
import chsi_mod as chsi


def template(x, nu, x_0=0.0, amp=1.0):
	dx = x - x_0
	return amp * (dx ** 1.5 * np.sin(2.0 * np.pi * nu * dx))


# settings
recompile = False
stan_constrain = False
n_samples = 10000
n_chains = 4
n_par = 2

# plotting settings
lw = 1.5
mp.rc('font', family = 'serif')
mp.rcParams['text.latex.preamble'] = [r'\boldmath']
mp.rcParams['axes.linewidth'] = lw
mp.rcParams['lines.linewidth'] = lw

# generate a template. simplest test case has h = 1 throughout
n_data = 50
nu = 2.0 / n_data
x = np.arange(n_data)
t = template(x, nu)
'''
x = np.arange(n_data) * 1.1 ** np.arange(n_data)
x = x / np.max(x) * n_data
y = x ** 1.5 * np.sin(2.0 * np.pi * nu * x)
'''
test = chsi.Interpolator(x, t)

# generate noisy observations with some lag
n_obs = n_data / 2
noise_sig = 5.0
amp = 2.0
lag = 0.25 / nu
x_obs = 15.0 + lag + np.linspace(0.0, 1.0 / nu, n_obs)
t_obs = template(x_obs, nu, lag, amp) + \
		npr.randn(n_obs) * noise_sig


'''
mp.plot(x, t)
mp.plot(x_obs, t_obs)
mp.show()

# now let's "sample" the correct lag and interpolate the template onto
# these points
x_int = x_obs - lag
t_int = amp * test.interpolate(x_int)
mp.plot(x, t)
mp.plot(x_obs, t_obs)
mp.scatter(x_int, t_int)
mp.show()


mp.plot(x_int, t_obs - t_int, 'b')
mp.axhline(noise_sig, color='grey')
mp.axhline(-noise_sig, color='grey')
mp.show()

n_lags = 20
lags = np.linspace(0.5, 1.5, n_lags) * lag
chisq = np.zeros(n_lags)
for i in range(n_lags):
	x_int = x_obs - lags[i]
	print np.min(x_int), np.max(x_int)
	t_int = amp * test.interpolate(x_int)
	chisq[i] = np.sum(((t_obs - t_int) / noise_sig) ** 2)
mp.plot(lags, chisq)
mp.show()
'''

# what do we actually need to do?
# 1) generate template with zero lag
# 2) observe without noise for now, with regular sampling, unit amplitude
# 3) generate dependent observations:
#     - noisy
#     - lagged
#     - linear function of template for now (amp * template(lag))
# 4) fit!
# 5) VERSIONS
#     - Python interpolant, Stan interpolate
#     - Python interpolant, C++ interpolate
#     - do these have different speeds? acceptance?


# compile Stan model
base = 'chsi'
if recompile:
    stan_model = ps.StanModel(base + '.stan')
    with open(base + '_model.pkl', 'wb') as f:
        pickle.dump(stan_model, f)
else:
    try:
        with open(base + '_model.pkl', 'rb') as f:
            stan_model = pickle.load(f)
    except EnvironmentError:
        print 'ERROR: pickled Stan model (' + base + \
              '_model.pkl) not found. Please set recompile = True'
        exit()

# set up stan inputs and sample
stan_data = {'n_obs': n_obs, 'x_obs': x_obs, 'y_obs': t_obs, \
             'noise_sig': noise_sig, 'n_knot': n_data, \
             'x_knot': test.x, 'y_knot': test.y, 'd_knot': test.d, \
             'c_knot': test.c, 'b_knot': test.b}
if stan_constrain:
    stan_seed = 23102014
else:
    stan_seed = None
fit = stan_model.sampling(data = stan_data, \
                          iter = n_samples, \
                          chains = n_chains, \
                          seed = stan_seed)
print fit

# plot fit
par_labels = [r"${\rm amp}$", r"${\rm lag}$"]
par_vals = [lag, amp]
raw_samples = fit.extract(permuted = False)
samples = np.zeros((n_chains * n_samples / 2, n_par))
for i in range(n_chains):
	samples[i * n_samples / 2: (i + 1) * n_samples / 2] = \
		raw_samples[:, i, 0: 2]
fig = co.corner(samples, labels=par_labels, \
				quantiles=[0.16, 0.84], show_titles=True, \
				title_kwargs={"fontsize": 12})
axes = fig.axes
for i in range(n_par):
	axes[i * n_par + i].axvline(par_vals[1 - i], color='red')
	for j in range(i):
		axes[i * n_par + j].axvline(par_vals[i], color='red')
		axes[i * n_par + j].axhline(par_vals[j], color='red')
mp.savefig(base + '_fit.pdf', bbox_inches='tight')
mp.show()
