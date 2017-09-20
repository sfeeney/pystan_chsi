functions {
    vector interpolate(int n_int, vector x_int, int n_knot, 
                       vector x_knot, vector y_knot, 
                       vector d_knot, vector c_knot, 
                       vector b_knot) {
        vector[n_int] y_int;
        for(j in 1: n_int) {
            int i;
            real s;
            i = 0;
            if(x_int[j] < x_knot[i+1] || x_int[j] > x_knot[n_knot]) {
                reject("Requested interpolation (at ", x_int[j], 
                       ") out of bounds (", x_knot[i+1], "<x<", 
                       x_knot[n_knot], ")");
            }
            while(x_int[j] > x_knot[i+1]) {
                i = i + 1;
            }
            s = x_int[j] - x_knot[i];
            y_int[j] = y_knot[i] +  
                       s * (d_knot[i] + 
                            s * (c_knot[i] + 
                                 s * b_knot[i]));
        }
        return y_int;
    }   
}
data {
    int<lower=0> n_obs;
    vector[n_obs] x_obs;
    vector[n_obs] y_obs;
    real noise_sig;
    int<lower=0> n_knot;
    vector[n_knot] x_knot;
    vector[n_knot] y_knot;
    vector[n_knot] d_knot;
    vector[n_knot-1] c_knot;
    vector[n_knot-1] b_knot;
}
parameters {
    real<lower=0,upper=5> amp;
    real<lower=3.0,upper=9.0> lag;
}
model {
    vector[n_obs] x_int;
    vector[n_obs] y_int;
    x_int = x_obs - lag;
    y_int = amp * interpolate(n_obs, x_int, n_knot, x_knot, y_knot, 
                              d_knot, c_knot, b_knot);
    y_obs ~ normal(y_int, noise_sig);
}