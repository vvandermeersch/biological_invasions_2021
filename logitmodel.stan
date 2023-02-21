data {
      int <lower=1> N;
      int <lower=1> P; //nombre de parcelles
      int <lower=1> S; //nombre de sols
      int <lower=1, upper=S> soil[N]; //soil
      vector[N] xedge; //covar edge distance
      vector[N] xsd; //dist seed source
      vector[N] fire; //covar fire
      vector[N] pba_remov; //ba removed after treatment
      int <lower=1, upper=P> plot[N]; //plot
      int <lower=0, upper=1> y[N]; //outcome
      }
      parameters {
      real theta_e; //edge effect
      real theta_sd; //seed source effect
      real theta_f; //1983 fire effect
      real theta_trmt; //treatment effect
      real theta_soil[S]; //soil effect
      //plot covariance
      real p[P]; 
      real <lower=0> sigma; 
      }
      model {
      vector[N] theta;
      for (n in 1:N)
        theta[n]=theta_e*xedge[n]+theta_sd*xsd[n]+theta_soil[soil[n]]+theta_f*fire[n]+theta_trmt*pba_remov[n]+p[plot[n]];
      y ~ bernoulli_logit(theta);
      p ~ normal(0, sigma);
      sigma ~ normal(0,1);
      }

