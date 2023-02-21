data {
      int <lower=1> Q; //nombre de quadrat
      int <lower=1> P; //nombre de parcelles
      vector[Q] p_ced; //proportion ced
      vector[Q] p_fire; //proportion fire
      vector[Q] edge_dist; //edge distance
      vector[Q] sd_dist; //seedsource distance
      vector[Q] pba_remov; //proportion of ba removed after treatment
      vector[Q] psoil1;
      vector[Q] psoil2;
      int <lower=1, upper=P> plot[Q]; //plot
      vector[Q] y; //simpson index
      }
      parameters {
      real theta_ced; //cedrela effect on diversity loss
      real theta_fire; //fire effect
      real theta_edge; //edge effect
      real theta_sd; //seedsource effect
      real theta_trmt; //treatment effect
      real theta_soil1;
      real theta_soil2;
      real <lower=0 > sigma;
      real <lower=0, upper=100> p[P];
      real <lower=0> sigmap;
      real mup;
      }
      model {
      real mu[Q];
      for (k in 1:Q)
        mu[k]=log((theta_fire*p_fire[k]+theta_edge*edge_dist[k]+theta_sd*sd_dist[k]+theta_trmt*pba_remov[k]+theta_soil1*psoil1[k]+theta_soil2*psoil2[k]+
              p[plot[k]])*exp(-(theta_ced)*p_ced[k]));
      y ~ lognormal(mu, sigma);
      p ~ lognormal(log(mup), sigmap);
      sigma ~ normal(0,1); 
      sigmap ~ normal(0,1); 
      }
