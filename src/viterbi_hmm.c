void argmax_geno(int n_ind, int n_pos, int n_gen, int *geno,
                 double *rf, double *rf2,
                 double error_prob, int *argmax,
                 double initf(int, int *),
                 double emitf(int, int, double, int *),
                 double stepf(int, int, double, double, int *))
{
    int i, j, v, v2;
    double s, t, *gamma, *tempgamma, *tempgamma2;
    int **Geno, **Argmax, **traceback;
    int cross_scheme[2];
    
    /* cross scheme hidden in argmax argument; used by hmm_bcsft */
        cross_scheme[0] = argmax[0];
    cross_scheme[1] = argmax[1];
    argmax[0] = geno[0];
    argmax[1] = geno[1];
    
    /* Read R's random seed */
    /* in the case of multiple "most likely" genotype sequences,
    we pick from them at random */
    GetRNGstate();
    
    /* allocate space and
    reorganize geno and argmax */
    reorg_geno(n_ind, n_pos, geno, &Geno);
    reorg_geno(n_ind, n_pos, argmax, &Argmax);
    allocate_imatrix(n_pos, n_gen, &traceback);
    allocate_double(n_gen, &gamma);
    allocate_double(n_gen, &tempgamma);
    allocate_double(n_gen, &tempgamma2);
    
    for(i=0; i<n_ind; i++) { /* i = individual */
    
    R_CheckUserInterrupt(); /* check for ^C */
    
    /* begin viterbi algorithm */
    if(n_pos > 1) { /* multiple markers */
    for(v=0; v<n_gen; v++)
    gamma[v] = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
    
    for(j=0; j<n_pos-1; j++) {
    for(v=0; v<n_gen; v++) {
    tempgamma[v] = s = gamma[0] + stepf(1, v+1, rf[j], rf2[j], cross_scheme);
    traceback[j][v] = 0;
    
    for(v2=1; v2<n_gen; v2++) {
    t = gamma[v2] + stepf(v2+1, v+1, rf[j], rf2[j], cross_scheme);
    if(t > s || (fabs(t-s) < TOL && unif_rand() < 0.5)) {
    tempgamma[v] = s = t;
    traceback[j][v] = v2;
    }
    }
    tempgamma2[v] = tempgamma[v] + emitf(Geno[j+1][i], v+1, error_prob, cross_scheme);
    }
    for(v=0; v<n_gen; v++) gamma[v] = tempgamma2[v];
    }
    
    /* finish off viterbi and then traceback to get most
    likely sequence of genotypes */
    Argmax[n_pos-1][i] = 0;
    s = gamma[0];
    for(v=1; v<n_gen; v++) {
    if(gamma[v] > s || (fabs(gamma[v]-s) < TOL &&
    unif_rand() < 0.5)) {
    s = gamma[v];
    Argmax[n_pos-1][i] = v;
    }
    }
    for(j=n_pos-2; j >= 0; j--)
    Argmax[j][i] = traceback[j][Argmax[j+1][i]];
    }
    else {  /* for exactly one marker */
    s = initf(1, cross_scheme) + emitf(Geno[0][i], 1, error_prob, cross_scheme);
    Argmax[0][i] = 0;
    for(v=1; v<n_gen; v++) {
    t = initf(v+1, cross_scheme) + emitf(Geno[0][i], v+1, error_prob, cross_scheme);
    if(t > s || (fabs(t-s) < TOL && unif_rand() < 0.5)) {
    s = t;
    Argmax[0][i] = v;
    }
    }
    }
    
    /* code genotypes as 1, 2, ... */
    for(j=0; j<n_pos; j++) Argmax[j][i]++;
    
    } /* loop over individuals */
    
    
    /* write R's random seed */
    PutRNGstate();
}
