using namespace Rcpp ;
using std::map ;
using std::pair ;
using std::vector ;
using std::make_pair ;
using std::sort ;
using std::cout ;
using std::endl ;

Function fun(f) ; // function called to calculate statistic
NumericMatrix data(d) ; // NxM data matrix
NumericVector obs(o) ; // N observed statistics
List pars(cntl) ; // R list with parameters
unsigned int i, j, k, N = obs.length(), M = data.ncol() ;
unsigned int B = as<unsigned int>(pars["B"]), 
	B0 = as<unsigned int>(pars["B0"]), 
	K = as<unsigned int>(pars["K"]), 
	tail = as<bool>(pars["lower"]) , // variables passed from List pars
	verbose = as<bool>(pars["verbose"]) ; // variables passed from List pars
double p0 = as<double>(pars["p0"]) ;
double *a, *p, *risk, *orig, *dest, *weights, pest, 
	stat, R, t, Tot = B * N, cs, rand ;
vector<unsigned int> numBoots = vector<unsigned int>(N) ; // bootstrap allocation
vector<double> bootStraps = vector<double>(M), 
	pval = vector<double>(N), vec,
	finalRisk = vector<double>(N) ;
vector< pair<double, int> > sorted ;
map<int, int> ntot, nstep ;
map<double, int> cumsum ;

// initialize GSL RNG
gsl_rng_env_setup() ;
gsl_rng * r = gsl_rng_alloc(gsl_rng_taus) ;

// initialize variables
a = new double[N] ; // # min(resamples > test statistic, resamples < test statistic)
p = new double[N] ; // p = a/B0
risk = new double[N] ; // r = IB(p0, a+1, n-a+1)
weights = new double[N] ; // r = IB(p0, a+1, n-a+1)
orig = new double[M] ;
dest = new double[M] ;
sorted = vector< pair<double, int> >(N) ;
for(i = 0 ; i < N ; i++)
	a[i] = p[i] = risk[i] = weights[i] = 0.0 ;

// burn-in (N * B0 resamples)
for(i = 0 ; i < N ; i++) {
	vec = as< vector<double> >(data(i,_)) ;
	for(j = 0 ; j < B0 ; j++) {
		for(k = 0 ; k < M ; k++) 
			orig[k] = vec[k] ;
		// resample data, and wrap it as SEXP
		gsl_ran_sample(r, dest, M, orig, M, sizeof(double)) ;
		for(k = 0 ; k < M ; k++) 
			bootStraps[k] = dest[k] ;
		// pass it to function, expecting vector of doubles
		stat = as<double>(fun(wrap(bootStraps))) ;
		if(stat >= obs[i]) a[i] += 1.0 ; 
	}
	ntot[i] = B0 ;
	pest = a[i]/((double) B0) ;
	switch(tail) {
		case 0:
			pval[i] = pest ;
		break ;
		case 1:
			pval[i] = pest > (1 - pest) ? 2 * (1 - pest) : 2 * pest ;
		break ;
		default: 
			exit(1) ;
	}
}
cout << "Completed." << endl ;
t = B0 * N ;
cout << "Beginning adaptive phase." << endl ;
//adaptive phase of algorithm
while(t < Tot) {
	//calculate risk for each bootstrap
	for(i = 0 ; i < N ; i++) {
		risk[i] = ntot[i] - a[i] > a[i] ? 1 - gsl_cdf_beta_P(p0, a[i] + 1, ntot[i] - a[i] + 1) :
			1 - gsl_cdf_beta_P(p0, ntot[i] - a[i] + 1, a[i] + 1) ;
		if(pval[i] >= p0) risk[i] = 1 - risk[i] ;
		R += risk[i] ;
	}
	for(i = 0 ; i < N ; i++) weights[i] = risk[i]/R ;
	//sort values by normalized 
	sorted.clear() ;
	for(i = 0 ; i < N ; i++) sorted.push_back(make_pair(weights[i], i)) ;
	sort(sorted.begin(), sorted.end()) ;
	cs = 0.0 ;
	for(i = 0 ; i < N ; i++)  {
		//cout << cs << " " << sorted[i].first << " " << sorted[i].second << endl ;
		cs += sorted[i].first ;
		cumsum[cs] = sorted[i].second ;
	}
	//allocate N*K samples proportional to weights
	nstep.clear() ;
	for(i = 0 ; i < (N * K) ; i++) {
		rand = gsl_rng_uniform(r) ;
		j = cumsum.upper_bound(rand)->second ;
		nstep[j]++ ;
		ntot[j]++ ;
	}
	t += i ;
	for(i = 0 ; i < N ; i++) {
		vec = as< vector<double> >(data(i,_)) ;
		for(j = 0 ; j < nstep[i] ; j++) {
			for(k = 0 ; k < M ; k++) orig[k] = vec[k] ;
			gsl_ran_sample(r, dest, M, orig, M, sizeof(double)) ;
			for(k = 0 ; k < M ; k++) bootStraps[k] = dest[k] ;
			stat = as<double>(fun(wrap(bootStraps))) ;
			if(stat >= obs[i]) a[i] += 1.0 ;
			//resample data, and wrap it as SEXP
			//pass it to function, expecting double
		}
		pest = a[i]/((double) ntot[i]) ;
		switch(tail) {
			case 0:
				pval[i] = pest ;
			break ;
			case 1:
				pval[i] = pest > (1 - pest) ? 2 * (1 - pest) : 2 * pest ;
			break ;
			default: 
				exit(1) ;
		}
	}
	cout << "Allocated " << t << " resamples." << endl ;
}
cout << "Finished." << endl ;
for(i = 0 ; i < N ; i++) {
	numBoots[i] = ntot[i] ;
	finalRisk[i] = weights[i] ;
}
return(wrap(List::create(Named("p", pval), Named("nboot", numBoots), 
	Named("risk", finalRisk)))) ;
