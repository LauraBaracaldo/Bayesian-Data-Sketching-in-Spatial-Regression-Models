// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>





using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
colvec rmvnormCPP(const colvec& mu, const mat& sigma){
    RNGScope scope;
    return mu + chol(sigma, "lower")*vec(rnorm(mu.n_rows));
}

// [[Rcpp::export]]
colvec uvnGausLogLik(const colvec& data, const colvec& mu, const colvec& sigma2){
	return -0.5*log(2*datum::pi*sigma2) - 0.5*square(data - mu)/sigma2;
}

// [[Rcpp::export]]
colvec Gibbsrho(const colvec& Z, const mat& PhiX, const colvec& Beta,  mat Hinv,  double logdetH,  double tau2,  double sig2, const double& atau, const double& btau, const double& asig, const double& bsig){
	colvec ZPhiXB = Z-PhiX*Beta;
	colvec ZHZ =  ZPhiXB.t()*Hinv*(ZPhiXB);
	return -0.5* logdetH  - 0.5*ZHZ -(atau+1)*log(tau2) - btau/tau2 -(asig+1)*log(sig2) -bsig/sig2 +log(tau2) + log(sig2);
	}

// [[Rcpp::export]]
mat expoRCPP(const mat& X){
	return exp(X);
}

// [[Rcpp::export]]
mat invRCPP(const mat& X){
	return inv_sympd(X);
}

// [[Rcpp::export]]
mat PseuinvRCPP(const mat& X){
	return pinv(X);
}

// [[Rcpp::export]]
double logdetRCPP(const mat& X){
	return sum(log_det_sympd(X));
}

// [[Rcpp::export]]
mat kronCPP(const mat& X, const mat& Y){
	return kron(X,Y);
}



// [[Rcpp::export]]
Col<int> sampleLabel(const colvec& data, const mat& forecastMean, const colvec& weight, const colvec& sigma2){
	int dataSize = data.n_rows, clusterNum = weight.n_rows;
	Col<int> labelValues = linspace<Col<int>>(1,clusterNum,clusterNum), labels(dataSize), labelValue;
	colvec distribution_i(clusterNum);
	for (int i=0; i<dataSize; ++i){
		colvec mu = trans(forecastMean.row(i));
		colvec repData_i(clusterNum);
		repData_i.fill(data(i));
		colvec loglik_i = uvnGausLogLik(repData_i, mu, sigma2);
		distribution_i = exp(log(weight) + loglik_i);
		distribution_i = distribution_i / sum(distribution_i);
		labelValue = Rcpp::RcppArmadillo::sample(labelValues, 1, 1, distribution_i);
		labels(i) = labelValue(0);
	}
	return labels;
}



// [[Rcpp::export]]
List sampleBeta(const colvec& data, const mat& forecast, const colvec& sigma2,
		            const Col<int>& labels, const Col<int>& label_geq2, const mat& priorMu, const colvec& priorSigma2){
	int clusterNum = sigma2.n_rows;
	mat beta(2, clusterNum);
	for (int j=0; j<clusterNum; ++j){
		if (any(label_geq2 == j+1)){
			uvec index = find(labels == j+1);
			colvec data_j = data.rows(index), forecast_j = forecast.col(j);
			mat X_j(data_j.n_rows, 2);
			X_j.col(0).fill(1);
			X_j.col(1) = forecast_j.rows(index);
			mat invSigma = inv_sympd(X_j.t()*X_j + eye(2,2));
			beta.col(j) = rmvnormCPP(invSigma*(X_j.t()*data_j+priorMu.col(j)), sigma2(j)*invSigma);
	    } else {
				mat I = zeros(2,2);
				I.diag() += priorSigma2(j);
		    beta.col(j) = rmvnormCPP(priorMu.col(j), I);
		  }
	}
    mat forecastMean = forecast.each_row() % beta.row(1);
	forecastMean.each_row() += beta.row(0); 
	return List::create(Named("beta")=beta, Named("mu")=forecastMean);
}


// [[Rcpp::export]]

colvec SampleBetaRcpp(const mat& PhiX, const colvec& Z, const mat& Hinv){
	int p = PhiX.n_cols;
	colvec beta(p);
	mat Sigma = inv_sympd(PhiX.t()*Hinv*PhiX + (1/1000)*eye(p,p));
	colvec mu = Sigma*PhiX.t()*Hinv*Z;
	beta = rmvnormCPP(mu, Sigma);
	return(beta);
}

// [[Rcpp::export]]

List trialfnRCPP(const mat& PhiX, const colvec& Z, colvec thetacur, const colvec& Beta,  mat Hinvcur,  double logdetHcur, const mat& M, 
		            const double& atau, const double& btau, const double& asig, const double& bsig){
	int m = Z.n_elem, p = thetacur.n_elem;
	double tau2cur = thetacur[1], sig2cur = thetacur[2];
		return List::create(Named("tau2cur")=tau2cur, Named("m")=m, Named("p")=p, Named("sig2cur")=sig2cur);
					}

// [[Rcpp::export]]
List Sampletau2sig2Rcpp(const mat& PhiX, const colvec& Z, colvec thetacur, const colvec& Beta,  mat Hinvcur,  double logdetHcur, const mat& M, 
		            const double& atau, const double& btau, const double& asig, const double& bsig){
	int m = Z.n_elem, p = thetacur.n_elem;
	double tau2cur = thetacur[1], sig2cur = thetacur[2];
	colvec thetapro = rmvnormCPP(log(thetacur), 0.1*eye(p,p)) ;
	double tau2pro = exp(thetapro[1]), sig2pro = exp(thetapro[2]);
	mat Hpro = sig2pro*(M) + tau2pro*eye(m,m);
	mat Hinvpro = inv_sympd(Hpro);
	double logdetHpro = log_det_sympd(Hpro);
	colvec logrhopro = Gibbsrho(Z,PhiX,  Beta,  Hinvpro,  logdetHpro, tau2pro, sig2pro,  atau,  btau,  asig,  bsig), logrhocur = Gibbsrho(Z,PhiX,  Beta,  Hinvcur,  logdetHcur, tau2cur, sig2cur,  atau,  btau,  asig,  bsig);
    double lr = sum(log(randu(1)))- sum(logrhopro - logrhocur);
	if( lr<0){
		thetacur = exp(thetapro);
		Hinvcur = Hinvpro;
		logdetHcur = logdetHpro;
		}
	return List::create(Named("thetacur")=thetacur, Named("Hinvcur")=Hinvcur, Named("logdetHcur")=logdetHcur);
}


