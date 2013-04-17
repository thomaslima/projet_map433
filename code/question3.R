# Paquets nécessaires
library(Rlab) #rbern

# Distributions exemple
gunif <- function (n, theta, mean = 0, sd = 1) {
	L = sqrt(3)*sd;
	return (runif(n, -L + theta + mean, L + theta + mean));
}

gnorm <- function (n, theta, mean = 0, sd = 1) {
	return (rnorm(n, mean + theta, sd))
}

gexp <- function (n, theta, mean = 0, sd = 1) {
	rate = 1/sd;
	return (rexp(n, rate) - 1/rate + mean + theta);
}

gbernoulli <- function (n, theta, p = 0.5, sd = 1) {
	factor = sd/sqrt(p*(1-p));
	return ((rbern(n, p) - p) * factor + theta);
}

# Algorithme Bootstrap
bootstrap <- function (n, theta, distr, M) {
	X_n = distr(n, theta);
	X_star_bar = array(dim = M);
	for (i in 1:M){
		X_star_bar[i] = mean(sample(X_n, M, replace=TRUE));
	}
	hist(X_star_bar);
	print(paste("Vn(boot)=", var(X_star_bar)));
	print(paste("Vn=", var(X_n)));
	print(paste("mean(X_star)=", mean(X_star_bar)));
	print(paste("mean(X_n)=", mean(X_n)));
}

bootstrap(100, 0.1, gbernoulli, 10000);