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

# Fonction Variance
variance <- function (x)
	return(sum((x-mean(x))^2)/length(x));

# Algorithme Bootstrap
bootstrap <- function (n, theta, distr, M) {
	X_n = distr(n, theta);
	X_star_bar = array(dim = M);
	for (i in 1:M){
		X_star_bar[i] = mean(sample(X_n, n, replace=TRUE));
	}
	hist(X_star_bar);
	print(paste("Vn(boot)=", var(X_star_bar)));
	print(paste("Vn/n=", var(X_n)/n));
	print(paste("mean(X_star)=", mean(X_star_bar)));
	print(paste("mean(X_n)=", mean(X_n)));
}

bootstrap(100, 0.1, gbernoulli, 10000);


# Question 4

bootstrap_q <- function(n, theta, distr, M) {
	X_n = distr(n, theta);
	X_star_bar = array(dim = M);
	for (i in 1:M){
		X_star_bar[i] = mean(sample(X_n, M, replace=TRUE));
	}
	x = seq(-1,1,length=500);
	F_star = array(dim=length(x));
	for (i in 1:length(x)){
		F_star[i] = sum((X_star_bar <= x[i]) * 1);
	}
	plot(x,F_star, type="n");
	lines(x,F_star,type="S");
}

bootstrap_q(100,0,gunif,100)
