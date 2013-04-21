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
bootstrap_v <- function (distr, theta, n, M, J=1) {
	result=matrix(nrow=J, ncol=6, 
			dimnames=list(NULL,c("Vn_boot","Vn","mean_X_star","mean_X_n","pvalue_boot", "pvalue_n")));
	for (j in 1:J) {
		X_n = distr(n, theta);
		X_star_bar = array(dim = M);
		for (i in 1:M){
			X_star_bar[i] = mean(sample(X_n, n, replace=TRUE));
		}
		result[j,]=cbind(variance(X_star_bar), 
					variance(X_n)/n, 
					mean(X_star_bar), 
					mean(X_n),
					2*pnorm(-abs(theta-mean(X_n)), sd=sqrt(variance(X_star_bar))),
					2*pnorm(-abs(theta-mean(X_n)), sd=sqrt(variance(X_n)/n)));
	}
	return(result);
}

# Question 3
# Plot longueur de l'intervalle (2*sqrt(V)*Phi(1-alpha/2)) sur n
n_list=seq(from=5,to=50,by=5);
alpha_list=seq(from=0.05,to=0.15,by=0.05);
distr=gunif;

longueur=matrix(nrow=0,ncol=4,dimnames=list(NULL,c("n", "alpha", "longueur_n", "longueur_boot")));
for (n in n_list) {
	boot=bootstrap_v(distr, theta=0.5, n=n, M=5*n, J=100);
	for (alpha in alpha_list) {
		if(sum((boot[,"pvalue_n"]>alpha)*(boot[,"pvalue_boot"]>alpha))>2){ # Considerer seulement si theta est dans l'intervalle.
			rows=cbind(n, alpha, sqrt(boot[boot[,"pvalue_n"]>alpha&boot[,"pvalue_boot"]>alpha,c("Vn","Vn_boot")])*2*qnorm(1-alpha/2));
			longueur=rbind(longueur,rows);
		}
	}
}
a=data.frame(longueur);
a$ratio=a$longueur_n/a$longueur_boot;
boxplot(data=a, subset=alpha==0.05, ratio~n, xlab="n", ylab="Ratio I_n/I_boot", main=paste("alpha =", 0.05)); abline(h=1,col="red",lwd=1.5);
boxplot(data=a, subset=alpha==0.10, ratio~n, xlab="n", ylab="Ratio I_n/I_boot", main=paste("alpha =", 0.10)); abline(h=1,col="red",lwd=1.5);
boxplot(data=a, subset=alpha==0.15, ratio~n, xlab="n", ylab="Ratio I_n/I_boot", main=paste("alpha =", 0.15)); abline(h=1,col="red",lwd=1.5);


# Question 4
# Fonction qui calcule les intervalles empiriques bootstrap et gaussien.
bootstrap_q <- function (distr, theta, n, M, alpha) {
	X_n = distr(n, theta);
	X_star_bar = array(dim = M);
	for (i in 1:M){
		X_star_bar[i] = mean(sample(X_n, n, replace=TRUE));
	}	

	F_star = ecdf(X_star_bar); # Fonction répartition empirique
	#F_norm <- function(x) {return(pnorm(x,mean=mean(X_n),sd=sqrt(variance(X_n)/length(X_n))))};

	F_star_inv <- function(alpha) {return (uniroot(function(x) F_star(x) - alpha, lower=min(X_star_bar)-0.5, upper=max(X_star_bar)+0.5)[1]$root) };
	I_star = F_star_inv (1-alpha/2)-F_star_inv(alpha/2);
	F_norm_inv <- function(x) qnorm(x, mean=mean(X_n), sd=sqrt(variance(X_n)/length(X_n)));
	I_norm = F_norm_inv(1-alpha/2)-F_norm_inv(alpha/2);
	

	#x=seq(min(X_star_bar),max(X_star_bar),length=max(length(X_star_bar),1000));
	#plot(x,F_star(x), type="n");
	#lines(x,pnorm(x,mean=mean(X_n),sd=sqrt(variance(X_n)/length(X_n))),type="S", col="blue");
	#lines(x,F_star(x), type="S", col=34)
	return(cbind(I_norm,I_star));
}

n_list=seq(from=5,to=35,by=3);
alpha_list=seq(from=0.05,to=0.15,by=0.05);
distr=gunif;
longueur=matrix(nrow=0,ncol=4,dimnames=list(NULL,c("n", "alpha", "longueur_n", "longueur_boot")));
for (n in n_list) {
	for (j in 1:100) {
		for (alpha in alpha_list) {
			boot=bootstrap_q(distr, theta=0.5, n=n, M=10*n, alpha);
			rows=cbind(n, alpha, boot);
			longueur=rbind(longueur,rows);
		}	
	}
}

a=data.frame(longueur);
a$ratio=a$longueur_n/a$longueur_boot;
boxplot(data=a, subset=alpha==0.05, ratio~n, xlab="n", ylab="Ratio I_n/I_boot", main=paste("alpha =", 0.05)); abline(h=1,col="red",lwd=1.5);
boxplot(data=a, subset=alpha==0.10, ratio~n, xlab="n", ylab="Ratio I_n/I_boot", main=paste("alpha =", 0.10)); abline(h=1,col="red",lwd=1.5);
boxplot(data=a, subset=alpha==0.15, ratio~n, xlab="n", ylab="Ratio I_n/I_boot", main=paste("alpha =", 0.15)); abline(h=1,col="red",lwd=1.5);

