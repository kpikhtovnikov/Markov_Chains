set.seed(42) # фиксируем генератор сл. чисел

# шаг цепи
mc_sim=function(s, pr_now, pr_mat){
  # s -- множество состояний,
  # pr_now -- текущий вектор распределений, pr_mat -- матрица перехода
  # выход: s -- множество состояний, s_next -- куда перешла цепь,
  # pr_next -- новое распределение, pr_mat -- матрица перехода
  pr_next=pr_now %*% pr_mat
  id=which(rmultinom(1, 1, pr_next) == 1)
  s_next=s[id]
  list(s=s,s_next=s_next,pr_next=pr_next,pr_mat=pr_mat)
}


s=c(1,2,3,4,5) # задаем множество состояний

pr_mat=matrix(c(0.1,0.3,0.1,0.5, 0,  
                0.4,0.2,0, 0.3, 0.1,    
                0.3,0, 0.1, 0.4, 0.2,  
                0.3,0.2, 0, 0.3,0.2,
                0,0.2,0.4,0.3,0.1), nrow=5, ncol=5)

pr_mat=t(pr_mat) # матрица перехода, транспонирую,т.к. ввел неправильно, вам не надо
# симуляция цепи
run_mc_sim=function(s,s_now, pr_now, pr_mat, n_iter){
  # n_iter -- число итераций, длинна цепи, pr_mat -- матрица перехода,
  # s -- множество состояний, s_now, pr_now -- начальные данные
  # выход: mc -- марковская цепь (траектория)
  # prob -- вектор вероятностей на последней итерации (при большом n_iter будет похож
  # на стационарное распределение, пригодится, если захотим продолжить итерации)
  pr=pr_now
  mc=rep(s_now,n_iter)
  pr_hist=matrix(NA, nrow=n_iter,ncol=length(pr_now))
  pr_hist[1,]=pr_now
  for (i in 2:n_iter) {
    c=mc_sim(s, pr, pr_mat)
    pr=c$pr_next
    mc[i]=c$s_next
    pr_hist[i,]=pr
  }
  list(mc=mc,prob=c$pr_next,pr_hist)
} 



# симуляция: запускаем цепи из каждого состояния
n_iter=10
m=length(s)

mc_m=matrix(NA, nrow=n_iter, ncol=m) # столбцы -- траектории цепей
freq_tab=matrix(NA,nrow=m, ncol=m) # частоты состояний для каждой траектории

for (i in 1:m) {
dum=rep(0,m)
dum[i]=1
mc1=run_mc_sim(s, i, dum, pr_mat, n_iter)  
mc_m[,i]=mc1$mc
c=table(mc1$mc)
idx=as.integer(row.names(c))
freq_tab[i,]=rep(0,m)
freq_tab[i,idx]=table(mc1$mc)/n_iter
}

# график
# анимация
#library(animation)
#oopt = ani.options(interval = 0.3, nmax = n_iter)
#for (j in 1:ani.options("nmax")) {
#  matplot(mc_m[1:j,], type='b',pch=16, lty=1, col=1:5, ylim=c(1,5), ylab='state', xlab='time')
#  ani.pause()
#}
# график

matplot(mc_m, type='b',pch=16, lty=1, col=1:5, ylim=c(1,5), ylab='state', xlab='time')
abline(h=1, lty=3)
abline(h=2, lty=3)
abline(h=3, lty=3)
abline(h=4, lty=3)
abline(h=5, lty=3)


# таблица частот содержится в freq_tab
# следим за поведением распределений
mc_sta=run_mc_sim(s, i, mc1$prob, pr_mat, n_iter) # что будет, когда стартуем
# из стационарного распределния
prob_hist_stat=as.matrix(mc_sta[[3]]) # здесь записано
# против того, что будет, например, когда стартуем с состояния 1
mc_sta_1=run_mc_sim(s, 1, c(1,0,0,0,0), pr_mat, n_iter)
prob_hist_1=as.matrix(mc_sta_1[[3]])
# смотрим на вектор финальных вероятностей из библиотеки markovchain
library(markovchain)
mc_ttt=new('markovchain', states=c('1','2','3','4','5'), transitionMatrix=pr_mat,name="dummy")
steadyStates(mc_ttt)
summary(mc_ttt)
is.irreducible(mc_ttt)
plot(mc_ttt)

st=steadyStates(mc_ttt)
mc_sta=run_mc_sim(s, 1,st, pr_mat, n_iter)
hist_1=mc_sta[[3]]


mc_sta2=run_mc_sim(s, 1, c(0,1,0,0,0) , pr_mat, n_iter)
hist_2=mc_sta2[[3]]

mc_sta3=run_mc_sim(s, 1, c(0,0,1,0,0), pr_mat, n_iter)
hist_3=mc_sta3[[3]]

mc_sta4=run_mc_sim(s, 1, c(0,0,0,1,0), pr_mat, n_iter)
hist_4=mc_sta4[[3]]

mc_sta5=run_mc_sim(s, 1, c(0,0,0,0,1), pr_mat, n_iter)
hist_5=mc_sta5[[3]]